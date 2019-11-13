/*
** Written by Oliver Frolovs, based on the examples by
** Copyright (c) 2014, David A.W. Barton
** (david.barton@bristol.ac.uk) All rights reserved.
**
** "A FOOLISH CONSISTENCY IS THE HOBGOBLIN OF LITTLE MINDS..." ~ R.W. Emerson
**
** Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**
** 1. Redistributions of source code must retain the above copyright
** notice, this list of conditions and the following disclaimer.
**
** 2. Redistributions in binary form must reproduce the above copyright
** notice, this list of conditions and the following disclaimer in the
** documentation and/or other materials provided with the distribution.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* ************************************************************************ */
/* * Includes, defines and global variables ******************************* */
/* ************************************************************************ */
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "rtc_data.h"
#include "rtc_main.h"
#include "rtc_util.h"
#include "rtc_user.h"
#include "neon_mathfun.h"

/* ************************************************************************ */
/* * Defines ************************************************************** */
/* ************************************************************************ */

// Number of lasers in the system
#define NUMBER_OF_LASERS 3

// Maths constants for convenience
#define M_2PI (2.0f * M_PI)

// Maxon motor controller 50/5 control voltage levels as per its datasheet
#define MAXON_ON   3.7f
#define MAXON_OFF  0.0f

/// Number of sin/cos values to precompute (???)
#define N_FOURIER_MODES 7  // N_FOURIER_MODES + 1 must be a multiple of 4
#define N_FOURIER_COEFF (2*N_FOURIER_MODES + 2)  // one will always be zero but it's easier/quicker than handling the special case

// Encoder data filter
#define INPUT_FILTER_ORDER 2
#define INPUT_FILTER_N_STATE (2*INPUT_FILTER_ORDER)

// Size of the speed history log. Used for enforcing the speed safety limit.
// Assumes the sample frequency is divisible by 20, which it should be.
#define SPEED_HISTORY_SIZE (TIMER_FREQ / 20)

// The motor cannot drive arbitrary slow speeds due to cogging and friction,
// so a realistic cut-off well below the speed range of interest is selected,
// for speed control purposes. Any value less than the cutoff value is considered
// to be "at rest".
#define RPM_AT_REST_CUTOFF 10

// Rig configuration and runtime errors are reported by setting corresponding bits in
// a 32-bit unsigned integer named status_flags. This is faster than using a whole integer
// for *each* flag and also allows to save them all easily. But this does require a certain
// support infrastructure to use from Matlab easily.
//
// !!! Never change the existing values as the data might have been saved using it.
//     Use the next available value for the new flags.
//
#define RIG_STATUS_OK                         0u
#define RIG_STATUS_UNKNOWN_ERROR              (0x1u << 0)
#define RIG_STATUS_SETUP_NOT_COMPLETED        (0x1u << 1)
#define RIG_STATUS_PID_NUMERIC_ERROR          (0x1u << 2)
#define RIG_STATUS_SPEED_SAFETY_LIMIT_REACHED (0x1u << 3)
#define RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MIN  (0x1u << 4)
#define RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MAX  (0x1u << 5)

/* ********************************************************************** */
/* * Types ************************************************************** */
/* ********************************************************************** */

struct biquad_filter_t {
	int order;
	float a[2];
	float b[3];
};

/* ********************************************************************** */
/* * Globals (internal) ************************************************* */
/* ********************************************************************** */

// Fourier calculations
//static float sinusoid_f[N_FOURIER_COEFF] MEM_ALIGN;  /* [sin(0*t), sin(1*t), sin(2*t), ..., cos(0*t), cos(1*t), cos(2*t), ...] */

// TIME-KEEPING is done both linearly (time) and (mod_2pi). Period start signals
// the start of a new period, that is it is related to (mod_2pi). Linear time
// just ticks along, increasing with each call to rtc_user_main() functions.
static uint32_t time;
static float time_mod_2pi ;
static float time_delta = M_2PI/TIMER_FREQ;
static uint32_t period_start = 1;  /* signal the start of a new period */

// Forcing frequency [Hz]. This parameter is accessible from Matlab
static float forcing_freq_hz = 10.0f;

// Forcing frequency (angular). This is computed by the _hz setter.
static float forcing_freq;

// The error and warning flags. Each bits of this variable corresponds to a single flag. 
// The list of available flags is above in definitions RIG_FLAG_... Zero bit value means
// the flag has not been set.
//
// Some variables such as motor voltage limits, PID coefficients, and possibly laser/encoder
// transfer function coefficients may not be set in firmware and require configuration from
// Matlab before the rig is run. Thus, on booting the controller, wait for the operator to
// explicitly declare that the rig has been set up and ready to go.
static uint32_t status_flags = RIG_STATUS_SETUP_NOT_COMPLETED;

// CONTROL TARGET: angular velocity to maintain, in rpm
static float rpm;

// CONTROL TARGET: angular velocity to maintain in rad/s. This value is computed
// by the update trigger of rpm variable.
static float target_speed;
 
// TODO digital signals to the motor become digital connections to free up analogue output:

// Motor control: 
//   * rotation sense and enable signal on one output channel; and 
//   * set (current) value on another output channel.
//
// This duplicates the settings in the ESCON controller.
static uint32_t output_channel_enable_motor    = 1; /* Output U1 (unipolar) */
static uint32_t output_channel_motor_set_level = 0; /* Output B0 (bipolar) */
static uint32_t output_channel_motor_direction = 3; /* Output U3 (unipolar) */

// Motor: set voltage level. This must be within [motor_min_voltage, motor_max_voltage]
static float motor_voltage_level;

// Motor: voltage limits set in the motor controller
static float motor_min_voltage;
static float motor_max_voltage;

// Motor PID controller coefficients
static float K_p;
static float K_i;
static float K_d;

// Motor PID controller error terms
static float pid_error_raw;
static float pid_error_filtered;
static float pid_error_before;
static float pid_error_derivative;
static float pid_error_integral;

// Lasers: distance to the object in mm
static float lasers_distance_mm[NUMBER_OF_LASERS];

// Lasers: distance to the reflector at rest. Every experiment should make note of this!
//         the default value is NAN for easy spotting of cases when the Matlab script
//         has not taken note of the distance at rest.
static float lasers_distance_at_rest[NUMBER_OF_LASERS] = {NAN};

// Lasers: voltage limits. For Omron ZX1-300 lasers with 250 Ohm resistor in 
// the signal conditioner, the signal varies from 1 to 5 volts. 
static float lasers_near_field_limit_volt = 5.0f;
static float lasers_far_field_limit_volt  = 1.0f;

// Lasers: transfer function gradient and intercept have been computed assuming:
//   * Omron ZX1-300 lasers: near field -150 mm (20 mA), far field 150 mm (4 mA)
//   * all (three) lasers are identical
//   * the resistor value in signal conditioner is 250 [Ohm]
//
// The calculations are in the `rig/Computing Lasers Transfer Function.pdf` file.
static float lasers_transfer_f_gradient  =   75.0f;
static float lasers_transfer_f_intercept = -225.0f;

// Rotary Encoder: voltage as read, and its interpretation as angular position,
// as the proportion of 2*pi. The encoder's analogue signal is connected to input zero.
static uint32_t encoder_input_channel;
static float encoder_voltage;
static float encoder_angular_position;
static float encoder_angular_position_prev;
static float encoder_speed_rpm;
static float encoder_speed;     // in rad/s

// Rotary Encoder: transfer function (voltage to SPEED) gradient and intercept
static float encoder_transfer_f_gradient;
static float encoder_transfer_f_intercept;

// Set to enable the speed control, unset if safety feature is not desired
static uint32_t enable_speed_safety_limit = 1;

// The motor will turn off if this speed is reached in either direction of rotation (rpm)
static float speed_safety_limit_rpm = 600;
 
// Indicates that the safety limit has been reached
//static uint32_t speed_safety_limit_reached_flag;

// Buffer to collect last N speed readings for averaging. Used for 
// enforcing the speed safety limit. The buffer wraps over at the end, using
// the index variable. The speed values are in rpm.
static float speed_history[SPEED_HISTORY_SIZE];
 
// Speed history buffer index
static uint32_t speed_history_idx;
 
// Since the buffer is circular, it's impossible to tell from the index alone,
// if it has been filled to full capacity. Therefore, there is a flag to indicate so.
static uint32_t speed_history_full;

// The mean spead will be calculated by averaging the full set of values from 
// the speed_history array.
static float mean_speed_rpm;

// The maximum instantaneous speed ever recorded
static float max_speed_rpm;

// The mean_speed_rpm value that caused activation of the safety stop
static float safety_triggered_speed_rpm;

/******************************************************************************
 *
 * OPTIONAL FEATURE: FILTER PID ERROR TERM
 */

static uint32_t enable_pid_error_filter;     // Motor PID filter for control signal enable flag
static float pid_error_filter_freq = 0.025f; // fraction of sampling frequency
static struct biquad_filter_t input_filter;
static float pid_error_filter_state[INPUT_FILTER_N_STATE];  // filtering encoder output

// Requests. Set to nonzero to activate. The triggers on update will do the job.
static uint32_t request_reset_status_flags;
static uint32_t request_reset_pid_error_terms;

/* ************************************************************************ */
/* * Internal prototypes ************************************************** */
/* ************************************************************************ */

// Reset the whole system
void bbb_reboot(void *new_value, void *trigger_data);

// Warning and error flag functions
static uint32_t get_status_flags(uint32_t idx);
static void set_status_flags(uint32_t idx);
static void clear_status_flags(uint32_t idx);

// Triggers for various exported variables...
void update_forcing_freq(void *new_value, void *trigger_data);   // Trigger
void update_target_speed(void *new_value, void *trigger_data);   // Trigger
void update_input_filter(void *new_value, void *trigger_data);   // Trigger
void reset_status_flags(void *new_value, void *trigger_data);    // Trigger
void reset_pid_error_terms(void *new_value, void *trigger_data); // Trigger

void reset_speed_history(void);

float biquad_filter(float x, const struct biquad_filter_t *filter, float state[]);
void update_filter(struct biquad_filter_t *filter, float freq);

void lasers_compute_distance_from_voltage(void);
void encoder_compute_speed_from_voltage(void);

void shutdown_motor(void);

float sum(float list[], uint32_t n);
float rad2rpm(float rad);
float rpm2rad(float rpm);


/* ************************************************************************ */
/* * Function definitions ************************************************* */
/* ************************************************************************ */
void rtc_user_init(void)
{
	// TODO better name, can this be a function in rtc() object instead?
	rtc_data_add_par("bbb_reboot", &time, RTC_TYPE_UINT32, sizeof(time), bbb_reboot, NULL);

	// Tracking of time, linear and periodic counters are provided
	rtc_data_add_par("time",  &time, RTC_TYPE_UINT32, sizeof(time), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("time_mod_2pi", &time_mod_2pi, RTC_TYPE_FLOAT, sizeof(time_mod_2pi), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("forcing_freq", &forcing_freq, RTC_TYPE_FLOAT, sizeof(forcing_freq), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("forcing_freq_hz", &forcing_freq_hz, RTC_TYPE_FLOAT, sizeof(forcing_freq_hz), update_forcing_freq, NULL);

	// Update forcing frequency (rad/s) value
	update_forcing_freq((void *)&forcing_freq_hz, NULL);

	// Rig flags available as a single variable
	rtc_data_add_par("status_flags",  &status_flags, RTC_TYPE_UINT32, sizeof(status_flags), NULL, NULL);
	
	// Desired speed of rotation in RPM
	rtc_data_add_par("rpm", &rpm, RTC_TYPE_FLOAT, sizeof(rpm), update_target_speed, NULL);

	// Motor: enable and set value channel numbers
	rtc_data_add_par("output_channel_enable_motor",  &output_channel_enable_motor, RTC_TYPE_UINT32, sizeof(output_channel_enable_motor),  NULL, NULL);
	rtc_data_add_par("output_channel_motor_set_level", &output_channel_motor_set_level, RTC_TYPE_UINT32, sizeof(output_channel_motor_set_level), NULL, NULL);
	rtc_data_add_par("output_channel_motor_direction", &output_channel_motor_direction, RTC_TYPE_UINT32, sizeof(output_channel_motor_direction), NULL, NULL);

	// Motor: voltage input translates into current drawn according to the motor controller settings
	rtc_data_add_par("motor_voltage_level", &motor_voltage_level, RTC_TYPE_FLOAT, sizeof(motor_voltage_level), NULL, NULL);

	// Motor: the limits on input voltage are set in the motor controller settings
	rtc_data_add_par("motor_min_voltage", &motor_min_voltage, RTC_TYPE_FLOAT, sizeof(motor_min_voltage), NULL, NULL);
	rtc_data_add_par("motor_max_voltage", &motor_max_voltage, RTC_TYPE_FLOAT, sizeof(motor_max_voltage), NULL, NULL);

	// Motor PID controller coefficients. All [RW]
	rtc_data_add_par("K_p", &K_p, RTC_TYPE_FLOAT, sizeof(K_p), NULL, NULL);
	rtc_data_add_par("K_i", &K_i, RTC_TYPE_FLOAT, sizeof(K_i), NULL, NULL);
	rtc_data_add_par("K_d", &K_d, RTC_TYPE_FLOAT, sizeof(K_d), NULL, NULL);

	// Motor PID controller propertional error filter, enable flag and value in ???
	rtc_data_add_par("enable_pid_error_filter", &enable_pid_error_filter, RTC_TYPE_UINT32, sizeof(enable_pid_error_filter), NULL, NULL);
	rtc_data_add_par("pid_error_filter_freq", &pid_error_filter_freq, RTC_TYPE_FLOAT, sizeof(pid_error_filter_freq), update_input_filter, NULL);

	// PID error terms. All [RW]
	rtc_data_add_par("pid_error_raw", &pid_error_raw, RTC_TYPE_FLOAT, sizeof(pid_error_raw), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("pid_error_filtered", &pid_error_filtered, RTC_TYPE_FLOAT, sizeof(pid_error_filtered), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("pid_error_integral", &pid_error_integral, RTC_TYPE_FLOAT, sizeof(pid_error_integral), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("pid_error_derivative", &pid_error_derivative, RTC_TYPE_FLOAT, sizeof(pid_error_derivative), rtc_data_trigger_read_only, NULL);
	
	// OPTIONLA: Speed safety limit in [rpm], and status flag (sticky). Both [RW].
	// The limit is a positive number, but the speed safety feature would not look at the sense of rotation, only at the actual speed.
	rtc_data_add_par("enable_speed_safety_limit",  &enable_speed_safety_limit,  RTC_TYPE_UINT32, sizeof(enable_speed_safety_limit), NULL, NULL);
	rtc_data_add_par("speed_safety_limit_rpm", &speed_safety_limit_rpm, RTC_TYPE_FLOAT, sizeof(speed_safety_limit_rpm), NULL, NULL);
	rtc_data_add_par("speed_history", &speed_history, RTC_TYPE_FLOAT, sizeof(speed_history), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("mean_speed_rpm", &mean_speed_rpm, RTC_TYPE_FLOAT, sizeof(mean_speed_rpm), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("max_speed_rpm", &max_speed_rpm, RTC_TYPE_FLOAT, sizeof(max_speed_rpm), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("safety_triggered_speed_rpm", &safety_triggered_speed_rpm, RTC_TYPE_FLOAT, sizeof(safety_triggered_speed_rpm), rtc_data_trigger_read_only, NULL);

	// The analogue voltage as read from the encoder and as a proportion of 2PI. This latter is computed from former.
	rtc_data_add_par("encoder_input_channel",  &encoder_input_channel, RTC_TYPE_UINT32, sizeof(encoder_input_channel),  NULL, NULL);
	rtc_data_add_par("encoder_voltage", &encoder_voltage, RTC_TYPE_FLOAT, sizeof(encoder_voltage), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("encoder_angular_position", &encoder_angular_position, RTC_TYPE_FLOAT, sizeof(encoder_angular_position), rtc_data_trigger_read_only, NULL);

	// Rotary Encoder: speed, computed or read, depending on the encoder type
	rtc_data_add_par("encoder_speed_rpm", &encoder_speed_rpm, RTC_TYPE_FLOAT, sizeof(encoder_speed_rpm), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("encoder_speed", &encoder_speed, RTC_TYPE_FLOAT, sizeof(encoder_speed), rtc_data_trigger_read_only, NULL);
	
	// Rotary Encoder: the linear transfer function (voltage to SPEED). These are computed by the triggers for the above four variables.
	rtc_data_add_par("encoder_transfer_f_gradient", &encoder_transfer_f_gradient, RTC_TYPE_FLOAT, sizeof(encoder_transfer_f_gradient), NULL, NULL);
	rtc_data_add_par("encoder_transfer_f_intercept", &encoder_transfer_f_intercept, RTC_TYPE_FLOAT, sizeof(encoder_transfer_f_intercept), NULL, NULL);

	// Actions. Set to nonzero value to activate and the triggers will do the rest.
	rtc_data_add_par("reset_status_flags", &request_reset_status_flags, RTC_TYPE_UINT32, sizeof(request_reset_status_flags), reset_status_flags, NULL);
	rtc_data_add_par("reset_pid_error_terms", &request_reset_pid_error_terms, RTC_TYPE_UINT32, sizeof(request_reset_pid_error_terms), reset_pid_error_terms, NULL);

	// Lasers:
	//   * raw output in volts. assumes lasers are connected by their index (laser 1 to input 1) [RO]
	//   * calculated position in mm [RO]
	//   * linear transfer function (volts to mm) gradient and intercept
	for (uint32_t i = 0; i < NUMBER_OF_LASERS; i++) {
		char laser_name_volt[] = "laser#_output_volts";
		char laser_name_distance[] = "laser#_distance_mm";
		
		uint8_t laser_no_char = '1' + (uint8_t) i;

		laser_name_volt[5]     = laser_no_char;
		laser_name_distance[5] = laser_no_char;

		// These two families of parameters index into their corresponding arrays, therefore sizeof is of a single element:
		rtc_data_add_par(laser_name_volt, &in_volt[i+1], RTC_TYPE_FLOAT, sizeof(float), rtc_data_trigger_read_only, NULL);
		rtc_data_add_par(laser_name_distance, &lasers_distance_mm[i], RTC_TYPE_FLOAT, sizeof(float), rtc_data_trigger_read_only, NULL);
	}
	rtc_data_add_par("lasers_distance_at_rest", &lasers_distance_at_rest, RTC_TYPE_FLOAT, sizeof(lasers_distance_at_rest), NULL, NULL);
	rtc_data_add_par("lasers_transfer_f_gradient", &lasers_transfer_f_gradient, RTC_TYPE_FLOAT, sizeof(lasers_transfer_f_gradient), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("lasers_transfer_f_intercept", &lasers_transfer_f_intercept, RTC_TYPE_FLOAT, sizeof(lasers_transfer_f_intercept), rtc_data_trigger_read_only, NULL);

	// Set up the filters
	update_filter(&input_filter, pid_error_filter_freq);
	input_filter.order = INPUT_FILTER_ORDER;
}

/* ************************************************************************ */
void rtc_user_main(void)
{
	static int32_t led = 0;

	// Will be set to filtered or unfiltered value, depending on the settings
	static float pid_error;
	
	/* ********************************************************************** */
	/*  Pre-calculations */
	/* ********************************************************************** */
	lasers_compute_distance_from_voltage();
	encoder_compute_speed_from_voltage();
	
	// In case a critical error has happened, don't do anything, sit and wait until 
	// the operator resets the corresponding flag. Critical conditions:
	//   * speed safety limit was reached
	//   * motor set level value by PID controller was not a valid fp number
	if (get_status_flags(RIG_STATUS_SPEED_SAFETY_LIMIT_REACHED | RIG_STATUS_PID_NUMERIC_ERROR))
		goto finalise;
	
	// If speed limit is enabled, compute the mean speed, and shut down the rig
	// if the limit has been reached.
	if (speed_safety_limit_rpm > 0) {
		
		// Update the max speed ever recorded, if the need be
		max_speed_rpm = (encoder_speed_rpm > max_speed_rpm) ? encoder_speed_rpm : max_speed_rpm;
		
		// Only limit the speed once the history buffer is full (small timescale anyway)
		if (speed_history_full) {
			mean_speed_rpm = sum(speed_history, SPEED_HISTORY_SIZE) / SPEED_HISTORY_SIZE;
			if (fabs(mean_speed_rpm) >= speed_safety_limit_rpm) {
				
				shutdown_motor();
				
				safety_triggered_speed_rpm = mean_speed_rpm;
				set_status_flags(RIG_STATUS_SPEED_SAFETY_LIMIT_REACHED);
				goto finalise;
			}
		}
	}
	
	// Sensing and PID output is computed at every tick, actuation is only done
	// at the beginning of the a period, defined by the forcing_freq_hz variable...
	
	// Both filtered and raw PID error is calculated. This is useful for comparison.
	pid_error_raw = target_speed - encoder_speed;
	pid_error_filtered = biquad_filter(pid_error_raw, &input_filter, pid_error_filter_state);
	
	// Decide whether to use filtered or raw error value.
	pid_error = enable_pid_error_filter ? pid_error_filtered : pid_error_raw;

	// Error integral computed as per its standard definition.
	// FIXME is there a risk of integral windup?
	pid_error_integral += (pid_error * TIMER_PERIOD_FLOAT);
	
	// Error derivate -- multiplication by frequency is the same as division by the period
	pid_error_derivative = (pid_error - pid_error_before) * TIMER_FREQ;
	
	// Compute the motor set level using PID principles.
	motor_voltage_level = K_p * pid_error + K_i * pid_error_integral + K_d * pid_error_derivative;
	
	// Prepare for the next time slice by rotating the error term
	pid_error_before = pid_error;
	
	// Save current speed in the history buffer (wraps over if full)
	if (speed_history_idx >= SPEED_HISTORY_SIZE) {
		speed_history_idx = 0;
		speed_history_full = 1;
	}
	speed_history[speed_history_idx++] = encoder_speed_rpm;
	
	// Save angular position
	// FIXME possibly do this AFTER control section
	encoder_angular_position_prev = encoder_angular_position;
	
	/* ********************************************************************** */
	/*  Real-time control */
	/* ********************************************************************** */
	
	// * The control is run at forcing_freq_hz rate. The time_mod_2pi is set up in 
	//   the way that a new period begins exactly at forcing_freq_hz rate.
	//
	// * This controller requires that ESCON is set to (?). Analogue input ? signal sets 
	//   the direction of rotation (using sign) and enables the motor. Analogue input ? 
	//   signal sets the magnitude of the current. The limits on the magnitude are set
	//   in the ESCON software and duplicated in this firmware in motor_{min|max}_voltage
	//   variables.
	//
	// * The PID controller above computes the voltage value not knowing that the sign 
	//   is handled by the motor_enable signal, so care must be taken to ignore the sign
	//   where it is not needed:
	//
	if (period_start) {		
		if (rpm != 0.0f) {
			
			// Shut down the system right now if the PID output value does not make sense!
			if (isfinite(motor_voltage_level) == 0) {
				set_status_flags(RIG_STATUS_PID_NUMERIC_ERROR);
				
				shutdown_motor();
				goto finalise;
			
			// The PID output value is numerically valid, carry on the actuation. Rather 
			// than fiddling with the value produced by the PID controller, split it into
			// the sense of rotation and the magnitude (set level) and send both signals
			// to the right outputs.
			} else {
				// The sense of rotation is given by the motor_enable signal sign:
				//   * CW is (?)
				//   * CCW is (?)
				
				// The voltage level is a magnitude, hence always non-negative
				float motor_set_level = fabs(motor_voltage_level);
				
				// Clip motor voltage levels, if necessary. Record the fact of clipping the output.
				// Rotation direction is set by the sign of voltage level, so ignore the sign.
				if (motor_set_level < motor_min_voltage) {
					motor_set_level = motor_min_voltage;
					set_status_flags(RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MIN);
				} else if (motor_set_level > motor_max_voltage) {
					motor_set_level = motor_max_voltage;
					set_status_flags(RIG_STATUS_MOTOR_VOLTAGE_CLIP_AT_MAX);
				}
				
				// Actuate!
				rtc_set_output(output_channel_motor_direction, 
					motor_voltage_level > 0 ? MAXON_OFF : MAXON_ON);
				rtc_set_output(output_channel_enable_motor, MAXON_ON);
				rtc_set_output(output_channel_motor_set_level, motor_set_level);
			}	
		}	
	}
	
	/* ********************************************************************** */
	/*  Update time */
	/* ********************************************************************** */

finalise:
	// Linear time just ticks along for now...
	time += 1;
	if (time >= SAMPLE_FREQ) {
		time = 0;
	}
	
	// Mod_2pi time drives the led and is periodic
  	period_start = 0;
  	time_mod_2pi += time_delta*forcing_freq;
  	if (time_mod_2pi > M_2PI) {
  		time_mod_2pi -= M_2PI;
  		rtc_led(2, led);
  		led = !led;
  		period_start = 1;
  	}
}


/*
 * AUX FUNCTIONS AHEAD...
 */
 
 
//	Calculate 1/tan(pi*x)
float inv_tan_pi(float x)
{
	v4sf t MEM_ALIGN, sine MEM_ALIGN, cosine MEM_ALIGN;
	float tt[4] MEM_ALIGN, sine_f[4] MEM_ALIGN, cosine_f[4] MEM_ALIGN;

	tt[0] = x; tt[1] = 0.0f; tt[2] = 0.0f; tt[3] = 0.0f;
	t = vld1q_f32(tt);
	t = vmulq_f32(t, vdupq_n_f32(M_PI));
	sincos_ps(t, &sine, &cosine);
	vst1q_f32(sine_f, sine);
	vst1q_f32(cosine_f, cosine);

	return (cosine_f[0]/sine_f[0]);
}

/* Trigger function.
   This allows to set how fast the controller runs. The user specifies the value
	 in times per second (Hz), but since internally the programme runs on 
	 time_mod_2pi, the value in Hz is also converted to rad/s and is used in 
	 subsequent calculations in this way.
	 
   The forcing frequency in Hz can be set by the user from Matlab.
   The forcing frequency in rad/s is recomputed when this happens. */
void update_forcing_freq(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_forcing_freq_hz = *((float *)new_value);
	
	/* Store both frequency and angular frequency */
	forcing_freq_hz = new_forcing_freq_hz;
	forcing_freq = M_2PI * new_forcing_freq_hz;
}

/* Trigger function.
   The target speed in rpm can be set by the user from Matlab.
   The target speed in rad/s is recomputed when this happens. */
void update_target_speed(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_rpm = *((float *)new_value);
	
	// The motor is stopped by setting the target to *EXACTLY* zero.
	// Since this is a well-defined specific value, we should not have any
	// floating point issues? As long as we *SET* it to *EXACTLY* zero lateral
	// constant rather than compute the value which "equals" zero. That value
	// may be *VERY* close, but it won't be exactly zero, so this won't work.
	// Always set the rpm value to exactly zero.
	if (fabs(new_rpm) == 0.0f) {
		shutdown_motor();
	// If new value is not zero, check whether the motor needs starting up,
	// in which case PID error terms and speed history need resetting. 
	} else if (rpm == 0.0f) {
		reset_pid_error_terms(NULL, NULL);
		reset_speed_history(); // TODO let speed history run forever, no problem with that
	}
	
	// Keep both rpm and rad/s values... FIXME why exactly?
	rpm = new_rpm;
	target_speed = rpm2rad(new_rpm);
}

/* **********************************************************************
	Calculate the Butterworth filter coefficients for a normalised
	frequency value (0 < freq < 1)
   ********************************************************************** */
void update_filter(struct biquad_filter_t *filter, float freq)
{
	float ita, q = 1.414213562f, b0;
	/* Check that the filter cut off is between 0 and 1 */
	if ((freq < 0.0f) || (freq >= 1.0f))
		freq = 0.01f; /* default value */
    ita = inv_tan_pi(freq);
    b0 = 1.0f / (1.0f + q*ita + ita*ita);
    filter->a[0] = 2.0f*(ita*ita - 1.0f)*b0;
    filter->a[1] = -(1.0f - q*ita + ita*ita)*b0;
    filter->b[0] = b0;
    filter->b[1] = 2.0f*b0;
    filter->b[2] = b0;
}

//	A biquad filter
float biquad_filter(float x, const struct biquad_filter_t *filter, float state[])
{
	float y = x;

	for (int i = 0; i < filter->order; i++) {
		float w = y + filter->a[0]*state[2*i+0] + filter->a[1]*state[2*i+1];
		y = filter->b[0]*w + filter->b[1]*state[2*i+0] + filter->b[2]*state[2*i+1];
		state[2*i+1] = state[2*i+0];
		state[2*i+0] = w;
	}

	return y;
}

//	Update the coefficients of the input filter
void update_input_filter(void *new_value, void *trigger_data)
{
	float freq = *(float *)new_value;
	if ((freq > 0.0f) && (freq <= 1.0f)) {
		pid_error_filter_freq = freq;
		update_filter(&input_filter, freq);
	}
}

// Shuts down the motor and updates the system's state variables to reflect that.
void shutdown_motor(void)
{
	rtc_set_output(output_channel_enable_motor, MAXON_OFF);
	rtc_set_output(output_channel_motor_set_level, 0.0f);
}

// Returns the value of the status bits indexed by a value idx
static uint32_t get_status_flags(uint32_t idx)
{
	return status_flags & idx;
}

// Sets to one status bits indexed by a value idx
static void set_status_flags(uint32_t idx)
{
	status_flags |= idx;
}

// Set to zero status bits indexed by a value idx
static void clear_status_flags(uint32_t idx)
{
	status_flags &= ~idx;
}

// Clear all status and error flags
void reset_status_flags(void *new_value, void *trigger_data)
{
	clear_status_flags(status_flags); // Whatever is set is going to be cleared
}

void reset_pid_error_terms(void *new_value, void *trigger_data)
{
	pid_error_before     = 0.0f;
	pid_error_raw        = 0.0f;
	pid_error_filtered   = 0.0f;
	pid_error_integral   = 0.0f;
	pid_error_derivative = 0.0f;
}

void reset_speed_history(void)
{
	speed_history_full = 0;
	speed_history_idx  = 0;
	mean_speed_rpm = 0.0f;
	max_speed_rpm  = 0.0f;
	safety_triggered_speed_rpm = 0.0f;
}

// Return the sum of an array of floating-point values.
float sum(float list[], uint32_t n)
{
	float Sigma = 0.0f;
	
	for (uint32_t i = 0; i < n; i++)
		Sigma += list[i];
	
	return Sigma;
}

// Input: speed in rad/s, Output: speed in revolutions per minute.
float rad2rpm(float rad)
{
	return (rad / M_2PI) * 60.0f;
}

// Input: speed in revolutions per minute, Output: speed in rad/s.
float rpm2rad(float rpm)
{
	return rpm * (M_2PI / 60.0f);
}

/*
 * Converts all entries from laser voltage output array to displacement in mm,
 * using the coefficients for the linear model. Assumes laser i is connected to
 * input i. Assumes all lasers are identical (both the laser and the signal
 * conditioner) so the same transfer function is applicable to all.
 */
void lasers_compute_distance_from_voltage(void)
{
	for (uint32_t i = 0; i < NUMBER_OF_LASERS; i++) {
		float voltage_read;
		
		voltage_read = in_volt[i+1];
		
		// Is the value outside of legit range of the lasers? To understand the comparison:
		// the laser outputs 4 to 20 mA, with the output being weaker at the distance,
		// therefore the limit of the near field is larger amount, even closer distance
		// would be even greater current; the far field is very weak signal, so going further
		// means even weaker signal.
		if (voltage_read < lasers_far_field_limit_volt || voltage_read > lasers_near_field_limit_volt)
			lasers_distance_mm[i] = NAN;
		else
			lasers_distance_mm[i] = lasers_transfer_f_gradient * voltage_read + lasers_transfer_f_intercept;
	}
}

// Maxon Motor Escon 50/5 encoder returns speed in rpm (as voltage)
void encoder_compute_speed_from_voltage(void)
{
	encoder_voltage = in_volt[encoder_input_channel];
	
	encoder_speed_rpm = encoder_voltage * encoder_transfer_f_gradient + encoder_transfer_f_intercept;
	encoder_speed = rpm2rad(encoder_speed_rpm);
}

void encoder_compute_speed_from_position(void)
{
	encoder_speed = (encoder_angular_position_prev - encoder_angular_position) / time_delta;
	encoder_speed_rpm = rad2rpm(encoder_speed);
}

// TODO this got to go to general RTC code, it's useful!
void bbb_reboot(void *new_value, void *trigger_data)
{
	HWREG(SOC_PRM_DEVICE_REGS) |= 0x00000001;
}