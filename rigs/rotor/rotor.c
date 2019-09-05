/*
** Written by Oliver Frolovs, based on the examples by
** Copyright (c) 2014, David A.W. Barton
** (david.barton@bristol.ac.uk) All rights reserved.
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
static float sinusoid_f[N_FOURIER_COEFF] MEM_ALIGN;  /* [sin(0*t), sin(1*t), sin(2*t), ..., cos(0*t), cos(1*t), cos(2*t), ...] */

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
static float forcing_freq = 1.0f;

// CONTROL TARGET: angular velocity to maintain, in rpm
static float rpm;

// CONTROL TARGET: angular velocity to maintain in rad/s. This value is computed
// by the update trigger of rpm variable.
static float target_angular_velocity;

// Motor control: two channels, "digital" to enable and set direction and analog to set value
static uint32_t output_channel_enable_motor    = 1; /* Output U1 (unipolar) */
static uint32_t output_channel_motor_set_level = 0; /* Output B0 (bipolar) */

// Motor: set voltage level. This must be within [motor_min_voltage, motor_max_voltage]
static float motor_voltage_level;

// Motor: voltage limits set in the motor controller
static float motor_min_voltage;
static float motor_max_voltage;

// Motor: min and max voltage attempted flag
static uint32_t motor_min_voltage_flag;
static uint32_t motor_max_voltage_flag;

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

// Lasers: displacement in mm
static float lasers_output_mm[NUMBER_OF_LASERS];

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
static float lasers_transfer_f_gradient = 75.0f;
static float lasers_transfer_f_intercept = -225.0f;

// Rotary Encoder: voltage as read, and its interpretation as angular position,
// as the proportion of 2*pi. The encoder's analogue signal is connected to input zero.
static uint32_t encoder_input_channel;
static float encoder_voltage;
static float encoder_angular_position;
static float encoder_angular_position_prev;
static float encoder_speed_rpm;
static float encoder_speed;     // in rad/s

// Rotary Encoder: Maxon Motor Controller Escon 50/5 built-in encoder outputs
// speed (12-bit resolution) as voltage. The linear range is set in ESCON
// studio software and the values should match those below for correct conversion.
static float encoder_low_voltage = 0.5f;
static float encoder_low_speed_rpm = 0;
static float encoder_high_voltage = 4.5f;
static float encoder_high_speed_rpm = 3000;

// Rotary Encoder: transfer function (voltage to SPEED) gradient and intercept
static float encoder_transfer_f_gradient;
static float encoder_transfer_f_intercept;

// Set to enable the speed control, unset if safety feature is not desired
static uint32_t enable_speed_safety_limit = 1;

// The motor will turn off if this speed is reached in either direction of rotation (rpm)
static float speed_safety_limit_rpm = 600;
 
// Indicates that the safety limit has been reached
static uint32_t speed_safety_limit_reached_flag;

// Buffer to collect last N speed readings for averaging. Used for 
// enforcing the speed safety limit. The buffer wraps over at the end, using
// the index variable.
static float speed_history[SPEED_HISTORY_SIZE];
 
// Speed history buffer index
static uint32_t speed_history_idx;
 
// Since the buffer is circular, it's impossible to tell from the index alone,
// if it has been filled to full capacity. Therefore, there is a flag to indicate so.
static uint32_t speed_history_full;

// The mean spead will be calculated by averaging the full set of values from 
// the speed_history array.
static float mean_speed;


/******************************************************************************
 *
 * OPTIONAL FEATURE: FILTER PID ERROR TERM
 */

static uint32_t enable_pid_error_filter;     // Motor PID filter for control signal enable flag
static float pid_error_filter_freq = 0.025f; // fraction of sampling frequency
static struct biquad_filter_t input_filter;
static float pid_error_filter_state[INPUT_FILTER_N_STATE];  // filtering encoder output


// Requests. Set to nonzero to activate. The triggers on update will do the job.
static uint32_t request_reset_error_and_status_flags;
static uint32_t request_reset_pid_error_terms;


/* ************************************************************************ */
/* * Internal prototypes ************************************************** */
/* ************************************************************************ */

// Triggers for various exported variables...
void update_forcing_freq(void *new_value, void *trigger_data);            // Trigger
void update_target_angular_velocity(void *new_value, void *trigger_data); // Trigger
void update_input_filter(void *new_value, void *trigger_data);            // Trigger
void reset_error_and_status_flags(void *new_value, void *trigger_data);   // Trigger
void reset_pid_error_terms(void *new_value, void *trigger_data);          // Trigger

// Triggers for the rotary encoder (voltage to SPEED conversion related variables)
void update_encoder_low_voltage(void *new_value, void *trigger_data);     // Trigger
void update_encoder_low_speed_rpm(void *new_value, void *trigger_data);   // Trigger
void update_encoder_high_voltage(void *new_value, void *trigger_data);    // Trigger
void update_encoder_high_speed_rpm(void *new_value, void *trigger_data);  // Trigger

void reset_speed_history(void);

// Encoder transfer function converts input voltage to speed in rpm. It is linear, so it requires gradient and intercept values to operate.
void update_encoder_transfer_f(float encoder_low_voltage, float encoder_low_speed_rpm, float encoder_high_voltage, float encoder_high_speed_rpm);

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
	char laser_name_volt[] = "laser#_output_volts";
	char laser_name_displacement[] = "laser#_output_mm";

	// Tracking of time, linear and periodic counters are provided
	rtc_data_add_par("time",  &time, RTC_TYPE_UINT32, sizeof(time), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("time_mod_2pi", &time_mod_2pi, RTC_TYPE_FLOAT, sizeof(time_mod_2pi), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("forcing_freq", &forcing_freq, RTC_TYPE_FLOAT, sizeof(forcing_freq), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("forcing_freq_hz", &forcing_freq_hz, RTC_TYPE_FLOAT, sizeof(forcing_freq_hz), update_forcing_freq, NULL);

	// Desired speed of rotation in RPM
	rtc_data_add_par("rpm", &rpm, RTC_TYPE_FLOAT, sizeof(rpm), update_target_angular_velocity, NULL);

	// Motor: enable and set value channel numbers
	rtc_data_add_par("output_channel_enable_motor",  &output_channel_enable_motor, RTC_TYPE_UINT32, sizeof(output_channel_enable_motor),  NULL, NULL);
	rtc_data_add_par("output_channel_motor_set_level", &output_channel_motor_set_level, RTC_TYPE_UINT32, sizeof(output_channel_motor_set_level), NULL, NULL);

	// Motor: voltage input translates into current drawn according to the motor controller settings
	rtc_data_add_par("motor_voltage_level", &motor_voltage_level, RTC_TYPE_FLOAT, sizeof(motor_voltage_level), NULL, NULL);

	// Motor: the limits on input voltage are set in the motor controller settings
	rtc_data_add_par("motor_min_voltage", &motor_min_voltage, RTC_TYPE_FLOAT, sizeof(motor_min_voltage), NULL, NULL);
	rtc_data_add_par("motor_max_voltage", &motor_max_voltage, RTC_TYPE_FLOAT, sizeof(motor_max_voltage), NULL, NULL);

	// Motor: min/max voltage status flag (sticky)
	rtc_data_add_par("motor_min_voltage_flag",  &motor_min_voltage_flag,  RTC_TYPE_UINT32, sizeof(motor_min_voltage_flag), NULL, NULL);
	rtc_data_add_par("motor_max_voltage_flag",  &motor_max_voltage_flag,  RTC_TYPE_UINT32, sizeof(motor_max_voltage_flag), NULL, NULL);

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
	rtc_data_add_par("speed_safety_limit_reached_flag",  &speed_safety_limit_reached_flag,  RTC_TYPE_UINT32, sizeof(speed_safety_limit_reached_flag), NULL, NULL);
	rtc_data_add_par("speed_history", &speed_history, RTC_TYPE_FLOAT, sizeof(speed_history), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("mean_speed", &mean_speed, RTC_TYPE_FLOAT, sizeof(mean_speed), rtc_data_trigger_read_only, NULL);

	// The analogue voltage as read from the encoder and as a proportion of 2PI. This latter is computed from former.
	rtc_data_add_par("encoder_input_channel",  &encoder_input_channel, RTC_TYPE_UINT32, sizeof(encoder_input_channel),  NULL, NULL);
	rtc_data_add_par("encoder_voltage", &encoder_voltage, RTC_TYPE_FLOAT, sizeof(encoder_voltage), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("encoder_angular_position", &encoder_angular_position, RTC_TYPE_FLOAT, sizeof(encoder_angular_position), rtc_data_trigger_read_only, NULL);

	// Rotary Encoder: the Maxon Motor Controller Escon 50/5 outputs a speed reading (as voltage), and the following values
	// are required to convert that to RPM value.
	rtc_data_add_par("encoder_low_voltage", &encoder_low_voltage, RTC_TYPE_FLOAT, sizeof(encoder_low_voltage), update_encoder_low_voltage, NULL);
	rtc_data_add_par("encoder_low_speed_rpm", &encoder_low_speed_rpm, RTC_TYPE_FLOAT, sizeof(encoder_low_speed_rpm), update_encoder_low_speed_rpm, NULL);
	rtc_data_add_par("encoder_high_voltage", &encoder_high_voltage, RTC_TYPE_FLOAT, sizeof(encoder_high_voltage), update_encoder_high_voltage, NULL);
	rtc_data_add_par("encoder_high_speed_rpm", &encoder_high_speed_rpm, RTC_TYPE_FLOAT, sizeof(encoder_high_speed_rpm), update_encoder_high_speed_rpm, NULL);

	// Rotary Encoder: the linear transfer function (voltage to SPEED). These are computed by the triggers for the above four variables.
	rtc_data_add_par("encoder_transfer_f_gradient", &encoder_transfer_f_gradient, RTC_TYPE_FLOAT, sizeof(encoder_transfer_f_gradient), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("encoder_transfer_f_intercept", &encoder_transfer_f_intercept, RTC_TYPE_FLOAT, sizeof(encoder_transfer_f_intercept), rtc_data_trigger_read_only, NULL);

	// Rotary Encoder: calculate the transfer function for the voltage to SPEED encoder
	update_encoder_transfer_f(encoder_low_voltage, encoder_low_speed_rpm, encoder_high_voltage, encoder_high_speed_rpm);

	// Actions. Set to nonzero value to activate and the triggers will do the rest.
	rtc_data_add_par("reset_error_and_status_flags", &request_reset_error_and_status_flags, RTC_TYPE_UINT32, sizeof(request_reset_error_and_status_flags), reset_error_and_status_flags, NULL);
	rtc_data_add_par("reset_pid_error_terms", &request_reset_pid_error_terms, RTC_TYPE_UINT32, sizeof(request_reset_pid_error_terms), reset_pid_error_terms, NULL);

	// Lasers:
	//   * raw output in volts. assumes lasers are connected by their index (laser 1 to input 1) [RO]
	//   * calculated position in mm [RO]
	//   * linear transfer function (volts to mm) gradient and intercept
	for (uint32_t i = 0; i < NUMBER_OF_LASERS; i++) {
		uint8_t laser_no_char = '1' + (uint8_t) i;

		laser_name_volt[5]          = laser_no_char;
		laser_name_displacement[5]  = laser_no_char;

		rtc_data_add_par(laser_name_volt, &in_volt[i+1], RTC_TYPE_FLOAT, sizeof(float), rtc_data_trigger_read_only, NULL);
		rtc_data_add_par(laser_name_displacement, &lasers_output_mm[i], RTC_TYPE_FLOAT, sizeof(float), rtc_data_trigger_read_only, NULL);
	}
	rtc_data_add_par("lasers_transfer_f_gradient", &lasers_transfer_f_gradient, RTC_TYPE_FLOAT, sizeof(lasers_transfer_f_gradient), NULL, NULL);
	rtc_data_add_par("lasers_transfer_f_intercept", &lasers_transfer_f_intercept, RTC_TYPE_FLOAT, sizeof(lasers_transfer_f_intercept), NULL, NULL);

	// Set up the filters
	update_filter(&input_filter, pid_error_filter_freq);
	input_filter.order = INPUT_FILTER_ORDER;
	
	// Forcing frequency in rad/s
	update_forcing_freq((void *)&forcing_freq_hz, NULL);	
}

/* ************************************************************************ */
void rtc_user_main(void)
{
	static int32_t led = 0;

	/* ********************************************************************** */
	/*  Pre-calculations */
	/* ********************************************************************** */
	lasers_compute_distance_from_voltage();
	/* ********************************************************************** */
	/*  Real-time control */
	/* ********************************************************************** */

	/* ********************************************************************** */
	/*  Update time */
	/* ********************************************************************** */

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

// Trigger function: on update, the transfer function for the encoder must be recomputed.
void update_encoder_low_voltage(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_encoder_low_voltage = *((float *)new_value);
	
	// TODO verify?
	encoder_low_voltage = new_encoder_low_voltage;
	
	/* Update the transfer function */
	update_encoder_transfer_f(encoder_low_voltage, encoder_low_speed_rpm, encoder_high_voltage, encoder_high_speed_rpm);
}

// Trigger function: on update, the transfer function for the encoder must be recomputed.
void update_encoder_low_speed_rpm(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_encoder_low_speed_rpm = *((float *)new_value);
	
	// TODO verify?
	encoder_low_speed_rpm = new_encoder_low_speed_rpm;
	
	/* Update the transfer function */
	update_encoder_transfer_f(encoder_low_voltage, encoder_low_speed_rpm, encoder_high_voltage, encoder_high_speed_rpm);
}

// Trigger function: on update, the transfer function for the encoder must be recomputed.
void update_encoder_high_voltage(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_encoder_high_voltage = *((float *)new_value);
	
	// TODO verify?
	encoder_high_voltage = new_encoder_high_voltage;
	
	/* Update the transfer function */
	update_encoder_transfer_f(encoder_low_voltage, encoder_low_speed_rpm, encoder_high_voltage, encoder_high_speed_rpm);
}

// Trigger function: on update, the transfer function for the encoder must be recomputed.
void update_encoder_high_speed_rpm(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_encoder_high_speed_rpm = *((float *)new_value);
	
	// TODO verify?
	encoder_high_speed_rpm = new_encoder_high_speed_rpm;
	
	/* Update the transfer function */
	update_encoder_transfer_f(encoder_low_voltage, encoder_low_speed_rpm, encoder_high_voltage, encoder_high_speed_rpm);
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
void update_target_angular_velocity(void *new_value, void *trigger_data)
{
	/* Extract the data from the pointers */
	float new_rpm = *((float *)new_value);
	
	// Full stop requested, stop the motor and reset the controller.
	if (fabs(new_rpm) == 0) {
		//shutdown_motor();
		reset_pid_error_terms(NULL, NULL);
		reset_speed_history();
	}
	
	// In any case, store both rpm and rad/s values.
	rpm = new_rpm;
	target_angular_velocity = rpm2rad(new_rpm);
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

// Stop the system, clear the error flags and controller error terms
void reset_error_and_status_flags(void *new_value, void *trigger_data)
{
	motor_min_voltage_flag = 0;
	motor_max_voltage_flag = 0;
	speed_safety_limit_reached_flag = 0;
}

void reset_pid_error_terms(void *new_value, void *trigger_data)
{
	pid_error_before = 0.0f;
	pid_error_raw    = 0.0f;
	pid_error_filtered = 0.0f;
	pid_error_integral = 0.0f;
	pid_error_derivative = 0.0f;
}

void reset_speed_history(void)
{
	speed_history_full = 0; // Flag to indicate that there is enough points to compute the mean value
	speed_history_idx = 0;  // Index to where the next data point is saved
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
			lasers_output_mm[i] = NAN;
		else
			lasers_output_mm[i] = lasers_transfer_f_gradient * voltage_read + lasers_transfer_f_intercept;
	}
}

// Maxon Motor Escon 50/5 encoder returns speed in rpm (as voltage)
void encoder_compute_speed_from_voltage(void)
{
	encoder_speed_rpm = encoder_voltage * encoder_transfer_f_gradient + encoder_transfer_f_intercept;
	encoder_speed = rpm2rad(encoder_speed);
}

void encoder_compute_speed_from_position(void)
{
	encoder_speed = (encoder_angular_position_prev - encoder_angular_position) / time_delta;
	encoder_speed_rpm = rad2rpm(encoder_speed);
}

// Solves a system of two linear equations in two variables in the most trivial way
void update_encoder_transfer_f(float encoder_low_voltage, float encoder_low_speed_rpm, float encoder_high_voltage, float encoder_high_speed_rpm)
{
	float x1 = encoder_low_voltage;
	float x2 = encoder_high_voltage;
	
	float y1 = encoder_low_speed_rpm;
	float y2 = encoder_high_speed_rpm;
	
	float gradient, intercept;
	
	if ((y1 - y2) == 0) return;
	
	gradient = (x1 - x2) / (y1 - y2);
	intercept = y1 - gradient * x1;
	
	encoder_transfer_f_gradient = gradient;
	encoder_transfer_f_intercept = intercept;
}
