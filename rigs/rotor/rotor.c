/*
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

// Rotor rig, Oliver Frolovs, 2018
// Use together with rigtest.mlx Matlab Live Script

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
#include "orbis.h"

/* ************************************************************************ */
/* * Defines ************************************************************** */
/* ************************************************************************ */

// Number of lasers in the system
#define LASER_N 3

// Maths constants for convenience
#define M_2PI (2.0f * M_PI)

// Maxon motor controller 50/5 control voltage levels as per its datasheet
#define MAXON_ON   3.5f
#define MAXON_OFF  0.0f

// Number of sin/cos values to precompute (???)
#define N_FOURIER_MODES 7  // N_FOURIER_MODES + 1 must be a multiple of 4
#define N_FOURIER_COEFF (2*N_FOURIER_MODES + 2)  // one will always be zero but it's easier/quicker than handling the special case

// Encoder data filter
#define INPUT_FILTER_ORDER 2
#define INPUT_FILTER_N_STATE (2*INPUT_FILTER_ORDER)

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

// Time keeping is done mod 2*pi (???)
static float time_mod_2pi;
static float time_delta = M_2PI/TIMER_FREQ;
static float period_start = 1.0f;

// Forcing frequency (???)
static float forcing_freq = 1.0f;

// Encoder data filter
static float pid_proportional_error_filter_freq = 0.025f; // fraction of sampling frequency
static struct biquad_filter_t input_filter;

// Angular velocity to maintain, in rpm
static float rpm;

// Motor control: two channels, "digital" to enable and set direction and analog to set value
static uint32_t chan_motor_enabled  = 1; /* Output U1 (unipolar) */
static uint32_t chan_motor_set = 0; /* Output B0 (bipolar) */

// Motor control: set to a positive value to turn on the motor
static uint32_t motor_enabled;

// Motor: set voltage level. This must be within [motor_set_min_voltage, motor_set_max_voltage]
static float motor_set_voltage;

// Motor: voltage limits set in the motor controller
static float motor_set_min_voltage;
static float motor_set_max_voltage;

// Motor: min and max voltage attempted flag
static uint32_t motor_min_voltage_flag;
static uint32_t motor_max_voltage_flag;

// Motor PID controller coefficients
static float forcing;
static float K_p;
static float K_i;
static float K_d;

// Backward whirl flag
static uint32_t backward_whirl_detected_flag;

// An accelerometer is attached to the rig frame. This is its input channel.
static uint32_t chan_frame_acceleration = 7;

// Accelerometer settings
static uint32_t frame_accelerometer_gain     = 1;
static uint32_t frame_accelerometer_filter   = 1000; // 1 kHz
static float frame_accelerometer_sensitivity = 98.4; // mV/g

// Frame acceleration
static float frame_acceleration_volts;
static float frame_acceleration_g;

// Frame accelerometer backward whirl threshold (in volts)
static float frame_acceleration_threshold_g;

// By how many rpm to backtrack when backward whirl is encountered?
static float backward_rpm_backtrack;

// Control backward whirl enable flag
static uint32_t backward_whirl_control_enabled;

// PID error terms
static float pid_proportional_error_raw;
static float pid_proportional_error;
static float pid_previous_proportional_error;
static float pid_derivative_error;
static float pid_integral_error;

// PID filter for control signal enable flag
static uint32_t pid_proportional_error_filter_enabled;

// Speed safety limit enable flag
static uint32_t speed_safety_limit_enabled = 1;

// The motor will turn off if this speed is reached in either direction of rotation (rpm)
static float speed_safety_limit_rpm;

// Indicates that the safety limit has been reached
static uint32_t speed_safety_limit_reached_flag;

// Fourier calculations
static float sinusoid_f[N_FOURIER_COEFF] MEM_ALIGN;  /* [sin(0*t), sin(1*t), sin(2*t), ..., cos(0*t), cos(1*t), cos(2*t), ...] */

// Filter state (???)
static float pid_proportional_error_filter_state[INPUT_FILTER_N_STATE];  // filtering encoder output

// Control parameter: angular velocity in rad/s
static float angular_velocity_target;

// Laser displacement in mm
static float lasers_mm[LASER_N];

// Laser coefficients [distance, mm] = ([output, volts] - b) / a
// The default values are for the lasers labelled {1,2,3} which I have
// They were computed by taking distance readings for two points for each laser,
// and computing the line through these two points.
// TODO lasers serial numbers for these values?
static float lasers_a[LASER_N] = {0.011617834394904f, 0.013365079365079f, 0.013323809523809f};
static float lasers_b[LASER_N] = {2.628232484076433f, 2.993146031746035f, 2.995842857142860f};

// Laser output at rest
static float lasers_initial_voltage[LASER_N];

// Laser output when out of range. Set to a very permissive value.
static float lasers_out_of_range_voltage = 6.0f;

// Requests. Set to nonzero to activate. Once the action is complete,
// the request flag is set to zero again.
static uint32_t request_reset_error_and_status_flags;
static uint32_t request_reset_pid_error_terms;
static uint32_t request_lasers_set_zero_point;
static float    request_rpm_to;

/* ************************************************************************ */
/* * Internal prototypes ************************************************** */
/* ************************************************************************ */
void calc_sinusoids(void);
float inv_tan_pi(float x);
float biquad_filter(float x, const struct biquad_filter_t *filter, float state[]);
void update_filter(struct biquad_filter_t *filter, float freq);
void update_input_filter(void *new_value, void *trigger_data);
float rpm2rad(float rpm);
void shutdown_motor(void);
void reset_error_and_status_flags(void);
void reset_pid_error_terms(void);
void lasers_set_zero_point(void);
void lasers_convert_voltage_to_distance(void);
float accelerometer_volts_to_g(float v, uint32_t gain, uint32_t sensitivity);

/* ************************************************************************ */
/* * Function definitions ************************************************* */
/* ************************************************************************ */
void rtc_user_init(void)
{
	char laser_name_volt[] = "laser#_volt";
	char laser_name_displacement[] = "laser#_displacement";
	char laser_name_coefficient_a[] = "laser#_a";
	char laser_name_coefficient_b[] = "laser#_b";

	// Tracking of time inside rtc_user_main is mod 2*pi
	rtc_data_add_par("time_mod_2pi", &time_mod_2pi, RTC_TYPE_FLOAT, sizeof(time_mod_2pi), NULL, NULL);
	rtc_data_add_par("forcing_freq", &forcing_freq, RTC_TYPE_FLOAT, sizeof(forcing_freq), NULL, NULL);

	// Desired speed of rotation in RPM
	rtc_data_add_par("rpm", &rpm, RTC_TYPE_FLOAT, sizeof(rpm), NULL, NULL);

	// Motor: enable and set value channel numbers
	rtc_data_add_par("chan_motor_enabled",  &chan_motor_enabled,  RTC_TYPE_UINT32, sizeof(chan_motor_enabled),  NULL, NULL);
	rtc_data_add_par("chan_motor_set", &chan_motor_set, RTC_TYPE_UINT32, sizeof(chan_motor_set), NULL, NULL);

	// Motor: the switch to turn the motor on or off
	rtc_data_add_par("motor_enabled", &motor_enabled, RTC_TYPE_UINT32, sizeof(motor_enabled), NULL, NULL);

	// Motor: voltage input translates into current drawn according to the motor controller settings
	rtc_data_add_par("motor_set_voltage", &motor_set_voltage, RTC_TYPE_FLOAT, sizeof(motor_set_voltage), NULL, NULL);

	// Motor: the limits on input voltage are set in the motor controller settings
	rtc_data_add_par("motor_set_min_voltage", &motor_set_min_voltage, RTC_TYPE_FLOAT, sizeof(motor_set_min_voltage), NULL, NULL);
	rtc_data_add_par("motor_set_max_voltage", &motor_set_max_voltage, RTC_TYPE_FLOAT, sizeof(motor_set_max_voltage), NULL, NULL);

	// Motor: min/max voltage status flag (sticky)
	rtc_data_add_par("motor_min_voltage_flag",  &motor_min_voltage_flag,  RTC_TYPE_UINT32, sizeof(motor_min_voltage_flag), NULL, NULL);
	rtc_data_add_par("motor_max_voltage_flag",  &motor_max_voltage_flag,  RTC_TYPE_UINT32, sizeof(motor_max_voltage_flag), NULL, NULL);

	// Frame accelerometer input channel
	rtc_data_add_par("chan_frame_acceleration", &chan_frame_acceleration, RTC_TYPE_UINT32, sizeof(chan_frame_acceleration), rtc_data_trigger_read_only, NULL);

	// Accelerometer settings
	rtc_data_add_par("frame_accelerometer_gain", &frame_accelerometer_gain, RTC_TYPE_UINT32, sizeof(frame_accelerometer_gain), NULL, NULL);
	rtc_data_add_par("frame_accelerometer_filter", &frame_accelerometer_filter, RTC_TYPE_UINT32, sizeof(frame_accelerometer_filter), NULL, NULL);
	rtc_data_add_par("frame_accelerometer_sensitivity", &frame_accelerometer_sensitivity, RTC_TYPE_FLOAT, sizeof(frame_accelerometer_sensitivity), NULL, NULL);

	// Acceleration of the frame
	rtc_data_add_par("frame_acceleration_volts", &frame_acceleration_volts, RTC_TYPE_FLOAT, sizeof(frame_acceleration_volts), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("frame_acceleration_g", &frame_acceleration_g, RTC_TYPE_FLOAT, sizeof(frame_acceleration_g), rtc_data_trigger_read_only, NULL);

	// Backward whirl flag and frame acceleration threshold
	rtc_data_add_par("backward_whirl_control_enabled", &backward_whirl_control_enabled, RTC_TYPE_UINT32, sizeof(backward_whirl_control_enabled), NULL, NULL);
	rtc_data_add_par("backward_whirl_detected_flag", &backward_whirl_detected_flag, RTC_TYPE_UINT32, sizeof(backward_whirl_detected_flag), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("frame_acceleration_threshold_g", &frame_acceleration_threshold_g, RTC_TYPE_FLOAT, sizeof(frame_acceleration_threshold_g), NULL, NULL);
	rtc_data_add_par("backward_rpm_backtrack", &backward_rpm_backtrack, RTC_TYPE_FLOAT, sizeof(backward_rpm_backtrack), NULL, NULL);

	// Motor PID controller coefficients. All [RW]
	rtc_data_add_par("forcing", &forcing, RTC_TYPE_FLOAT, sizeof(forcing), NULL, NULL);
	rtc_data_add_par("K_p", &K_p, RTC_TYPE_FLOAT, sizeof(K_p), NULL, NULL);
	rtc_data_add_par("K_i", &K_i, RTC_TYPE_FLOAT, sizeof(K_i), NULL, NULL);
	rtc_data_add_par("K_d", &K_d, RTC_TYPE_FLOAT, sizeof(K_d), NULL, NULL);

	// Motor PID controller propertional error filter, enable flag and value in ???
	rtc_data_add_par("pid_proportional_error_filter_enabled", &pid_proportional_error_filter_enabled, RTC_TYPE_UINT32, sizeof(pid_proportional_error_filter_enabled), NULL, NULL);
	rtc_data_add_par("pid_proportional_error_filter_freq", &pid_proportional_error_filter_freq, RTC_TYPE_FLOAT, sizeof(pid_proportional_error_filter_freq), update_input_filter, NULL);

	// PID error terms. All [RW]
	rtc_data_add_par("pid_proportional_error_raw", &pid_proportional_error_raw, RTC_TYPE_FLOAT, sizeof(pid_proportional_error_raw), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("pid_proportional_error", &pid_proportional_error, RTC_TYPE_FLOAT, sizeof(pid_proportional_error), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("pid_integral_error", &pid_integral_error, RTC_TYPE_FLOAT, sizeof(pid_integral_error), rtc_data_trigger_read_only, NULL);
	rtc_data_add_par("pid_derivative_error", &pid_derivative_error, RTC_TYPE_FLOAT, sizeof(pid_derivative_error), rtc_data_trigger_read_only, NULL);

	// Speed safety limit enable flag, value [rpm], and status flag (sticky). All [RW]
	rtc_data_add_par("speed_safety_limit_enabled",  &speed_safety_limit_enabled,  RTC_TYPE_UINT32, sizeof(speed_safety_limit_enabled), NULL, NULL);
	rtc_data_add_par("speed_safety_limit_rpm", &speed_safety_limit_rpm, RTC_TYPE_FLOAT, sizeof(speed_safety_limit_rpm), NULL, NULL);
	rtc_data_add_par("speed_safety_limit_reached_flag",  &speed_safety_limit_reached_flag,  RTC_TYPE_UINT32, sizeof(speed_safety_limit_reached_flag), NULL, NULL);

	// Requests. All [RW]
	rtc_data_add_par("lasers_set_zero_point", &request_lasers_set_zero_point, RTC_TYPE_UINT32, sizeof(request_lasers_set_zero_point), NULL, NULL);
	rtc_data_add_par("reset_error_and_status_flags", &request_reset_error_and_status_flags, RTC_TYPE_UINT32, sizeof(request_reset_error_and_status_flags), NULL, NULL);
	rtc_data_add_par("reset_pid_error_terms", &request_reset_pid_error_terms, RTC_TYPE_UINT32, sizeof(request_reset_pid_error_terms), NULL, NULL);
	rtc_data_add_par("request_rpm_to", &request_rpm_to, RTC_TYPE_FLOAT, sizeof(request_rpm_to), NULL, NULL);

	// Lasers:
	//   * raw output in volts. assumes lasers are connected by their index (laser 1 to input 1) [RO]
	//   * calculated displacement in mm [RO]
	//   * linear coefficients to convert output voltage to displacement [RW]
	for (uint32_t i = 0; i < LASER_N; i++) {
		uint8_t laser_no_char = '1' + (uint8_t) i;

		laser_name_volt[5]          = laser_no_char;
		laser_name_displacement[5]  = laser_no_char;
		laser_name_coefficient_a[5] = laser_no_char;
		laser_name_coefficient_b[5] = laser_no_char;

		rtc_data_add_par(laser_name_volt, &in_volt[i+1], RTC_TYPE_FLOAT, sizeof(float), rtc_data_trigger_read_only, NULL);
		rtc_data_add_par(laser_name_displacement, &lasers_mm[i], RTC_TYPE_FLOAT, sizeof(float), rtc_data_trigger_read_only, NULL);
		rtc_data_add_par(laser_name_coefficient_a, &lasers_a[i], RTC_TYPE_FLOAT, sizeof(float), NULL, NULL);
		rtc_data_add_par(laser_name_coefficient_b, &lasers_b[i], RTC_TYPE_FLOAT, sizeof(float), NULL, NULL);
	}
	rtc_data_add_par("lasers_out_of_range_voltage", &lasers_out_of_range_voltage, RTC_TYPE_FLOAT, sizeof(lasers_out_of_range_voltage), NULL, NULL);

	// Set up the filters
	update_filter(&input_filter, pid_proportional_error_filter_freq);
	input_filter.order = INPUT_FILTER_ORDER;
}

/* ************************************************************************ */
void rtc_user_main(void)
{
	static int led;

	// Will be set to filtered or unfiltered value, depending on the settings
	static float pid_error;

	// Backward whirl is detected at present time
	static uint32_t backward_whirl_detected;

	// Calculate sinusoids for use in the Fourier calculations (???)
	calc_sinusoids();

	// Prepare laser output for reading or using in control signal calculations
	lasers_convert_voltage_to_distance();

	// Convert desired angular velocity from rpm to rad/s
	angular_velocity_target = rpm2rad(rpm);

	// Acceleration of the frame [in volts]
	frame_acceleration_volts = in_volt[chan_frame_acceleration];
	frame_acceleration_g = accelerometer_volts_to_g(frame_acceleration_volts,
		frame_accelerometer_gain, frame_accelerometer_sensitivity);

	/*
	 * Safety...
	 */
	if (orbisInTrouble()) shutdown_motor();

	// On by default: if the speed exceeds the safety limit, shut down the motor
	if (speed_safety_limit_enabled) {
		if (fabs(orbisAngularVelocity) > fabs(rpm2rad(speed_safety_limit_rpm))) {
			shutdown_motor();
			speed_safety_limit_reached_flag = 1;
		}
	}

	/*
	 * Processing user requests...
	 */

	if (request_reset_error_and_status_flags) {
		reset_error_and_status_flags();
		request_reset_error_and_status_flags = 0;
	}

	if (request_reset_pid_error_terms) {
		reset_pid_error_terms();
		request_reset_pid_error_terms = 0;
	}

	if (request_lasers_set_zero_point) {
		lasers_set_zero_point();
		request_lasers_set_zero_point = 0;
	}

	// Keep slowly increasing the speed if requested. Currently ~5 rpm/sec
	if (request_rpm_to > 0.0f) {
		if (request_rpm_to > rpm) {
			rpm += 0.005f;
		} else {
			request_rpm_to = 0.0f;
		}
	}

	/*
	 * Detecting and handling backward whirl, if requested
	 */

	// Backward whirl is detected via abnormal acceleration of the frame
	backward_whirl_detected = (fabs(frame_acceleration_g) >= fabs(frame_acceleration_threshold_g)) ? 1 : 0;

	// The global flag is sticky
	backward_whirl_detected_flag |= backward_whirl_detected;

	// Reduce angular velocity target if backward whirl control is enabled
	if (backward_whirl_control_enabled && backward_whirl_detected) {
		rpm -= backward_rpm_backtrack;
		// Don't change the sense of rotation
		if (rpm < 0) {
			rpm = 1.0f;
		}
		angular_velocity_target = rpm2rad(rpm);
	}

	/*
	 *  PID controller error terms...
	 */

	pid_proportional_error_raw = angular_velocity_target - orbisAngularVelocity;
	pid_proportional_error = biquad_filter(pid_proportional_error_raw, &input_filter, pid_proportional_error_filter_state);

	// Use filtered error signal or raw error signal in PID control
	if (pid_proportional_error_filter_enabled)
		pid_error = pid_proportional_error;
	else
		pid_error = pid_proportional_error_raw;

	// THAT case when comparing fp to zero is actually correct
	// we SET it to zero when we want the error terms reset
	if (rpm == 0.0f) {
		reset_pid_error_terms();
	} else {
		pid_integral_error += (pid_error * TIMER_PERIOD_FLOAT);
		pid_derivative_error = (pid_error - pid_previous_proportional_error) * TIMER_FREQ;
		pid_previous_proportional_error = pid_error;
	}

	// FIXME handle potential integral windup

	// Compute the motor set level using PID principles. Forcing is (optional) constant input.
	motor_set_voltage = forcing
				+ K_p * pid_error
				+ K_i * pid_integral_error
				+ K_d * pid_derivative_error;

	// Limit the motor_set_voltage to levels accepted by the controller and set the corresponding flag
	if (motor_set_voltage < motor_set_min_voltage) {
		motor_set_voltage = motor_set_min_voltage;
		motor_min_voltage_flag = 1;
	} else if (motor_set_voltage > motor_set_max_voltage) {
		motor_set_voltage = motor_set_max_voltage;
		motor_max_voltage_flag = 1;
	}

	// Enable or disable the motor, as requested and set voltage level
	if (motor_enabled) {
		rtc_set_output(chan_motor_enabled, MAXON_ON);
		rtc_set_output(chan_motor_set, motor_set_voltage);
	} else {
		shutdown_motor();
	}

	//  Update time and blink the LED on the board
	period_start = 0.0f;
	time_mod_2pi += time_delta*forcing_freq;
	if (time_mod_2pi > M_2PI) {
		time_mod_2pi -= M_2PI;
		rtc_led(2, led);
		led = !led;
		period_start = 1.0f;
	}
}

//	Calculate sinusoids
void calc_sinusoids()
{
	int i;
	v4sf t MEM_ALIGN, sine MEM_ALIGN, cosine MEM_ALIGN;
	float tt[4] MEM_ALIGN = {0.0f, 1.0f, 2.0f, 3.0f};

	t = vld1q_f32(tt);
	t = vmulq_f32(t, vdupq_n_f32(time_mod_2pi));
	sincos_ps(t, &sine, &cosine);
	vst1q_f32(&sinusoid_f[0], sine);
	vst1q_f32(&sinusoid_f[N_FOURIER_COEFF/2], cosine);

	for (i = 1; i < N_FOURIER_COEFF/8; i++) {
		t = vaddq_f32(t, vdupq_n_f32(4.0f*time_mod_2pi));
		sincos_ps(t, &sine, &cosine);
		vst1q_f32(&sinusoid_f[4*i], sine);
		vst1q_f32(&sinusoid_f[4*i + N_FOURIER_COEFF/2], cosine);
	}
}

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

//	Update the coefficients of the input filter
void update_input_filter(void *new_value, void *trigger_data)
{
	float freq = *(float *)new_value;
	if ((freq > 0.0f) && (freq <= 1.0f)) {
		pid_proportional_error_filter_freq = freq;
		update_filter(&input_filter, freq);
	}
}

// Convert given revolutions per minute value to radians per second
float rpm2rad(float rpm)
{
	return rpm * (M_2PI / 60.0f);
}

// Shuts down the motor and updates the system's state variables to reflect that
void shutdown_motor(void)
{
	rtc_set_output(chan_motor_enabled, MAXON_OFF);
	rtc_set_output(chan_motor_set, 0.0f);
	motor_enabled = 0;
}

// Stop the system, clear the error flags and controller error terms
void reset_error_and_status_flags(void)
{
	motor_min_voltage_flag = 0;
	motor_max_voltage_flag = 0;
	speed_safety_limit_reached_flag = 0;
	backward_whirl_detected_flag = 0;

	orbisReset();
}

void reset_pid_error_terms(void)
{
	pid_previous_proportional_error = 0.0f;
	pid_proportional_error_raw      = 0.0f;
	pid_proportional_error          = 0.0f;
	pid_integral_error              = 0.0f;
	pid_derivative_error            = 0.0f;
}

void lasers_set_zero_point(void)
{
	// FXIME would rather take an average, but not sure how to within a single call to rtc_user_main?
	// FIXME does this mean that taking an averaged zero point is an operation best done on Matlab side?
	// FIXME before anything else happens?
	for (uint32_t i = 0; i < LASER_N; i++) {
		lasers_initial_voltage[i] = in_volt[i+1];
	}
}

/*
 * Converts all entries from laser voltage output array to displacement in mm,
 * using the coefficients for the linear model. Out of range value is NAN,
 * but the default threshold is very permissive, so unless adjusted, the out
 * of range voltage values will be converted to (incorrect) distance.
 */
void lasers_convert_voltage_to_distance(void)
{
	for (uint32_t i = 0; i < LASER_N; i++) {
		float value = in_volt[i+1];

		if (value < lasers_out_of_range_voltage) {
			lasers_mm[i] = (in_volt[i+1] - lasers_b[i]) / lasers_a[i];
		} else {
			lasers_mm[i] = NAN;
		}
	}
}

// Convert accelerometer output value in volts to value in g,
// using gain and sensitivity parameters
float accelerometer_volts_to_g(float v, uint32_t gain, uint32_t sensitivity)
{
	// the gain is provided in mV, so convert voltage to mV first
	return ((v * 1000.0f * gain) / sensitivity);
}
