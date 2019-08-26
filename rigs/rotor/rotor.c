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
#include "rtc_data.h"
#include "rtc_main.h"
#include "rtc_util.h"
#include "rtc_user.h"
#include "neon_mathfun.h"
#include <string.h>

/* ************************************************************************ */
/* * Defines ************************************************************** */
/* ************************************************************************ */
#define M_2PI  6.283185307f  /* 2*pi */
#define MOTOR_ON  4.2
#define MOTOR_OFF 0.0
/* ********************************************************************** */
/* * Types ************************************************************** */
/* ********************************************************************** */

/* ********************************************************************** */
/* * Globals (internal) ************************************************* */
/* ********************************************************************** */
static float time_mod_2pi;
static float time_delta;
static float period_start;  /* signal the start of a new period */

static unsigned int motor_enable_chan = 1;
static unsigned int motor_set_value_chan = 0;


static unsigned int motor_enable;
static float motor_set_value;

/* ************************************************************************ */
/* * Function definitions ************************************************* */
/* ************************************************************************ */
void rtc_user_init(void)
{
	/* initialise time and other fundamental stuff */
	rtc_data_add_par("time_mod_2pi", &time_mod_2pi, RTC_TYPE_FLOAT, sizeof(time_mod_2pi), NULL, NULL);
	time_delta = M_2PI/TIMER_FREQ;
	period_start = 1;
	
	/* How the motor is connected to mcu */
	rtc_data_add_par("motor_enable_chan", &motor_enable_chan, RTC_TYPE_UINT32, sizeof(motor_enable_chan), NULL, NULL);
	rtc_data_add_par("motor_set_value_chan", &motor_set_value_chan, RTC_TYPE_UINT32, sizeof(motor_set_value_chan), NULL, NULL);
	
	/* How the motor is actually operated */
	rtc_data_add_par("motor_enable", &motor_enable, RTC_TYPE_UINT32, sizeof(motor_enable), NULL, NULL);
	rtc_data_add_par("motor_set_value", &motor_set_value, RTC_TYPE_FLOAT, sizeof(motor_set_value), NULL, NULL);
}

/* ************************************************************************ */
void rtc_user_main(void)
{
	static int led = 0;
	static unsigned int time = 0;

	/* ********************************************************************** */
	/*  Pre-calculations */
	/* ********************************************************************** */

	/* ********************************************************************** */
	/*  Real-time control */
	/* ********************************************************************** */
	rtc_set_output(motor_enable_chan, motor_enable);
	rtc_set_output(motor_set_value_chan, motor_set_value);

	/* ********************************************************************** */
	/*  Update time */
	/* ********************************************************************** */

	time += 1;
	time_mod_2pi += time_delta;
	if (time >= SAMPLE_FREQ) {
		time = 0;
		rtc_led(2, led);
		led = !led;
	}

}
