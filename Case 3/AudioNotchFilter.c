/*********************************************************************************

Copyright(c) 2015 Analog Devices, Inc. All Rights Reserved.

This software is proprietary and confidential.  By using this software you agree
to the terms of the associated Analog Devices License Agreement.

Modified by Kim Bjerge, Aarhus University, 27-06-2017

*********************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "AudioCallback.h"

# define BITS 14
# define BIT_SHIFT_FACTOR (1<<BITS)
#define SAMPLES_PER_CHAN   (NUM_AUDIO_SAMPLES/2)

// Input samples, 24 LSBs
int32_t inLeft[SAMPLES_PER_CHAN];
int32_t inRight[SAMPLES_PER_CHAN];

// Output samples, 24 LSBs
int32_t outLeft[SAMPLES_PER_CHAN];
int32_t outRight[SAMPLES_PER_CHAN];

// Processing mode BYPASS or FILTER active
FILTER_MODE currentMode = PASS_THROUGH;

// Filter coefficients (parameters + Q14 storage).
double r = 0.99;                 // Pole radius
double f0 = 785.0;               // Notch frequency in Hz
double fs = 48000.0;             // Sampling frequency in Hz

/* compute coeffs at runtime in FilterInit - do NOT call cos() at file scope */
short b0 = 0;
short b1 = 0;
short b2 = 0;
short a1 = 0;
short a2 = 0;

// ...existing code...

/* initialize the IIR filter, called once when button PB1 pressed */
void FilterInit(FILTER_MODE mode)
{
    switch (mode)
    {
        case PASS_THROUGH:
            currentMode = PASS_THROUGH;
            printf("Pass through\n");
            break;
        case IIR_FILTER_ACTIVE:
            currentMode = IIR_FILTER_ACTIVE;
            printf("Filter on\n");

            /* Minimal, low-level coefficient conversion to Q14 (no lround, no extra libs).
               Just cast the double*BIT_SHIFT_FACTOR to short (truncation). */
            {
                double omega = 2.0 * 3.14159 * (f0 / fs);
                double c = cos(omega);

                b0 = (short)(1.0 * (double)BIT_SHIFT_FACTOR);
                b1 = (short)(-2.0 * c * (double)BIT_SHIFT_FACTOR);
                b2 = (short)(1.0 * (double)BIT_SHIFT_FACTOR);

                a1 = (short)(-2.0 * r * c * (double)BIT_SHIFT_FACTOR);
                a2 = (short)(r * r * (double)BIT_SHIFT_FACTOR);
            }
            break;
    }
}

// Modify and insert your notch filter here!!!!
short myNotchFilter(short x)
{
    // Static variables to store previous input and output samples
    static short x1 = 0, x2 = 0;  // Previous input samples
    static short y1 = 0, y2 = 0;  // Previous output samples
    
    int32_t y;  // changed from int to int32_t to reduce overflow risk

    // IIR filter 
    y = b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2;
    y = y >> BITS; // Adjust for fixed-point scaling

    // Update previous samples

    x2 = x1;
    x1 = x;
    y2 = y1;
    y1 = y;

    return (short)y;
}

/* Compute filter response or bypass  */
void AudioNotchFilter(const int32_t dataIn[], int32_t dataOut[])
{
	int n, i;
	short xn, yn;

	// Copy dataIn to left/right input buffers
	i = 0;
	for (n=0; n<SAMPLES_PER_CHAN; n++){
		inLeft[n] = dataIn[i++];
		inRight[n] = dataIn[i++];
	}

	// Perform pass through or filtering
	switch (currentMode)
	{
		case PASS_THROUGH:
			// Copy input to output
			for (n=0; n<SAMPLES_PER_CHAN; n++)
			{
				outLeft[n] = inLeft[n];
				outRight[n] = inRight[n];
			}
			break;

		case IIR_FILTER_ACTIVE:
			// IIR filter active
			for (n=0; n<SAMPLES_PER_CHAN; n++)
			{
				// Filter in left channel only
				xn = inLeft[n]>>8; // 24 bits samples right aligned, remove 8 LSBs
				yn = myNotchFilter(xn);
				outLeft[n] = yn << 8;

				// Mute right channel
				xn = 0; //inRight[n]>>8; // 24 bits samples right aligned, remove 8 LSBs
				outRight[n] = xn << 8;
			}
			break;

	}

	i = 0;
	// Copy left/right output buffers to dataOut
	for (n=0; n<SAMPLES_PER_CHAN; n++){
		dataOut[i++] = outLeft[n];
		dataOut[i++] = outRight[n];
	}

}
