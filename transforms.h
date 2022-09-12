#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "./kissfft/kiss_fft.h"
// #include "DTW_cpp/include/DTW.hpp"

/* Useful conversion macros. */
#define HZ_TO_MEL( hz ) ( 1127 * log( 1 + ( hz ) / 700.0 ) )            /* log() is natural log. */
#define MEL_TO_HZ( mel ) ( 700 * ( pow( M_E, ( mel ) / 1127.0 ) - 1 ) ) /* M_E is eulers number. */
#define MS_TO_FRAMES( ms, sampleRate ) ( ( uint32_t ) sampleRate / 1000 ) * ms

/* Window functions. */
#define TRIANGLE_WINDOW( n, N ) ( 1 - fabs( ( 2.0 * n - N ) / N ) )     /* N is window width. */

#define BUFFER_SIZE         8192

#define MFCC_COEFFICEINTS_   13

#ifdef MFCC_FIRST_DERIVATIVE
    #define MFCC_COEFFICEINTS   2 * MFCC_COEFFICEINTS_
#endif

#ifdef MFCC_SECOND_DERIVATIVE
    #define MFCC_COEFFICEINTS   3 * MFCC_COEFFICEINTS_
#endif

#ifndef MFCC_COEFFICEINTS
    #define MFCC_COEFFICEINTS   MFCC_COEFFICEINTS_
#endif


void vConvertToMono( short* psBuffer, uint8_t uiChannels );
short sCalculateNoiseLevel( short* psBuffer, uint32_t uiNoiseWindow );
void vFFT( std::vector<int> &xVector, uint16_t usWindowSize );
void vSpectralMelCoefficients( std::vector<double> &xMelCoefficients, std::vector<int> &xPowerSpectrum, uint16_t usSampleRate );
void vCosineTransform( std::vector<double> &xMelCepstrum, std::vector<double> &xMelCoefficients );

#endif