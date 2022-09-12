#include "transforms.h"

void vConvertToMono( short* psBuffer, uint8_t uiChannels )
{
    if( uiChannels == 1 ) return;

    for( uint32_t i = 0; i < BUFFER_SIZE; i++ )
    {
        int iAverage = 0;
        for( uint8_t c = 0; c < uiChannels; c++ )
        {
            iAverage += psBuffer[ i * uiChannels + c ];
        }
        iAverage /= uiChannels;
        psBuffer[i] = iAverage;
    }
}

short sCalculateNoiseLevel( short* psBuffer, uint32_t uiNoiseWindow )
{
    long long llMean = 0;
    long long llStandardDeviation = 0;
    
    for( uint32_t i = 0; i < uiNoiseWindow; i++ )
    {
        llMean += psBuffer[i];
    }

    llMean = ( float ) llMean / uiNoiseWindow;

    for( uint32_t i = 0; i < uiNoiseWindow; i++ )
    {
        llStandardDeviation +=  ( llMean - psBuffer[i] ) * ( llMean - psBuffer[i] );
    }
    llStandardDeviation = sqrt( llStandardDeviation / uiNoiseWindow );

    // printf( "Mean: %lld, STD: %lld\n", llMean, llStandardDeviation );

    return llMean + 2 * llStandardDeviation;
}

void vFFT( std::vector<int> &xVector, uint16_t usWindowSize )
{
     /* We need to allocate these only once. */
    static kiss_fft_cfg xFftCfg;
    static kiss_fft_cpx* pxComplexInput ;
    static kiss_fft_cpx* pxSpectrum;
    if( pxComplexInput == NULL )
    {
        xFftCfg = kiss_fft_alloc( usWindowSize, 0, NULL, NULL );
        pxComplexInput    = ( kiss_fft_cpx* ) malloc( usWindowSize * sizeof( kiss_fft_cpx ) );
        pxSpectrum        = ( kiss_fft_cpx* ) malloc( usWindowSize * sizeof( kiss_fft_cpx ) );
    }

    for( uint32_t i = 0; i < usWindowSize; i++ )
    {
        pxComplexInput[i].r = xVector[i];
        pxComplexInput[i].i = 0;
    }

    kiss_fft( xFftCfg, pxComplexInput, pxSpectrum );

    /* Absolute value of a complex number is just a magnitude of a vector. */
    for( uint32_t i = 0; i < usWindowSize; i++ )
    {
        xVector[i] = sqrt( pxSpectrum[i].r * pxSpectrum[i].r + pxSpectrum[i].i * pxSpectrum[i].i );
    }
}

void vSpectralMelCoefficients( std::vector<double> &xMelCoefficients, std::vector<int> &xPowerSpectrum, uint16_t usSampleRate )
{
    uint16_t usDeltaMel         = HZ_TO_MEL( usSampleRate / 2 ) / xMelCoefficients.size(); /* freqMax / nFilters. */
    uint16_t usFftResolution    = usSampleRate / xPowerSpectrum.size();
    uint16_t usPreviousCenter   = 0;

    for( uint8_t i = 0; i < xMelCoefficients.size(); i++ )
    {
        uint16_t usFilterCenter = round( MEL_TO_HZ( usDeltaMel * ( i + 1 ) ) / usFftResolution ); /* Absolute index in Power Spectrum. */
        uint16_t usFilterEnd    = 2 * usFilterCenter - usPreviousCenter;
        uint16_t usWidth        = usFilterEnd - usPreviousCenter;
        double dCoeff = 0.0;

        for( uint16_t j = usPreviousCenter; j < usFilterEnd; j++ )
        {
            uint16_t n = j - usPreviousCenter;
            dCoeff += TRIANGLE_WINDOW( n, usWidth ) * xPowerSpectrum[j];
        }

        xMelCoefficients[i] = dCoeff;
        usPreviousCenter = usFilterCenter;
    }
}

void vCosineTransform( std::vector<double> &xMelCepstrum, std::vector<double> &xMelCoefficients )
{
    for( uint8_t i = 0; i < MFCC_COEFFICEINTS_; i++ )
    {   
        double dMFCC = 0.0;
        for( uint8_t j = 0; j < xMelCoefficients.size(); j++ )
        {
            dMFCC += log( xMelCoefficients[j] ) * cos( M_PI * i * ( 2 * j + 1 ) / ( 2.0 * xMelCoefficients.size() ) );
        }
        xMelCepstrum[i] = dMFCC;
    }

    #if defined( MFCC_FIRST_DERIVATIVE ) || defined( MFCC_SECOND_DERIVATIVE )
    {
        xMelCepstrum[ MFCC_COEFFICEINTS_ ] = xMelCepstrum[2];
        xMelCepstrum[ MFCC_COEFFICEINTS_ + 1 ] = xMelCepstrum[3];
        xMelCepstrum[ MFCC_COEFFICEINTS_ * 2 - 1 ] = xMelCepstrum[ MFCC_COEFFICEINTS_ - 3 ];
        xMelCepstrum[ MFCC_COEFFICEINTS_ * 2 - 2 ] = xMelCepstrum[ MFCC_COEFFICEINTS_ - 4 ];

        for( uint8_t i = 2; i < MFCC_COEFFICEINTS_ - 2; i++ )
        {
            xMelCepstrum[ MFCC_COEFFICEINTS_ + i ] = xMelCepstrum[ i + 2 ] - xMelCepstrum[ i - 2 ];
        } 
    }
    #endif
    #ifdef MFCC_SECOND_DERIVATIVE
    {
        xMelCepstrum[ MFCC_COEFFICEINTS_ * 2 ] = xMelCepstrum[2];
        xMelCepstrum[ MFCC_COEFFICEINTS_ * 2 + 1 ] = xMelCepstrum[3];
        xMelCepstrum[ MFCC_COEFFICEINTS_ * 3 - 1 ] = xMelCepstrum[ MFCC_COEFFICEINTS_ - 3 ];
        xMelCepstrum[ MFCC_COEFFICEINTS_ * 3 - 2 ] = xMelCepstrum[ MFCC_COEFFICEINTS_ - 4 ];

        for( uint8_t i = 2; i < MFCC_COEFFICEINTS_ - 2; i++ )
        {
            xMelCepstrum[ MFCC_COEFFICEINTS_ * 2 + i ] = xMelCepstrum[ MFCC_COEFFICEINTS_ + i + 2 ] - xMelCepstrum[ MFCC_COEFFICEINTS_ + i - 2 ];
        } 
    }
    #endif

    // for( uint8_t i = 0; i < MFCC_COEFFICEINTS; i++ )
    // {
    //     printf( "%f, ", xMelCepstrum[i] );
    // }
    // printf( "\n" );
    // exit( 0 );
}