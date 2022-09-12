#include <sndfile.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <list>
#include <Python.h>
#include "./kissfft/kiss_fft.h"
#include "./matplotlib-cpp/matplotlibcpp.h"
#include "DTW_cpp/include/DTW.hpp"


#include "transforms.h"

namespace plt = matplotlibcpp;

#define PLOT_NOISE          0
#define PLOT_SPECTROGRAM    1
#define PLOT_SONOGRAM       2
#define CLASSIFIER          3

typedef struct {
    char pcClassName[15];
    std::vector< std::vector<double> > xPrototypeMatrix;
} xClass;

uint8_t ucHandleArgs( int argc, char** argv )
{
    uint8_t ucProgramMode = PLOT_NOISE;
    if( argc < 3 )
    {
        printf( "Usage: program filename.wav [ noise | spectrogram | sonogram | classifier database.txt ]\n" );
        exit( 1 );
    }
    else if( !strcmp( argv[2], "spectrogram" ) )
    {
        ucProgramMode = PLOT_SPECTROGRAM;
    }
    else if( !strcmp( argv[2], "sonogram" ) )
    {
        ucProgramMode = PLOT_SONOGRAM;
    }
    else if( !strcmp( argv[2], "classifier" ) )
    {
        ucProgramMode = CLASSIFIER;
    }

    return ucProgramMode;
}

void vExtractMfccVectors( std::vector< std::vector< double > > &xMfccVectors, char* pcFileName )
{  
    SF_INFO* pxSfInfo = ( SF_INFO* ) malloc( sizeof( SF_INFO ) );
    SNDFILE* pxSoundFile = sf_open( pcFileName, SFM_READ, pxSfInfo );


    short* psBuffer = ( short* ) malloc( BUFFER_SIZE * pxSfInfo->channels * sizeof( short ) );

    uint32_t uiNoiseFrames = MS_TO_FRAMES( 150, pxSfInfo->samplerate );
    uint16_t usFftWindowLength = MS_TO_FRAMES( 20, pxSfInfo->samplerate );

    sf_readf_short( pxSoundFile, psBuffer, uiNoiseFrames );
    vConvertToMono( psBuffer, pxSfInfo->channels );

    short sNoiseLimit = sCalculateNoiseLevel( psBuffer, uiNoiseFrames );

    /* Checking for begining of the word/sound. */
    uint32_t uiFirst = 0, uiLast = 0;

    /* Reading file in chunks of BUFFER_SIZE. */
    const uint16_t usOverflow = uiNoiseFrames + ( ( pxSfInfo->frames - uiNoiseFrames ) / BUFFER_SIZE + 1 ) * BUFFER_SIZE - pxSfInfo->frames;
    for( uint64_t uiOffset = uiNoiseFrames; uiOffset < pxSfInfo->frames; uiOffset += BUFFER_SIZE )
    {
        sf_readf_short( pxSoundFile, psBuffer, BUFFER_SIZE );   /* Reads N * number of channels. */
        vConvertToMono( psBuffer, pxSfInfo->channels );

        /* Zero padding the buffer after EOF. */
        if( pxSfInfo->frames - uiOffset < BUFFER_SIZE )
        {
            for( uint16_t i = 1; i < usOverflow; i++ )
            {
                psBuffer[ BUFFER_SIZE - i ] = 0;
            }
        }
        else if( uiFirst != 0 )
        {
            /* Calculating FFT for the window of size usFftWindowLength */
            for( uint32_t i = 0; i < BUFFER_SIZE; i += usFftWindowLength / 2 )
            {
                std::vector<int> xSpectVect( psBuffer + i, psBuffer + i + usFftWindowLength ); //&xPltVect[i], &xPltVect[i] + usFftWindowLength );
                vFFT( xSpectVect , usFftWindowLength );

                std::vector<double> xMelCoefficients( 22 );
                vSpectralMelCoefficients( xMelCoefficients, xSpectVect, pxSfInfo->samplerate );

                std::vector<double> xMelCepstrum( MFCC_COEFFICEINTS );
                vCosineTransform( xMelCepstrum, xMelCoefficients );

                xMfccVectors.push_back( xMelCepstrum );
            }
        }

        /* Checking for begining and end of the word/sound. */
        for( uint16_t i = 0; i < BUFFER_SIZE; i++ )
        {
            if( abs( psBuffer[i] ) > sNoiseLimit && uiFirst == 0 )
            {
                uiFirst = uiOffset + i ;
            }
            else if( abs( psBuffer[i] ) > sNoiseLimit )
            {
                uiLast = uiOffset + i;
            }
        }
    }

    free( pxSfInfo );
    free( psBuffer );
}

std::vector< std::vector<double> >  xMakePrototype(
                    std::vector< std::vector< double > > &xMaleMfcc, 
                    std::vector< std::vector< double > > &xFemaleMfcc )
{
    DTW::DTW xDtw( xMaleMfcc, xFemaleMfcc, 1.0 );
    std::vector< std::vector< int > > xPath = xDtw.path();
    std::vector< std::vector<double> > xPrototypeMfcc;
    
    for ( uint16_t i = 0; i < xPath.size(); i++ )
    {
        std::vector<double> xMfccVector( MFCC_COEFFICEINTS );
        for( uint8_t j = 0; j < MFCC_COEFFICEINTS; j++ )
        {
            xMfccVector[j] = ( xMaleMfcc[ xPath[i][0] ][j] + xFemaleMfcc[ xPath[i][1] ][j] ) / 2.0;
        }
        xPrototypeMfcc.push_back( xMfccVector );
        // printf( "( %d, %d ), ", xPath[i][0], xPath[i][1] );

    }

    return xPrototypeMfcc;
}

void vGenerateMfccDatabase( std::list<xClass> &xDatabase, char* pcFileName )
{
    FILE* pxFile = fopen( pcFileName , "r" );
    char * pcLine = NULL;
    size_t lLength = 0;
    ssize_t lRead;

    if ( pxFile == NULL ) exit( EXIT_FAILURE );

    char pcAudioName[30];
    char pcLineCopy[10];

    while ( ( lRead = getline( &pcLine, &lLength, pxFile ) ) != -1 )
    {
        // xClass* pxNewClass = ( xClass* ) malloc( sizeof( xClass ) );
        xClass xNewClass;

        strcpy( pcAudioName, "samples/" );
        strcpy( pcLineCopy, pcLine );

        if( pcLineCopy[ lRead - 1 ] == '\n' ) pcLineCopy[ lRead - 1 ] = '\0';
        
        strcat( pcAudioName, pcLineCopy );
        strcat( pcAudioName, "-m.wav" );
        strcpy( xNewClass.pcClassName, pcLineCopy );

        std::vector< std::vector<double> > xMaleMfcc;
        vExtractMfccVectors( xMaleMfcc, pcAudioName );

        strcpy( pcAudioName, "samples/" );
        strcpy( pcLineCopy, pcLine );
        if( pcLineCopy[ lRead - 1 ] == '\n' ) pcLineCopy[ lRead - 1 ] = '\0';
        strcat( pcAudioName, pcLineCopy );
        strcat( pcAudioName, "-f.wav" );

        std::vector< std::vector<double> > xFemaleMfcc;
        vExtractMfccVectors( xFemaleMfcc, pcAudioName );
        
        xNewClass.xPrototypeMatrix = xMakePrototype( xMaleMfcc, xFemaleMfcc );

        xDatabase.push_back( xNewClass );
       
    }

    fclose( pxFile );
        
}

int main( int argc, char** argv )
{
    uint8_t ucProgramMode =  ucHandleArgs( argc, argv );

    if( ucProgramMode == CLASSIFIER )
    {
        std::list<xClass> xDatabase;
        vGenerateMfccDatabase( xDatabase, argv[3] );

        std::vector< std::vector<double> > xQueryFeatures;
        vExtractMfccVectors( xQueryFeatures, argv[1] );

        for( auto &xClass: xDatabase )
        {
            // DTW::DTW xDtw( xQueryFeatures, xClass, 2.0 );
            double dDistance = DTW::dtw_distance_only( xQueryFeatures, xClass.xPrototypeMatrix, 1.0  );
            printf( "%s - %f\n", xClass.pcClassName, dDistance );
        }

        exit( 0 );
    }
    
    SF_INFO* pxSfInfo = ( SF_INFO* ) malloc( sizeof( SF_INFO ) );
    SNDFILE* pxSoundFile = sf_open( argv[1], SFM_READ, pxSfInfo );


    short* psBuffer = ( short* ) malloc( BUFFER_SIZE * pxSfInfo->channels * sizeof( short ) );

    uint32_t uiNoiseFrames = MS_TO_FRAMES( 150, pxSfInfo->samplerate );
    uint16_t usFftWindowLength = MS_TO_FRAMES( 20, pxSfInfo->samplerate );

    sf_readf_short( pxSoundFile, psBuffer, uiNoiseFrames );
    vConvertToMono( psBuffer, pxSfInfo->channels );

    short sNoiseLimit = sCalculateNoiseLevel( psBuffer, uiNoiseFrames );

    /*
    There was a bug where matplotlib overflows values of the X-axis. That's why <int> is used.
    The error occures in matplotlibcpp.h around line 1832 in the function plot( y ); 
    std::vector<short> xPltVect( psBuffer, psBuffer + uiNoiseFrames );
    */
    std::vector<int> xPltVect( psBuffer, psBuffer + uiNoiseFrames );
    std::vector< std::vector<double> >xMFCC;
    std::vector< std::vector<int> > xSonogramData;

    /* Checking for begining of the word/sound. */
    uint32_t uiFirst = 0, uiLast = 0;
    // for( uint16_t i = 0; i < uiNoiseFrames; i++ )
    // {
    //     if( abs( psBuffer[i] ) > sNoiseLimit && uiFirst == 0 )
    //     {
    //         uiFirst = i;
    //     }
    // }

    /* Reading file in chunks of BUFFER_SIZE. */
    const uint16_t usOverflow = uiNoiseFrames + ( ( pxSfInfo->frames - uiNoiseFrames ) / BUFFER_SIZE + 1 ) * BUFFER_SIZE - pxSfInfo->frames;
    for( uint64_t uiOffset = uiNoiseFrames; uiOffset < pxSfInfo->frames; uiOffset += BUFFER_SIZE )
    {
        sf_readf_short( pxSoundFile, psBuffer, BUFFER_SIZE );   /* Reads N * number of channels. */
        vConvertToMono( psBuffer, pxSfInfo->channels );

        /* Zero padding the buffer after EOF. */
        if( pxSfInfo->frames - uiOffset < BUFFER_SIZE )
        {
            for( uint16_t i = 1; i < usOverflow; i++ )
            {
                psBuffer[ BUFFER_SIZE - i ] = 0;
            }
        }
        else if( ucProgramMode != PLOT_NOISE && uiFirst != 0 )
        {
            /* Calculating FFT for the window of size usFftWindowLength */
            for( uint32_t i = 0; i < BUFFER_SIZE; i += usFftWindowLength / 2 )
            {
                std::vector<int> xSpectVect( psBuffer + i, psBuffer + i + usFftWindowLength ); //&xPltVect[i], &xPltVect[i] + usFftWindowLength );
                vFFT( xSpectVect , usFftWindowLength );

                // std::vector<double> xMelCoefficients( 22 );
                // vSpectralMelCoefficients( xMelCoefficients, xSpectVect, pxSfInfo->samplerate );

                // printf( "FRAME: %llu ", uiOffset + i );
                //std::vector<int> xMelCepstrum( 13 );
                //vCosineTransform( xMelCepstrum, xMelCoefficients );
                
                // exit( 0 );

                xSonogramData.push_back( xSpectVect );
            }
        }
            
        std::vector<int> xBufferVect( psBuffer, psBuffer + BUFFER_SIZE );
        xPltVect.insert( std::end( xPltVect ), std::begin( xBufferVect ), std::end( xBufferVect ) );

        /* Checking for begining and end of the word/sound. */
        for( uint16_t i = 0; i < BUFFER_SIZE; i++ )
        {
            if( abs( psBuffer[i] ) > sNoiseLimit && uiFirst == 0 )
            {
                uiFirst = uiOffset + i ;
            }
            else if( abs( psBuffer[i] ) > sNoiseLimit )
            {
                uiLast = uiOffset + i;
            }
        }
    }

    std::map<std::string, std::string> xMapArgs = { { "linewidth", "5.0" } };
    if( ucProgramMode == PLOT_NOISE )
    {
        if( uiFirst == 0 || uiLast == 0 )
        {
            printf( "Error: No words found, only noise.\n" );
            return 0;
        }

        plt::plot( xPltVect );// xMapArgs );
        plt::axhline( sNoiseLimit );
        plt::axhline( -sNoiseLimit );
        plt::axvline( uiFirst );
        plt::axvline( uiLast );
    } 
    else if( ucProgramMode == PLOT_SPECTROGRAM )
    { 
        plt::xlabel( "Frequency ( Hz )" );
        plt::ylabel( "Amplitude" );
        plt::plot( xSonogramData[ xSonogramData.size() / 2 ], pxSfInfo->samplerate / usFftWindowLength, xMapArgs ); /* Converting to Hz. */
    }
    else if( ucProgramMode == PLOT_SONOGRAM )
    {
        plt::xlabel( "Time 10ms" );
        plt::ylabel( "Frequency ( Hz )" );
        for( uint32_t uiTime = 0; uiTime < xSonogramData.size(); uiTime++ )
        {
            std::vector<int> xFreqVect( xSonogramData[uiTime].begin(), xSonogramData[uiTime].begin() + usFftWindowLength / 2 );
            std::vector<int> a{ 0 };
            std::vector<int> b{ 0 };

            char pcOpacity[20];
            uint16_t usMin = 0xffff, usMax = 0;

            for( int i = 0; i < xFreqVect.size(); i++ )
            {
                if( xFreqVect[i] > usMax ) usMax = xFreqVect[i];
                if( xFreqVect[i] < usMin ) usMin = xFreqVect[i];
            }

            for( uint16_t usFreq = 0; usFreq < xFreqVect.size(); usFreq++ )
            {
                a[0] = uiTime;
                b[0] = xFreqVect[usFreq];
                
                if( b[0] / 500 > 50 ) /* No need to plot really low amplitude frequencies, saves a lot of time. */
                {
                    double dOpacity = ( ( double )( b[0] - usMin ) ) / ( usMax - usMin );   /* Normalizing amplitudes. */
                    sprintf( pcOpacity, "%lf", dOpacity );
                    std::string str( pcOpacity );
                    
                    std::map<std::string, std::string> xMapOpacity = { { "alpha", str } };
                    b[0] = usFreq * ( pxSfInfo->samplerate / usFftWindowLength );

                    plt::scatter( a, b, 2.0, xMapOpacity );
                }
            }
        }
    }
    

    plt::show();

    return 0;   
}