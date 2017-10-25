#ifdef DEBUG
    #include <cstdio>
    #include <cstdlib>
#endif // DEBUG
#include <cstring>
#include <cmath>

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

#include "minmax.h"
#include "rawimgtk_tmodrago03.h"

////////////////////////////////////////////////////////////////////////////////

static inline double biasFunction( const double b, const double x )
{
    // pow(x, log(bias)/log(0.5)
	return pow(x, b);
}

static inline double pade_log( const double x )
{
	if( x < 1.0 )
    {
		return ( x * ( 6.0 + x ) / ( 6.0 + 4.0 * x ) );
	}
	else
    if( x < 2.0 )
    {
		return ( x * ( 6.0 + 0.7662 * x ) / ( 5.9897 + 3.7658 * x ) );
	}

	return log( x + 1.0 );
}

static bool toneMappingDrago03( float* src, unsigned srcsz, const float maxLum, const float avgLum, float biasParam, const float exposure)
{
    if ( src == NULL )
    {
        return false;
    }

	const float LOG05 = -0.693147F;	/// == log(0.5)

	double Lmax;
	double divider;
	double interpol;
	double biasP;
	double L;

	// arbitrary Bias Parameter
	if( biasParam == 0.0 )
    {
        biasParam = 0.85f;
    }

	// normalize maximum luminance by average luminance
	Lmax = maxLum / avgLum;

	divider = log10(Lmax+1);
	biasP   = log(biasParam)/LOG05;

	#pragma omp parallel for
	for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        double Yw   = (double)src[cnt] / (double)avgLum;
        Yw         *= exposure;
        interpol    = log( 2.0 + biasFunction( biasP, Yw / Lmax ) * 8.0 );
        L           = pade_log(Yw);// log(Yw + 1)
        src[cnt]    = (float)( ( L / interpol) / divider );
	}

	return true;
}

static bool getLuminance( float* src, unsigned srcsz, float* maxLum, float* minLum, float* worldLum )
{
    if ( src == NULL )
    {
        return false;
    }

	float max_lum = 0;
	float min_lum = 0;
	double sum = 0;

	#pragma omp pararelle for
	for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        const float fp = MAX( 0.0f, src[cnt] );
        max_lum = ( max_lum < fp ) ? fp : max_lum;	/// max Luminance in the scene
        min_lum = ( min_lum < fp ) ? min_lum : fp;	/// min Luminance in the scene
        sum += log( 2.3e-5f + fp );				    /// contrast constant in Tumblin paper
	}

	// maximum luminance
	*maxLum = max_lum;
	// minimum luminance
	*minLum = min_lum;

	// average log luminance
	double avgLogLum = (double)sum / (double)srcsz ;

	// world adaptation luminance
	*worldLum = (float)exp( avgLogLum );

	return true;
}

static bool rec709GammaCorrection( unsigned short* src, unsigned srcsz, const float gammaval )
{
    if ( ( src == NULL ) || ( srcsz == 0 ) )
    {
        return false;
    }

	float slope = 4.5f;
	float start = 0.018f;

	const float fgamma = (float)( ( 0.45f / gammaval ) * 2.0f );

	if( gammaval >= 2.1f )
    {
		start = (float)( 0.018f / ( ( gammaval - 2.0f ) * 7.5f ) );
		slope = (float)( 4.5f   * ( ( gammaval - 2.0f ) * 7.5f ) );
	}
	else
    if( gammaval <= 1.9f )
    {
		start = (float)( 0.018f * ( ( 2.0f - gammaval ) * 7.5f ) );
		slope = (float)( 4.5f   / ( ( 2.0f - gammaval ) * 7.5f ) );
	}

    #pragma omp parallel for
	for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
		float pixel = src[ cnt ];
        pixel = ( pixel <= start ) ? pixel * slope : ( 1.099f * pow( pixel, fgamma ) - 0.099f );
        src[ cnt ] = (unsigned short)pixel;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////

bool RAWImageToolKit::tmoDrago03( unsigned short *src, unsigned srcsz, float maxf, float normalf, double gamma, double exposure)
{
    if ( ( src == NULL ) || ( srcsz == 0 ) )
        return false;

    // Convert unsigned short to float ...
    float* convf = new float[ srcsz ];
    if ( convf == NULL )
        return false;

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        convf[ cnt ] = (float) src[ cnt ] / maxf;
    }

	// default algorithm parameters
	const float biasParam = 0.85f;
	const float expoParam = (float)pow( 2.0f, exposure ); //default exposure is 1, 2^0

	float maxLum = 0.0f;
	float minLum = 0.0f;
	float avgLum = 0.0f;

    if ( getLuminance( convf, srcsz, &maxLum, &minLum, &avgLum ) == true )
    {
        toneMappingDrago03( convf, srcsz, maxLum, avgLum, biasParam, expoParam );

        if ( gamma != 1.0f )
        {
            rec709GammaCorrection( src, srcsz, gamma );
        }
    }
    else
    {
        delete[] convf;

        return false;
    }

    // Check overflown
    float recoverf = 1.0f;

    for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        if ( convf[ cnt ] > recoverf )
        {
            recoverf = convf[ cnt ];
        }
    }

    if ( recoverf > 1.0f )
    {
        float divf = 1.0f / recoverf;
        #pragma omp parallel for
        for( unsigned cnt=0; cnt<srcsz; cnt++ )
        {
            convf[ cnt ] *= divf;
        }
    }

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        //src[cnt] = MIN( (unsigned short)( convf[ cnt ] * normalf ), normalf );

        //Fixing clipped over exposures.
        float tmpfv = convf[ cnt ] * normalf;

        if ( tmpfv > 65535.f )
        {
            tmpfv = 65535.f;
        }

        src[ cnt ] = (unsigned short)tmpfv;
    }

    delete[] convf;

    return true;
}
