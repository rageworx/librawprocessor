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
#include "rawimgtk_tmoreinhard05.h"

////////////////////////////////////////////////////////////////////////////////

static unsigned convertUStoF( const unsigned short* src, unsigned srcsz, float** out )
{
    if ( ( src != NULL ) && ( srcsz > 0 ) )
    {
        float* outf = new float[ srcsz ];
        if ( outf != NULL )
        {
            #pragma omp parallel for
            for( unsigned cnt=0; cnt<srcsz; cnt++ )
            {
                outf[ cnt ] = (float)src[ cnt ];
            }
        }
        *out = outf;

        return srcsz;
    }

    return 0;
}

static bool luminanceFromY( float* src, unsigned srcsz,
                            float *maxLum, float *minLum, float *Lav, float *Llav)
{
    if ( ( src == NULL ) || ( srcsz == 0 ) )
    {
        return false;
    }

	float max_lum = -1e20f;
	float min_lum = 1e20f;
	double sumLum = 0;
	double sumLogLum = 0;

	#pragma omp parallel for
    for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        const float Y = src[ cnt ];
        max_lum = (max_lum < Y) ? Y : max_lum;				// max Luminance in the scene
        min_lum = ((Y > 0) && (min_lum < Y)) ? min_lum : Y;	// min Luminance in the scene
        sumLum += Y;										// average luminance
        sumLogLum += log(2.3e-5F + Y);						// contrast constant in Tumblin paper
    }

	// maximum luminance
	*maxLum = max_lum;
	// minimum luminance
	*minLum = min_lum;
	// average luminance
	*Lav = (float)(sumLum / srcsz);
	// average log luminance, a.k.a. world adaptation luminance
	*Llav = (float)exp(sumLogLum / srcsz);

	return true;
}

static bool toneMappingReinhard05( float* src, float* lumiY, unsigned srcsz,
                                   float ins, float cont, float adapt, float cc = 0.0f )
{
	double Cav = 0;		// channel average
	float Lav = 0;		// average luminance
	float Llav = 0;		// log average luminance
	float minLum = 1;	// min luminance
	float maxLum = 1;	// max luminance

	float L;		    // pixel luminance
	float I_g, I_l;     // global and local light adaptation
	float I_a;		    // interpolated pixel light adaptation
	float keyf;	        // key (low-key means overall dark image, high-key means overall light image)

	// Check intensity range ...
	if( ins < -8.0f )
        ins = -8.0f;

    if( ins > 8.0f )
        ins = 8.0f;
    // ..........................

    // Check contrast range
    if( cont < 0.3 )
        cont = 0.3f;

    if( cont >= 1.0 )
        cont = 0.99999f;
    // ..........................

    // Check adaptation range
    if( adapt < 0.0f )
        adapt = 0.0f;

    if( adapt > 2.0f )
        adapt = 2.0f;
    // ..........................

    // Check color_correction range
    if( cc < 0.0f )
        cc = 0.0f;

    if( cc > 1.0f )
        cc = 1.0f;
    // .........................

	ins = exp(-ins);

	if( ( cont == 0 ) || ( ( adapt != 1 ) && ( cc != 1 ) ) )
    {
		luminanceFromY( lumiY, srcsz, &maxLum, &minLum, &Lav, &Llav );
		keyf = (log(maxLum) - Llav) / (log(maxLum) - log(minLum));
		if( keyf < 0.0 )
		{
			keyf = ( log( maxLum ) - log( Llav ) ) /
			       ( log( maxLum ) - log( minLum ) );

			if( keyf < 0.0f ) cont = 0.3f;
		}
	}

	cont = ( cont > 0.0 ) ? cont : (float)( 0.3f + 0.7f * pow( keyf, 1.4f ) );

	float max_color = -1e6F;
	float min_color = +1e6F;

	// tone map image
	if( ( adapt == 1.0f ) && ( cc == 0.0f ) )
    {
		// when using default values, use a fastest code

		for( unsigned cnt=0; cnt<srcsz; cnt++ )
        {
            I_a = lumiY[ cnt ];	// luminance(x, y)

            src[ cnt ] /= ( src[ cnt ] + pow( ins * I_a, cont ) );

            max_color = (src[ cnt ] > max_color) ? src[ cnt ] : max_color;
            min_color = (src[ cnt ] < min_color) ? src[ cnt ] : min_color;
		}
	}
	else
    {
		// complete algorithm

		// Average of color/Level
		if( ( adapt != 1.0f ) && ( cc != 0.0f ) )
        {
			for(unsigned cnt=0; cnt<srcsz; cnt++ )
			{
				Cav += src[ cnt ];
			}

            Cav /= srcsz;
		}

		// perform tone mapping
		for( unsigned cnt=0; cnt<srcsz; cnt++ )
        {
            L = lumiY[ cnt ];	// luminance(x, y)

            I_l = cc * src[ cnt ] + ( 1 - cc ) * L;
            I_g = cc * Cav + ( 1 - cc ) * Lav;
            I_a = adapt * I_l + ( 1 - adapt ) * I_g;
            src[ cnt ] /= ( src[ cnt ] + pow( ins * I_a, cont ) );

            max_color = ( src[ cnt ] > max_color) ? src[ cnt ] : max_color;
            min_color = ( src[ cnt ] < min_color) ? src[ cnt ] : min_color;
		}
	}

	// normalize intensities

	if( max_color != min_color )
    {
		const float range = max_color - min_color;

		#pragma omp parallel for
		for( unsigned cnt=0; cnt<srcsz; cnt++ )
        {
            src[ cnt ] = ( src[ cnt ] - min_color ) / range;
		}
	}

	return true;
}

////////////////////////////////////////////////////////////////////////////////

bool RAWImageToolKit::tmoReinhard2005( unsigned short* src, unsigned srcsz,
                                       float maxf, float normalf,
                                       float intensity, float contrast, float adaptation,
                                       float colorcorrection )
{
    if ( ( src == NULL ) || ( srcsz == 0 ) )
        return false;

    float* convf = new float[ srcsz ];
    if ( convf == NULL )
        return false;

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        convf[ cnt ] = (float) src[ cnt ] / maxf;
    }

    float* convy = NULL;

    if ( convertUStoF( src, srcsz, &convy ) > 0 )
    {

        toneMappingReinhard05( convf,
                               convy,
                               srcsz,
                               intensity,
                               contrast,
                               adaptation,
                               colorcorrection );
        delete[] convy;
    }
    else
    {
        delete[] convf;

        return false;
    }

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<srcsz; cnt++ )
    {
        //src[cnt] = MIN( (unsigned short)( convf[ cnt ] * (float)normalf ), normalf );
        float tmpfv = convf[ cnt ] * (float)normalf;

        if ( tmpfv > 65535.f )
        {
            tmpfv = 65535.f;
        }

        src[cnt] = (unsigned short)tmpfv;
    }

    delete[] convf;

    return true;
}
