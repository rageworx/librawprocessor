#include "rawfilters.h"

#include <cmath>
#ifdef USE_INTERNAL_FFT
#include <complex>
#include <iostream>
#include <valarray>
#endif /// of USE_INTERNAL_FFT
#include <cstring>
#include <string>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;

#ifdef USE_INTERNAL_FFT
const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
#endif /// of USE_INTERNAL_FFT

const float matrixdata_blur[] =
{
    0.0, 0.2,  0.0,
    0.2, 0.4,  0.2,
    0.0, 0.2,  0.0
};

const float matrixdata_blurmore[] =
{
    0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 0.0,
    1.0, 1.0, 4.0, 1.0, 1.0,
    0.0, 1.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0
};

const float matrixdata_gaussianblur[] =
{
    1.0,  4.0,  7.0,  4.0,  1.0,
    4.0, 16.0, 26.0, 16.0,  4.0,
    7.0, 26.0, 41.0, 26.0,  7.0,
    4.0, 16.0, 26.0, 16.0,  4.0,
    1.0,  4.0,  7.0,  4.0,  1.0
};

const float matrixdata_edgeonly[] =
{
    0.0, 0.0, -1.0, 0.0, 0.0,
    0.0, 0.0, -1.0, 0.0, 0.0,
    0.0, 0.0,  2.0, 0.0, 0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,
    0.0, 0.0,  0.0, 0.0, 0.0
};

const float matrixdata_edgemore[] =
{
   -1.0,  0.0,  0.0,  0.0,  0.0,
    0.0, -2.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  6.0,  0.0,  0.0,
    0.0,  0.0,  0.0, -2.0,  0.0,
    0.0,  0.0,  0.0,  0.0, -1.0
};

const float matrixdata_edgestrong[] =
{
   -1.0, -1.0, -1.0,
   -1.0,  8.0, -1.0,
   -1.0, -1.0, -1.0
};

const float matrixdata_sharpen[] =
{
    0.0, -1.0,  0.0,
   -1.0,  5.0, -1.0,
    0.0, -1.0,  0.0
};

const float matrixdata_sharpenhard[] =
{
   -1.0, -1.0, -1.0,
   -1.0,  9.0, -1.0,
   -1.0, -1.0, -1.0
};

const float matrixdata_sharpensubtle[] =
{
  -1.0, -1.0, -1.0, -1.0, -1.0,
  -1.0,  2.0,  2.0,  2.0, -1.0,
  -1.0,  2.0,  8.0,  2.0, -1.0,
  -1.0,  2.0,  2.0,  2.0, -1.0,
  -1.0, -1.0, -1.0, -1.0, -1.0
};

const float matrixdata_unsharpenmask[] =
{
  0.0, 0.2, 0.0,
  0.2, 0.2, 0.2,
  0.0, 0.2, 0.0
};


////////////////////////////////////////////////////////////////////////////////

void stdvector_f_put_array( vector< float > &v, const float* a )
{
    if ( a != NULL )
    {
        for( unsigned cnt=0; cnt<v.size(); cnt++ )
        {
            v[cnt] = a[cnt];
        }
    }
}

#ifdef USE_INTERNAL_FFT

// FFT and IFFT codes from
// https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
//
// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
vool ct_fft(CArray& x)
{
    const size_t N = x.size();

    if (N <= 1)
        return false;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }

    return true;
}

// inverse fft (in-place)
bool ct_ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);

    // forward fft
    if ( ct_fft( x ) == false )
        return false;

    // conjugate the complex numbers again
    x = x.apply(std::conj);

    // scale the numbers
    x /= x.size();

    return true;
}
#endif /// of USE_INTERNAL_FFT

////////////////////////////////////////////////////////////////////////////////

RAWProcessor::FilterConfig* RAWImageFilterKit::GetPresetFilter( const char* name )
{
    string filtername = name;

    RAWProcessor::FilterConfig* newfilter = new RAWProcessor::FilterConfig;

    if ( filtername == RAWFILTER_PRESET_BLUR )
    {
        newfilter->width  = 3;
        newfilter->height = 3;
        newfilter->factor = 1.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 9 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_blur );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_BLURMORE )
    {
        newfilter->width  = 5;
        newfilter->height = 5;
        newfilter->factor = 1.0 / 13.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 25 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_blurmore );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_GAUSSIAN )
    {
        newfilter->width  = 5;
        newfilter->height = 5;
        newfilter->factor = 1.0 / 13.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 25 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_gaussianblur );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_EDGE )
    {
        newfilter->width  = 5;
        newfilter->height = 5;
        newfilter->factor = 1.0 / 13.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 25 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_edgeonly );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_EDGE2 )
    {
        newfilter->width  = 5;
        newfilter->height = 5;
        newfilter->factor = 1.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 25 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_edgemore );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_EDGE3 )
    {
        newfilter->width  = 3;
        newfilter->height = 3;
        newfilter->factor = 1.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 9 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_edgestrong );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_SHARPEN )
    {
        newfilter->width  = 3;
        newfilter->height = 3;
        newfilter->factor = 1.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 9 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_sharpen );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_SHARPEN2 )
    {
        newfilter->width  = 3;
        newfilter->height = 3;
        newfilter->factor = 1.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 9 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_sharpenhard );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_SHARPENSB )
    {
        newfilter->width  = 5;
        newfilter->height = 5;
        newfilter->factor = 1.0 / 8.0;
        newfilter->bias   = 0.0;
        newfilter->matrix.resize( 9 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_sharpensubtle );

        return newfilter;
    }

    if ( filtername == RAWFILTER_PRESET_UNSHARPENM )
    {
        newfilter->width  = 3;
        newfilter->height = 3;
        newfilter->factor = 1.0;
        newfilter->bias   = 2.0;
        newfilter->matrix.resize( 9 );
        stdvector_f_put_array( newfilter->matrix, matrixdata_unsharpenmask );

        return newfilter;
    }

    delete newfilter;

    return NULL;
}

void RAWImageFilterKit::RemoveFilter( RAWProcessor::FilterConfig* fp )
{
    if ( fp != NULL )
    {
        fp->matrix.clear();

        delete fp;
    }
}

bool RAWImageFilterKit::ApplyLowFreqFilter( unsigned short* ptr, unsigned w, unsigned h, unsigned fsz, unsigned iter )
{
    if ( ptr == NULL )
        return false;

    unsigned imgsz = w * h;

    unsigned short* imgSRC = new unsigned short[ imgsz ];

    if ( imgSRC == NULL )
        return false;

    unsigned short* imgLF = new unsigned short[ imgsz ];

    if ( imgLF == NULL )
    {
        delete[] imgSRC;

        return false;
    }

    memcpy( imgSRC, ptr, imgsz * sizeof( unsigned short ) );
    memset( imgLF, 0 , imgsz * sizeof( unsigned short ) );

    long     sfsz  = (long)fsz;
    unsigned itcnt = 0;
    unsigned cntw  = 0;
    unsigned cnth  = 0;
    int      maxw  = (int)w;
    int      maxh  = (int)h;

    for( itcnt=0; itcnt<iter; itcnt++ )
    {
        #pragma omp parallel for private( cntw )
        for( cnth=0; cnth<h; cnth++ )
        {
            for( cntw=0; cntw<w; cntw++ )
            {
                unsigned long long sum = 0;
                unsigned sumsz = 0;

                for( int lp1=(-sfsz/2); lp1<=(sfsz/2); lp1++ )
                {
                    for( int lp2=(-sfsz/2); lp2<=(sfsz/2); lp2++ )
                    {
                        long long chkh = cnth + lp1;
                        long long chkv = cntw + lp2;

                        if ( ( chkv >= 0 ) && ( chkh >= 0 ) &&
                             ( chkv < maxw ) && ( chkh < maxh ) )
                        {
                            sum += imgSRC[ chkh * w + chkv ];
                            sumsz++;
                        }
                    }
                }

                imgLF[ cnth * w + cntw ] = sum / sumsz;
            }
        }

        if ( itcnt == ( iter - 1 ) )
        {
            memcpy( ptr, imgLF, imgsz * sizeof( unsigned short ) );
        }
        else
        {
            memcpy( imgSRC, imgLF, imgsz * sizeof( unsigned short ) );
        }
    }

    delete[] imgSRC;
    delete[] imgLF;

    return true;
}

bool RAWImageFilterKit::ApplyEdgeLowFreqFilter( unsigned short* ptr, unsigned w, unsigned h, unsigned fsz )
{
    if ( ptr == NULL )
        return false;

    unsigned imgsz = w * h;

    unsigned short* imgSRC = new unsigned short[ imgsz ];

    if ( imgSRC == NULL )
        return false;

    unsigned short* imgLF = new unsigned short[ imgsz ];

    if ( imgLF == NULL )
    {
        delete[] imgSRC;

        return false;
    }

    memcpy( imgSRC, ptr, imgsz * sizeof( unsigned short ) );
    memset( imgLF, 0 , imgsz * sizeof( unsigned short ) );

    long halfsz = (int)fsz / 2;

    if ( halfsz == 0 )
    {
        delete[] imgSRC;
        delete[] imgLF;

        return false;
    }

    unsigned cntw = 0;
    unsigned cnth = 0;
    unsigned powfsz = ( fsz * fsz );

    #pragma omp parallel for private( cntw )
    for( cnth=halfsz; cnth<(h-halfsz); cnth++)
    {
        for( cntw=halfsz; cntw<(w-halfsz); cntw++ )
        {
            unsigned long long sum = 0;

            for( int lp1=(-halfsz); lp1<(halfsz+1); lp1++ )
            {
                for( int lp2=(-halfsz); lp2<(halfsz+1); lp2++ )
                {
                    sum += imgSRC[ ( lp1 + cnth ) * w + ( lp2 + cntw ) ];
                }
            }

            imgLF[ cnth * w + cntw ] = sum / powfsz;
        }
    }

    memcpy( ptr, imgLF, imgsz * sizeof( unsigned short ) );

    delete[] imgSRC;
    delete[] imgLF;

    return true;
}

bool RAWImageFilterKit::ApplyAnisotropicFilter( unsigned short* ptr, unsigned w, unsigned h, unsigned fstr, unsigned fparam )
{
    if ( ptr == NULL )
        return false;

    if ( fparam == 0 )
        return false;

    unsigned imgsz = w * h;

    const float diamf  = pow( (float)( fparam / 4 ), 2.0f );
    const float lambda = 0.005f;

    float* imgSRC = new float[ imgsz ];

    if ( imgSRC == NULL )
        return false;

    float* imgAF = new float[ imgsz ];

    if ( imgAF == NULL )
    {
        delete[] imgSRC;

        return false;
    }

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        // Normalize to vector.
        imgSRC[ cnt ] = (float)ptr[ cnt ] / 65536.0f;
    }

    memset( imgAF, 0 , imgsz * sizeof( float ) );

    for( unsigned cnt=0; cnt<fstr; cnt++ )
    {
        float grad[8] = {0.0f};
        float gcol[8] = {0.0f};

        unsigned cnth = 0;
        unsigned cntw = 0;

        #pragma omp parallel for private( cntw, grad, gcol )
        for( cnth=0; cnth<h; cnth++ )
        {
            for( cntw=0; cntw<w; cntw++ )
            {
                unsigned pos = cnth * w + cntw;

                if ( ( cntw > 0 ) && ( cntw < ( w - 1 ) ) &&
                     ( cnth > 0 ) && ( cnth < ( h - 1 ) ) )
                {
                    float decf   = imgSRC[ pos ];

                    grad[0] = imgSRC[ ( cnth - 1 ) * w + ( cntw - 1 ) ] - decf;
                    grad[1] = imgSRC[ ( cnth - 1 ) * w + ( cntw + 0 ) ] - decf;
                    grad[2] = imgSRC[ ( cnth - 1 ) * w + ( cntw + 1 ) ] - decf;
                    grad[3] = imgSRC[ ( cnth + 0 ) * w + ( cntw - 1 ) ] - decf;
                    grad[4] = imgSRC[ ( cnth + 0 ) * w + ( cntw + 1 ) ] - decf;
                    grad[5] = imgSRC[ ( cnth + 1 ) * w + ( cntw - 1 ) ] - decf;
                    grad[6] = imgSRC[ ( cnth + 1 ) * w + ( cntw + 0 ) ] - decf;
                    grad[7] = imgSRC[ ( cnth + 1 ) * w + ( cntw + 1 ) ] - decf;

                    for( unsigned gcnt=0; gcnt<8; gcnt ++ )
                    {
                        if ( diamf >= ( pow ( grad[ gcnt ], 2.0f ) / 9.0f ) )
                        {
                            gcol[ gcnt ] = 1.0f - ( pow( grad[ gcnt ] , 2.0f ) / ( 9.0f * diamf ) );
                        }
                        else
                        {
                            gcol[ gcnt ] = 1.0f - ( 1.0f + ( pow( grad[ gcnt ], 2.0f ) / diamf ) );
                        }
                    }

                    float sumf = 0.0f;

                    for( unsigned gcnt=0; gcnt<8; gcnt ++ )
                    {
                        sumf += gcol[ gcnt ] * grad[ gcnt ];
                    }

                    imgAF[ pos ] = imgSRC[ pos ] + lambda * sumf;
                }
                else
                {
                    imgAF[ pos ] = imgSRC[ pos ];
                }
            }
        }

        memcpy( imgSRC, imgAF, imgsz * sizeof( float ) );
    }

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<imgsz; cnt++ )
    {
        // Restore .
        ptr[ cnt ] = (unsigned short)( imgSRC[ cnt ] * 65536.0f );
    }

    delete[] imgSRC;
    delete[] imgAF;

    return true;
}
