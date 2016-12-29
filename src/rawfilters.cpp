#include "rawfilters.h"
#include <string>

/*
        typedef struct
        {
            unsigned char   width;
            unsigned char   height;
            float           factor;
            float           bias;
            std::vector< float > matrix;
        }FilterConfig;
*/

using namespace std;

const float matrixdata_blur[] =
{
   0.0, 0.2,  0.0,
   0.2, 0.2,  0.2,
   0.0, 0.2,  0.0
};

const float matrixdata_blurmore[] =
{
  0.0, 0.0, 1.0, 0.0, 0.0,
  0.0, 1.0, 1.0, 1.0, 0.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
  0.0, 1.0, 1.0, 1.0, 0.0,
  0.0, 0.0, 1.0, 0.0, 0.0
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
