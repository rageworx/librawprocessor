#ifdef __APPLE__
    #include <sys/uio.h>
#endif // __APPLE__

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <string>
#include <cstring>

#include <list>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

#include "rawprocessor.h"
#include "rawscale.h"
#include "rawimgtk.h"
#include "rawfilters.h"
#include "rawimgtk_tmodrago03.h"
#include "rawimgtk_tmoreinhard05.h"
#include "minmax.h"
#include "stdunicode.h"

////////////////////////////////////////////////////////////////////////////////

using namespace std;

////////////////////////////////////////////////////////////////////////////////

#ifdef __APPLE__
    #ifdef UNICODE
        #undef UNICODE
    #endif
#endif // __APPLE__

#ifdef WCHAR_SUPPORTED
    #define _TSTRING        wstring
    #define _TCM2W( _x_ )   convertM2W( (const char*)_x_ )
    #define _TCW2M( _x_ )   convertW2M( (const wchar_t*)_x_ )
    #ifndef _T
        #define _T( _x_ )       L##_x_
    #endif /// of _T
    #define DEF_MEMORY_LOADED   L"//MEMORY_LOAD//"
#else
    #define _TSTRING        string
    #define _TCM2W( _x_ )   _x_
    #define _TCW2M( _x_ )   _x_
    #ifndef _T
        #define _T( _x_ )       _x_
    #endif /// of _T
    #define DEF_MEMORY_LOADED   "//MEMORY_LOAD//"
#endif //// of WCHAR_SUPPORTED

#define DEF_RAW_I_HEIGHT    ( 1024 )
#define DEF_window    ( 1.0f )

#define DEF_PIXEL32_MAX     ( 0xFFFFFFFFF )
#define DEF_PIXEL16_MAX     ( 0xFFFF )
#define DEF_PIXEL8_MAX      ( 0xFF )

#define DEF_CALC_F_DMAX     ( 4294967295.f )
#define DEF_CALC_I_DMAX     ( 4294967295 )
#define DEF_CALC_F_WMAX     ( 65535.0f )
#define DEF_CALC_I_WMAX     ( 65535 )
#define DEF_CALC_F_BMAX     ( 255.0f )
#define DEF_CALC_I_BMAX     ( 255 )

#define DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY     0,9,60,160

////////////////////////////////////////////////////////////////////////////////

#define rawimgtk            RAWImageToolKit
#define rawimgfk            RAWImageFilterKit

////////////////////////////////////////////////////////////////////////////////


#define BYTE_SWAP_16(n) ((((n) >> 8) & 0xffu) | (((n) & 0xffu) << 8))
#define BYTE_SWAP_32(x)  \
     ((((x) & 0xff000000u) >> 24) | (((x) & 0x00ff0000u) >>  8) |  \
      (((x) & 0x0000ff00u) <<  8) | (((x) & 0x000000ffu) << 24))

////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    uint32_t min_v;
    uint32_t max_v;
}minmaxpair;

typedef struct
{
    float min_v;
    float max_v;
}minmaxfpair;

////////////////////////////////////////////////////////////////////////////////

bool RAWProcessor_sortcondition( int i,int j )
{
    return ( i < j );
}

int RAWProcessor_uscompare( const void * a, const void * b )
{
    if( *(const float*)a < *(const float*)b )
        return -1;

    return *(const float*)a > *(const float*)b;
}

int RAWProcessor_fcompare( const void * a, const void * b )
{
    if( *(const float*)a < *(const float*)b )
        return -1;

    return *(const float*)a > *(const float*)b;
}

inline void _bswap2( void* ptr )
{
    uint8_t* pcast = (uint8_t*)ptr;
    uint8_t  tmpstr[2] = { pcast[0], pcast[1] };
    pcast[0] = tmpstr[1];
    pcast[1] = tmpstr[2];
}

inline void _bswap4( void* ptr )
{
    uint8_t* pcast = (uint8_t*)ptr;
    uint8_t  tmpstr[4] = { pcast[0], pcast[1], pcast[2], pcast[3] };
    pcast[0] = tmpstr[3];
    pcast[1] = tmpstr[2];
    pcast[2] = tmpstr[1];
    pcast[3] = tmpstr[0];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

RAWProcessor::RAWProcessor()
 : raw_loaded(false),
   pixel_arrays(NULL),
   pixel_arrays_paddedsz(0),
   window_max(0),
   pixel_bpp(32),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(0),
   img_width(0)
{
    resetWindow();
}

RAWProcessor::RAWProcessor( const char* raw_file, uint32_t height )
 : raw_loaded(false),
   pixel_arrays(NULL),
   pixel_arrays_paddedsz(0),
   window_max( 0 ),
   pixel_bpp(10),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(height),
   img_width(0)
{
    resetWindow();

    if ( raw_file != NULL )
    {
        Load( raw_file, TRANSFORM_NONE, height );
    }
}

#ifdef WCHAR_SUPPORTED
RAWProcessor::RAWProcessor( const wchar_t* raw_file, uint32_t height )
 : raw_loaded(false),
   pixel_arrays(NULL),
   pixel_arrays_paddedsz(0),
   window_max( 0 ),
   pixel_bpp(10),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(height),
   img_width(0)
{
    resetWindow();

    if ( raw_file != NULL )
    {
        Load( raw_file, TRANSFORM_NONE, height );
    }
}
#endif /// of WCHAR_SUPPORTED

RAWProcessor::~RAWProcessor()
{
    Unload();
    resetWindow();
}

void RAWProcessor::CutoffLevels( float minv, float maxv )
{
    #pragma omp parallel for
    for( size_t cnt=0; cnt<pixel_arrays_paddedsz; cnt++ )
    {
        if ( pixel_arrays[cnt] > maxv )
        {
            pixel_arrays[cnt] = maxv;
        }
        else
        if ( pixel_arrays[cnt] < minv )
        {
            pixel_arrays[cnt] = minv;
        }
    }
}

void RAWProcessor::CutoffLevelsRanged( float minv, float maxv, float valmin, float valmax )
{
    #pragma omp parallel for
    for( size_t cnt=0; cnt<pixel_arrays_paddedsz; cnt++ )
    {
        if ( pixel_arrays[cnt] > maxv )
        {
            pixel_arrays[cnt] = valmax;
        }
        else
        if ( pixel_arrays[cnt] < minv )
        {
            pixel_arrays[cnt] = valmin;
        }
    }
}

void RAWProcessor::Version( char** retverstr )
{
    if ( retverstr == NULL )
    {
        return;
    }

    int32_t retia[4] = { DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY };
    char retvstr[32] = {0};

    sprintf( retvstr, "%d.%d.%d.%d", retia[0], retia[1], retia[2], retia[3] );

    int32_t retvstrsz = strlen( retvstr );

    *retverstr = (char*)malloc( retvstrsz + 1 );
    if ( *retverstr != NULL )
    {
        memset( *retverstr, 0, retvstrsz + 1 );
        strcpy( *retverstr, retvstr );
    }
}

void RAWProcessor::Version( int* retverints )
{
    if ( retverints == NULL )
    {
        return;
    }

    int32_t retia[4] = { DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY };

    for( size_t cnt=0; cnt<4; cnt++ )
    {
        retverints[ cnt ] = retia[ cnt ];
    }
}

#ifdef WCHAR_SUPPORTED
bool RAWProcessor::Load( const wchar_t* raw_file, uint32_t trnsfm, size_t height, uint32_t dtype, bool byteswap )
{
    if ( height <= 0 )
        return false;

    const char* pconv = convertW2M( raw_file );

    if ( pconv != NULL )
    {
        string fname = pconv;

        return Load( fname.c_str(), trnsfm, height );
    }

    return false;
}
#endif // of WCHAR_SUPPORTED

bool RAWProcessor::Load( const char* raw_file, uint32_t trnsfm, size_t height, uint32_t dtype, bool byteswap )
{
    if ( height == 0 )
        return false;

    string fname = raw_file;

    FILE* fp = fopen( fname.c_str(), "rb" );

    if ( fp != NULL )
    {
        fseek( fp, 0L, SEEK_END );
        size_t fsize = ftell( fp );
        rewind( fp );

        if ( fsize == 0 )
        {
            fclose( fp );
            raw_loaded = false;

            return false;
        }

        pixel_min_level = 0.f;
        pixel_med_level = 0.f;
        pixel_max_level = 0.f;
        img_height = height;
        size_t readsz = 1;
        size_t arraysz = 0;

        switch( dtype )
        {
            case DATATYPE_UINT8:
                arraysz = fsize;
                break;

            case DATATYPE_UINT16:
                readsz = sizeof( uint16_t );
                arraysz = fsize / readsz;
                break;

            case DATATYPE_UINT32:
                readsz = sizeof( uint32_t );
                arraysz = fsize / readsz;
                break;

            case DATATYPE_FLOAT:
                readsz = sizeof( float );
                arraysz = fsize / readsz;
                break;

            default:
                return false;
        }

        size_t blancsz = 0;

        if (  arraysz > 0 )
        {
            // set temporary width ...
            img_width = arraysz / img_height;

            // Check omit pixels ...
            if ( arraysz < ( img_width * img_height ) )
            {
                blancsz = ( img_width * img_height ) - arraysz;
            }
        }

        if ( pixel_arrays != NULL )
        {
            delete[] pixel_arrays;
            pixel_arrays = NULL;
        }

        pixel_arrays_paddedsz = arraysz + blancsz;
        pixel_arrays = new float[pixel_arrays_paddedsz];

        resetWindow();

        // To save memory, it doesn't using direct load to memory.
        for( size_t cnt=0; cnt<arraysz; cnt++)
        {
            char  chardata[4] = {0,0,0,0};
            float convdata = 0.f;

            memset( chardata, 0, 4 );
            fread( chardata, readsz, 1, fp );

            // convert each data type to floating point 4 byte.
            switch( dtype )
            {
                case DATATYPE_UINT8:
                    convdata = (float)chardata[0] / DEF_CALC_F_BMAX;
                    break;

                case DATATYPE_UINT16:
                    {
                        uint16_t usdata = *(uint16_t*)chardata;

                        if ( byteswap == false )
                        {
                            convdata = (float)usdata / DEF_CALC_F_WMAX;
                        }
                        else
                        {
                            _bswap2( &usdata );
                            convdata = (float)usdata / DEF_CALC_F_WMAX;
                        }
                    }
                    break;

                case DATATYPE_UINT32:
                    {
                        uint32_t usdata = *(uint16_t*)chardata;

                        if ( byteswap == false )
                        {
                            convdata = (float)usdata / DEF_CALC_F_DMAX;
                        }
                        else
                        {
                            _bswap4( &usdata );
                            convdata = (float)usdata / DEF_CALC_F_DMAX;
                        }
                    }
                    break;

                case DATATYPE_FLOAT:
                    {
                        convdata = *(float*)chardata;

                        if ( byteswap == true )
                        {
                            _bswap4( &convdata );
                        }
                    }
                    break;
            }

            pixel_arrays[cnt] = convdata;

            if ( window_max < convdata )
            {
                window_max = convdata;
            }
        }

        // actual pixel count, not paaded.
        pixel_arrays_srcsz = arraysz;
        pixel_med_level = ( pixel_max_level + pixel_min_level ) / 2.f;

        fclose( fp );

        calcWindow();
        analyse();

        raw_loaded = true;

        ApplyTransform( trnsfm );

        return raw_loaded;
    }

    return false;
}

bool RAWProcessor::LoadFromMemory( void* buffer, size_t bufferlen, uint32_t trnsfm, size_t height, uint32_t dtype, bool byteswap )
{
    if ( height <= 0 )
        return false;

    if ( ( buffer != NULL ) && ( bufferlen > 0 ) )
    {
        const uint8_t* pbuff = (const uint8_t*)buffer;
        pixel_min_level = 0;
        pixel_med_level = 0;
        pixel_max_level = 0;
        img_height = height;
        size_t readsz = 1;
        size_t arraysz = 0;

        switch( dtype )
        {
            case DATATYPE_UINT8:
                arraysz = bufferlen;
                break;

            case DATATYPE_UINT16:
                readsz = sizeof( uint16_t );
                arraysz = bufferlen / readsz;
                break;

            case DATATYPE_UINT32:
                readsz = sizeof( uint32_t );
                arraysz = bufferlen / readsz;
                break;

            case DATATYPE_FLOAT:
                readsz = sizeof( float );
                arraysz = bufferlen/ readsz;
                break;

            default:
                return false;
        }

#ifdef DEBUG
        printf( "RAWProcessor::LoadFromMemory, arraysz = %zu * readsz = %zu.\n", arraysz, readsz );
        fflush( stdout );
#endif /// of DEBUG

        size_t blancsz   = 0;

        if (  arraysz > 0 )
        {
            // set temporary width ...
            img_width = arraysz / img_height;

            // Check omit pixels ...
            if ( arraysz < ( img_width * img_height ) )
            {
                blancsz = ( img_width * img_height ) - arraysz;
            }
        }

        if ( pixel_arrays != NULL )
        {
            delete[] pixel_arrays;
            pixel_arrays = NULL;
        }

        pixel_arrays_paddedsz = arraysz + blancsz;
        pixel_arrays = new float[pixel_arrays_paddedsz];

        resetWindow();

        // To save memory, it doesn't using direct load to memory.
        //#pragma omp parallel for private(window_max)
        for( size_t cnt=0; cnt<arraysz; cnt++)
        {
            float convdata = 0;

            switch( dtype )
            {
                case DATATYPE_UINT8:
                    convdata = (float)pbuff[cnt] / DEF_CALC_F_BMAX;
                    break;

                case DATATYPE_UINT16:
                    {
                        uint16_t* ptr = (uint16_t*)pbuff;
                        uint16_t  usdata = ptr[cnt];

                        if ( byteswap == false )
                        {
                            convdata = (float)usdata / DEF_CALC_F_WMAX;
                        }
                        else
                        {
                            _bswap2( &usdata );
                            convdata = (float)usdata / DEF_CALC_F_WMAX;
                        }
                    }
                    break;

                case DATATYPE_UINT32:
                    {
                        uint32_t* ptr = (uint32_t*)pbuff;
                        uint32_t  usdata = ptr[cnt];

                        if ( byteswap == false )
                        {
                            convdata = (float)usdata / DEF_CALC_F_DMAX;
                        }
                        else
                        {
                            _bswap4( &usdata );
                            convdata = (float)usdata / DEF_CALC_F_DMAX;
                        }
                    }
                    break;

                case DATATYPE_FLOAT:
                    {
                        float* ptr = (float*)pbuff;
                        convdata = ptr[cnt];

                        if ( byteswap == true )
                        {
                            _bswap4( &convdata );
                        }
                    }
                    break;
            }

            pixel_arrays[ cnt ] = convdata;

            if ( window_max < convdata )
            {
                window_max = convdata;
            }
        }

#ifdef DEBUG
        printf( "load from memory #1\n" );
        fflush( stdout );
#endif
        calcWindow();
        analyse();

#ifdef DEBUG
        printf( "load from memory #2\n" );
        fflush( stdout );
#endif

        raw_loaded = true;

        ApplyTransform( trnsfm );

        raw_file_name = DEF_MEMORY_LOADED;

        return true;
    }
#ifdef DEBUG
    else
    {
        printf( "ERROR: Nothing to load from memory, p=%p, s=%zu\n", buffer,  bufferlen );
        fflush( stdout );
    }
#endif

    return false;
}

#ifdef WCHAR_SUPPORTED
bool RAWProcessor::Reload( const wchar_t* raw_file, uint32_t trnsfm, uint32_t height )
{
    return Load( raw_file, trnsfm, height );
}
#endif /// of WCHAR_SUPPORTED

bool RAWProcessor::Reload( const char* raw_file, uint32_t trnsfm, uint32_t height )
{
    return Load( raw_file, trnsfm, height );
}

bool RAWProcessor::Reload()
{
    if ( raw_file_name == DEF_MEMORY_LOADED )
        return false;

#ifndef WCHAR_SUPPORTED
    return Reload( raw_file_name.c_str() );
#else
    return Reload( _TCM2W(raw_file_name.c_str()) );
#endif /// of WCHAR_SUPPORTED
}

void RAWProcessor::Unload()
{
    if ( raw_loaded == true )
    {
        raw_loaded = false;
    }

    if ( pixel_arrays != NULL )
    {
        delete[] pixel_arrays;
        pixel_arrays = NULL;
    }
    pixel_arrays_paddedsz = 0;

    resetWindow();
}

bool RAWProcessor::ApplyTransform( uint32_t trnsfm )
{
    if ( ( raw_loaded == true ) && ( trnsfm > 0 ) )
    {
        bool retb = false;
        float* ptrd = pixel_arrays;

        if ( ( trnsfm & TRANSFORM_SWAP ) > 0 )
        {
            if ( ( trnsfm & TRANSFORM_PARAM_SWAP_H ) > 0 )
            {
                retb = rawimgtk::FlipHorizontal( ptrd, img_width, img_height );
            }

            if ( ( trnsfm & TRANSFORM_PARAM_SWAP_V ) > 0 )
            {
                retb = rawimgtk::FlipVertical( ptrd, img_width, img_height );
            }
        }

        if ( ( trnsfm & TRANSFORM_ROTATE ) > 0 )
        {
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_C90 ) > 0 )
            {
                retb = rawimgtk::Rotate90( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_C180 ) > 0 )
            {
                retb = rawimgtk::Rotate180( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_C270 ) > 0 )
            {
                retb = rawimgtk::Rotate270( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_RC90 ) > 0 )
            {
                retb = rawimgtk::Rotate270( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_RC180 ) > 0 )
            {
                retb = rawimgtk::Rotate180( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_RC270 ) > 0 )
            {
                retb = rawimgtk::Rotate90( ptrd, &img_width, &img_height );
            }
        }

        if ( retb == true )
        {
            current_transform = trnsfm;
        }

        return retb;
    }

    return false;
}

bool RAWProcessor::Invert()
{
    return InvertAuto();
}

bool RAWProcessor::InvertAuto()
{
    size_t dlen = pixel_arrays_paddedsz;

    if ( dlen == 0 )
        return false;

    #pragma omp parallel for
    for( size_t cnt=0; cnt<dlen; cnt++ )
    {
        pixel_arrays[ cnt ] = pixel_max_level - pixel_arrays[ cnt ];
    }

    return true;
}

void RAWProcessor::ChangeHeight( uint32_t h )
{
    if ( raw_loaded == false )
        return;

    if ( ( h > 0 ) && ( img_height != h ) )
    {
        img_height = h;
        analyse();
    }
}

bool RAWProcessor::Get8bitDownscaled( vector<uint8_t>& b_arrays, DownscaleType dntype, bool reversed )
{
    if ( raw_loaded == false )
        return false;

    size_t arrsz = pixel_arrays_paddedsz;

    if ( arrsz == 0 )
        return false;

    b_arrays.clear();
    b_arrays.reserve( arrsz );
    b_arrays.resize( arrsz );

    const float* ref_pixel_arrays = pixel_arrays;

    if ( dntype == DNSCALE_NORMAL )
    {
        float dscale_ratio = DEF_CALC_F_BMAX / float( pixel_max_level );
#ifdef DEBUG
        printf( "dscale_ratio = %f\n", dscale_ratio );
        fflush( stdout );
#endif /// of DEBUG
        #pragma omp parallel for
        for( size_t cnt=0; cnt<arrsz; cnt++ )
        {
            float fdspixel = float(ref_pixel_arrays[cnt]) * dscale_ratio;

            if ( fdspixel > DEF_CALC_F_BMAX )
                fdspixel = DEF_CALC_F_BMAX;
            else
            if ( fdspixel < 0.f )
                fdspixel = 0.f;

            uint8_t dspixel = (uint8_t)MIN(fdspixel,DEF_CALC_F_BMAX);

            if ( reversed == true )
            {
                dspixel = DEF_PIXEL8_MAX - dspixel;
            }

            b_arrays[cnt] = dspixel;
        }

    }
    else
    if ( dntype == DNSCALE_FULLDOWN )
    {
#ifdef DEBUG
        printf( "scaling ratio = %f\n", DEF_CALC_F_BMAX );
        fflush( stdout );
#endif /// of DEBUG

        #pragma omp parallel for
        for( uint32_t cnt=0; cnt<arrsz; cnt++ )
        {
            float fdspixel = ref_pixel_arrays[cnt] * DEF_CALC_F_BMAX;

            if ( fdspixel > DEF_CALC_F_BMAX )
                fdspixel = DEF_CALC_F_BMAX;
            else
            if ( fdspixel < 0.f )
                fdspixel = 0.f;

            uint8_t dspixel = (uint8_t)MIN(fdspixel,DEF_CALC_F_BMAX);

            if ( reversed == true )
            {
                dspixel = DEF_PIXEL8_MAX - dspixel;
            }

            b_arrays[cnt] = dspixel;
        }
    }

    return true;
}

bool RAWProcessor::Get16bitScaledImage( vector<uint16_t>& w_arrays, bool reversed )
{
    return f2dexport( 16, &w_arrays, reversed );
}

bool RAWProcessor::GetRawImage( std::vector<float>& f_arrays, bool reversed )
{
    if ( raw_loaded == false )
        return false;

    size_t arrsz = pixel_arrays_paddedsz;

    if ( arrsz == 0 )
        return false;

    f_arrays.clear();
    f_arrays.reserve( arrsz );
    f_arrays.resize( arrsz );

    float* ref_pixel_arrays = pixel_arrays;

    if ( reversed == true )
    {
        #pragma omp parallel for
        for ( size_t cnt=0; cnt<arrsz; cnt++ )
        {
            f_arrays[cnt] = ( pixel_max_level - pixel_arrays[cnt] );
        }
    }
    else
    {
        #pragma omp parallel for
        for ( size_t cnt=0; cnt<arrsz; cnt++ )
        {
            f_arrays[cnt] = pixel_arrays[cnt];
        }
    }

    return true;
}

bool RAWProcessor::GetAnalysisReport( WindowAnalysisReport& report, bool start_minlevel_zero )
{
    if ( ( window_max == 0 ) || ( pixel_arrays_paddedsz == 0 ) )
        return false;

    // # phase 01
    // get current time stamp.
    time_t curtime;
    report.timestamp = (uint32_t)time(&curtime);
    report.wide_min = window_min;
    report.wide_max = window_max;
    report.wide_max_amount = window_wide;

    if ( start_minlevel_zero == true )
    {
        report.wide_min = 0;
        report.wide_max_amount += window_min;
    }

    return true;
}

// floating point doesn't have trhesholding method, just make uint16_t type image.
bool RAWProcessor::GetWindowedImage( WindowAnalysisReport& report, std::vector<uint32_t>& d_arrays, bool reversed )
{
    return f2dthldexport( 32, &d_arrays, reversed, &report );
}

bool RAWProcessor::GetWindowedImage( WindowAnalysisReport &report,  vector<uint16_t>& w_arrays, bool reversed )
{
    return f2dthldexport( 16, &w_arrays, reversed, &report );
}

bool RAWProcessor::GetWindowedImage( WindowAnalysisReport &report, std::vector<uint8_t>& b_arrays, bool reversed )
{
    return f2dthldexport( 8, &b_arrays, reversed, &report );
}

bool RAWProcessor::GetPixel( uint32_t x, uint32_t y, float &px )
{
    if ( pixel_arrays_paddedsz == 0 )
        return false;

    if ( ( x > img_width ) || ( y > img_height ) )
        return false;

    size_t pixpos = ( y * img_height ) + x;

    px = pixel_arrays[ pixpos ];

    return true;
}

bool RAWProcessor::GetHistogram( std::vector<size_t> &d_histo, size_t& max_lvl )
{
    if ( pixel_arrays_paddedsz == 0 )
        return false;

    // histogram amounts always follow window_max; 0 to window_max
    d_histo.clear();
    d_histo.resize( window_max );

    size_t kinds_amount = 0;
    max_lvl = window_max;

    #pragma omp parallel for
    for( size_t cnt=0; cnt<pixel_arrays_paddedsz; cnt++ )
    {
        size_t pos32 = pixel_arrays[cnt] * window_max;
        if ( pos32 < window_max )
            d_histo[pos32]++;
    }

    return true;
}

bool RAWProcessor::SaveToFile( const char* path )
{
    if ( pixel_arrays_paddedsz == 0 )
        return false;

    if ( path != NULL )
    {
        if ( access( path, F_OK ) == 0 )
        {
            if ( unlink( path ) != 0 )
                return false;
        }

        FILE* fp = fopen( path, "wb" );
        if ( fp != NULL )
        {
            fwrite( pixel_arrays, 2, pixel_arrays_paddedsz, fp );
            fclose( fp );

            return true;
        }
    }

    return false;
}

#ifdef WCHAR_SUPPORTED
bool RAWProcessor::SaveToFile( const wchar_t* path )
{
    if ( pixel_arrays_paddedsz == 0 )
        return false;

    if ( path != NULL )
    {
        if ( _waccess( path, F_OK ) == 0 )
        {
            if ( _wunlink( path ) != 0 )
                return false;
        }

        FILE* fp = _wfopen( path, L"wb" );
        if ( fp != NULL )
        {
            fwrite( pixel_arrays, 2, pixel_arrays_paddedsz, fp );
            fclose( fp );

            return true;
        }
    }

    return false;
}
#endif /// of WCHAR_SUPPORTED

bool RAWProcessor::RotateFree( float degree )
{
    float* src = pixel_arrays;
    float* dst = NULL;

    uint32_t dst_w = img_width;
    uint32_t dst_h = img_height;

    float bgl = pixel_min_level;

    bool retb = false;

    if ( rawimgtk::RotateFree( src, &dst_w, &dst_h, &dst, degree, bgl ) == true )
    {
        float* dstcrop = NULL;

        if ( rawimgtk::CropCenter( dst, dst_w, dst_h, &dstcrop, img_width, img_height ) == true )
        {
            memcpy( src, dstcrop, pixel_arrays_paddedsz * sizeof( float ) );

            retb = true;

            delete[] dstcrop;
        }

        delete[] dst;
    }

    return retb;
}

RAWProcessor* RAWProcessor::RotateFree( float degree, float background )
{
    float* src = pixel_arrays;
    float* dst = NULL;

    uint32_t dst_w = img_width;
    uint32_t dst_h = img_height;

    float bgl = background;

    bool retb = false;

    if ( rawimgtk::RotateFree( src, &dst_w, &dst_h, &dst, degree, bgl ) == true )
    {
        RAWProcessor* newrp = new RAWProcessor();
        if ( newrp != NULL )
        {
            const char* ptr     = (const char*)dst;
            uint32_t ptrsz =  dst_w * dst_h * sizeof(float);

            retb = newrp->LoadFromMemory( (void*)ptr, ptrsz, TRANSFORM_NONE, dst_h );

            delete[] dst;

            if ( retb == true )
            {
                return newrp;
            }

            delete newrp;
        }
    }

    return NULL;
}

RAWProcessor* RAWProcessor::Rescale( uint32_t w, uint32_t h, RescaleType st )
{
    if ( ( pixel_arrays_paddedsz > 0 ) && ( w > 0 ) && ( h > 0 ) )
    {
        RAWGenericFilter* afilter = NULL;

        switch( st )
        {
            default:
            case RESCALE_NEAREST:
                afilter = new RAWBoxFilter();
                break;

            case RESCALE_BILINEAR:
                afilter = new RAWBilinearFilter();
                break;

            case RESCALE_BICUBIC:
                afilter = new RAWBicubicFilter();
                break;

            case RESCALE_BSPLINE:
                afilter = new RAWBSplineFilter();
                break;

            case RESCALE_LANZCOS3:
                afilter = new RAWLanczos3Filter();
                break;
        }

        if ( afilter != NULL )
        {
            RAWResizeEngine rawrse( afilter );

            const float* refsrc = pixel_arrays;
            float* dst = NULL;
            size_t retsz = rawrse.scale( refsrc, img_width, img_height, w, h, &dst );
            delete afilter;

            if ( retsz > 0 )
            {
                RAWProcessor* newme = new RAWProcessor();
                if ( newme != NULL )
                {
                    bool retb = newme->LoadFromMemory( (void*)dst, retsz,
                                                       TRANSFORM_NONE, h );
                    if ( retb == true )
                    {
                        return newme;
                    }

                    delete newme;
                }
            }
        }
    }

    return NULL;
}

RAWProcessor* RAWProcessor::Clone()
{
    if ( ( pixel_arrays_paddedsz > 0 ) && ( img_width > 0 ) && ( img_height > 0 ) )
    {
        RAWProcessor* newone = new RAWProcessor();
        if ( newone != NULL )
        {
            bool retb = newone->LoadFromMemory( (void*)data(),
                                                datasize() * sizeof( float ),
                                                TRANSFORM_NONE,
                                                img_height );
            if( retb == true )
            {
                return newone;
            }

            delete newone;
        }
    }

    return NULL;
}

void RAWProcessor::GetLinearPixels( uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2, vector<float>* pixels )
{
    if ( pixels == NULL )
    {
        return;
    }

    int _x0   = x1;
    int _x1   = x2;
    int _y0   = y1;
    int _y1   = y2;

    int inc1  = 0;
    int inc2  = 0;
    int cnt   = 0;
    int y_adj = 0;
    int dy    = 0;
    int dx    = 0;
    int x_adj = 0;

    if ( _x0 == _x1 )
    {
        if ( _y0 > _y1 )
        {
            swap( _y0, _y1 );
        }

        int cnt = _y1 - _y0 + 1;

        while( cnt-- )
        {
            addpixelarray( pixels, _x0, _y0 + cnt );
        }
    }
    else
    {
        if ( _y0 == _y1 )
        {
            if ( _x0 > _x1 )
            {
                swap( _x0, _x1 );
            }

            dx = _x1 - _x0 + 1;

            for( int cnt=0; cnt<dx; cnt++ )
            {
                addpixelarray( pixels, _x0 + cnt, _y0 );
            }
        }
        else
        {
            dy = _y1 - _y0;
            dx = _x1 - _x0;

            if ( abs( dy ) < abs( dx ) )
            {
                if ( _x0 > _x1 )
                {
                    swap( _x0, _x1 );
                    swap( _y0, _y1 );
                }

                dy = _y1 - _y0;
                dx = _x1 - _x0;

                if ( dy < 0 )
                {
                    dy    = -dy;
                    y_adj = -1;
                }
                else
                    y_adj = 1;

                inc1 = dy << 1;
                inc2 = ( dy - dx ) << 1;
                cnt  = ( dy << 1 ) - dx;

                dx++;
                int py = y_adj;

                while ( dx-- )
                {
                    addpixelarray( pixels, _x0 + dx, _y1 - py );

                    if ( cnt >= 0 )
                    {
                        cnt += inc2;
                        py  += y_adj;
                    }
                    else
                    {
                        cnt += inc1;
                    }
                }
            }
            else
            {
                if ( _y0 > _y1 )
                {
                    swap( _x0, _x1 );
                    swap( _y0, _y1 );
                }

                dy = _y1 - _y0;
                dx = _x1 - _x0;

                if ( dx < 0)
                {
                    dx    = -dx;
                    x_adj = -1;
                }
                else
                {
                    x_adj = 1;
                }

                inc1 = dx << 1;
                inc2 = ( dx - dy ) << 1;
                cnt  = ( dx << 1 ) - dy;

                dy++;
                int px = x_adj;

                while ( dy-- )
                {
                     addpixelarray( pixels, _x0 + px, _y1 - dy );

                    if ( cnt >= 0 )
                    {
                        cnt += inc2;
                        px  += x_adj;
                    }
                    else
                    {
                        cnt += inc1;
                    }
                }
            }
        }
    }
}

void RAWProcessor::GetRectPixels( uint32_t x, uint32_t y, uint32_t w, uint32_t h, vector<float>* pixels)
{
    if ( pixels == NULL )
    {
        return;
    }

    uint32_t rw = w;
    uint32_t rh = h;

    if ( ( rw + x ) > img_width )
    {
        rw = img_width - x;
    }

    if ( ( rh + y ) > img_height )
    {
        rh = img_height - y;
    }

    for( uint32_t cnty=y; cnty<rh; cnty++ )
    {
        for( uint32_t cntx=x;cntx<rw; cntx++ )
        {
            uint32_t mpos = img_width * cnty + cntx;

            if ( pixel_arrays_paddedsz > mpos )
            {
                pixels->push_back( pixel_arrays[ mpos ] );
            }
        }
    }
}

void RAWProcessor::GetPolygonPixels( vector<polygoncoord>* coords, vector<float>* pixels)
{
    if ( ( coords == NULL ) || ( pixels == NULL ) )
    {
        return;
    }

    reordercoords( coords );

    uint32_t ptsz = coords->size();

    if ( ptsz < 2 )
    {
        // need it starts triangle or more.
        return;
    }

    const uint32_t max_y = img_height;
    const uint32_t max_x = img_width;

    vector< double > node_x;

    for( uint32_t cur_y=0; cur_y<max_y; cur_y++ )
    {
        uint32_t    rcnt        = ptsz - 1;

        for( uint32_t cnt=0; cnt<ptsz; cnt++ )
        {
            if ( ( ( (double)coords->at(cnt).y   < (double)cur_y ) &&
                   ( (double)coords->at(rcnt).y >= (double)cur_y ) ) ||
                 ( ( (double)coords->at(rcnt).y  < (double)cur_y ) &&
                   ( (double)coords->at(cnt).y  >= (double)cur_y ) ) )
            {
                double newv = (double)coords->at(cnt).x +
                              ( (double)cur_y              - (double)coords->at(cnt).y ) /
                              ( (double)coords->at(rcnt).y - (double)coords->at(cnt).y ) *
                              ( (double)coords->at(rcnt).x - (double)coords->at(cnt).x );

                node_x.push_back( newv );
            }

            rcnt = cnt;
        }

        uint32_t node_x_sz = node_x.size();

        // sort nodes ..
        if ( node_x_sz > 1 )
        {
            sort( node_x.begin(),
                  node_x.begin() + node_x_sz,
                  RAWProcessor_sortcondition );
        }

        for( uint32_t dcnt=0; dcnt<node_x_sz; dcnt+=2 )
        {
            if ( node_x[dcnt] >= max_x )
                break;

            if ( node_x[dcnt+1] > 0 )
            {
                if ( node_x[dcnt] < 0 )
                {
                    node_x[dcnt] = 0;
                }

                if ( node_x[dcnt+1] > max_x )
                {
                    node_x[dcnt+1] = max_x;
                }

                for( int cur_x=node_x[dcnt]; cur_x<=node_x[dcnt+1]; cur_x++ )
                {
                    addpixelarray( pixels, cur_x, cur_y );
                }
            }
        }

        node_x.clear();
    }
}

void RAWProcessor::GetAnalysisFromPixels( vector<float>& pixels, SimpleAnalysisInfo& info )
{
    float local_p_min_l = _MAX_F_;
    float local_p_max_l = _MIN_F_;
    float local_p_med_l = 0.0f;
    size_t local_w_w = 0;

    info.variance = 0.0;
    info.deviation = 0.0;
    double sd = 0.0;

    // find min/max pixel level ...
    #pragma omp for reduction(+:local_p_max_l) reduction(-:local_p_min_l) private(sd)
    for( size_t cnt=0; cnt<pixels.size(); cnt++ )
    {
        float pixelf = pixels[cnt];

        if ( ( pixelf > 0.f ) && ( local_w_w == 0 ) )
        {
            local_w_w = cnt;
        }

        if ( pixelf > pixel_max_level )
        {
            local_p_max_l = pixels[cnt];
        }
        else
        if ( pixelf < pixel_min_level )
        {
            local_p_min_l = pixels[cnt];
        }

        double dpixel = (double)pixels[cnt];
        sd += pow( dpixel - info.average , 2 );
    }

    // build historgram amounts as compressed into pixel counts.
    const float compressf = (float)pixels.size();
    vector< size_t > bHist;
    bHist.resize( (size_t)compressf  );

    #pragma omp parallel for
    for( size_t cnt=0; cnt<pixels.size(); cnt++ )
    {
        size_t castpos = pixels[cnt] * compressf;
        if ( castpos < bHist.size() )
            bHist[castpos]++;
    }

    size_t local_w_min = 0;
    size_t local_w_max = 0;
    window_max = pixel_arrays_paddedsz - 1;

    // find zero amounts of bHist.
    for ( size_t cnt=0; bHist.size(); cnt++ )
    {
        if ( bHist[cnt] > 0 )
        {
            local_w_min = cnt;
            break;
        }
    }

    for( size_t cnt=bHist.size()-1; cnt!=0; cnt-- )
    {
        if ( bHist[cnt] > 0 )
        {
            local_w_max = cnt;
            break;
        }
    }

    bHist.clear();

    local_p_med_l = ( pixel_max_level + pixel_min_level ) / 2.f;
    local_w_w = local_w_max - local_w_min;

    // Make information report.
    info.minLevel = local_p_min_l;
    info.maxLevel = local_p_max_l;
    info.average  = ( local_p_max_l - local_p_min_l ) / 2;

    info.variance  = sd / (double)pixels.size();
    info.deviation = sqrt( info.variance );
}

bool RAWProcessor::ApplyFilter( FilterConfig* fconfig )
{
    if ( fconfig != NULL )
    {
        if ( ( fconfig->width > 0 ) && ( fconfig->height > 0 ) &&
             ( pixel_arrays_paddedsz > 0 ) )
        {
            float* copy_arrays = new float[ pixel_arrays_paddedsz ];

            if ( copy_arrays == NULL )
            {
                return false;
            }

            memset( copy_arrays, 0, pixel_arrays_paddedsz * sizeof( float ) );

            #pragma omp parallel for
            for( uint32_t cntx=0; cntx<img_width; cntx++ )
            {
                for( uint32_t cnty=0; cnty<img_height; cnty++ )
                {
                    double adjustedp = 0.0;

                    // -- applying matrix ---
                    for( uint32_t fcntx=0; fcntx<fconfig->width; fcntx++ )
                    {
                        for( uint32_t fcnty=0; fcnty<fconfig->height; fcnty++ )
                        {
                            uint32_t posX = ( cntx - fconfig->width / 2 + fcntx + img_width )
                                            % img_width;
                            uint32_t posY = ( cnty - fconfig->height / 2 + fcnty + img_height )
                                            % img_height;

                            uint32_t posM = posY * img_width + posX;

                            if ( posM < pixel_arrays_paddedsz )
                            {
                                adjustedp += (double)pixel_arrays[ posM ] *
                                             (double)fconfig->matrix[ fcnty * fconfig->width + fcntx ];
                            }
                        }
                    }
                    // -- applying matrix ---

                    float rpixel = MIN( MAX( fconfig->factor * adjustedp + fconfig->bias, 0) \
                                        , 65535 );

                    copy_arrays[ cnty * img_width + cntx ] = rpixel;
                }
            }

            memcpy( pixel_arrays, copy_arrays, pixel_arrays_paddedsz * sizeof( float ) );

            delete[] copy_arrays;

            return true;
        }
    }

    return false;
}

bool RAWProcessor::ApplyMedianFilter()
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        float* copy_arrays = new float[ pixel_arrays_paddedsz ];

        if ( copy_arrays == NULL )
        {
            return false;
        }

        memset( copy_arrays, 0, pixel_arrays_paddedsz * sizeof( float ) );

        #pragma omp parallel for
        for( uint32_t cntx=0; cntx<img_width; cntx++ )
        {
            for( uint32_t cnty=0; cnty<img_height; cnty++ )
            {
                float  medimatrix[9] = {0};
                size_t medimatrixsz  = 0;

                // -- applying matrix ---
                for( uint32_t fcntx=0; fcntx<3; fcntx++ )
                {
                    for( uint32_t fcnty=0; fcnty<3; fcnty++ )
                    {
                        uint32_t posX = ( cntx - 3 / 2 + fcntx + img_width )
                                        % img_width;
                        uint32_t posY = ( cnty - 3 / 2 + fcnty + img_height )
                                        % img_height;

                        uint32_t posM = posY * img_width + posX;

                        if ( posM < pixel_arrays_paddedsz )
                        {
                             //medimatrix.push_back ( pixel_arrays[ posM ] );
                             medimatrix[ medimatrixsz ] = pixel_arrays[ posM ];
                             medimatrixsz++;
                        }
                    }
                }
                // -- applying matrix ---

                // sort it !
                qsort( medimatrix, medimatrixsz, sizeof(float), RAWProcessor_uscompare );

                copy_arrays[ cnty * img_width + cntx ] = medimatrix[ medimatrixsz / 2 ];
            }
        }

        memcpy( pixel_arrays, copy_arrays, pixel_arrays_paddedsz * sizeof( float ) );

        delete[] copy_arrays;
        return true;
    }

    return false;
}

RAWProcessor* RAWProcessor::CloneWithFilter( FilterConfig* fconfig )
{
    if ( fconfig != NULL )
    {
        RAWProcessor* newRP = Clone();
        if ( newRP != NULL )
        {
            if ( newRP->ApplyFilter( fconfig ) == true )
            {
                return newRP;
            }

            delete newRP;
        }
    }

    return NULL;
}

RAWProcessor::FilterConfig* RAWProcessor::GetPresetFilter( uint32_t fnum )
{
    FilterConfig* retfp = NULL;

    switch( fnum )
    {
        case PRESET_FILTER_BLUR:
            retfp = RAWImageFilterKit::GetPresetFilter( RAWFILTER_PRESET_BLUR );
            break;

        case PRESET_FILTER_BLURMORE:
            retfp = RAWImageFilterKit::GetPresetFilter( RAWFILTER_PRESET_BLURMORE );
            break;

        case PRESET_FILTER_SHARPEN:
            retfp = RAWImageFilterKit::GetPresetFilter( RAWFILTER_PRESET_SHARPEN );
            break;

        case PRESET_FILTER_UNSHARPEN:
            retfp = RAWImageFilterKit::GetPresetFilter( RAWFILTER_PRESET_UNSHARPENM );
            break;
    }

    return retfp;
}

void RAWProcessor::DiscardFilter( FilterConfig* fp )
{
    RAWImageFilterKit::RemoveFilter( fp );
}

bool RAWProcessor::AdjustGamma( float gamma )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        return rawimgtk::AdjustGamma( pixel_arrays,
                                      pixel_arrays_paddedsz,
                                      gamma );
    }

    return false;
}

bool RAWProcessor::AdjustBrightness( float percent )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        return rawimgtk::AdjustBrightness( pixel_arrays,
                                           pixel_arrays_paddedsz,
                                           percent );
    }

    return false;
}

bool RAWProcessor::AdjustContrast( float percent )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        return rawimgtk::AdjustContrast( pixel_arrays,
                                         pixel_arrays_paddedsz,
                                         percent );
    }

    return false;
}

bool RAWProcessor::ApplyToneMapping( uint32_t ttype, float p1, float p2, float p3, float p4 )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        switch ( ttype )
        {
            case TONEMAP_TYPE_DRAGO:
                {
                    float gamma = p1;
                    float exposure = p2;

                    return rawimgtk::tmoDrago03( pixel_arrays,
                                                 pixel_arrays_paddedsz,
                                                 pixel_max_level,
                                                 pixel_max_level,
                                                 gamma,
                                                 exposure );
                }
                break;

            case TONEMAP_TYPE_REINHARD:
                {
                    float intensity = p1;
                    float contrast = p2;
                    float adaptation = p3;
                    float colorcorrection = p4;

                    return rawimgtk::tmoReinhard2005( pixel_arrays,
                                                      pixel_arrays_paddedsz,
                                                      pixel_max_level,
                                                      pixel_max_level,
                                                      intensity,
                                                      contrast,
                                                      adaptation,
                                                      colorcorrection );
                }
                break;

        }
    }

    return false;
}

bool RAWProcessor::ApplyCLAHE( WindowAnalysisReport &report, uint32_t applysz, uint32_t bins, float slope )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        uint32_t minv = report.wide_min;
        uint32_t maxv = report.wide_max;
        float* ptr = pixel_arrays;

        return rawimgtk::ApplyCLAHE( ptr, img_width, img_height, minv, maxv,
                                     applysz, applysz, 0, slope );
    }

    return false;
}

bool RAWProcessor::ApplyLowFrequency( uint32_t filtersz, uint32_t repeat )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        float* ptr = pixel_arrays;

        return rawimgfk::ApplyLowFreqFilter( ptr,
                                             img_width, img_height,
                                             filtersz, repeat );
    }

    return false;
}

bool RAWProcessor::ApplyEdgeEnhance( uint32_t fszh, uint32_t fszv, uint32_t edgesz, uint32_t margin )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        float* ptr    = pixel_arrays;
        size_t imgsz  = pixel_arrays_paddedsz;
        float* imgEH1 = new float[ imgsz ];
        if ( imgEH1 == NULL )
            return false;

        float* imgEH2 = new float[ imgsz ];
        if ( imgEH2 == NULL )
        {
            delete[] imgEH1;

            return false;
        }

        memcpy( imgEH1, ptr, imgsz * sizeof( float ) );
        memcpy( imgEH2, ptr, imgsz * sizeof( float ) );

        bool retb = false;

        if ( fszh < 2 )
            fszh = 2;

        retb = RAWImageFilterKit::ApplyEdgeLowFreqFilter( imgEH1,
                                                          img_width, img_height,
                                                          fszh );

        if ( retb == false )
        {
            delete[] imgEH1;
            delete[] imgEH2;

            return false;
        }

        if ( fszv < 2 )
            fszv = 2;

        retb = RAWImageFilterKit::ApplyEdgeLowFreqFilter( imgEH2,
                                                          img_width, img_height,
                                                          fszv );

        if ( retb == false )
        {
            delete[] imgEH1;
            delete[] imgEH2;

            return false;
        }

        if ( ( img_width < (margin * 2) ) || ( img_height < (margin * 2) ) )
        {
            margin = 0;
        }

        uint32_t cnth = 0;
        uint32_t cntw = 0;
        uint32_t mgnx = margin;
        uint32_t mgny = margin;
        uint32_t mgnw = img_width - (margin * 2);
        uint32_t mgnh = img_height - (margin * 2);

        float fedgev = (float)edgesz / 8.0f;

        #pragma omp parallel for private( cntw )
        for( cnth=mgny; cnth<mgnh; cnth++ )
        {
            for( cntw=mgnx; cntw<mgnw; cntw++ )
            {
                uint32_t pos = cnth * img_width + cntw;
                double pixelv = abs( (float)ptr[pos] + (float)imgEH1[pos] + (float)imgEH2[pos] )
                                / fedgev;
                if ( pixelv < 0.0 )
                {
                    pixelv = 0.0;
                }
                else
                if ( pixelv > DEF_CALC_F_WMAX )
                {
                    pixelv = DEF_CALC_F_WMAX;
                }

                ptr[ pos ] = (float) pixelv;
            }
        }

        delete[] imgEH1;
        delete[] imgEH2;

        return true;
    }

    return false;
}

bool RAWProcessor::ApplyAnisotropicFilter( uint32_t strength, uint32_t param )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        float* ptr = pixel_arrays;

        return RAWImageFilterKit::ApplyAnisotropicFilter( ptr,
                                                          img_width, img_height,
                                                          strength, param );
    }

    return false;

}

const size_t RAWProcessor::datasize()
{
    return pixel_arrays_paddedsz;
}

const float* RAWProcessor::data()
{
    if ( pixel_arrays_paddedsz == 0 )
        return NULL;

    return pixel_arrays;
}

void RAWProcessor::analyse()
{
    //if ( pixel_counts > 0 )
    if ( pixel_arrays_paddedsz > 0 )
    {
        //img_width = pixel_counts / img_height;
        img_width = pixel_arrays_paddedsz / img_height;
    }

    // floating point version don't calculating bits width .
    // always 32bit.
    pixel_bpp = 32;
}

void RAWProcessor::resetWindow()
{
    window_min = 0;
    window_max = pixel_arrays_paddedsz;
}

void RAWProcessor::calcWindow()
{
    pixel_min_level = _MAX_F_;
    pixel_max_level = _MIN_F_;
    pixel_med_level = 0.0f;

    resetWindow();

    // find min/max pixel level ...
    #pragma omp for reduction(+:pixel_max_level) reduction(-:pixel_min_level)
    for( size_t cnt=0; cnt<pixel_arrays_paddedsz; cnt++ )
    {
        float pixelf = pixel_arrays[cnt];

        if ( ( pixelf > 0.f ) && ( window_min == 0 ) )
        {
            window_min = cnt;
        }

        if ( pixelf > pixel_max_level )
        {
            pixel_max_level = pixel_arrays[cnt];
        }
        else
        if ( pixelf < pixel_min_level )
        {
            pixel_min_level = pixel_arrays[cnt];
        }
    }

    // build historgram amounts as compressed into pixel counts.
    const float compressf = (float)pixel_arrays_paddedsz;
    vector< size_t > bHist;
    bHist.resize( (size_t)compressf  );

    #pragma omp parallel for
    for( size_t cnt=0; cnt<pixel_arrays_paddedsz; cnt++ )
    {
        size_t castpos = pixel_arrays[cnt] * compressf;
        if ( castpos < bHist.size() )
            bHist[castpos]++;
    }

    window_min = 0;
    window_max = pixel_arrays_paddedsz - 1;

    // find zero amounts of bHist.
    for ( size_t cnt=0; bHist.size(); cnt++ )
    {
        if ( bHist[cnt] > 0 )
        {
            window_min = cnt;
            break;
        }
    }

    for( size_t cnt=bHist.size()-1; cnt!=0; cnt-- )
    {
        if ( bHist[cnt] > 0 )
        {
            window_max = cnt;
            break;
        }
    }

    bHist.clear();

    pixel_med_level = ( pixel_max_level + pixel_min_level ) / 2.f;
    window_wide = window_max - window_min;

#ifdef DEBUG
    printf( "calcWindow() result min,max,med = %.5f %.5f, %.5f | wide = %zu, min = %zu, max = %zu\n",
            pixel_min_level, pixel_max_level, pixel_med_level,
            window_wide, window_min, window_max );
#endif // DEBUG
}

void RAWProcessor::findWideness(float& minf, float& maxf )
{
    if ( pixel_arrays_paddedsz > 0 )
    {
        if ( window_min < pixel_arrays_paddedsz )
            minf = pixel_arrays[ window_min ];

        if ( window_max < pixel_arrays_paddedsz )
            maxf = pixel_arrays[ window_max ];
    }
}

void RAWProcessor::addpixelarray( std::vector<float>* outpixels, uint32_t x, uint32_t y )
{
    if ( outpixels == NULL )
        return;

    if ( ( x < img_width ) && ( y < img_height ) )
    {
        size_t mpos = img_width * y + x;
        if ( pixel_arrays_paddedsz > mpos )
        {
            outpixels->push_back( pixel_arrays[ mpos ] );
        }
    }
}

void RAWProcessor::reordercoords( std::vector<polygoncoord>* coords )
{
    if ( coords == NULL )
        return;

    size_t ptsz = coords->size();

    if ( ptsz > 2 )
    {
        size_t idxFirst = -1;
        size_t minX = img_width;
        size_t minY = img_height;

        for( size_t cnt=0; cnt<(ptsz - 1); cnt++ )
        {
            if ( coords->at(cnt).y < minY )
            {
                minX = coords->at(cnt).x;
                minY = coords->at(cnt).y;
                idxFirst = cnt;
            }
        }

        // X is next.
        for( size_t cnt=0; cnt<(ptsz - 1); cnt++ )
        {
            if ( ( coords->at(cnt).x < minX ) && ( coords->at(cnt).y < minY ) )
            {
                minX = coords->at(cnt).x;
                minY = coords->at(cnt).y;
                idxFirst = cnt;
                break;
            }
        }

        if ( idxFirst > 0 )
        {
            vector< polygoncoord > copydummy;
            copydummy.resize( ptsz );

            // copy all to dummy.
            for( size_t cpcnt=0; cpcnt<ptsz-1; cpcnt++ )
            {
                copydummy[cpcnt].x = coords->at(cpcnt).x;
                copydummy[cpcnt].y = coords->at(cpcnt).y;
            }

            size_t lastQ = 0;

            // copy back all except last coordination
            for( size_t cnt=idxFirst; cnt<(ptsz - 1); cnt++ )
            {
                coords->at(cnt-idxFirst).x = copydummy[cnt].x;
                coords->at(cnt-idxFirst).y = copydummy[cnt].y;

                lastQ = cnt - idxFirst;
            }

            // Last end point will be re-assign for new vector's first.

            for( size_t cnt=0; cnt<idxFirst; cnt++ )
            {
                coords->at( cnt + lastQ + 1 ).x = copydummy[ cnt ].x;
                coords->at( cnt + lastQ + 1 ).y = copydummy[ cnt ].y;
            }

            // make perfect ring vector. start pos = end pos.
            coords->at( ptsz - 1 ).x = coords->at( 0 ).x;
            coords->at( ptsz - 1 ).y = coords->at( 0 ).y;

            copydummy.clear();
        }
    }
}

bool RAWProcessor::f2dthldexport( uint8_t tp, void* pd, bool rvs, WindowAnalysisReport* report )
{
    if ( report == NULL )
        return false;

    if ( report->timestamp == 0 )
        return false;

    size_t array_sz = pixel_arrays_paddedsz;

    if ( array_sz == 0 )
        return false;

    vector< uint8_t >*  pvt8  = NULL;
    vector< uint16_t >* pvt16 = NULL;
    vector< uint32_t >* pvt32 = NULL;
    uint8_t* pdt_src = NULL;
    size_t   pdt_sz = 0;
    size_t   pdt_elem_sz = 0;

    switch( tp )
    {
        case 8: /// == uint8_t
            pvt8 = (vector< uint8_t >*)pd;
            pvt8->reserve( array_sz );
            pvt8->resize( array_sz, 0 );
            pdt_src = (uint8_t*)pvt8->data();
            pdt_sz  = pvt8->size();
            pdt_elem_sz = 1;
            break;

        case 10: /// == uint16_t
        case 12: /// == uint16_t
        case 16: /// == uint16_t
            pvt16 = (vector< uint16_t >*)pd;
            pvt16->reserve( array_sz );
            pvt16->resize( array_sz, 0 );
            pdt_src = (uint8_t*)pvt16->data();
            pdt_sz  = pvt16->size();
            pdt_elem_sz = 2;
            break;

        case 24: /// == uint32_t
        case 32: /// == uint32_t
            pvt32 = (vector< uint32_t >*)pd;
            pvt32->reserve( array_sz );
            pvt32->resize( array_sz, 0 );
            pdt_src = (uint8_t*)pvt32->data();
            pdt_sz  = pvt32->size();
            pdt_elem_sz = 4;
            break;

        default: /// unsupported.
            return false;
    }

    if ( pdt_sz == 0 )
        return false;

    float thld_min = report->wide_min;
    float thld_max = report->wide_max;
    float maxbf    = thld_max - thld_min;
    float normf    = 1.f / maxbf;
    float multf    = (float)pow(2, tp) - 1.f;

#ifdef DEBUG
    printf( "[debug]f2dthldexport()\n"
            "thld_min, thld_max = %f, %f\n"
            "maxbf = %f, normf = %f, multf = %f\n",
            thld_min, thld_max, maxbf, normf, multf );
    fflush( stdout );
#endif /// of DEBUG

    #pragma omp parallel for private
    for( size_t cnt=0; cnt<array_sz; cnt++ )
    {
        float apixel = MIN(MAX(pixel_arrays[cnt],thld_min),thld_max) - thld_min;

        if ( rvs == true )
        {
            //apixel = 1.0f - apixel;
            //maxbf  = 1.0f - maxbf;
        }

        switch( tp )
        {
            case 8: /// == uint8_t
                {
                    uint8_t* cast8 = &pdt_src[pdt_elem_sz * cnt];
                    *cast8 = (uint8_t)MIN((apixel * multf * normf),DEF_CALC_F_BMAX);
                }
                break;

            case 10: /// == uint16_t
            case 12: /// == uint16_t
            case 16: /// == uint16_t
                {
                    uint16_t* cast16 = (uint16_t*)&pdt_src[pdt_elem_sz * cnt];
                    *cast16 = (uint16_t)MIN((apixel * multf * normf),DEF_CALC_F_WMAX);
                }
                break;

            case 24: /// == uint32_t
            case 32: /// == uint32_t
                {
                    uint32_t* cast32 = (uint32_t*)&pdt_src[pdt_elem_sz * cnt];
                    *cast32 = (uint32_t)MIN((apixel * multf * normf ),DEF_CALC_F_DMAX);
                }
                break;
        }
    }

    return true;
}

bool RAWProcessor::f2dexport( uint8_t tp, void* pd, bool rvs )
{
    size_t array_sz = pixel_arrays_paddedsz;

    if ( array_sz == 0 )
        return false;

    vector< uint8_t >*  pvt8  = NULL;
    vector< uint16_t >* pvt16 = NULL;
    vector< uint32_t >* pvt32 = NULL;
    uint8_t* pdt_src = NULL;
    size_t   pdt_sz = 0;
    size_t   pdt_elem_sz = 0;

    switch( tp )
    {
        case 8: /// == uint8_t
            pvt8 = (vector< uint8_t >*)pd;
            pvt8->reserve( array_sz );
            pvt8->resize( array_sz );
            pdt_src = (uint8_t*)pvt8->data();
            pdt_sz  = pvt8->size();
            pdt_elem_sz = 1;
            break;

        case 10: /// == uint16_t
        case 12: /// == uint16_t
        case 16: /// == uint16_t
            pvt16 = (vector< uint16_t >*)pd;
            pvt16->reserve( array_sz );
            pvt16->resize( array_sz );
            pdt_src = (uint8_t*)pvt16->data();
            pdt_sz  = pvt16->size();
            pdt_elem_sz = 2;

            break;

        case 24: /// == uint32_t
        case 32: /// == uint32_t
            pvt32 = (vector< uint32_t >*)pd;
            pvt32->reserve( array_sz );
            pvt32->resize( array_sz );
            pdt_src = (uint8_t*)pvt32->data();
            pdt_sz  = pvt32->size();
            pdt_elem_sz = 4;
            break;

        default: /// unsupported.
            return false;
    }

    if ( pdt_sz == 0 )
        return false;

    float multf    = (float)pow(2, tp) - 1.f;
    float thld_min = 0.f;
    float thld_max = pow( 2, tp );
    float dividerf = 1.f / pixel_max_level;

    #pragma omp parallel for
    for( size_t cnt=0; cnt<array_sz; cnt++ )
    {
        /*
        float apixel = pixel_arrays[cnt] * dividerf * multf;

        // cut off threhold pixel value.
        if ( apixel > thld_max )
        {
            apixel = thld_max;
        }
        else
        if ( apixel < thld_min )
        {
            apixel = thld_min;
        }
        */

        float apixel = MIN(MAX(pixel_arrays[cnt] * dividerf * multf,thld_min),thld_max);

        if ( rvs == true )
        {
            apixel = ( pixel_max_level * multf ) - apixel;
        }

        switch( tp )
        {
            case 8: /// == uint8_t
                {
                    uint8_t* cast8 = &pdt_src[pdt_elem_sz * cnt];
                    *cast8 = (uint8_t)MIN(apixel,thld_max);
                }
                break;

            case 10: /// == uint16_t
            case 12: /// == uint16_t
            case 16: /// == uint16_t
                {
                    uint16_t* cast16 = (uint16_t*)&pdt_src[pdt_elem_sz * cnt];
                    *cast16 = (uint16_t)MIN(apixel,thld_max);
                }
                break;

            case 24: /// == uint32_t
            case 32: /// == uint32_t
                {
                    uint32_t* cast32 = (uint32_t*)&pdt_src[pdt_elem_sz * cnt];
                    *cast32 = (uint32_t)MIN(apixel,thld_max);
                }
                break;
        }
    }

    return true;
}

// all pixels makes to normalized as 0.0 to 1.0 range --
void RAWProcessor::normalize()
{
    size_t array_sz = pixel_arrays_paddedsz;

    if ( array_sz == 0 )
        return;

    float normal_maxf = 0.f;

    if ( pixel_max_level > 0.f )
    {
        normal_maxf = pixel_max_level / 1.f;

        #pragma omp parallel for
        for( size_t cnt=0; cnt<array_sz; cnt++ )
        {
            pixel_arrays[cnt] *= normal_maxf;
        }

        pixel_max_level = 1.f;
    }

#ifdef DEBUG
    printf( "noramlize() result normalize maxf = %f\n", normal_maxf );
#endif // DEBUG
}
