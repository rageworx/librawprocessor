#ifdef __APPLE__
    #include <sys/uio.h>
#endif // __APPLE__

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <ctime>
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
    #undef UNICODE
#endif // __APPLE__

#ifdef UNICODE
    #define _TSTRING        wstring
    #define _TCM2W( _x_ )   convertM2W( (const char*)_x_ )
    #define _TCW2M( _x_ )   convertW2M( (const wchar_t*)_x_ )
    #ifndef _T
        #define _T( _x_ )       L##_x_
    #endif
	#define DEF_MEMORY_LOADED	_T("//MEMORY_LOAD//")
#else
    #define _TSTRING        string
    #define _TCM2W( _x_ )   _x_
    #define _TCW2M( _x_ )   _x_
    #ifndef _T
        #define _T( _x_ )       _x_
    #endif
	#define DEF_MEMORY_LOADED	"//MEMORY_LOAD//"
#endif //// of UNICODE

#define DEF_RAW_I_HEIGHT    1024
#define DEF_PIXEL_WEIGHTS   65535

#define DEF_PIXEL16_MAX     65536
#define DEF_PIXEL8_MAX      256

#define DEF_CALC_F_WMAX     65535.0f
#define DEF_CALC_I_WMAX     65535
#define DEF_CALC_F_BMAX     255.0f
#define DEF_CALC_I_BMAX     255

#define DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY     0,9,53,150

////////////////////////////////////////////////////////////////////////////////

#define rawimgtk            RAWImageToolKit
#define rawimgfk            RAWImageFilterKit

////////////////////////////////////////////////////////////////////////////////

typedef struct
{
    unsigned min_v;
    unsigned max_v;
}minmaxpair;

////////////////////////////////////////////////////////////////////////////////

bool RAWProcessor_sortcondition( int i,int j )
{
    return ( i < j );
}

int RAWProcessor_uscompare (const void * a, const void * b)
{
    if( *(const unsigned short*)a < *(const unsigned short*)b )
        return -1;

    return *(const unsigned short*)a > *(const unsigned short*)b;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

RAWProcessor::RAWProcessor()
 : raw_loaded(false),
   pixel_swaping(false),
   pixel_arrays_realsz(0),
   pixel_weights(NULL),
   pixel_weights_max(0),
   pixel_bpp(10),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(0),
   img_width(0),
   userscaler(NULL)
{
    pixel_weights = new unsigned int[ DEF_PIXEL_WEIGHTS + 1 ];
    resetWeights();
}

RAWProcessor::RAWProcessor( const char* raw_file, unsigned int height )
 : raw_loaded(false),
   pixel_arrays_realsz(0),
   pixel_weights( NULL ),
   pixel_weights_max( 0 ),
   pixel_bpp(10),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(height),
   img_width(0)
{
    pixel_weights = new unsigned int[ DEF_PIXEL_WEIGHTS + 1];
    resetWeights();

    if ( raw_file != NULL )
    {
        Load( raw_file, TRANSFORM_NONE, height );
    }
}

#ifdef WCHAR_SUPPORTED
RAWProcessor::RAWProcessor( const wchar_t* raw_file, unsigned int height )
 : raw_loaded(false),
   pixel_arrays_realsz(0),
   pixel_weights( NULL ),
   pixel_weights_max( 0 ),
   pixel_bpp(10),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(height),
   img_width(0)
{
    pixel_weights = new unsigned int[ DEF_PIXEL_WEIGHTS + 1];
    resetWeights();

    if ( raw_file != NULL )
    {
        Load( raw_file, TRANSFORM_NONE, height );
    }
}
#endif /// of WCHAR_SUPPORTED

RAWProcessor::~RAWProcessor()
{
    Unload();

    if  ( pixel_weights != NULL )
    {
        delete[] pixel_weights;
        pixel_weights = NULL;
    }
}

void RAWProcessor::Version( char** retverstr )
{
    if ( retverstr == NULL )
    {
        return;
    }

    int retia[4] = { DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY };

    char retvstr[32] = {0};
    snprintf( retvstr, 32,
              "%d.%d.%d.%d", retia[0], retia[1], retia[2], retia[3] );

    int retvstrsz = strlen( retvstr );

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

    int retia[4] = { DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY };
    for( unsigned cnt=0; cnt<4; cnt++ )
    {
        retverints[ cnt ] = retia[ cnt ];
    }
}

#ifdef WCHAR_SUPPORTED
bool RAWProcessor::Load( const wchar_t* raw_file, unsigned int trnsfm, unsigned height )
{
    if ( height <= 0 )
        return false;

    string fname = convertW2M( raw_file );

    return Load( fname.c_str(), trnsfm, height );
}
#endif // of WCHAR_SUPPORTED

bool RAWProcessor::Load( const char* raw_file, unsigned int trnsfm, unsigned height )
{
    if ( height <= 0 )
        return false;

    fstream rfstrm;
    string  fname = raw_file;

    rfstrm.open( fname.c_str(), fstream::app | fstream::binary | fstream::in );

    if ( rfstrm.is_open() == true )
    {
        rfstrm.seekg( 0, ios::end );
        unsigned int fsize = rfstrm.tellg();
        rfstrm.seekg( 0, ios::beg );
        rfstrm.clear();

        if ( fsize == 0 )
        {
            rfstrm.close();
            raw_loaded = false;

            return false;
        }

        pixel_min_level = 0;
        pixel_med_level = 0;
        pixel_max_level = 0;
        img_height = height;

        pixel_arrays_srcsz   = fsize / sizeof( unsigned short );
        unsigned blancsz = 0;

        if (  pixel_arrays_srcsz > 0 )
        {
            // set temporary width ...
            img_width = pixel_arrays_realsz / img_height;

            // Check omit pixels ...
            if ( pixel_arrays_srcsz < ( img_width * img_height ) )
            {
                blancsz = ( img_width * img_height ) - pixel_arrays_srcsz;
            }
        }

        pixel_arrays_realsz = pixel_arrays_srcsz + blancsz;
        pixel_arrays.clear();
        pixel_arrays.reserve( pixel_arrays_realsz );
        pixel_arrays.resize( pixel_arrays_realsz );

        resetWeights();

        // To save memory, it doesn't using direct load to memory.
        for( unsigned cnt=0; cnt<fsize / sizeof(unsigned short); cnt++)
        {
            char           chardata[2] = {0};
            unsigned short worddata = 0;

            rfstrm.read( chardata, sizeof(unsigned short) );
            if ( pixel_swaping == false )
            {
                worddata = ( chardata[1] << 8 ) + ( chardata[0] & 0x00FF );
            }
            else
            {
                worddata = ( chardata[0] << 8 ) + ( chardata[1] & 0x00FF );
            }
            //pixel_arrays.push_back( worddata );
            pixel_arrays[ cnt ] = worddata;
            pixel_weights[ worddata ] += 1;

            if ( pixel_weights_max < worddata )
            {
                pixel_weights_max = worddata;
            }
        }

        pixel_med_level = ( pixel_max_level + pixel_min_level ) / 2;

        rfstrm.clear();
        rfstrm.close();

        calcWeights();
        analyse();

        raw_loaded = true;

        ApplyTransform( trnsfm );

        return raw_loaded;
    }

    return false;
}

bool RAWProcessor::LoadFromMemory( const char* buffer, unsigned long bufferlen, unsigned int trnsfm, unsigned height )
{
    if ( height <= 0 )
        return false;

    if ( ( buffer != NULL ) && ( bufferlen > 0 ) )
    {
        pixel_min_level = 0;
        pixel_med_level = 0;
        pixel_max_level = 0;
        img_height = height;

        pixel_arrays_srcsz = bufferlen / sizeof( unsigned short );
        unsigned blancsz   = 0;

        if (  pixel_arrays_srcsz > 0 )
        {
            // set temporary width ...
            img_width = pixel_arrays_srcsz / img_height;

            // Check omit pixels ...
            if ( pixel_arrays_srcsz < ( img_width * img_height ) )
            {
                blancsz = ( img_width * img_height ) - pixel_arrays_srcsz;
            }
        }

        pixel_arrays_realsz = pixel_arrays_srcsz + blancsz;

        pixel_arrays.clear();
        pixel_arrays.reserve( pixel_arrays_realsz );
        pixel_arrays.resize( pixel_arrays_realsz );

        resetWeights();

        // To save memory, it doesn't using direct load to memory.
        for( unsigned cnt=0; cnt<bufferlen / sizeof(unsigned short); cnt++)
        {
            char           chardata[2] = {0};
            unsigned short worddata = 0;

            chardata[0] = buffer[cnt*2];
            chardata[1] = buffer[cnt*2+1];
            worddata = ( chardata[1] << 8 ) + ( chardata[0] & 0x00FF );
            //pixel_arrays.push_back( worddata );
            pixel_arrays[ cnt ] = worddata;
            pixel_weights[ worddata ] += 1;

            if ( pixel_weights_max < worddata )
            {
                pixel_weights_max = worddata;
            }
        }

        calcWeights();
        analyse();

        raw_loaded = true;

        ApplyTransform( trnsfm );

        raw_file_name = DEF_MEMORY_LOADED;

        return true;
    }

    return false;
}

#ifdef WCHAR_SUPPORTED
bool RAWProcessor::Reload( const wchar_t* raw_file, unsigned int trnsfm, unsigned height )
{
    return Load( raw_file, trnsfm, height );
}
#endif /// of WCHAR_SUPPORTED

bool RAWProcessor::Reload( const char* raw_file, unsigned int trnsfm, unsigned height )
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

    pixel_arrays.clear();
    pixel_arrays.reserve( 0 );
    pixel_arrays.resize( 0 );

    pixel_arrays_realsz = 0;

    resetWeights();
}

bool RAWProcessor::ApplyTransform( unsigned int trnsfm )
{
    if ( ( raw_loaded == true ) && ( trnsfm > 0 ) )
    {
        bool retb = false;
        unsigned short* ptrd = pixel_arrays.data();

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

void RAWProcessor::SetUserScale( RAWUserScaleIF* ptr )
{
    userscaler = ptr;
}

bool RAWProcessor::Invert( unsigned char maxbits )
{
    if ( maxbits < 8 )
        return false;

    unsigned short maxval = pow( 2, maxbits );

    unsigned dlen = pixel_arrays.size();

    if ( dlen == 0 )
        return false;

    for( unsigned cnt=0; cnt<dlen; cnt++ )
    {
        pixel_arrays[ cnt ] = maxval - pixel_arrays[ cnt ];
    }

    return true;
}

bool RAWProcessor::InvertAuto()
{
    unsigned dlen = pixel_arrays.size();

    if ( dlen == 0 )
        return false;

    for( unsigned cnt=0; cnt<dlen; cnt++ )
    {
        pixel_arrays[ cnt ] = pixel_max_level - pixel_arrays[ cnt ];
    }

    return true;
}

void RAWProcessor::ChangeHeight( unsigned h )
{
    if ( raw_loaded == false )
        return;

    if ( ( h > 0 ) && ( img_height != h ) )
    {
        img_height = h;
        analyse();
    }
}

bool RAWProcessor::Get8bitDownscaled( vector<unsigned char>* byte_arrays, DownscaleType dntype, bool reversed )
{
    if ( raw_loaded == false )
        return false;

    unsigned arrsz = pixel_arrays.size();

    if ( arrsz == 0 )
        return false;

    byte_arrays->clear();
    byte_arrays->reserve( arrsz );
    byte_arrays->resize( arrsz );

    unsigned short* ref_pixel_arrays = pixel_arrays.data();

    if ( dntype == DNSCALE_NORMAL )
    {
        float dscale_ratio = DEF_CALC_F_BMAX / float( pixel_max_level );

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<arrsz; cnt++ )
        {
            float fdspixel = float(ref_pixel_arrays[cnt]) * dscale_ratio;
            if ( fdspixel > 255.f )
                fdspixel = 255.f;
            unsigned char dspixel = (unsigned char)fdspixel;

            if ( reversed == true )
            {
                dspixel = DEF_PIXEL8_MAX - dspixel - 1;
            }

            byte_arrays->at( cnt ) = dspixel;
        }

    }
    else
    if ( dntype == DNSCALE_FULLDOWN )
    {
        float uscale_ratio = DEF_CALC_F_WMAX / float( pixel_max_level );
        float dscale_ratio = DEF_CALC_F_BMAX / DEF_CALC_F_WMAX;

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<arrsz; cnt++ )
        {
            float fuspixel = float(ref_pixel_arrays[cnt]) * uscale_ratio;
            float fdspixel = fuspixel * dscale_ratio;
            if ( fdspixel > 255.f )
                fdspixel = 255.f;
            unsigned char dspixel = (unsigned char)fdspixel;

            if ( reversed == true )
            {
                dspixel = DEF_PIXEL8_MAX - dspixel - 1;
            }

            byte_arrays->at( cnt ) = dspixel;
        }
    }
    else
    if ( userscaler != NULL )
    {
        return userscaler->processUserScale( pixel_arrays, byte_arrays );
    }

    return true;
}

bool RAWProcessor::Get16bitRawImage( std::vector<unsigned short>* word_arrays, bool reversed )
{
    if ( raw_loaded == false )
        return false;

    unsigned arrsz = pixel_arrays.size();

    if ( arrsz == 0 )
        return false;

    word_arrays->clear();
    word_arrays->reserve( arrsz );
    word_arrays->resize( arrsz );

    unsigned short* ref_pixel_arrays = pixel_arrays.data();

    if ( reversed == true )
    {
        for ( unsigned cnt=0; cnt<arrsz; cnt++ )
        {
            word_arrays->at( cnt ) = DEF_PIXEL_WEIGHTS - ref_pixel_arrays[cnt];
        }
    }
    else
    {
        std::copy( pixel_arrays.begin(), pixel_arrays.end(), word_arrays->begin() );
    }

    return true;
}

bool RAWProcessor::GetWeights( std::vector<unsigned int>* weight_arrays )
{
    if ( raw_loaded == false )
        return false;

    int arrsz = pixel_weights_max;

    if ( arrsz == 0 )
        return false;

    weight_arrays->clear();
    weight_arrays->reserve( arrsz );
    weight_arrays->resize( arrsz );

    //weight_arrays.assign( pixel_weights.begin(), pixel_weights.end() );
    //std::copy( pixel_weights.begin(), pixel_weights.end(), weight_arrays.begin() );
    for( int cnt=0; cnt<arrsz; cnt++ )
    {
        weight_arrays->at( cnt ) = pixel_weights[cnt];
    }

    return true;
}

bool RAWProcessor::GetAnalysisReport( WeightAnalysisReport &report, bool start_minlevel_zero )
{
    if ( pixel_weights_max == 0 )
        return false;

    // # phase 01
    // get current time stamp.
    time_t curtime;
    report.timestamp = (unsigned long)time(&curtime);

    // # phase 02
    // find highest value in pixels ... ( 50% of maximum level )
    unsigned identify_min_level = (unsigned)( float(pixel_max_level) * 0.2f );
    if ( start_minlevel_zero == true )
    {
        identify_min_level = 0;
    }

    unsigned index_center_thld  = 0;

    for( unsigned cnt=identify_min_level; cnt<pixel_max_level; cnt++ )
    {
        if ( pixel_weights[ cnt ] > index_center_thld )
        {
            //report.base_threshold_index = pixel_weights[ cnt ];
            //index_center_thld = pixel_weights[ cnt ];
            index_center_thld = cnt;
        }
    }

    // # phase 03
    // find change pixel count fall into min level.
    unsigned short min_weight_wide = 0;
    unsigned short max_weight_wide = 0;

    vector< minmaxpair > mmpairs;

    unsigned avr_l = 0;
    const unsigned min_l = 100;   /// Minimal amount of pixel counts.
    unsigned min_q = 0;
    bool raised = false;

    // Makes min-max pairs.
    for( unsigned cnt=0; cnt<pixel_max_level; cnt++ )
    {
        avr_l += pixel_weights[ cnt ];
        avr_l /= 2;

        if ( ( avr_l > min_l ) && ( raised == false ) )
        {
            min_q = cnt;
            raised = true;
        }
        else
        if ( ( avr_l < min_l ) && ( raised == true ) )
        {
            minmaxpair ppair = {0,0};
            ppair.min_v = min_q;
            ppair.max_v = cnt;
            mmpairs.push_back( ppair );
            min_q = cnt;
            raised = false;
        }
    }

    // Find longest pairs.
    int longest_idx = -1;
    int next_idx = -1;
    unsigned test_long = 0;
    unsigned scnd_long = 0;

    for( unsigned cnt=0; cnt<mmpairs.size(); cnt++ )
    {
        unsigned tmp_long = mmpairs[ cnt ].max_v - mmpairs[ cnt ].min_v;

        if ( tmp_long > test_long )
        {
            scnd_long = test_long;
            test_long = tmp_long;
            longest_idx = cnt;
        }
        else
        if ( tmp_long > scnd_long )
        {
            scnd_long = tmp_long;
            next_idx = cnt;
        }
    }

    bool retb = false;

    if ( longest_idx >= 0 )
    {
        min_weight_wide = (float)mmpairs[ longest_idx ].min_v;

        if ( ( next_idx >= 0 ) && ( next_idx > longest_idx ) )
        {
            max_weight_wide = (float)mmpairs[ next_idx ].max_v;
        }
        else
        {
            max_weight_wide = (float)mmpairs[ longest_idx ].max_v;
        }

        retb = true;
    }
    else
    {
        min_weight_wide = 0;
        max_weight_wide = DEF_PIXEL_WEIGHTS;
    }

    mmpairs.clear();

    // # phase 04
    // get wide count.
    report.threshold_wide_min = min_weight_wide;
    report.threshold_wide_max = max_weight_wide;

    return retb;
}

bool RAWProcessor::Get16bitThresholdedImage( WeightAnalysisReport &report,  vector<unsigned short>* word_arrays, bool reversed )
{
    if ( report.timestamp == 0 )
        return false;

    float thld_min = (float)report.threshold_wide_min;
    float thld_max = (float)report.threshold_wide_max;
    float maxbf    = thld_max - thld_min;
    float normf    = maxbf / maxbf;

    word_arrays->clear();

    int array_max = pixel_arrays.size();

    word_arrays->reserve( array_max );
    word_arrays->resize( array_max );

    #pragma omp parallel for
    for( int cnt=0; cnt<array_max; cnt++ )
    {
        //unsigned short apixel = pixel_arrays[cnt] - thld_min - 1;
        unsigned short apixel = pixel_arrays[cnt];
        float          fpixel = 0.0f;

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

        apixel -= thld_min;

        // Rescaling pixel  !
        fpixel = (float)(apixel) * normf;
        apixel = (unsigned short)( fpixel );

        if ( reversed == true )
        {
            apixel = DEF_PIXEL16_MAX - apixel - 1;
        }

        word_arrays->at(cnt) = apixel;
    }

    return true;
}

bool RAWProcessor::Get8bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned char>* byte_arrays, bool reversed )
{
    if ( report.timestamp == 0 )
        return false;

    float thld_min = report.threshold_wide_min;
    float thld_max = (float)report.threshold_wide_max;
    float normfbpp = thld_max - thld_min;
    float normf    = DEF_CALC_F_BMAX / normfbpp;

    byte_arrays->clear();

    int array_max = pixel_arrays.size();

    byte_arrays->reserve( array_max );
    byte_arrays->resize( array_max );

    #pragma omp parallel for
    for( int cnt=0; cnt<array_max; cnt++ )
    {
        unsigned short apixel = pixel_arrays[cnt];
        float          fpixel = 0.0f;
        unsigned char  bpixel = 0;

        // cut off threhold pixel value.
        if ( apixel > thld_max )
        {
            //apixel = DEF_CALC_F_WMAX;
            apixel = thld_max;
        }
        else
        if ( apixel < thld_min )
        {
            apixel = thld_min;
        }

        apixel -= thld_min;

        // Rescaling pixel  !
        fpixel = (float)(apixel) * normf;

        // over flow check.
        if ( fpixel >= DEF_CALC_F_BMAX )
        {
            // if it overflows , set to 255.0f
            bpixel = DEF_PIXEL8_MAX - 1;
        }
        else
        {
            bpixel = (unsigned char)( fpixel );
        }

        if ( reversed == true )
        {
            bpixel = DEF_PIXEL8_MAX - bpixel - 1;
        }

        //byte_arrays.push_back( bpixel );
        byte_arrays->at( cnt ) = bpixel;
    }

    return true;
}

bool RAWProcessor::Get16bitPixel( unsigned x, unsigned y, unsigned short &px )
{
    if ( pixel_arrays.size() == 0 )
        return false;

    if ( ( x > img_width ) || ( y > img_height ) )
        return false;

    int pixpos = ( y * img_height ) + x;

    px = pixel_arrays[ pixpos ];

    return true;
}

bool RAWProcessor::SaveToFile( const char* path )
{
    if ( pixel_arrays.size() == 0 )
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
            fwrite( pixel_arrays.data(), 2, pixel_arrays.size(), fp );
            fclose( fp );

            return true;
        }
    }

    return false;
}

#ifdef WCHAR_SUPPORTED
bool RAWProcessor::SaveToFile( const wchar_t* path )
{
    if ( pixel_arrays.size() == 0 )
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
            fwrite( pixel_arrays.data(), 2, pixel_arrays.size(), fp );
            fclose( fp );

            return true;
        }
    }

    return false;
}
#endif /// of WCHAR_SUPPORTED

bool RAWProcessor::RotateFree( float degree )
{
    unsigned short* src = (unsigned short*)pixel_arrays.data();
    unsigned short* dst = NULL;

    unsigned dst_w = img_width;
    unsigned dst_h = img_height;

    unsigned short bgl = pixel_min_level;

    bool retb = false;

    if ( rawimgtk::RotateFree( src, &dst_w, &dst_h, &dst, degree, bgl ) == true )
    {
        unsigned short* dstcrop = NULL;

        if ( rawimgtk::CropCenter( dst, dst_w, dst_h, &dstcrop, img_width, img_height ) == true )
        {
            memcpy( src, dstcrop, pixel_arrays_srcsz * sizeof( unsigned short ) );

            retb = true;

            delete[] dstcrop;
        }

        delete[] dst;
    }

    return retb;
}

RAWProcessor* RAWProcessor::RotateFree( float degree, unsigned int background )
{
    unsigned short* src = (unsigned short*)pixel_arrays.data();
    unsigned short* dst = NULL;

    unsigned dst_w = img_width;
    unsigned dst_h = img_height;

    unsigned short bgl = background;

    bool retb = false;

    if ( rawimgtk::RotateFree( src, &dst_w, &dst_h, &dst, degree, bgl ) == true )
    {
        RAWProcessor* newrp = new RAWProcessor();
        if ( newrp != NULL )
        {
            const char* ptr     = (const char*)dst;
            unsigned long ptrsz =  dst_w * dst_h * sizeof(unsigned short);

            retb = newrp->LoadFromMemory( ptr, ptrsz, TRANSFORM_NONE, dst_h );

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

RAWProcessor* RAWProcessor::Rescale( unsigned w, unsigned h, RescaleType st )
{
    if ( ( pixel_arrays.size() > 0 ) && ( w > 0 ) && ( h > 0 ) )
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

            const unsigned short* refsrc = pixel_arrays.data();
            unsigned short* dst = NULL;
            unsigned long retsz = rawrse.scale( refsrc, img_width, img_height, w, h, &dst );
            delete afilter;

            if ( retsz > 0 )
            {
                RAWProcessor* newme = new RAWProcessor();
                if ( newme != NULL )
                {
                    bool retb = newme->LoadFromMemory( (const char*)dst, retsz,
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
    if ( ( pixel_arrays.size() > 0 ) && ( img_width > 0 ) && ( img_height > 0 ) )
    {
        RAWProcessor* newone = new RAWProcessor();
        if ( newone != NULL )
        {
            bool retb = newone->LoadFromMemory( (const char*)data(),
                                                datasize() * sizeof( unsigned short ),
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

void RAWProcessor::GetLinearPixels( unsigned x1, unsigned y1, unsigned x2, unsigned y2, vector<unsigned short>* pixels )
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

void RAWProcessor::GetRectPixels( unsigned x, unsigned y, unsigned w, unsigned h, vector<unsigned short>* pixels)
{
    if ( pixels == NULL )
    {
        return;
    }

    unsigned rw = w;
    unsigned rh = h;

    if ( ( rw + x ) > img_width )
    {
        rw = img_width - x;
    }

    if ( ( rh + y ) > img_height )
    {
        rh = img_height - y;
    }

    for( unsigned cnty=y; cnty<rh; cnty++ )
    {
        for( unsigned cntx=x;cntx<rw; cntx++ )
        {
            unsigned mpos = img_width * cnty + cntx;

            if ( pixel_arrays_realsz > mpos )
            {
                pixels->push_back( pixel_arrays[ mpos ] );
            }
        }
    }
}

void RAWProcessor::GetPolygonPixels( vector<polygoncoord>* coords, vector<unsigned short>* pixels)
{
    if ( ( coords == NULL ) || ( pixels == NULL ) )
    {
        return;
    }

    reordercoords( coords );

    unsigned ptsz = coords->size();

    if ( ptsz < 2 )
    {
        // need it starts triangle or more.
        return;
    }

    const unsigned max_y = img_height;
    const unsigned max_x = img_width;

    vector< double > node_x;

    for( unsigned cur_y=0; cur_y<max_y; cur_y++ )
    {
        unsigned    rcnt        = ptsz - 1;

        for( unsigned cnt=0; cnt<ptsz; cnt++ )
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

        unsigned node_x_sz = node_x.size();

        // sort nodes ..
        if ( node_x_sz > 1 )
        {
            sort( node_x.begin(),
                  node_x.begin() + node_x_sz,
                  RAWProcessor_sortcondition );
        }

        for( unsigned dcnt=0; dcnt<node_x_sz; dcnt+=2 )
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

void RAWProcessor::GetAnalysisFromPixels( std::vector<unsigned short>* pixels, std::vector<unsigned int>* weights, SimpleAnalysisInfo* info )
{
    if ( ( pixels == NULL ) || ( weights == NULL ) || ( info == NULL ) )
    {
        return;
    }

    // Make weights ...
    weights->resize( DEF_PIXEL_WEIGHTS + 1 );

    for( unsigned cnt=0; cnt<(DEF_PIXEL_WEIGHTS + 1); cnt++ )
    {
        weights->at( cnt ) = 0;
    }

    unsigned       pxsz      = pixels->size();
    unsigned short max_level = 0;
    unsigned short min_level = DEF_PIXEL_WEIGHTS;
    double         summ      = 0.0;

    for( unsigned cnt=0; cnt<pxsz; cnt++ )
    {
        unsigned short apixel = pixels->at( cnt );

        weights->at( apixel ) ++;

        if ( apixel > max_level )
        {
            max_level = apixel;
        }
        else
        if ( apixel < min_level )
        {
            min_level = apixel;
        }

        summ += (double)apixel;
    }

    // Make information report.
    info->minLevel = min_level;
    info->maxLevel = max_level;
    info->average  = summ / pxsz;


    info->variance = 0.0;
    info->deviation = 0.0;

    double sd = 0.0;
    // getting variance & deviation ...
    for( unsigned cnt=0; cnt<pxsz; cnt++ )
    {
        //info->deviation
        double dpixel = (double)pixels->at( cnt );
        sd += pow( dpixel - info->average , 2 );
    }

    info->variance  = sd / (double)pxsz;
    info->deviation = sqrt( info->variance );
}

bool RAWProcessor::ApplyFilter( FilterConfig* fconfig )
{
    if ( fconfig != NULL )
    {
        if ( ( fconfig->width > 0 ) && ( fconfig->height > 0 ) &&
             ( pixel_arrays_realsz > 0 ) )
        {
            unsigned short* copy_arrays = new unsigned short[ pixel_arrays_realsz ];

            if ( copy_arrays == NULL )
            {
                return false;
            }

            memset( copy_arrays, 0, pixel_arrays_realsz * sizeof( unsigned short ) );

            #pragma omp parallel for
            for( unsigned cntx=0; cntx<img_width; cntx++ )
            {
                for( unsigned cnty=0; cnty<img_height; cnty++ )
                {
                    double adjustedp = 0.0;

                    // -- applying matrix ---
                    for( unsigned fcntx=0; fcntx<fconfig->width; fcntx++ )
                    {
                        for( unsigned fcnty=0; fcnty<fconfig->height; fcnty++ )
                        {
                            unsigned posX = ( cntx - fconfig->width / 2 + fcntx + img_width )
                                            % img_width;
                            unsigned posY = ( cnty - fconfig->height / 2 + fcnty + img_height )
                                            % img_height;

                            unsigned posM = posY * img_width + posX;

                            if ( posM < pixel_arrays_realsz )
                            {
                                adjustedp += (double)pixel_arrays[ posM ] *
                                             (double)fconfig->matrix[ fcnty * fconfig->width + fcntx ];
                            }
                        }
                    }
                    // -- applying matrix ---

                    unsigned short rpixel = MIN( \
                                                MAX( fconfig->factor * adjustedp + fconfig->bias \
                                                    , 0) \
                                            , 65535 );

                    copy_arrays[ cnty * img_width + cntx ] = rpixel;
                }
            }

            memcpy( pixel_arrays.data(), copy_arrays, pixel_arrays_realsz * sizeof( unsigned short ) );

            delete[] copy_arrays;

            return true;
        }
    }

    return false;
}

bool RAWProcessor::ApplyMedianFilter()
{
    if ( pixel_arrays_realsz > 0 )
    {
        unsigned short* copy_arrays = new unsigned short[ pixel_arrays_realsz ];

        if ( copy_arrays == NULL )
        {
            return false;
        }

        memset( copy_arrays, 0, pixel_arrays_realsz * sizeof( unsigned short ) );

        #pragma omp parallel for
        for( unsigned cntx=0; cntx<img_width; cntx++ )
        {
            for( unsigned cnty=0; cnty<img_height; cnty++ )
            {
                unsigned short medimatrix[9] = {0};
                unsigned char  medimatrixsz  = 0;

                // -- applying matrix ---
                for( unsigned fcntx=0; fcntx<3; fcntx++ )
                {
                    for( unsigned fcnty=0; fcnty<3; fcnty++ )
                    {
                        unsigned posX = ( cntx - 3 / 2 + fcntx + img_width )
                                        % img_width;
                        unsigned posY = ( cnty - 3 / 2 + fcnty + img_height )
                                        % img_height;

                        unsigned posM = posY * img_width + posX;

                        if ( posM < pixel_arrays_realsz )
                        {
                             //medimatrix.push_back ( pixel_arrays[ posM ] );
                             medimatrix[ medimatrixsz ] = pixel_arrays[ posM ];
                             medimatrixsz++;
                        }
                    }
                }
                // -- applying matrix ---

                // sort it !
                qsort( medimatrix, medimatrixsz, sizeof(unsigned short), RAWProcessor_uscompare );

                copy_arrays[ cnty * img_width + cntx ] = medimatrix[ medimatrixsz / 2 ];
            }
        }

        memcpy( pixel_arrays.data(), copy_arrays, pixel_arrays_realsz * sizeof( unsigned short ) );

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

RAWProcessor::FilterConfig* RAWProcessor::GetPresetFilter( unsigned fnum )
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
    if ( pixel_arrays_realsz > 0 )
    {
        return rawimgtk::AdjustGamma( pixel_arrays.data(),
                                      pixel_arrays.size(),
                                      gamma );
    }

    return false;
}

bool RAWProcessor::AdjustBrightness( float percent )
{
    if ( pixel_arrays_realsz > 0 )
    {
        return rawimgtk::AdjustBrightness( pixel_arrays.data(),
                                           pixel_arrays.size(),
                                           percent );
    }

    return false;
}

bool RAWProcessor::AdjustContrast( float percent )
{
    if ( pixel_arrays_realsz > 0 )
    {
        return rawimgtk::AdjustContrast( pixel_arrays.data(),
                                         pixel_arrays.size(),
                                         percent );
    }

    return false;
}

bool RAWProcessor::AdjustToneMapping( unsigned ttype, float p1, float p2, float p3, float p4 )
{
    if ( pixel_arrays_realsz > 0 )
    {
        switch ( ttype )
        {
            case TONEMAP_TYPE_DRAGO:
                {
                    float gamma = p1;
                    float exposure = p2;

                    return rawimgtk::tmoDrago03( pixel_arrays.data(),
                                                 pixel_arrays.size(),
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

                    return rawimgtk::tmoReinhard2005( pixel_arrays.data(),
                                                      pixel_arrays.size(),
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

bool RAWProcessor::ApplyCLAHE( WeightAnalysisReport &report, unsigned applysz, unsigned bins, float slope )
{
    if ( pixel_arrays_realsz > 0 )
    {
        unsigned minv = report.threshold_wide_min;
        unsigned maxv = report.threshold_wide_max;
        unsigned short* ptr = pixel_arrays.data();

        return rawimgtk::ApplyCLAHE( ptr, img_width, img_height, minv, maxv,
                                     applysz, applysz, 0, slope );
    }

    return false;
}

bool RAWProcessor::ApplyLowFrequency( unsigned filtersz, unsigned repeat )
{
    if ( pixel_arrays_realsz > 0 )
    {
        unsigned short* ptr = pixel_arrays.data();

        return rawimgfk::ApplyLowFreqFilter( ptr,
                                             img_width, img_height,
                                             filtersz, repeat );
    }

    return false;
}

bool RAWProcessor::ApplyEdgeEnhance( unsigned fszh, unsigned fszv, unsigned edgesz, unsigned margin )
{
    if ( pixel_arrays_realsz > 0 )
    {
        unsigned short* ptr = pixel_arrays.data();
        unsigned        imgsz = pixel_arrays.size();

        unsigned short* imgEH1 = new unsigned short[ imgsz ];

        if ( imgEH1 == NULL )
            return false;

        unsigned short* imgEH2 = new unsigned short[ imgsz ];

        if ( imgEH2 == NULL )
        {
            delete[] imgEH1;

            return false;
        }

        memcpy( imgEH1, ptr, imgsz * sizeof( unsigned short ) );
        memcpy( imgEH2, ptr, imgsz * sizeof( unsigned short ) );

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

        unsigned cnth = 0;
        unsigned cntw = 0;
		unsigned mgnx = margin;
		unsigned mgny = margin;
		unsigned mgnw = img_width - (margin * 2);
		unsigned mgnh = img_height - (margin * 2);

        float fedgev = (float)edgesz / 8.0f;

        #pragma omp parallel for private( cntw )
        for( cnth=mgny; cnth<mgnh; cnth++ )
        {
            for( cntw=mgnx; cntw<mgnw; cntw++ )
            {
				unsigned pos = cnth * img_width + cntw;
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

				ptr[ pos ] = (unsigned short) pixelv;
            }
        }

        delete[] imgEH1;
        delete[] imgEH2;

        return true;
    }

    return false;
}

bool RAWProcessor::ApplyAnisotropicFilter( unsigned strength, unsigned param )
{
    if ( pixel_arrays_realsz > 0 )
    {
        unsigned short* ptr = pixel_arrays.data();

        return RAWImageFilterKit::ApplyAnisotropicFilter( ptr,
                                                          img_width, img_height,
                                                          strength, param );
    }

    return false;

}

const unsigned long RAWProcessor::datasize()
{
    //return pixel_arrays.size();
    return pixel_arrays_realsz;
}

const unsigned short* RAWProcessor::data()
{
    if ( pixel_arrays.size() == 0 )
        return NULL;

    return pixel_arrays.data();
}

void RAWProcessor::analyse()
{
    //if ( pixel_arrays.size() == 0 )
    if ( pixel_arrays_srcsz == 0 )
        return;

    //int pixel_counts = pixel_arrays.size();
    //if ( pixel_counts > 0 )
    if ( pixel_arrays_srcsz > 0 )
    {
        //img_width = pixel_counts / img_height;
        img_width = pixel_arrays_srcsz / img_height;
    }

    // measure Bits Per Pixel.
    // it starts from 10 bits.
    pixel_bpp = 10;
    if ( pixel_max_level > 0x07FF )
    {
        pixel_bpp = 12;

        if ( pixel_max_level > 0x1FFF )
        {
            pixel_bpp = 14;

            if ( pixel_max_level > 0x4FFF )
            {
                pixel_bpp = 16;
            }
        }
    }
}

void RAWProcessor::resetWeights()
{
#ifdef DEBUG
    if ( pixel_weights == NULL )
    {
        printf("????? why pixel_weights is NULL ???\n");
        pixel_weights = new unsigned int[ DEF_PIXEL_WEIGHTS + 1 ];
    }
#endif // DEBUG
    if ( pixel_weights != NULL )
    {
        memset( pixel_weights, 0 , ( DEF_PIXEL_WEIGHTS  + 1 ) * sizeof( unsigned int ) );
    }

    pixel_weights_max = 0;
}

void RAWProcessor::calcWeights()
{
    pixel_min_level = DEF_PIXEL16_MAX - 1;
    pixel_max_level = 0;
    pixel_med_level = DEF_PIXEL16_MAX / 2;

    for( unsigned cnt=0; cnt<pixel_arrays_srcsz; cnt++ )
    {
        if ( pixel_arrays[cnt] > pixel_max_level )
        {
            pixel_max_level = pixel_arrays[cnt];
        }
        else
        if ( pixel_arrays[cnt] < pixel_min_level )
        {
            pixel_min_level = pixel_arrays[cnt];
        }
    }

    pixel_med_level = ( pixel_max_level + pixel_min_level ) / 2;

#ifdef DEBUG
    printf( "calcWeights() result min,max,med = %d, %d, %d\n",
            pixel_min_level, pixel_max_level, pixel_med_level );
#endif // DEBUG
}

void RAWProcessor::addpixelarray( std::vector<unsigned short>* outpixels, unsigned x, unsigned y )
{
    if ( outpixels == NULL )
    {
        return;
    }

    if ( ( x < img_width ) && ( y < img_height ) )
    {
        unsigned mpos = img_width * y + x;
        if ( pixel_arrays_realsz > mpos )
        {
            outpixels->push_back( pixel_arrays[ mpos ] );
        }
    }
}

void RAWProcessor::reordercoords( std::vector<polygoncoord>* coords )
{
    if ( coords == NULL )
    {
        return;
    }

    unsigned ptsz = coords->size();

    if ( ptsz > 2 )
    {
        unsigned idxFirst = -1;
        unsigned minX = img_width;
        unsigned minY = img_height;

        for( unsigned cnt=0; cnt<(ptsz - 1); cnt++ )
        {
            if ( coords->at(cnt).y < minY )
            {
                minX = coords->at(cnt).x;
                minY = coords->at(cnt).y;
                idxFirst = cnt;
            }
        }

        // X is next.
        for( unsigned cnt=0; cnt<(ptsz - 1); cnt++ )
        {
            if ( ( coords->at(cnt).x < minX ) && ( coords->at(cnt).y < minY ) )
            {
                minX = coords->at(cnt).x;
                minY = coords->at(cnt).y;
                idxFirst = cnt;
            }
        }

        if ( idxFirst > 0 )
        {
            vector< polygoncoord > copydummy;
            copydummy.resize( ptsz );

            // copy all to dummy.
            for( unsigned cpcnt=0; cpcnt<ptsz-1; cpcnt++ )
            {
                copydummy[cpcnt].x = coords->at(cpcnt).x;
                copydummy[cpcnt].y = coords->at(cpcnt).y;
            }

            unsigned lastQ = 0;

            // copy back all except last coordination
            for( unsigned cnt=idxFirst; cnt<(ptsz - 1); cnt++ )
            {
                coords->at(cnt-idxFirst).x = copydummy[cnt].x;
                coords->at(cnt-idxFirst).y = copydummy[cnt].y;

                lastQ = cnt - idxFirst;
            }

            // Last end point will be re-assign for new vector's first.

            for( unsigned cnt=0; cnt<idxFirst; cnt++ )
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

void RAWProcessor::CutoffLevels( unsigned short minv, unsigned short maxv )
{
    #pragma omp parallel for
    for( unsigned cnt=0; cnt<pixel_arrays_srcsz; cnt++ )
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

void RAWProcessor::CutoffLevelsRanged( unsigned short minv, unsigned short maxv, unsigned short valmin, unsigned short valmax )
{
    #pragma omp parallel for
    for( unsigned cnt=0; cnt<pixel_arrays_srcsz; cnt++ )
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
