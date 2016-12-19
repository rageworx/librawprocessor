#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cstring>
#include <list>
#include <iostream>
#include <fstream>

#ifdef RAWPROCESSOR_USE_LOCALTCHAR
	#include "tchar.h"
#else
	#ifdef _WIN32
		#include <tchar.h>
	#else
		#include "tchar.h"
	#endif
#endif

#include "rawprocessor.h"
#include "rawscale.h"
#include "rawimgtk.h"
#include "stdunicode.h"

////////////////////////////////////////////////////////////////////////////////

using namespace std;

////////////////////////////////////////////////////////////////////////////////

#ifdef UNICODE
    #define _TSTRING    wstring
    #define _TCM2W( _x_ )  convertM2W( (const char*)_x_ )
    #define _TCW2M( _x_ )  convertW2M( (const wchar_t*)_x_ )
#else
    #define _TSTRING    string
    #define _TCM2W( _x_ )  _x_
    #define _TCW2M( _x_ )  _x_
#endif

#define DEF_MEMORY_LOADED   "//MEMORY_LOAD//"
#define DEF_RAW_I_HEIGHT    1024
#define DEF_PIXEL_WEIGHTS   65535

#define DEF_PIXEL16_MAX     65536
#define DEF_PIXEL8_MAX      256

#define DEF_CALC_F_WMAX     65535.0f
#define DEF_CALC_I_WMAX     65535
#define DEF_CALC_F_BMAX     255.0f
#define DEF_CALC_I_BMAX     255

#define DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY     0,9,10,54

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

RAWProcessor::RAWProcessor()
 : raw_loaded(false),
   pixel_min_level(0),
   pixel_max_level(0),
   pixel_arrays_realsz( 0 ),
   pixel_weights( NULL ),
   pixel_weights_max( 0 ),
   img_height(0),
   img_width(0),
   userscaler(NULL)
{
    pixel_weights = new unsigned int[ DEF_PIXEL_WEIGHTS + 1 ];
    resetWeights();
}

RAWProcessor::RAWProcessor( const char* raw_file, int height )
 : raw_loaded(false),
   pixel_min_level(0),
   pixel_max_level(0),
   pixel_weights( NULL ),
   pixel_weights_max( 0 ),
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
        return;

    int retia[4] = { DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY };

    char retvstr[32] = {0};
    sprintf( retvstr, "%d.%d.%d.%d", retia[0], retia[1], retia[2], retia[3] );

    int retvstrsz = strlen( retvstr );

    *retverstr = (char*)malloc( retvstrsz + 1 );
    if ( *retverstr != NULL )
    {
        memset( *retverstr, 0, retvstrsz + 1 );
        strcpy( *retverstr, retvstr );
    }
}

void RAWProcessor::Version( int** retverints )
{
    if ( retverints == NULL )
        return;

    int retia[4] = { DEF_LIBRAWPROCESSOR_VERSION_I_ARRAY };
    memcpy( *retverints, retia, sizeof(int)*4 );
}

bool RAWProcessor::Load( const wchar_t* raw_file, unsigned int trnsfm, int height )
{
    if ( height <= 0 )
        return false;

    string fname = convertW2M( raw_file );

    return Load( fname.c_str(), trnsfm, height );
}

bool RAWProcessor::Load( const char* raw_file, unsigned int trnsfm, int height )
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

        int pxlsz   = fsize / sizeof( unsigned short );
        int blancsz = 0;

        if (  pxlsz > 0 )
        {
            // set temporary width ...
            img_width = pixel_arrays_realsz / img_height;

            // Check omit pixels ...
            if ( pxlsz < ( img_width * img_height ) )
            {
                blancsz = ( img_width * img_height ) - pxlsz;
            }
        }

        pixel_arrays_realsz = pxlsz + blancsz;
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
            worddata = ( chardata[1] << 8 ) + ( chardata[0] & 0x00FF );
            //pixel_arrays.push_back( worddata );
            pixel_arrays[ cnt ] = worddata;
            pixel_weights[ worddata ] += 1;

            if ( pixel_weights_max < worddata )
            {
                pixel_weights_max = worddata;
            }

            if ( worddata < pixel_min_level )
            {
                pixel_min_level = worddata;
            }
            else
            if ( worddata > pixel_max_level )
            {
                pixel_max_level = worddata;
                index_max_pixel = cnt;
            }

        }

        pixel_med_level = ( pixel_max_level + pixel_min_level ) / 2;

        rfstrm.clear();
        rfstrm.close();

        analyse();

        raw_loaded = true;

        ApplyTransform( trnsfm );

        return raw_loaded;
    }

    return false;
}

bool RAWProcessor::LoadFromMemory( const char* buffer, unsigned long bufferlen, unsigned int trnsfm, int height )
{
    if ( height <= 0 )
        return false;

    if ( ( buffer != NULL ) && ( bufferlen > 0 ) )
    {
        pixel_min_level = 0;
        pixel_med_level = 0;
        pixel_max_level = 0;
        img_height = height;

        int pxlsz   = bufferlen / sizeof( unsigned short );
        int blancsz = 0;

        if (  pxlsz > 0 )
        {
            // set temporary width ...
            img_width = pixel_arrays_realsz / img_height;

            // Check omit pixels ...
            if ( pxlsz < ( img_width * img_height ) )
            {
                blancsz = ( img_width * img_height ) - pxlsz;
            }
        }

        pixel_arrays_realsz = pxlsz + blancsz;

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

            if ( worddata < pixel_min_level )
            {
                pixel_min_level = worddata;
            }
            else
            if ( worddata > pixel_max_level )
            {
                pixel_max_level = worddata;
                index_max_pixel = cnt;
            }

        }

        pixel_med_level = ( pixel_max_level + pixel_min_level ) / 2;

        analyse();

        raw_loaded = true;

        ApplyTransform( trnsfm );

        raw_file_name = _TEXT(DEF_MEMORY_LOADED);

        return true;
    }

    return false;
}

bool RAWProcessor::Reload( const wchar_t* raw_file, unsigned int trnsfm, int height )
{
    return Load( raw_file, trnsfm, height );
}

bool RAWProcessor::Reload( const char* raw_file, unsigned int trnsfm, int height )
{
    return Load( raw_file, trnsfm, height );
}

bool RAWProcessor::Reload()
{
    if ( raw_file_name == _TEXT(DEF_MEMORY_LOADED) )
        return false;

    return Reload( _TCM2W(raw_file_name.c_str()) );
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
                retb = RAWImageToolKit::FlipHorizontal( ptrd, img_width, img_height );
            }

            if ( ( trnsfm & TRANSFORM_PARAM_SWAP_V ) > 0 )
            {
                retb = RAWImageToolKit::FlipVertical( ptrd, img_width, img_height );
            }
        }

        if ( ( trnsfm & TRANSFORM_ROTATE ) > 0 )
        {
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_C90 ) > 0 )
            {
                retb = RAWImageToolKit::Rotate90( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_C180 ) > 0 )
            {
                retb = RAWImageToolKit::Rotate180( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_C270 ) > 0 )
            {
                retb = RAWImageToolKit::Rotate270( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_RC90 ) > 0 )
            {
                retb = RAWImageToolKit::Rotate270( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_RC180 ) > 0 )
            {
                retb = RAWImageToolKit::Rotate180( ptrd, &img_width, &img_height );
            }
            else
            if ( ( trnsfm & TRANSFORM_PARAM_ROT_RC270 ) > 0 )
            {
                retb = RAWImageToolKit::Rotate90( ptrd, &img_width, &img_height );
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

bool RAWProcessor::Reverse( unsigned char maxbits )
{
    if ( maxbits < 8 )
        return false;

    unsigned short maxval = pow( 2, maxbits );

    unsigned dlen = pixel_arrays.size();

    for( unsigned cnt=0; cnt<dlen; cnt++ )
    {
        pixel_arrays[ cnt ] = maxval - pixel_arrays[ cnt ];
    }

    return true;
}

void RAWProcessor::ChangeHeight( int h )
{
    if ( raw_loaded == false )
        return;

    if ( ( h > 0 ) && ( img_height != h ) )
    {
        img_height = h;
        analyse();
    }
}

bool RAWProcessor::Get8bitDownscaled( vector<unsigned char> *byte_arrays, DownscaleType dntype, bool reversed )
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

        for( unsigned cnt=0; cnt<arrsz; cnt++ )
        {
            float fdspixel = float(ref_pixel_arrays[cnt]) * dscale_ratio;
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

        for( unsigned cnt=0; cnt<arrsz; cnt++ )
        {
            float fuspixel = float(ref_pixel_arrays[cnt]) * uscale_ratio;
            float fdspixel = fuspixel * dscale_ratio;
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

bool RAWProcessor::Get16bitRawImage( std::vector<unsigned short> *word_arrays, bool reversed )
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

bool RAWProcessor::GetWeights( std::vector<unsigned int> *weight_arrays )
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
    int identify_min_level = int( float(pixel_max_level) * 0.2f );
    if ( start_minlevel_zero == true )
    {
        identify_min_level = 0;
    }

    unsigned index_center_thld  = 0;

    for( int cnt=identify_min_level; cnt<pixel_max_level; cnt++ )
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

    if ( index_center_thld > 1 )
    {
        report.base_threshold_index = index_center_thld;
        int            weight_idx   = index_center_thld;
        unsigned int   curweight    = pixel_weights[ weight_idx ];

        // # phase 03-01
        // under counting ...
        while( true )
        {
            weight_idx--;

            curweight = pixel_weights[ weight_idx ];

            if ( ( curweight <= (unsigned)identify_min_level ) ||
                 ( weight_idx == identify_min_level ) )
            {
                min_weight_wide = weight_idx;
                break;
            }

            if ( curweight > report.threshold_max_amount )
            {
                report.threshold_max_amount = curweight;
            }
        }

        // #phase 03-02
        // over counting ...
        weight_idx = index_center_thld;
        curweight  = pixel_weights[ weight_idx ];
        while( true )
        {
            weight_idx++;
            curweight = pixel_weights[ weight_idx ];

            if ( ( curweight <= (unsigned)identify_min_level ) ||
                 ( weight_idx == ( DEF_PIXEL_WEIGHTS - 1 ) ) )
            {
                max_weight_wide = weight_idx;
                break;
            }

            if ( curweight > report.threshold_max_amount )
            {
                report.threshold_max_amount = curweight;
            }
        }
    }
    else
    {
        return false;
    }

    // # phase 04
    // get wide count.
    report.threshold_wide_min = min_weight_wide;
    report.threshold_wide_max = max_weight_wide;

    return true;
}

bool RAWProcessor::Get16bitThresholdedImage( WeightAnalysisReport &report,  vector<unsigned short> *word_arrays, bool reversed )
{
    if ( report.timestamp == 0 )
        return false;

    int thld_min = report.threshold_wide_min;
    int thld_max = report.threshold_wide_max;
    int thld_wid = report.threshold_wide_max - report.threshold_wide_min;

    float normf = float( DEF_PIXEL_WEIGHTS ) / float( thld_wid );

    word_arrays->clear();

    int array_max = pixel_arrays.size();

    word_arrays->reserve( array_max );
    word_arrays->resize( array_max );

    for( int cnt=0; cnt<array_max; cnt++ )
    {
        unsigned short apixel = pixel_arrays[cnt] - thld_min - 1;
        float          fpixel = 0.0f;

        // cut off threhold pixel value.
        if ( apixel > thld_max )
        {
            //apixel = thld_max;
            apixel = DEF_CALC_F_WMAX;
        }
        else
        if ( apixel < thld_min )
        {
            apixel = 0;
        }

        apixel -= thld_min;

        // Rescaling pixel  !
        fpixel = float( apixel ) * normf;

        apixel = (unsigned short)( fpixel );

        if ( reversed == true )
        {
            apixel = DEF_PIXEL_WEIGHTS - apixel;
        }

        //word_arrays.push_back( apixel );
        word_arrays->at(cnt) = apixel;
    }

    return true;
}

bool RAWProcessor::Get8bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned char> *byte_arrays, bool reversed )
{
    if ( report.timestamp == 0 )
        return false;

    int thld_min = report.threshold_wide_min;
    int thld_max = report.threshold_wide_max;

    float normf  = float( DEF_CALC_F_BMAX ) / float( thld_max - thld_min );

    byte_arrays->clear();

    int array_max = pixel_arrays.size();

    byte_arrays->reserve( array_max );
    byte_arrays->resize( array_max );

    for( int cnt=0; cnt<array_max; cnt++ )
    {
        unsigned short apixel = pixel_arrays[cnt];
        float          fpixel = 0.0f;
        unsigned char  bpixel = 0;

        // cut off threhold pixel value.
        if ( apixel > thld_max )
        {
            //apixel = thld_max;
            apixel = DEF_CALC_F_WMAX;
        }
        else
        if ( apixel < thld_min )
        {
            apixel = 0;
        }

        // Rescaling pixel  !
        fpixel = float( apixel ) * normf;

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

bool RAWProcessor::Get16bitPixel( int x, int y, unsigned short &px )
{
    if ( pixel_arrays.size() == 0 )
        return false;

    if ( ( x < 0 ) || ( y < 0 ) || ( x > img_width ) || ( y > img_height ) )
        return false;

    int pixpos = ( y * img_height ) + x;

    px = pixel_arrays[ pixpos ];

    return true;
}

RAWProcessor* RAWProcessor::rescale( int w, int h, RescaleType st )
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
                                                       current_transform, h );
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
    if ( pixel_arrays_realsz == 0 )
        return;

    //int pixel_counts = pixel_arrays.size();
    //if ( pixel_counts > 0 )
    if ( pixel_arrays_realsz > 0 )
    {
        //img_width = pixel_counts / img_height;
        img_width = pixel_arrays_realsz / img_height;
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
