#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <list>
#include <iostream>
#include <fstream>

#if defined(_DEBUG)&&defined(_USEFLTK)
    #include <FL/Fl_Ask.H>
#endif

#include "stdunicode.h"
#include "rawprocessor.h"

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

#define DEF_RAW_I_HEIGHT    1024
#define DEF_PIXEL_WEIGHTS   65535

#define DEF_CALC_F_WMAX     65536.0f
#define DEF_CALC_I_WMAX     65536
#define DEF_CALC_F_BMAX     256.0f
#define DEF_CALC_I_BMAX     256

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

RAWProcessor::RAWProcessor()
 : raw_loaded(false),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(DEF_RAW_I_HEIGHT),
   img_width(0),
   userscaler(NULL)
{
    resetWeights();
}

RAWProcessor::RAWProcessor( const char* raw_file, int height )
 : raw_loaded(false),
   pixel_min_level(0),
   pixel_max_level(0),
   img_height(height),
   img_width(0)
{
    if ( raw_file != NULL )
    {
        Load( raw_file );
    }
}

RAWProcessor::~RAWProcessor()
{
    if ( raw_loaded == true )
    {
        Unload();
    }
}

bool RAWProcessor::Load( const char* raw_file, int height )
{
    fstream rfstrm;
    //string  fname = _TCW2M( raw_file );
	string  fname = raw_file;

    rfstrm.open( fname.c_str(), fstream::app | fstream::binary | fstream::in );

    if ( rfstrm.is_open() == true )
    {
        rfstrm.seekg( 0, ios::end );
        unsigned int fsize = rfstrm.tellg();
        rfstrm.seekg( 0, ios::beg );
        rfstrm.clear();

        pixel_min_level = 0;
        pixel_med_level = 0;
        pixel_max_level = 0;
        img_height = height;

        pixel_arrays.clear();
        resetWeights();

        // To save memory, it doesn't using direct load to memory.
        for( unsigned cnt=0; cnt<fsize / sizeof(unsigned short); cnt++)
        {
            char           chardata[2] = {0};
            unsigned short worddata = 0;

            rfstrm.read( chardata, sizeof(unsigned short) );
            worddata = ( chardata[1] << 8 ) + ( chardata[0] & 0x00FF );
            pixel_arrays.push_back( worddata );
            pixel_weights[ worddata ] ++;

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

        if ( pixel_arrays.size() > 0 )
        {
            raw_loaded = true;
        }

        return raw_loaded;
    }
#if defined(_DEBUG)&&defined(_USEFLTK)
    fl_alert("DEBUG.ERROR: can not open : %s !!", fname.c_str() );
#endif
    return false;
}

bool RAWProcessor::Reload( const char* raw_file, int height )
{
    Unload();
    return Load( raw_file, height );
}

bool RAWProcessor::Reload()
{
    //return Reload( _TCM2W(raw_file_name.c_str()) );
    return Reload( raw_file_name.c_str() );
}

void RAWProcessor::Unload()
{
    raw_loaded = false;

    pixel_arrays.clear();
    raw_file_name.clear();
    pixel_weights.clear();
}

void RAWProcessor::SetUserScale( RAWUserScaleIF* ptr )
{
    userscaler = ptr;
}

bool RAWProcessor::Get8bitDownscaled( vector<unsigned char> &byte_arrays, DownscaleType dntype )
{
    if ( raw_loaded == false )
        return false;

    byte_arrays.clear();

    if ( dntype == DNSCALE_NORMAL )
    {
        float dscale_ratio = DEF_CALC_F_BMAX / float( pixel_max_level );

        for( unsigned cnt=0; cnt<pixel_arrays.size(); cnt++ )
        {
            float fdspixel = float(pixel_arrays[cnt]) * dscale_ratio;
            unsigned char dspixel = (unsigned char)fdspixel;
            byte_arrays.push_back( dspixel );
        }

    }
    else
    if ( dntype == DNSCALE_FULLDOWN )
    {
        float uscale_ratio = DEF_CALC_F_WMAX / float( pixel_max_level );
        float dscale_ratio = DEF_CALC_F_BMAX / DEF_CALC_F_WMAX;

        for( unsigned cnt=0; cnt<pixel_arrays.size(); cnt++ )
        {
            float fuspixel = float(pixel_arrays[cnt]) * uscale_ratio;
            float fdspixel = fuspixel * dscale_ratio;
            unsigned char dspixel = (unsigned char)fdspixel;
            byte_arrays.push_back( dspixel );
        }
    }
    else
    if ( userscaler != NULL )
    {
        return userscaler->processUserScale( pixel_arrays, byte_arrays );
    }

    return true;
}

bool RAWProcessor::Get16bitRawImage( std::vector<unsigned short> &word_arrays )
{
    if ( raw_loaded == false )
        return false;

    word_arrays.clear();
    word_arrays.assign( pixel_arrays.begin(), pixel_arrays.end() );

    return true;
}

bool RAWProcessor::GetWeights( std::vector<int> &weight_arrays )
{
    if ( raw_loaded == false )
        return false;

    weight_arrays.clear();
    weight_arrays.assign( pixel_weights.begin(), pixel_weights.end() );

    return true;
}

bool RAWProcessor::GetAnalysisReport( WeightAnalysisReport &report )
{
    if ( pixel_weights.size() == 0 )
        return false;

    // # phase 01
    // get current time stamp.
    time_t curtime;
    report.timestamp = (unsigned long)time(&curtime);

    // # phase 02
    // find highest value in pixels ... ( 20% of maximum level )
    int identify_min_level = int( float(pixel_max_level) * 0.1f );
    int index_center_thld  = 0;

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

            if ( curweight > report.htreshold_max_amount )
            {
                report.htreshold_max_amount = curweight;
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

            if ( curweight > report.htreshold_max_amount )
            {
                report.htreshold_max_amount = curweight;
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

bool RAWProcessor::Get16bitThresholdedImage( WeightAnalysisReport &report,  vector<unsigned short> &word_arrays )
{
    if ( report.timestamp == 0 )
        return false;

    int thld_min = report.threshold_wide_min;
    int thld_max = report.threshold_wide_max;
    int thld_wid = report.threshold_wide_max - report.threshold_wide_min;

    float normf = DEF_PIXEL_WEIGHTS / float( thld_wid );

    word_arrays.clear();

    int array_max = pixel_arrays.size();

    for( int cnt=0; cnt<array_max; cnt++ )
    {
        unsigned short apixel = pixel_arrays[cnt] - thld_min - 1;
        float          fpixel = 0.0f;

        // cut off threhold pixel value.
        if ( apixel > thld_max )
        {
            //apixel = thld_max;
            apixel = DEF_CALC_I_WMAX - 1;
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

        word_arrays.push_back( apixel );
    }

    return true;
}

bool RAWProcessor::Get8bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned char> &byte_arrays )
{
    if ( report.timestamp == 0 )
        return false;

    int thld_min = report.threshold_wide_min;
    int thld_max = report.threshold_wide_max;

    float normf  = DEF_CALC_F_BMAX / float( report.threshold_wide_max );

    byte_arrays.clear();

    int array_max = pixel_arrays.size();

    for( int cnt=0; cnt<array_max; cnt++ )
    {
        unsigned short apixel = pixel_arrays[cnt];
        float          fpixel = 0.0f;
        unsigned char  bpixel = 0;

        // cut off threhold pixel value.
        if ( apixel > thld_max )
        {
            //apixel = thld_max;
            apixel = DEF_CALC_I_WMAX - 1;
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
            bpixel = 255;
        }
        else
        {
            bpixel = (unsigned char)( fpixel );
        }

        byte_arrays.push_back( bpixel );
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

const unsigned long RAWProcessor::datasize()
{
    return pixel_arrays.size();
}

const unsigned short* RAWProcessor::data()
{
    if ( pixel_arrays.size() == 0 )
        return NULL;

    return pixel_arrays.data();
}

void RAWProcessor::analyse()
{
    if ( pixel_arrays.size() == 0 )
        return;

    int pixel_counts = pixel_arrays.size();
    img_width = pixel_counts / img_height;
}

void RAWProcessor::resetWeights()
{
    pixel_weights.clear();
    pixel_weights.reserve(DEF_PIXEL_WEIGHTS);

    for (int cnt=0; cnt<DEF_PIXEL_WEIGHTS; cnt++)
    {
        pixel_weights.push_back(0);
    }
}
