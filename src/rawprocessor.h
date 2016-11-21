#ifndef __RAWPROCESSOR_H__
#define __RAWPROCESSOR_H__
////////////////////////////////////////////////////////////////////////////////
//
//  16bit gray scale RAW image processor for stdc++ w/ FLTK
//  =========================================================================
//  (C)Copyright 2013, Raphael Kim
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

typedef enum
{
    DNSCALE_NORMAL = 0,
    DNSCALE_FULLDOWN,
    DNSCALE_USER,
    DNSCALE_MAX
}DownscaleType;

struct WeightAnalysisReport
{
    unsigned long       timestamp;
    unsigned short      base_threshold_index;
    unsigned short      threshold_wide_min;
    unsigned short      threshold_wide_max;
    unsigned int        htreshold_max_amount;
};

class RAWUserScaleIF
{
    public:
        virtual bool processUserScale( std::vector<unsigned short> &src,
                                       std::vector<unsigned char> &dst ) = 0;
};

class RAWProcessor
{
    public:
        RAWProcessor();
        RAWProcessor( const char* raw_file, int height = 1024 );
        virtual~RAWProcessor();

    public:
        bool            isLoaded()          { return raw_loaded; }
        int             getPixelCount()     { return pixel_arrays.size(); }
        unsigned short  getMinimumLevel()   { return pixel_min_level; }
        unsigned short  getMaximumLevel()   { return pixel_max_level; }
        unsigned short  getMediumLevel()    { return pixel_med_level; }
        unsigned short* refPixelPointer()   { return pixel_arrays.data(); }
        int             getWidth()          { return img_width; }
        int             getHeight()         { return img_height; }
        int             getWeightsCount()   { return pixel_weights.size(); }

    public:
        bool Load( const char* raw_file, int height = 1024 );
        bool Reload( const char* raw_file, int height = 1024 );
        bool Reload();
        void Unload();
        void SetUserScale( RAWUserScaleIF* ptr );
        bool Get8bitDownscaled( std::vector<unsigned char> &byte_arrays, DownscaleType dntype = DNSCALE_NORMAL );
        bool Get16bitRawImage( std::vector<unsigned short> &word_arrays );
        bool GetWeights( std::vector<int> &weight_arrays );
        bool GetAnalysisReport( WeightAnalysisReport &report );
        bool Get16bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned short> &word_arrays );
        bool Get8bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned char> &byte_arrays );
        bool Get16bitPixel( int x, int y, unsigned short &px );

    public:
        const unsigned long         datasize();
        const unsigned short*       data();

    protected:
        void analyse();
        void resetWeights();

    protected:
        bool                        raw_loaded;
        std::vector<unsigned short> pixel_arrays;
        std::vector<int>            pixel_weights;
        unsigned short              pixel_min_level;
        unsigned short              pixel_max_level;
        unsigned short              pixel_med_level;
        unsigned short              index_max_pixel;
#ifdef UNICODE
        std::wstring                raw_file_name;
#else
        std::string                 raw_file_name;
#endif
        int                         img_height;
        int                         img_width;
        RAWUserScaleIF*             userscaler;

};


#endif /// of __RAWPROCESSOR_H__
