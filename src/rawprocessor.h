#ifndef __RAWPROCESSOR_H__
#define __RAWPROCESSOR_H__
////////////////////////////////////////////////////////////////////////////////
//
//  16bit gray scale RAW image processor for stdc++ w/ FLTK or other Graphics.
//  =========================================================================
//  (C)Copyright 2013~, Raphael Kim (rageworx@gmail.com)
//
//  * Processing read pixels to down-scaling, threshold cutting.
//
//  _ updates _
//
//  2016-12-19 :
//       - Added some graphical help methods : AdjustGamma, AdjustBrightness,
//                                             AdjustContrast, AdjustContrast,
//                                             AdjustCurve
//       - Added Version check method.
//
//  2016-11-25 :
//       - set vector reserve() and resize() for allocating more faster !
//       - some optimized code for accelerated by AVX instruction.
//
//  2016-12-09
//       - Moved some enum and struct to inside of class.
//       - RAW resize included, source code refer to FreeImage open source.
//         FreeImage project is place at http://freeimage.sourceforge.net/
//
// 2016-12-23
//       - Changed some method name in public, removed "get*()" scheme to "*()".
//       - Renamed rescale() to Rescale()
//       - Added SaveToFile() methods.
//
// 2016-12-29
//       - Added Filter processing.
//         referenced to http://lodev.org/cgtutor/filtering.html
//         basic filter algorithm by Lode Vandevenne
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

class RAWUserScaleIF
{
    public:
        virtual bool processUserScale( std::vector<unsigned short> &src,
                                       std::vector<unsigned char> *dst ) = 0;
};

class RAWProcessor
{
    public:
        typedef enum
        {
            DNSCALE_NORMAL = 0,
            DNSCALE_FULLDOWN,
            DNSCALE_USER,
            DNSCALE_MAX
        }DownscaleType;

        typedef enum
        {
            RESCALE_NEAREST = 0,      /// == BOX
            RESCALE_BILINEAR,
            RESCALE_BICUBIC,
            RESCALE_BSPLINE,
            RESCALE_LANZCOS3
        }RescaleType;

        #define TRANSFORM_NONE              0x00000000

        #define TRANSFORM_SWAP              0x00000010
        #define TRANSFORM_PARAM_SWAP_H      0x00000002
        #define TRANSFORM_PARAM_SWAP_V      0x00000004

        #define TRANSFORM_ROTATE            0x00010000
        #define TRANSFORM_PARAM_ROT_C90     0x00000200
        #define TRANSFORM_PARAM_ROT_C180    0x00000400
        #define TRANSFORM_PARAM_ROT_C270    0x00000800
        #define TRANSFORM_PARAM_ROT_RC90    0x00002000
        #define TRANSFORM_PARAM_ROT_RC180   0x00004000
        #define TRANSFORM_PARAM_ROT_RC270   0x00008000

    public:
        struct WeightAnalysisReport
        {
            unsigned long       timestamp;
            unsigned short      base_threshold_index;
            unsigned short      threshold_wide_min;
            unsigned short      threshold_wide_max;
            unsigned int        threshold_max_amount;
        };

        struct SimpleAnalysisInfo
        {
            unsigned short      minLevel;
            unsigned short      maxLevel;
            double              average;
            double              variance;
            double              deviation;
        };

    public:
        RAWProcessor();
        RAWProcessor( const char* raw_file, unsigned int height = 0 );
        virtual~RAWProcessor();

    public:
        bool Loaded()                   { return raw_loaded; }
        unsigned PixelCount()           { return datasize(); }
        unsigned short MinimumLevel()   { return pixel_min_level; }
        unsigned short MaximumLevel()   { return pixel_max_level; }
        unsigned short MediumLevel()    { return pixel_med_level; }
        unsigned Width()                { return img_width; }
        unsigned Height()               { return img_height; }
        unsigned WeightsCount()         { return pixel_weights_max; }
        unsigned char BPP()             { return pixel_bpp; }

    public:
        void Version( char** retverstr ); /// put NULL initialized char* array.
        void Version( int** retverints ); /// put int[4] array.
        bool Load( const char* raw_file, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
        bool Load( const wchar_t* raw_file, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
        bool LoadFromMemory( const char* buffer, unsigned long bufferlen, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
        bool Reload( const char* raw_file, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
        bool Reload( const wchar_t* raw_file, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
        bool Reload();
        void Unload();
        bool ApplyTransform( unsigned int trnsfm = TRANSFORM_NONE );
        void ChangeHeight( unsigned h );
        void SetUserScale( RAWUserScaleIF* ptr = NULL );
        bool Reverse( unsigned char maxbits = 16 );
        bool Get8bitDownscaled( std::vector<unsigned char>* byte_arrays, DownscaleType dntype = DNSCALE_NORMAL, bool reversed = false );
        bool Get16bitRawImage( std::vector<unsigned short>* word_arrays, bool reversed = false );
        bool GetWeights( std::vector<unsigned int>* weight_arrays );
        bool GetAnalysisReport( WeightAnalysisReport &report, bool start_minlevel_zero = false );
        bool Get16bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned short>* word_arrays, bool reversed = false );
        bool Get8bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned char>* byte_arrays, bool reversed = false );
        bool Get16bitPixel( unsigned x, unsigned y, unsigned short &px );

    // Some additional tools here ...
    public:
        bool SaveToFile( const char* path );
        bool SaveToFile( const wchar_t* path );

    public:
        RAWProcessor* Rescale( unsigned w, unsigned h, RescaleType st = RESCALE_NEAREST );
        RAWProcessor* Clone();

    public:
        typedef struct { unsigned x; unsigned y; } polygoncoord;
        // Weights == Histogram.
        void GetLinearPixels( unsigned x1, unsigned y1, unsigned x2, unsigned y2, std::vector<unsigned short>* pixels );
        void GetRectPixels( unsigned x, unsigned y, unsigned w, unsigned h, std::vector<unsigned short>* pixels);
        void GetPolygonPixels( std::vector<polygoncoord>* coords, std::vector<unsigned short>* pixels);
        void GetAnalysisFromPixels( std::vector<unsigned short>* pixels, std::vector<unsigned int>* weights, SimpleAnalysisInfo* info );

    public:
        typedef struct
        {
            unsigned char   width;
            unsigned char   height;
            float           factor;
            float           bias;
            std::vector< float > matrix;
        }FilterConfig;

        #define PRESET_FILTER_NONE      0
        #define PRESET_FILTER_BLUR      1
        #define PRESET_FILTER_BLURMORE  2
        #define PRESET_FILTER_SHARPEN   3
        #define PRESET_FILTER_UNSHARPEN 4

        bool          ApplyFilter( FilterConfig* fconfig );
        RAWProcessor* CloneWithFilter( FilterConfig* fconfig );
        bool          ApplyMedianFilter();

        FilterConfig* GetPresetFilter( unsigned fnum );
        void          DiscardFilter( FilterConfig* fp );

    public:
        bool AdjustGamma( float gamma );
        bool AdjustBrightness( float percent );
        bool AdjustContrast( float percent );

    public:
        const unsigned long         datasize();
        const unsigned short*       data();

    protected:
        void analyse();
        void resetWeights();

    protected:
        void addpixelarray( std::vector<unsigned short>* outpixels, unsigned x, unsigned y );
        void reordercoords( std::vector<polygoncoord>* coords );

    protected:
        bool                        bigendian;
        bool                        raw_loaded;
        std::vector<unsigned short> pixel_arrays;
        unsigned long               pixel_arrays_srcsz;
        unsigned long               pixel_arrays_realsz;
        unsigned int*               pixel_weights;
        unsigned short              pixel_weights_max;
        unsigned char               pixel_bpp;
        unsigned short              pixel_min_level;
        unsigned short              pixel_max_level;
        unsigned short              pixel_med_level;
        unsigned short              index_max_pixel;
        unsigned int                current_transform;
#ifdef UNICODE
        std::wstring                raw_file_name;
#else
        std::string                 raw_file_name;
#endif
        unsigned                    img_height;
        unsigned                    img_width;
        RAWUserScaleIF*             userscaler;

};


#endif /// of __RAWPROCESSOR_H__
