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
//  2016-11-25 :
//       - set vector reserve() and resize() for allocating more faster !
//       - some optimized code for accelerated by AVX instruction.
//
//  2016-12-0
//       - Moved some enum and struct to inside of class.
//       - RAW resize included, source code refer to FreeImage open source.
//         FreeImage project is place at http://freeimage.sourceforge.net/
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

    public:
        RAWProcessor();
        RAWProcessor( const char* raw_file, int height = 1024 );
        virtual~RAWProcessor();

    public:
        bool isLoaded() { return raw_loaded; }
        int  getPixelCount() { return pixel_arrays.size(); }
        unsigned short getMinimumLevel() { return pixel_min_level; }
        unsigned short getMaximumLevel() { return pixel_max_level; }
        unsigned short getMediumLevel()  { return pixel_med_level; }
        unsigned short* refPixelPointer() { return pixel_arrays.data(); }
        int  getWidth() { return img_width; }
        int  getHeight() { return img_height; }
        int  getWeightsCount() { return pixel_weights_max; }

    public:
        bool Load( const char* raw_file, unsigned int trnsfm = TRANSFORM_NONE, int height = 0 );
        bool Load( const wchar_t* raw_file, unsigned int trnsfm = TRANSFORM_NONE, int height = 0 );
        bool LoadFromMemory( const char* buffer, unsigned long bufferlen, unsigned int trnsfm = TRANSFORM_NONE, int height = 0 );
        bool Reload( const char* raw_file, unsigned int trnsfm = TRANSFORM_NONE, int height = 0 );
        bool Reload( const wchar_t* raw_file, unsigned int trnsfm = TRANSFORM_NONE, int height = 0 );
        bool Reload();
        void Unload();
        bool ApplyTransform( unsigned int trnsfm = TRANSFORM_NONE );
        void ChangeHeight( int h );
        void SetUserScale( RAWUserScaleIF* ptr );
        bool Reverse( unsigned char maxbits = 16 );
        bool Get8bitDownscaled( std::vector<unsigned char> *byte_arrays, DownscaleType dntype = DNSCALE_NORMAL, bool reversed = false );
        bool Get16bitRawImage( std::vector<unsigned short> *word_arrays, bool reversed = false );
        bool GetWeights( std::vector<unsigned int> *weight_arrays );
        bool GetAnalysisReport( WeightAnalysisReport &report, bool start_minlevel_zero = false );
        bool Get16bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned short> *word_arrays, bool reversed = false );
        bool Get8bitThresholdedImage( WeightAnalysisReport &report, std::vector<unsigned char> *byte_arrays, bool reversed = false );
        bool Get16bitPixel( int x, int y, unsigned short &px );

    public:
        RAWProcessor* rescale( int w, int h, RescaleType st = RESCALE_NEAREST );

    public:
        const unsigned long         datasize();
        const unsigned short*       data();

    protected:
        void analyse();
        void resetWeights();

    protected:
        bool                        raw_loaded;
        std::vector<unsigned short> pixel_arrays;
        unsigned long               pixel_arrays_realsz;
        unsigned int*               pixel_weights;
        unsigned short              pixel_weights_max;
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
        int                         img_height;
        int                         img_width;
        RAWUserScaleIF*             userscaler;

};


#endif /// of __RAWPROCESSOR_H__
