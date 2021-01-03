#ifndef __RAWPROCESSOR_H__
#define __RAWPROCESSOR_H__
////////////////////////////////////////////////////////////////////////////////
//
//  32bit precise gray scale RAW image processor for stdc++.
//  ============================================================================
//  (C)Copyright 2013~, Raphael Kim (rageworx@gmail.com)
//
//  * Processing read pixels to down-scaling, threshold cutting.
//
//  _ updates _
//
//  See CHANGES file.
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

class RAWUserScaleIF
{
    public:
        virtual bool processUserScale( std::vector<float> &src,
                                       std::vector<float> *dst ) = 0;
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

        #define DATATYPE_BYTE               0x00000000
        #define DATATYPE_USHORT             0x00000001
        #define DATATYPE_FLOAT              0x00000002

    public:
        struct WeightAnalysisReport
        {
            unsigned long   timestamp;
            float           base_threshold_index;
            float           threshold_wide_min;
            float           threshold_wide_max;
            unsigned long   threshold_max_amount;
        };

        struct SimpleAnalysisInfo
        {
            float           minLevel;
            float           maxLevel;
            float           average;
            float           variance;
            float           deviation;
        };

    public:
        RAWProcessor();
        RAWProcessor( const char* raw_file, unsigned int height = 0 );
#ifdef WCHAR_SUPPORTED
        RAWProcessor( const wchar_t* raw_file, unsigned int height = 0 );
#endif // WCHAR_SUPPORTED
        virtual~RAWProcessor();

    public:
        bool Loaded()                   { return raw_loaded; }
        unsigned PixelCount()           { return datasize(); }
        float MinimumLevel()            { return pixel_min_level; }
        float MaximumLevel()            { return pixel_max_level; }
        float MediumLevel()             { return pixel_med_level; }
        unsigned Width()                { return img_width; }
        unsigned Height()               { return img_height; }
        unsigned WeightsCount()         { return pixel_weights_max; }
        unsigned char BPP()             { return pixel_bpp; }
        void RecalcLevels()             { calcWeights(); }

    public:
        void CutoffLevels( float minv, float maxv );
        void CutoffLevelsRanged( float minv, float maxv, float valmin = 0, float valmax = 1.f );

    public:
        // Get Versions in string or, integer array.
        void Version( char** retverstr ); /// put NULL initialized char* array.
        void Version( int* retverints ); /// put int[4] array.

        // Load from file.
        bool Load( const char* raw_file, unsigned int trnsfm = TRANSFORM_NONE, size_t height= 0, unsigned int dtype = DATATYPE_USHORT, bool byteswap = false );
#ifdef WCHAR_SUPPORTED
        bool Load( const wchar_t* raw_file, unsigned int trnsfm = TRANSFORM_NONE, size_t height= 0, unsigned int dtype = DATATYPE_USHORT, bool byteswap = false );
#endif // WCHAR_SUPPORTED

        // Load from memory.
        bool LoadFromMemory( void* buffer, size_t bufferlen, unsigned int trnsfm = TRANSFORM_NONE, size_t height = 0, unsigned int dtype = DATATYPE_FLOAT, bool byteswap = false );

        // Reload from file, actually same as like Load.
        bool Reload( const char* raw_file, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
#ifdef WCHAR_SUPPORTED
		bool Reload( const wchar_t* raw_file, unsigned int trnsfm = TRANSFORM_NONE, unsigned height = 0 );
#endif // WCHAR_SUPPORTED
		bool Reload();

		// Unload clears internal memory.
        void Unload();

        // Transform related...
        bool ApplyTransform( unsigned int trnsfm = TRANSFORM_NONE );
        // rotated image may cropped to same size of previous in center.
        // recommended degrees : 0.0 ~ 359.99
        bool RotateFree( float degree );
        void ChangeHeight( unsigned h );
        bool Get8bitDownscaled( std::vector<unsigned char>* byte_arrays, DownscaleType dntype, bool reversed );
        bool Get16bitRawImage( std::vector<unsigned short>& word_arrays, bool reversed );
        void SetUserScale( RAWUserScaleIF* ptr = NULL );
        bool Invert( unsigned char maxbits = 16 );
        bool InvertAuto();

        // GetXXXX methods --
        bool GetDownscaled( std::vector<unsigned char>* byte_arrays, DownscaleType dntype = DNSCALE_NORMAL, bool reversed = false );
        bool GetRawImage( std::vector<float>* word_arrays, bool reversed = false );
        bool GetRawImage( std::vector<float>* word_arrays );
        bool GetAnalysisReport( WeightAnalysisReport& report, bool start_minlevel_zero = false );
        bool GetThresholdedImage( WeightAnalysisReport& report, std::vector<float>* word_arrays, bool reversed = false );
        bool GetThresholdedImage( WeightAnalysisReport& report, std::vector<unsigned char>* byte_arrays, bool reversed = false );
        bool GetThresholdedImage( WeightAnalysisReport& report, std::vector<float>* byte_arrays );
        bool GetPixel( unsigned x, unsigned y, float &px );

    // Some additional tools here ...
    public:
        bool SaveToFile( const char* path );
#ifdef WCHAR_SUPPORTED
        bool SaveToFile( const wchar_t* path );
#endif // WCHAR_SUPPORTED

    public:
        // image may enlarged.
        RAWProcessor* RotateFree( float degree, float background = 0 );
        RAWProcessor* Rescale( unsigned w, unsigned h, RescaleType st = RESCALE_NEAREST );
        RAWProcessor* Clone();

    public:
        typedef struct { unsigned x; unsigned y; } polygoncoord;
        // Weights == Histogram.
        void GetLinearPixels( unsigned x1, unsigned y1, unsigned x2, unsigned y2, std::vector<float>* pixels );
        void GetRectPixels( unsigned x, unsigned y, unsigned w, unsigned h, std::vector<float>* pixels);
        void GetPolygonPixels( std::vector<polygoncoord>* coords, std::vector<float>* pixels);
        void GetAnalysisFromPixels( std::vector<float>& pixels, SimpleAnalysisInfo& info );

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
        #define TONEMAP_TYPE_DRAGO      0   /// Adaptive Logarithmic Mapping for Displaying High Contrast Scenes
        #define TONEMAP_TYPE_REINHARD   1   /// Erik Reinhard and Kate Devlin, 'Dynamic Range Reduction Inspired by Photographer Physiology'

        bool AdjustGamma( float gamma );
        bool AdjustBrightness( float percent );
        bool AdjustContrast( float percent );
        bool AdjustToneMapping( unsigned ttype, float p1, float p2, float p3, float p4 );

    public:
        // CLAHE ( Contrast Limited Adaptive Histogram Equalization )
        bool ApplyCLAHE( WeightAnalysisReport &report, unsigned applysz, unsigned bins, float slope );

    public:
        bool ApplyLowFrequency( unsigned filtersz = 3, unsigned repeat = 1 );
        bool ApplyEdgeEnhance( unsigned fszh = 5, unsigned fszv = 5, unsigned edgesz = 3, unsigned margin = 0 );
        bool ApplyAnisotropicFilter( unsigned strength, unsigned param );

    public:
        const size_t datasize();
        const float* data();

    protected:
        void analyse();
        void resetWeights();
        void calcWeights();
        void findWideness(float& minf, float& maxf );

    protected:
        void addpixelarray( std::vector<float>* outpixels, unsigned x, unsigned y );
        void reordercoords( std::vector<polygoncoord>* coords );

    protected:
        bool               raw_loaded;
        std::vector<float> pixel_arrays;
        unsigned long      pixel_arrays_srcsz;
        unsigned long      pixel_arrays_realsz;
        float              pixel_weights_max;
        unsigned char      pixel_bpp;
        float              pixel_min_level;
        float              pixel_max_level;
        float              pixel_med_level;
        unsigned long      index_max_pixel;
        unsigned int       current_transform;
#ifdef WCHAR_SUPPORTED
        std::wstring       raw_file_name;
#else
        std::string        raw_file_name;
#endif
        unsigned           img_height;
        unsigned           img_width;
        RAWUserScaleIF*    userscaler;
};
#endif /// of __RAWPROCESSOR_H__
