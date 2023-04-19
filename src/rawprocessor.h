#ifndef __RAWPROCESSOR_H__
#define __RAWPROCESSOR_H__
////////////////////////////////////////////////////////////////////////////////
//
//  32bit floating point precise gray scale RAW image processor for stdc++.
//  ============================================================================
//  (C)Copyright 2013-2023, Raphael Kim (rageworx@gmail.com)
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
#include <cstdint>

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

        #define TRANSFORM_NONE              ( 0x00000000 )

        #define TRANSFORM_SWAP              ( 0x00000010 )
        #define TRANSFORM_PARAM_SWAP_H      ( 0x00000002 )
        #define TRANSFORM_PARAM_SWAP_V      ( 0x00000004 )

        #define TRANSFORM_ROTATE            ( 0x00010000 )
        #define TRANSFORM_PARAM_ROT_C90     ( 0x00000200 )
        #define TRANSFORM_PARAM_ROT_C180    ( 0x00000400 )
        #define TRANSFORM_PARAM_ROT_C270    ( 0x00000800 )
        #define TRANSFORM_PARAM_ROT_RC90    ( 0x00002000 )
        #define TRANSFORM_PARAM_ROT_RC180   ( 0x00004000 )
        #define TRANSFORM_PARAM_ROT_RC270   ( 0x00008000 )

        #define DATATYPE_BYTE               ( 0x00000000 )
        #define DATATYPE_USHORT             ( 0x00000001 )
        #define DATATYPE_FLOAT              ( 0x00000002 )

    public:
        struct WindowAnalysisReport
        {
            uint32_t    timestamp;
            float       base_threshold_index;
            float       threshold_wide_min;
            float       threshold_wide_max;
            uint32_t    threshold_max_amount;
        };

        struct SimpleAnalysisInfo
        {
            float   minLevel;
            float   maxLevel;
            float   average;
            float   variance;
            float   deviation;
        };

    public:
        RAWProcessor();
        RAWProcessor( const char* raw_file, uint32_t height = 0 );
#ifdef WCHAR_SUPPORTED
        RAWProcessor( const wchar_t* raw_file, uint32_t height = 0 );
#endif // WCHAR_SUPPORTED
        virtual~RAWProcessor();

    public:
        bool Loaded()                   { return raw_loaded; }
        size_t PixelCount()             { return datasize(); }
        float MinimumLevel()            { return pixel_min_level; }
        float MaximumLevel()            { return pixel_max_level; }
        float MediumLevel()             { return pixel_med_level; }
        uint32_t Width()                { return img_width; }
        uint32_t Height()               { return img_height; }
        uint32_t WindowCount()          { return pixel_window_max; }
        uint32_t BPP()                  { return pixel_bpp; }
        void RecalcLevels()             { calcWindow(); }

    public:
        void CutoffLevels( float minv, float maxv );
        void CutoffLevelsRanged( float minv, float maxv, float valmin = 0, float valmax = 1.f );

    public:
        // Get Versions in string or, integer array.
        void Version( char** retverstr ); /// put NULL initialized char* array.
        void Version( int* retverints ); /// put int[4] array.

        // Load from file.
        bool Load( const char* raw_file, uint32_t trnsfm = TRANSFORM_NONE, size_t height= 0, uint32_t dtype = DATATYPE_USHORT, bool byteswap = false );
#ifdef WCHAR_SUPPORTED
        bool Load( const wchar_t* raw_file, uint32_t trnsfm = TRANSFORM_NONE, size_t height= 0, uint32_t dtype = DATATYPE_USHORT, bool byteswap = false );
#endif // WCHAR_SUPPORTED

        // Load from memory.
        bool LoadFromMemory( void* buffer, size_t bufferlen, uint32_t trnsfm = TRANSFORM_NONE, size_t height = 0, uint32_t dtype = DATATYPE_FLOAT, bool byteswap = false );

        // Reload from file, actually same as like Load.
        bool Reload( const char* raw_file, uint32_t trnsfm = TRANSFORM_NONE, uint32_t height = 0 );
#ifdef WCHAR_SUPPORTED
        bool Reload( const wchar_t* raw_file, uint32_t trnsfm = TRANSFORM_NONE, uint32_t height = 0 );
#endif // WCHAR_SUPPORTED
        bool Reload();

        // Unload clears internal memory.
        void Unload();

        // Transform related...
        bool ApplyTransform( uint32_t trnsfm = TRANSFORM_NONE );
        // rotated image may cropped to same size of previous in center.
        // recommended degrees : 0.0 ~ 359.99
        bool RotateFree( float degree );
        void ChangeHeight( uint32_t h );
        void SetUserScale( RAWUserScaleIF* ptr = NULL );
        bool Invert(); /// same as InvertAuto() in float version.
        bool InvertAuto();

        // GetXXXX methods --
        bool Get8bitDownscaled( std::vector<uint8_t>* byte_arrays, DownscaleType dntype, bool reversed );
        bool Get16bitRawImage( std::vector<uint16_t>& word_arrays, bool reversed );
        bool GetDownscaled( std::vector<uint8_t>* byte_arrays, DownscaleType dntype = DNSCALE_NORMAL, bool reversed = false );
        bool GetRawImage( std::vector<float>* f_arrays, bool reversed = false );
        bool GetRawImage( std::vector<float>* f_arrays );
        bool GetAnalysisReport( WindowAnalysisReport& report, bool start_minlevel_zero = false );
        bool GetThresholdedImage( WindowAnalysisReport& report, std::vector<uint32_t>* d_arrays, bool reversed = false );
        bool GetThresholdedImage( WindowAnalysisReport& report, std::vector<uint16_t>* w_arrays, bool reversed = false );
        bool GetThresholdedImage( WindowAnalysisReport& report, std::vector<uint8_t>* b_arrays, bool reversed = false );
        bool GetPixel( uint32_t x, uint32_t y, float &px );
        bool GetHistography( std::vector<float>* f_histo );

    // Some additional tools here ...
    public:
        bool SaveToFile( const char* path );
#ifdef WCHAR_SUPPORTED
        bool SaveToFile( const wchar_t* path );
#endif // WCHAR_SUPPORTED

    public:
        // image may enlarged.
        RAWProcessor* RotateFree( float degree, float background = 0 );
        RAWProcessor* Rescale( uint32_t w, uint32_t h, RescaleType st = RESCALE_NEAREST );
        RAWProcessor* Clone();

    public:
        typedef struct { uint32_t x; uint32_t y; } polygoncoord;
        void GetLinearPixels( uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2, std::vector<float>* pixels );
        void GetRectPixels( uint32_t x, uint32_t y, uint32_t w, uint32_t h, std::vector<float>* pixels);
        void GetPolygonPixels( std::vector<polygoncoord>* coords, std::vector<float>* pixels);
        void GetAnalysisFromPixels( std::vector<float>& pixels, SimpleAnalysisInfo& info );

    public:
        typedef struct
        {
            uint32_t width;
            uint32_t height;
            float    factor;
            float    bias;
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

        FilterConfig* GetPresetFilter( uint32_t fnum );
        void          DiscardFilter( FilterConfig* fp );

    public:
        #define TONEMAP_TYPE_DRAGO      0   /// Adaptive Logarithmic Mapping for Displaying High Contrast Scenes
        #define TONEMAP_TYPE_REINHARD   1   /// Erik Reinhard and Kate Devlin, 'Dynamic Range Reduction Inspired by Photographer Physiology'

        bool AdjustGamma( float gamma );
        bool AdjustBrightness( float percent );
        bool AdjustContrast( float percent );
        bool AdjustToneMapping( uint32_t ttype, float p1, float p2, float p3, float p4 );

    public:
        // CLAHE ( Contrast Limited Adaptive Histogram Equalization )
        bool ApplyCLAHE( WindowAnalysisReport &report, uint32_t applysz, uint32_t bins, float slope );

    public:
        bool ApplyLowFrequency( uint32_t filtersz = 3, uint32_t repeat = 1 );
        bool ApplyEdgeEnhance( uint32_t fszh = 5, uint32_t fszv = 5, uint32_t edgesz = 3, uint32_t margin = 0 );
        bool ApplyAnisotropicFilter( uint32_t strength, uint32_t param );

    public:
        const size_t datasize();
        const float* data();

    protected:
        void analyse();
        void resetWindow();
        void calcWindow();
        void findWideness(float& minf, float& maxf );

    protected:
        void addpixelarray( std::vector<float>* outpixels, uint32_t x, uint32_t y );
        void reordercoords( std::vector<polygoncoord>* coords );
        bool f2dthldexport( uint8_t tp = 0, void* pd = NULL, bool rvs = false, WindowAnalysisReport* report = NULL );
        bool f2dexport( uint8_t tp = 0, void* pd = NULL, bool rvs = false );

    protected:
        bool               raw_loaded;
        std::vector<float> pixel_arrays;
        size_t             pixel_arrays_srcsz;
        size_t             pixel_arrays_realsz;
        float              pixel_window_max;
        uint32_t           pixel_bpp;
        float              pixel_min_level;
        float              pixel_max_level;
        float              pixel_med_level;
        size_t             index_max_pixel;
        uint32_t           current_transform;
#ifdef WCHAR_SUPPORTED
        std::wstring       raw_file_name;
#else
        std::string        raw_file_name;
#endif
        uint32_t           img_height;
        uint32_t           img_width;
        RAWUserScaleIF*    userscaler;
};
#endif /// of __RAWPROCESSOR_H__
