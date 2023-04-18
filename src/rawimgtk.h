#ifndef __RAWIMGTOOLKIT_H__
#define __RAWIMGTOOLKIT_H__

////////////////////////////////////////////////////////////////////////////////

#include <cstdint>

////////////////////////////////////////////////////////////////////////////////

namespace RAWImageToolKit
{
    bool FlipHorizontal( float* ptr,  uint32_t w,  uint32_t h );
    bool FlipVertical( float* ptr,  uint32_t w,  uint32_t h );
    bool Rotate90( float* ptr, uint32_t* w, uint32_t* h );
    bool Rotate180( float* ptr, uint32_t* w, uint32_t* h );
    bool Rotate270( float* ptr, uint32_t* w, uint32_t* h );

    bool RotateFree( float* ptr, uint32_t* w, uint32_t* h, \
                     float** newptr, float degree,         \
                     float background );

    bool Crop( float* ptr,  uint32_t w,  uint32_t h,       \
               float** newptr,                             \
               uint32_t cx,  uint32_t cy,  uint32_t cw,  uint32_t ch );

    bool CropCenter( float* ptr,  uint32_t w,  uint32_t h, \
                     float** newptr,  uint32_t cw,  uint32_t ch );

    /***
    *
    * these following methods: AdjustGamma, AdjustBrightness, AdjustContrast
    *
    * are motivated from FreeImageToolkit open source.
    * --  http://freeimage.sourceforge.net/
    * ,by Following FreeImage License.
    *
    * reprogrammed by Raph.K. (rageworx@gmail.com) for librawprocessor.
    *
    * -- parameter notice ------------------------------------------------------
    * ** gamma == 0.0 to 0.1f.
    * ** perc == percent, 0.0 to 100.0f.
    ***/
    bool AdjustGamma( float* ptr, size_t arraysz, float gamma );
    bool AdjustBrightness( float* ptr, size_t arraysz, float perc );
    bool AdjustContrast( float* ptr, size_t arraysz, float perc );

    /***
    * CLAHE ( Contrast Limited Adaptive Histogram Equalization )
    * ----------------------------------------------------------
    * Motivated on
    * https://en.wikipedia.org/wiki/Adaptive_histogram_equalization#Contrast_Limited_AHE
    *
    * Source refered to
    * Karel Zuiderveld, karel@cv.ruu.nl of "Contrast Limited Adaptive Histogram Equalization"
    * http://www.realtimerendering.com/resources/GraphicsGems/gemsiv/clahe.c
    ***/
    bool ApplyCLAHE( float* ptr,  uint32_t w,  uint32_t h,
                     float minlvl, float maxlvl,
                     uint32_t blkw,  uint32_t blkh, uint32_t bins, float slope );
}

#endif /// of __RAWIMGTOOLKIT_H__
