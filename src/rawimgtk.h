#ifndef __RAWIMGTOOLKIT_H__
#define __RAWIMGTOOLKIT_H__

////////////////////////////////////////////////////////////////////////////////

namespace RAWImageToolKit
{
    bool FlipHorizontal( float* ptr, unsigned w, unsigned h );
    bool FlipVertical( float* ptr, unsigned w, unsigned h );
    bool Rotate90( float* ptr, unsigned* w, unsigned* h );
    bool Rotate180( float* ptr, unsigned* w, unsigned* h );
    bool Rotate270( float* ptr, unsigned* w, unsigned* h );

    bool RotateFree( float* ptr, unsigned* w, unsigned* h,
                     float** newptr, float degree,
                     float background );

    bool Crop( float* ptr, unsigned w, unsigned h,
               float** newptr,
               unsigned cx, unsigned cy, unsigned cw, unsigned ch );

    bool CropCenter( float* ptr, unsigned w, unsigned h,
                     float** newptr, unsigned cw, unsigned ch );


    /***
    *
    * these following methods: AdjustGamma, AdjustBrightness, AdjustContrast,
    *                          AdjustContrast, AdjustCurve
    *
    * are motivated from FreeImageToolkit open source.
    * --  http://freeimage.sourceforge.net/
    * ,by Following FreeImage License.
    *
    * reprogrammed by Raph.K. (rageworx@gmail.com) for librawprocessor.
    *
    ***/
    bool AdjustGamma( float* ptr, unsigned arraysz, double gamma );
    bool AdjustBrightness( float* ptr, unsigned arraysz, double perc );
    bool AdjustContrast( float* ptr, unsigned arraysz, double perc );

    bool AdjustCurve( float* ptr, unsigned arraysz, float* LUT );

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
    bool ApplyCLAHE( float* ptr, unsigned w, unsigned h,
                     float minlvl, float maxlvl,
                     unsigned blkw, unsigned blkh, unsigned bins, float slope );
}

#endif /// of __RAWIMGTOOLKIT_H__
