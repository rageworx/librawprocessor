#ifndef __RAWIMGTOOLKIT_H__
#define __RAWIMGTOOLKIT_H__

namespace RAWImageToolKit
{
    bool FlipHorizontal( unsigned short* ptr, unsigned w, unsigned h );
    bool FlipVertical( unsigned short* ptr, unsigned w, unsigned h );
    bool Rotate90( unsigned short* ptr, unsigned* w, unsigned* h );
    bool Rotate180( unsigned short* ptr, unsigned* w, unsigned* h );
    bool Rotate270( unsigned short* ptr, unsigned* w, unsigned* h );

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
    bool AdjustGamma( unsigned short* ptr, unsigned arraysz, double gamma );
    bool AdjustBrightness( unsigned short* ptr, unsigned arraysz, double perc );
    bool AdjustContrast( unsigned short* ptr, unsigned arraysz, double perc );

    bool AdjustCurve( unsigned short* ptr, unsigned arraysz, unsigned short* LUT );

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
    bool ApplyCLAHE( unsigned short* ptr, unsigned w, unsigned h,
                     unsigned short minlvl, unsigned short maxlvl,
                     unsigned blkw, unsigned blkh, unsigned bins, float slope );
}

#endif /// of __RAWIMGTOOLKIT_H__
