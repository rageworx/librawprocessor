#ifndef __RAWIMGTOOLKIT_H__
#define __RAWIMGTOOLKIT_H__

namespace RAWImageToolKit
{
    bool FlipHorizontal( unsigned short* ptr, int w, int h );
    bool FlipVertical( unsigned short* ptr, int w, int h );
    bool Rotate90( unsigned short* ptr, int* w, int* h );
    bool Rotate180( unsigned short* ptr, int* w, int* h );
    bool Rotate270( unsigned short* ptr, int* w, int* h );

    /***
    *
    * these following methods: ApplyGamma, AdjustBrightness, AdjustContrast,
    *                          AdjustContrast, AdjustCurve
    *
    * are motivated from FreeImageToolkit open source.
    * --  http://freeimage.sourceforge.net/
    * ,by Following FreeImage License.
    *
    * reprogrammed by Raph.K. (rageworx@gmail.com) for librawprocessor.
    *
    ***/
    bool ApplyGamma( unsigned short* ptr, unsigned arraysz, double gamma );
    bool AdjustBrightness( unsigned short* ptr, unsigned arraysz, double perc );
    bool AdjustContrast( unsigned short* ptr, unsigned arraysz, double perc );

    bool AdjustCurve( unsigned short* ptr, unsigned arraysz, unsigned short* LUT );
}

#endif /// of __RAWIMGTOOLKIT_H__
