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
    * ToneMapping reprogrammed from FreeImage 3 library, ToneMapping.cpp
    * for gray scaled image.
    *
    * tonemaptypes : 0 = Adaptive logarithmic mapping (F. Drago, 2003)
    *                1 = Dynamic range reduction inspired by
    *                    photoreceptor phhysiology (E. Reinhard, 2005)
    ***/
    bool ToneMapping( unsigned short* ptr, unsigned sz, unsigned maxlvl, unsigned tmtype, float p1, float p2, float p3, float p4 );
}

#endif /// of __RAWIMGTOOLKIT_H__
