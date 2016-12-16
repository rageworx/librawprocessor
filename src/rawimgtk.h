#ifndef __RAWIMGTOOLKIT_H__
#define __RAWIMGTOOLKIT_H__

namespace RAWImageToolKit
{
    bool FlipHorizontal( unsigned short* ptr, int w, int h );
    bool FlipVertical( unsigned short* ptr, int w, int h );
    bool Rotate90( unsigned short* ptr, int* w, int* h );
    bool Rotate180( unsigned short* ptr, int* w, int* h );
    bool Rotate270( unsigned short* ptr, int* w, int* h );

    bool AdjustCurve( unsigned short* ptr, unsigned arraysz, unsigned short* LUT, unsigned LUTsz );
    bool ApplyGamma( unsigned short* ptr, unsigned arraysz, double gamma );
    bool AdjustBrightness( unsigned short* ptr, unsigned arraysz, double perc );
    bool AdjustContrast( unsigned short* ptr, unsigned arraysz, double perc );
}

#endif /// of __RAWIMGTOOLKIT_H__
