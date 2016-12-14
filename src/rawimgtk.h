#ifndef __RAWIMGTOOLKIT_H__
#define __RAWIMGTOOLKIT_H__

namespace RAWImageToolKit
{
    bool FlipHorizontal( unsigned short* ptr, int w, int h );
    bool FlipVertical( unsigned short* ptr, int w, int h );
    bool Rotate90( unsigned short* ptr, int* w, int* h );
    bool Rotate180( unsigned short* ptr, int* w, int* h );
    bool Rotate270( unsigned short* ptr, int* w, int* h );
}

#endif /// of __RAWIMGTOOLKIT_H__
