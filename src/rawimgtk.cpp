#ifdef DEBUG
#include <cstdio>
#include <cstdlib>
#endif // DEBUG
#include <cstring>
#include <math.h>
#include "rawimgtk.h"

////////////////////////////////////////////////////////////////////////////////

#ifndef MIN
    #define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef MAX
    #define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#define RAWImageToolKitSwapUS( _a_, _b_ )   unsigned short t=_a_; _a_=_b_; _b_=t;

////////////////////////////////////////////////////////////////////////////////

bool RAWImageToolKit::FlipHorizontal( unsigned short* ptr, int w, int h )
{
#ifdef DEBUG
    printf("RAWImageToolKit::FlipHorizontal()\n");
#endif // DEBUG

    if ( ( w > 0 ) && ( h > 0 ) )
    {
        int hcenter = h/2;

        for( int cnth=0; cnth<hcenter; cnth++ )
        {
            for( int cntw=0; cntw<w; cntw++ )
            {
                RAWImageToolKitSwapUS( ptr[ w * ( h - cnth ) + cntw ],
                                       ptr[ w * cnth + cntw ] );
            }
        }

        return true;
    }

    return false;
}

bool RAWImageToolKit::FlipVertical( unsigned short* ptr, int w, int h )
{
#ifdef DEBUG
    printf("RAWImageToolKit::FlipVertical()\n");
#endif // DEBUG
    if ( ( w > 0 ) && ( h > 0 ) )
    {
        int wcenter = w/2;

        for( int cntw=0; cntw<wcenter; cntw++ )
        {
            for( int cnth=0; cnth<h; cnth++ )
            {
                RAWImageToolKitSwapUS( ptr[ w * cnth + ( w - cntw ) ],
                                       ptr[ w * cnth + cntw ] );
            }
        }

        return true;
    }

    return false;
}

bool RAWImageToolKit::Rotate90( unsigned short* ptr, int* w, int* h )
{
#ifdef DEBUG
    printf("RAWImageToolKit::Rotate90()\n");
#endif // DEBUG
    int cur_w = *w;
    int cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        int src_x = 0;
        int src_y = 0;
        int new_w = cur_h;
        int new_h = cur_w;

        unsigned short* tempbuff = new unsigned short[ new_w * new_h ];
        if ( tempbuff == NULL )
            return false;

        memset( tempbuff, 0, new_w * new_h * sizeof( unsigned short ) );

        for( int cntw=new_w-1; cntw>=0; cntw-- )
        {
            for( int cnth=0; cnth<new_h; cnth++ )
            {
                tempbuff[ new_w * cnth + cntw ] = ptr[ cur_w * src_y + src_x ];

                src_x++;
                if ( src_x >= cur_w )
                {
                    src_x = 0;
                    src_y++;
                }
            }
        }

        memcpy( ptr, tempbuff,  new_w * new_h * sizeof( unsigned short ) );

        delete[] tempbuff;

        *w = new_w;
        *h = new_h;

        return true;
    }

    return false;
}

bool RAWImageToolKit::Rotate180( unsigned short* ptr, int* w, int* h )
{
#ifdef DEBUG
    printf("RAWImageToolKit::Rotate180()\n");
#endif // DEBUG
    int cur_w = *w;
    int cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        for( int cntw=0; cntw<cur_w; cntw++ )
        {
            for( int cnth=0; cnth<cur_h; cnth++ )
            {
                RAWImageToolKitSwapUS( ptr[ cur_w * cnth + cntw ],
                                       ptr[ cur_w * ( cur_h - cnth ) + ( cur_w - cntw ) ] );
            }
        }

        return true;
    }

    return false;
}

bool RAWImageToolKit::Rotate270( unsigned short* ptr, int* w, int* h )
{
#ifdef DEBUG
    printf("RAWImageToolKit::Rotate270()\n");
#endif // DEBUG
    int cur_w = *w;
    int cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        int src_x = 0;
        int src_y = 0;
        int new_w = cur_h;
        int new_h = cur_w;

        unsigned short* tempbuff = new unsigned short[ new_w * new_h ];
        if ( tempbuff == NULL )
            return false;

        memset( tempbuff, 0,  new_w * new_h * sizeof( unsigned short ) );

        for( int cntw=0; cntw<new_w; cntw++ )
        {
            for( int cnth=new_h-1; cnth>=0; cnth-- )
            {
                tempbuff[ new_w * cnth + cntw ] = ptr[ cur_w * src_y + src_x ];

                src_x++;
                if ( src_x >= cur_w )
                {
                    src_x = 0;
                    src_y++;
                }
            }
        }

        memcpy( ptr, tempbuff,  new_w * new_h * sizeof( unsigned short ) );

        delete[] tempbuff;

        *w = new_w;
        *h = new_h;

        return true;
    }

    return false;
}

bool RAWImageToolKit::AdjustCurve( unsigned short* ptr, unsigned arraysz, unsigned short* LUT, unsigned LUTsz )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) || ( LUT == NULL ) || ( LUTsz < 65535 ) )
        return false;

    for( unsigned cnt=0; cnt<arraysz; cnt++ )
    {
        ptr[ cnt ] = LUT[ ptr[cnt] ];
    }

    return true;
}

bool RAWImageToolKit::ApplyGamma( unsigned short* ptr, unsigned arraysz, double gamma )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    unsigned short LUT[65536] = {0};

    double expn = 1.0 / gamma;
    double newv = 65535.0 * (double)pow((double)65535, -expn );

    for( unsigned cnt=0; cnt<arraysz; cnt++ )
    {
        double newlvl = (double)pow((double)cnt, expn) * newv;
        if( newlvl > 65535.0 )
        {
            newlvl = 65535;
        }
        LUT[ cnt ] = (unsigned short)floor( newlvl + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT, 65535 );
}

bool RAWImageToolKit::AdjustBrightness( unsigned short* ptr, unsigned arraysz, double perc )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    unsigned short LUT[65536] = {0};
    const double bscaled = ( 100.0 + perc ) / 100.0;
    double bval = 0.0;

    for( unsigned cnt=0; cnt<65536; cnt++ )
    {
        bval = (double)cnt * bscaled;
        bval = MAX( 0.0, MIN( bval, 65535.0 ) );
        LUT[ cnt ] = (unsigned short)floor( bval + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT, 65535 );
}

bool RAWImageToolKit::AdjustContrast( unsigned short* ptr, unsigned arraysz, double perc )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    unsigned short LUT[65536] = {0};
    const double bscaled = ( 100.0 + perc ) / 100.0;
    double bval = 0.0;

    for( unsigned cnt=0; cnt<65536; cnt++ )
    {
        bval = (double)32768.0 + ( (double)cnt - 32768 ) * bscaled;
        bval = MAX( 0.0, MIN( bval, 65535.0 ) );
        LUT[ cnt ] = (unsigned short)floor( bval + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT, 65535 );
}
