#ifdef DEBUG
#include <cstdio>
#include <cstdlib>
#endif // DEBUG
#include <cstring>
#include <math.h>
#include "rawimgtk.h"

////////////////////////////////////////////////////////////////////////////////

#define RAWIMGTK_MAX_D_ARRAY_SZ     65536
#define RAWIMGTK_MAX_F_VAL          65535.0
#define RAWIMGTK_HALF_F_VAL         32768.0

#ifndef MIN
    #define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef MAX
    #define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#define RAWImageToolKitSwapUS( _a_, _b_ )   unsigned short t=_a_; _a_=_b_; _b_=t;

////////////////////////////////////////////////////////////////////////////////

bool RAWImageToolKit::FlipHorizontal( unsigned short* ptr, unsigned w, unsigned h )
{
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

bool RAWImageToolKit::FlipVertical( unsigned short* ptr, unsigned w, unsigned h )
{
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

bool RAWImageToolKit::Rotate90( unsigned short* ptr, unsigned* w, unsigned* h )
{
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

bool RAWImageToolKit::Rotate180( unsigned short* ptr, unsigned* w, unsigned* h )
{
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

bool RAWImageToolKit::Rotate270( unsigned short* ptr, unsigned* w, unsigned* h )
{
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

bool RAWImageToolKit::AdjustGamma( unsigned short* ptr, unsigned arraysz, double gamma )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    unsigned short LUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};
    double expn     = 1.0 / gamma;
    double newv     = RAWIMGTK_MAX_F_VAL * (double)pow( (double)RAWIMGTK_MAX_F_VAL, -expn );
    double newlvl   = 0.0;

    for( unsigned cnt=0; cnt<arraysz; cnt++ )
    {
        newlvl = (double)pow((double)cnt, expn) * newv;
        if( newlvl > RAWIMGTK_MAX_F_VAL )
        {
            newlvl = RAWIMGTK_MAX_F_VAL;
        }
        LUT[ cnt ] = (unsigned short)floor( newlvl + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT );
}

bool RAWImageToolKit::AdjustBrightness( unsigned short* ptr, unsigned arraysz, double perc )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    unsigned short LUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};
    const \
    double bscaled  = ( 100.0 + perc ) / 100.0;
    double bval     = 0.0;

    for( unsigned cnt=0; cnt<RAWIMGTK_MAX_D_ARRAY_SZ; cnt++ )
    {
        bval = (double)cnt * bscaled;
        bval = MAX( 0.0, MIN( bval, RAWIMGTK_MAX_F_VAL ) );
        LUT[ cnt ] = (unsigned short)floor( bval + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT );
}

bool RAWImageToolKit::AdjustContrast( unsigned short* ptr, unsigned arraysz, double perc )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    unsigned short LUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};
    const \
    double bscaled  = ( 100.0 + perc ) / 100.0;
    double bval     = 0.0;

    for( unsigned cnt=0; cnt<RAWIMGTK_MAX_D_ARRAY_SZ; cnt++ )
    {
        bval = (double)RAWIMGTK_HALF_F_VAL + ( (double)cnt - RAWIMGTK_HALF_F_VAL ) * bscaled;
        bval = MAX( 0.0, MIN( bval, RAWIMGTK_MAX_F_VAL ) );
        LUT[ cnt ] = (unsigned short)floor( bval + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT );
}

bool RAWImageToolKit::AdjustCurve( unsigned short* ptr, unsigned arraysz, unsigned short* LUT )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) || ( LUT == NULL ) )
        return false;

    for( unsigned cnt=0; cnt<arraysz; cnt++ )
    {
        ptr[ cnt ] = LUT[ ptr[ cnt ] ];
    }

    return true;
}
