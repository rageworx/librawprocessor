#ifdef DEBUG
#include <cstdio>
#include <cstdlib>
#endif // DEBUG
#include <cstring>
#include <cmath>

#ifdef USE_OMP
#include <omp.h>
#endif // USE_OMP

#include "rawimgtk.h"
#include "minmax.h"
#include "clahe.h"

////////////////////////////////////////////////////////////////////////////////

#define RAWIMGTK_MAX_D_ARRAY_SZ     65536
#define RAWIMGTK_MAX_F_VAL          65535.0
#define RAWIMGTK_HALF_F_VAL         32768.0

#define CLAHE_MAX_REG_W             16
#define CLAHE_MAX_REG_H             16

#define RAWImageToolKitSwapUS( _a_, _b_ )   unsigned short t=_a_; _a_=_b_; _b_=t;

////////////////////////////////////////////////////////////////////////////////

bool RAWImageToolKit::FlipHorizontal( unsigned short* ptr, unsigned w, unsigned h )
{
    if ( ( w > 0 ) && ( h > 0 ) )
    {
        unsigned hcenter = h/2;
        unsigned cnth = 0;
        unsigned cntw = 0;

        #pragma omp parallel for private(cntw)
        for( cnth=0; cnth<hcenter; cnth++ )
        {
            for( cntw=0; cntw<w; cntw++ )
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
        unsigned wcenter = w/2;
        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<wcenter; cntw++ )
        {
            for( cnth=0; cnth<h; cnth++ )
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
    unsigned cur_w = *w;
    unsigned cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        unsigned src_x = 0;
        unsigned src_y = 0;
        unsigned new_w = cur_h;
        unsigned new_h = cur_w;

        unsigned short* tempbuff = new unsigned short[ new_w * new_h ];

        if ( tempbuff == NULL )
            return false;

        memset( tempbuff, 0, new_w * new_h * sizeof( unsigned short ) );

        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=new_w-1; cntw>=0; cntw-- )
        {
            for( cnth=0; cnth<new_h; cnth++ )
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
    unsigned cur_w = *w;
    unsigned cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<cur_w; cntw++ )
        {
            for( cnth=0; cnth<cur_h; cnth++ )
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
    unsigned cur_w = *w;
    unsigned cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        unsigned src_x = 0;
        unsigned src_y = 0;
        unsigned new_w = cur_h;
        unsigned new_h = cur_w;

        unsigned short* tempbuff = new unsigned short[ new_w * new_h ];

        if ( tempbuff == NULL )
            return false;

        memset( tempbuff, 0,  new_w * new_h * sizeof( unsigned short ) );

        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<new_w; cntw++ )
        {
            for( cnth=new_h-1; cnth>=0; cnth-- )
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

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<RAWIMGTK_MAX_D_ARRAY_SZ; cnt++ )
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

    #pragma omp parallel for
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

    #pragma omp parallel for
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

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<arraysz; cnt++ )
    {
        ptr[ cnt ] = LUT[ ptr[ cnt ] ];
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////
// Some CLAHE depends functions here:

#ifdef CLAHE_NATIVE_FEATURE

void CLAHE_ClipHistogram( unsigned long* pHisto, unsigned greylevels, unsigned cliplimit)
/* This function performs clipping of the histogram and redistribution of bins.
 * The histogram is clipped and the number of excess pixels is counted. Afterwards
 * the excess pixels are equally redistributed across the whole histogram (providing
 * the bin count is smaller than the cliplimit).
 */
{
    unsigned long* pBinPtr = NULL;
    unsigned long* pEndPtr = NULL;
    unsigned long* pHistoPtr = NULL;

    unsigned long nrexcess = 0;
    unsigned long upper;
    unsigned long bininc;
    unsigned long stepsz;

    long     binexcess;

    pBinPtr = pHisto;

    for ( unsigned cnt=0; cnt<greylevels; cnt++ )
    { /* calculate total number of excess pixels */
        binexcess = (long)pBinPtr[cnt] - (long)cliplimit;

        if ( binexcess > 0 )
        {
            /* excess in current bin */
            nrexcess += binexcess;
        }
    }

    /* Second part: clip histogram and redistribute excess pixels in each bin */
    bininc = nrexcess / greylevels;     /// average binincrement
    upper  = cliplimit - bininc;	    /// Bins larger than upper set to cliplimit

    for ( unsigned cnt=0; cnt<greylevels; cnt++ )
    {
        if ( pHisto[cnt] > cliplimit )
        {
            pHisto[cnt] = cliplimit; /* clip bin */
        }
        else
        {
            if ( pHisto[cnt] > upper )
            {
                /* high bin count */
                nrexcess   -= pHisto[cnt] - upper;
                pHisto[cnt] = cliplimit;
            }
            else
            {
                /* low bin count */
                nrexcess    -= bininc;
                pHisto[cnt] += bininc;
            }
        }
    }

    while ( nrexcess > 0 )
    {   /* Redistribute remaining excess  */
        pEndPtr   = &pHisto[greylevels];
        pHistoPtr = pHisto;

        while ( ( nrexcess > 0 ) && ( pHistoPtr < pEndPtr ) )
        {
            stepsz = greylevels / nrexcess;

            if ( stepsz < 1 )
                stepsz = 1;		  /* stepsize at least 1 */

            for ( pBinPtr = pHistoPtr;
                  ( pBinPtr < pEndPtr ) && ( nrexcess > 0 );
                  pBinPtr += stepsz )
            {
                if ( *pBinPtr < cliplimit )
                {
                    (*pBinPtr)++;
                    /* reduce excess */
                    nrexcess--;
                }
            }

            /* restart redistributing on other bin location */
            pHistoPtr++;
        }
    }
}

void CLAHE_MakeHistogram ( unsigned short* pImage,
                           unsigned wsz,
                           unsigned sizex, unsigned sizey,
                           unsigned long* pHisto,
                           unsigned greylevels, unsigned short* pLUT )
/* This function classifies the greylevels present in the array image into
 * a greylevel histogram. The pLookupTable specifies the relationship
 * between the greyvalue of the pixel (typically between 0 and 4095) and
 * the corresponding bin in the histogram (usually containing only 128 bins).
 */
{
    unsigned short* pImagePointer = NULL;

    memset( pHisto, 0, greylevels * sizeof( unsigned ) );

    for ( unsigned cnt=0; cnt<sizey; cnt++ )
    {
        pImagePointer = &pImage[ sizex ];

        while ( pImage < pImagePointer )
        {
            pHisto[ pLUT[ *pImage++ ] ]++;
        }

        pImagePointer += wsz;

        pImage = &pImagePointer[-sizex];
    }
}

void CLAHE_MapHistogram ( unsigned long* pHisto,
                          unsigned short minlvl, unsigned short maxlvl,
                          unsigned greylevels, unsigned sz )
/* This function calculates the equalized lookup table (mapping) by
 * cumulating the input histogram. Note: lookup table is rescaled in range [Min..Max].
 */
{
    const float fScale       = ((float)(maxlvl - minlvl)) / sz;
    const unsigned long uMin = (unsigned long) minlvl;

    unsigned long sum = 0;

    for ( unsigned cnt=0; cnt<greylevels; cnt++ )
    {
        sum        += pHisto[cnt];
        pHisto[cnt] = (unsigned long)( uMin + sum * fScale );

        if ( pHisto[ cnt ] > maxlvl )
        {
            pHisto[ cnt ] = maxlvl;
        }
    }
}

void CLAHE_MakeLut( unsigned short* pLUT, unsigned short minlvl, unsigned short maxlvl, unsigned int bins )
/* To speed up histogram clipping, the input image [Min,Max] is scaled down to
 * [0,uiNrBins-1]. This function calculates the LUT.
 */
{
    const unsigned long binsz = 1 + ( maxlvl - minlvl ) / bins;

    for ( unsigned cnt=minlvl; cnt<=maxlvl; cnt++ )
    {
        pLUT[ cnt ] = ( cnt - minlvl ) / binsz;
    }
}

/***
 * ptr         - pointer to input/output image
 * wsz         - resolution of image in x-direction
 * mapXXX      - mappings of grey levels from histograms
 * w           - Width of image submatrix
 * h           - Height of image submatrix
 * pLUT	       - lookup table containing mapping greyvalues to bins
 * This function calculates the new greylevel assignments of pixels within a submatrix
 * of the image with size Width and Height. This is done by a bilinear interpolation
 * between four different mappings in order to eliminate boundary artifacts.
 * It uses a division; since division is often an expensive operation, I added code to
 * perform a logical shift instead when feasible.
***/
void CLAHE_Interpolate( unsigned short* pImg, unsigned wsz,
                        unsigned long* mapLU, unsigned long* mapRU,
                        unsigned long* mapLB, unsigned long* mapRB,
                        unsigned subw, unsigned subh, unsigned short* pLUT )
{
    /* Pointer increment after processing row */
    const unsigned incVal = wsz - subw;

    unsigned short greyVal = 0;

    /* Normalization factor */
    unsigned int counts   = subw * subh;

    unsigned int XCoef    = 0;
    unsigned int YCoef    = 0;
    unsigned int XInvCoef = 0;
    unsigned int YInvCoef = 0;
    unsigned int Shifts   = 0;

    /* If Counts is not a power of two, use division */
    if ( ( counts > 0 ) & ( ( counts - 1 ) > 0 ) )
    {
        for ( YCoef = 0, YInvCoef = subh;
              YCoef < subh;
              YCoef++, YInvCoef--, pImg+=incVal )
        {
            for ( XCoef = 0, XInvCoef = subw;
                  XCoef < subw;
                  XCoef++, XInvCoef-- )
            {
                /* get histogram bin value */
                greyVal = pLUT[ *pImg ];

                *pImg++ = (unsigned short)
                          ( YInvCoef * (XInvCoef*mapLU[greyVal]
                            + XCoef * mapRU[greyVal])
                            + YCoef * (XInvCoef * mapLB[greyVal]
                            + XCoef * mapRB[greyVal]) ) / counts;
            }
        }
    }
    else
    {   /* avoid the division and use a right shift instead */
        while ( counts >>= 1 )
        {
            /* Calculate 2log of Counts */
            Shifts++;
        }

        for ( YCoef = 0, YInvCoef = subh;
              YCoef < subh;
              YCoef++, YInvCoef--, pImg+=incVal)
        {
            for ( XCoef = 0, XInvCoef = subw;
                  XCoef < subw;
                  XCoef++, XInvCoef-- )
            {
                greyVal = pLUT[*pImg];	  /* get histogram bin value */
                *pImg++ = (unsigned short)
                          ( YInvCoef* (XInvCoef * mapLU[greyVal]
                            + XCoef * mapRU[greyVal])
                            + YCoef * (XInvCoef * mapLB[greyVal]
                            + XCoef * mapRB[greyVal]) ) >> Shifts;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

/*   ptr    - Pointer to the input/output image
 *   w      - Image resolution in the X direction
 *   h      - Image resolution in the Y direction
 *   minlvl - Minimum greyvalue of input image (also becomes minimum of output image)
 *   maxlvl - Maximum greyvalue of input image (also becomes maximum of output image)
 *   blkw   - Number of contextial regions in the X direction (min 2, max uiMAX_REG_X)
 *   blkh   - Number of contextial regions in the Y direction (min 2, max uiMAX_REG_Y)
 *   bins   - Number of greybins for histogram ("dynamic range")
 *   slope  - Normalized cliplimit (higher values give more contrast)
 * The number of "effective" greylevels in the output image is set by uiNrBins; selecting
 * a small value (eg. 128) speeds up processing and still produce an output image of
 * good quality. The output image will have the same minimum and maximum value as the input
 * image. A clip limit smaller than 1 results in standard (non-contrast limited) AHE.
 */

bool RAWImageToolKit::ApplyCLAHE( unsigned short* ptr, unsigned w, unsigned h,
                                  unsigned short minlvl, unsigned short maxlvl,
                                  unsigned blkw, unsigned blkh, unsigned bins, float slope )
{
    /* counters */
    unsigned cntX;
    unsigned cntY;

    /* size of context. reg. and subimages */
    unsigned subW = 0;
    unsigned subH = 0;
    unsigned subX = 0;
    unsigned subY = 0;

    /* auxiliary variables interpolation routine */
    unsigned inpolXL;
    unsigned inpolXR;
    unsigned inpolYU;
    unsigned inpolYB;

    /* clip limit and region pixel count */
    unsigned long cliplimit;
    unsigned long nrpixels;

    /* pointer to image */
    unsigned short* pImPointer = NULL;

    /* lookup table used for scaling of input image */
    unsigned short aLUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};

    /* pointer to histogram and mappings*/
    unsigned long* pulHist = NULL;
    unsigned long* pulMapArray = NULL;

    /* auxiliary pointers interpolation */
    unsigned long* pulLU = NULL;
    unsigned long* pulLB = NULL;
    unsigned long* pulRU = NULL;
    unsigned long* pulRB = NULL;

    if ( blkw > CLAHE_MAX_REG_W )
        blkw = CLAHE_MAX_REG_W;

    if ( blkh > CLAHE_MAX_REG_H )
        blkh = CLAHE_MAX_REG_H;

    if ( w % blkw )
        return false;   /// x-resolution no multiple of blkw

    if ( h & blkh )
        return false;  /// y-resolution no multiple of blkh

    if ( maxlvl >= RAWIMGTK_MAX_D_ARRAY_SZ )
        return false;   /// maximum too large

    if ( minlvl >= maxlvl )
        return false;   /// minimum equal or larger than maximum

    if ( ( blkw< 2 ) || ( blkh< 2 ) )
        return false;   /// at least 4 contextual regions required

    if ( slope == 1.0f )
        return true;    /// is OK, immediately returns original image.

    if ( bins == 0 )
        bins = 128;     /// default value when not specified

    pulMapArray = new unsigned long[ blkw * blkh * bins ];

    if ( pulMapArray == NULL )
        return false;   /// Not enough memory! (try reducing bins)

    /* Actual size of contextual regions */
    subW     = w / blkw;
    subH     = h / blkh;
    nrpixels = (unsigned long)subW * (unsigned long)subH;

    if ( slope < 1.0f )
        slope = 1.0f;   /// if it going under 1.0f, there's a bug in CLAHE_ClipHistogram() !

    /* Calculate actual cliplimit */
    cliplimit = (unsigned long)( slope * ( subW * subH ) / bins );
    cliplimit = (cliplimit < 1UL) ? 1UL : cliplimit;

    /* Make lookup table for mapping of greyvalues */
    CLAHE_MakeLut( aLUT, minlvl, maxlvl, bins );

    /* Calculate greylevel mappings for each contextual region */
    for ( cntY=0, pImPointer=ptr; cntY<blkh; cntY++ )
    {
        for ( cntX=0; cntX<blkw; cntX++, pImPointer += subW )
        {
            pulHist = &pulMapArray[bins * (cntY * blkw + cntX)];

            CLAHE_MakeHistogram( pImPointer, w, subW, subH, pulHist, bins, aLUT );
            CLAHE_ClipHistogram( pulHist, bins, cliplimit );
            CLAHE_MapHistogram( pulHist, minlvl, maxlvl, bins, nrpixels );
        }

        /* skip lines, set pointer */
        pImPointer += ( subH - 1 ) * w;
    }

    /* Interpolate greylevel mappings to get CLAHE image */
    for ( pImPointer = ptr, cntY = 0; cntY<=blkh; cntY++ )
    {
        if (cntY == 0)
        {
            /* special case: top row */
            subY    = subH >> 1;
            inpolYU = 0;
            inpolYB = 0;
        }
        else
        {
            if (cntY == blkh)
            {
                /* special case: bottom row */
                subY    = subH >> 1;
                inpolYU = blkh-1;
                inpolYB = inpolYU;
            }
            else
            {
                /* default values */
                subY    = subH;
                inpolYU = cntY - 1;
                inpolYB = inpolYU + 1;
            }
        }

        for ( cntX = 0; cntX <= blkw; cntX++ )
        {
            if ( cntX == 0 )
            {
                /* special case: left column */
                subX    = subW >> 1;
                inpolXL = 0;
                inpolXR = 0;
            }
            else
            {
                if ( cntX == blkw )
                {
                    /* special case: right column */
                    subX    = subW >> 1;
                    inpolXL = blkw - 1;
                    inpolXR = inpolXL;
                }
                else
                {
                    /* default values */
                    subX    = subW;
                    inpolXL = cntX - 1;
                    inpolXR = inpolXL + 1;
                }
            }

            pulLU = &pulMapArray[ bins * ( inpolYU * blkw + inpolXL ) ];
            pulRU = &pulMapArray[ bins * ( inpolYU * blkw + inpolXR ) ];
            pulLB = &pulMapArray[ bins * ( inpolYB * blkw + inpolXL ) ];
            pulRB = &pulMapArray[ bins * ( inpolYB * blkw + inpolXR ) ];

            CLAHE_Interpolate( pImPointer, w,
                               pulLU, pulRU, pulLB, pulRB,
                               subX, subY, aLUT );

            /* set pointer on next matrix */
            pImPointer += subX;
        }

        pImPointer += (subY - 1) * w;
    }

    /* free space for histograms */
    delete[] pulMapArray;

    return true;
}

#else
bool RAWImageToolKit::ApplyCLAHE( unsigned short* ptr, unsigned w, unsigned h,
                                  unsigned short minlvl, unsigned short maxlvl,
                                  unsigned blkw, unsigned blkh, unsigned bins, float slope )
{
    int ret = CLAHE( ptr, w, h, minlvl, maxlvl, blkw, blkh, bins, slope );

    if ( ret == 0 )
        return true;

    return false;
}
#endif /// of CLAHE_NATIVE_FEATURE
