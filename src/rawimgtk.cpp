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

/* CLAHE_ClipHistogram() :
 * This function performs clipping of the histogram and redistribution of bins.
 * The histogram is clipped and the number of excess pixels is counted. Afterwards
 * the excess pixels are equally redistributed across the whole histogram (providing
 * the bin count is smaller than the cliplimit).
 */
void CLAHE_ClipHistogram( unsigned long* pHisto, unsigned int greyLvl, unsigned long clipLimit )
{
    unsigned long* pRangePtr = pHisto;
	unsigned long* pEndPtr = NULL;
	unsigned long* pHistPtr = NULL;

    unsigned long excessSz = 0;
	unsigned long upperSz;
	unsigned long rangeInc;
	unsigned long stepsz;
	unsigned long cnt;

    long lBinExcess;

	/* calculate total number of excess pixels */
    for ( cnt = 0; cnt < greyLvl; cnt++)
	{
		lBinExcess = (long) pRangePtr[cnt] - (long) clipLimit;

		/* excess in current bin */
		if ( lBinExcess > 0 )
		{
			excessSz += lBinExcess;
		}
    }

    /* Second part: clip histogram and redistribute excess pixels in each bin */
    rangeInc = excessSz / greyLvl;		 /// average binincrement
    upperSz  = clipLimit - rangeInc;	 /// Bins larger than upperSz set to cliplimit

    for ( cnt=0; cnt<greyLvl; cnt++ )
	{
		if (pHisto[cnt] > clipLimit)
		{
			pHisto[cnt] = clipLimit;    /// clip bin
		}
		else
		{
			if (pHisto[cnt] > upperSz)
			{	/* high bin count */
				excessSz    -= pHisto[cnt] - upperSz;
				pHisto[cnt]  = clipLimit;
			}
			else
			{	/* low bin count */
				excessSz    -= rangeInc;
				pHisto[cnt] += rangeInc;
			}
		}
    }

    while ( excessSz > 0 )
	{   /* Redistribute remaining excess  */
		pEndPtr = &pHisto[greyLvl]; pHistPtr = pHisto;

		while ( ( excessSz > 0 ) && ( pHistPtr < pEndPtr )  )
		{
			stepsz = greyLvl / excessSz;

			if ( stepsz < 1 )
			{
				/* stepsize at least 1 */
				stepsz = 1;
			}

			for ( pRangePtr=pHistPtr;
			      pRangePtr < pEndPtr && excessSz;
				  pRangePtr += stepsz)
			{
				if (*pRangePtr < clipLimit)
				{
					/* reduce excess */
					(*pRangePtr)++;
					excessSz--;
				}
			}

			/* restart redistributing on other bin location */
			pHistPtr++;
		}
    }
}

/* CLAHE_MakeHistogram() :
 * This function classifies the greylevels present in the array image into
 * a greylevel histogram. The pLookupTable specifies the relationship
 * between the greyvalue of the pixel (typically between 0 and 4095) and
 * the corresponding bin in the histogram (usually containing only 128 bins).
 */
void CLAHE_MakeHistogram ( unsigned short* pImage,
                           unsigned int imgWidth,
                           unsigned int rgnszW, unsigned int rgnszH,
                           unsigned long* pHisto,
                           unsigned int greyLvl, unsigned short* pLUT )
{
    unsigned short* pImgPtr = NULL;
    unsigned int cnt;

	/* clear histogram */
    for ( cnt=0; cnt<greyLvl; cnt++ )
	{
		pHisto[cnt] = 0L;
	}

    for ( cnt=0; cnt< rgnszH; cnt++ )
	{
		pImgPtr = &pImage[rgnszW];

		while ( pImage < pImgPtr )
		{
			pHisto[ pLUT[ *pImage++ ] ]++;
		}

		pImgPtr += imgWidth;
		pImage = &pImgPtr[-rgnszW];
    }
}

/* CLAHE_MapHistogram() :
 * This function calculates the equalized lookup table (mapping) by
 * cumulating the input histogram. Note: lookup table is rescaled in range [Min..Max].
 */
void CLAHE_MapHistogram ( unsigned long* pHisto,
                          unsigned short Min, unsigned short Max,
                          unsigned int greyLvl, unsigned long pixelsz )
{
    unsigned int  cnt = 0;
	unsigned long sum = 0;

    const float         fScale = ( (float)(Max - Min) ) / pixelsz;
    const unsigned long ulMin  = (unsigned long) Min;

    for ( cnt=0; cnt<greyLvl; cnt++)
	{
		sum += pHisto[cnt];
		pHisto[cnt] = (unsigned long)( ulMin + sum * fScale );

		if ( pHisto[cnt] > Max )
		{
			pHisto[cnt] = Max;
		}
    }
}

/* CLAHE_MakeLut() :
 * To speed up histogram clipping, the input image [Min,Max] is scaled down to
 * [0,uiNrBins-1]. This function calculates the LUT.
 */
void CLAHE_MakeLut( unsigned short * pLUT, unsigned short Min, unsigned short Max, unsigned int ranges )
{
    const unsigned short BinSize = (unsigned short) (1 + (Max - Min) / ranges);

    for ( int cnt=Min; cnt<=Max; cnt++)
    {
        pLUT[cnt] = (cnt - Min) / BinSize;
    }
}

/* CLAHE_Interpolate() :
 * pImage      - pointer to input/output image
 * imgWidth      - resolution of image in x-direction
 * pMap*     - mappings of greylevels from histograms
 * subszW     - subszW of image submatrix
 * subszH     - subszH of image submatrix
 * pLUT	       - lookup table containing mapping greyvalues to bins
 * This function calculates the new greylevel assignments of pixels within a submatrix
 * of the image with size subszW and subszH. This is done by a bilinear interpolation
 * between four different mappings in order to eliminate boundary artifacts.
 * It uses a division; since division is often an expensive operation, I added code to
 * perform a logical shift instead when feasible.
 */
void CLAHE_Interpolate( unsigned short * pImage,
                        int imgWidth, unsigned long * pMapLU,
                        unsigned long * pMapRU, unsigned long * pMapLB,  unsigned long * pMapRB,
                        unsigned int subszW, unsigned int subszH, unsigned short * pLUT)
{
    const unsigned int incSz = imgWidth-subszW; /* Pointer increment after processing row */
    unsigned short GreyValue;
	unsigned int normFactor = subszW * subszH; /* Normalization factor */

    unsigned int coefW = 0;
	unsigned int coefH = 0;
	unsigned int invcoefW = 0;
	unsigned int invcoefH = 0;
	unsigned int shifts = 0;

	/* If normFactor is not a power of two, use division */
    if ( normFactor & (normFactor - 1) )
	{
		for ( coefH=0, invcoefH = subszH;
		      coefH < subszH;
		      coefH++, invcoefH--,pImage+=incSz )
		{
			for ( coefW=0, invcoefW = subszW;
			      coefW < subszW;
			      coefW++, invcoefW--)
			{
				/* get histogram bin value */
				GreyValue = pLUT[*pImage];

				*pImage++ = (unsigned short )( ( invcoefH *
												 ( invcoefW*pMapLU[GreyValue] +  coefW * pMapRU[GreyValue] ) +
						                         coefH * (invcoefW * pMapLB[GreyValue] +
												 coefW * pMapRB[GreyValue])) / normFactor );
			}
		}
	}
    else
	{	/* avoid the division and use a right shift instead */
		while ( normFactor >>= 1 )
		{
			/* Calculate 2log of normFactor */
			shifts++;
		}

		for ( coefH = 0, invcoefH = subszH;
		      coefH < subszH;
			  coefH++, invcoefH--,pImage+=incSz )
		{
			 for ( coefW = 0, invcoefW = subszW;
			       coefW < subszW;
			       coefW++, invcoefW-- )
			{
				/* get histogram bin value */
				GreyValue = pLUT[*pImage];

				*pImage++ = (unsigned short)( ( invcoefH *
				                                ( invcoefW * pMapLU[GreyValue] + coefW * pMapRU[GreyValue] ) +
												coefH * (invcoefW * pMapLB[GreyValue] +
												coefW * pMapRB[GreyValue])) >> shifts );
			}
		}
    }
}

////////////////////////////////////////////////////////////////////////////////

/* RAWImageToolKit::ApplyCLAHE() :
 *   pImage - Pointer to the input/output image
 *   imgWidth - Image resolution in the X direction
 *   imgHeight - Image resolution in the Y direction
 *   Min - Minimum greyvalue of input image (also becomes minimum of output image)
 *   Max - Maximum greyvalue of input image (also becomes maximum of output image)
 *   rgnWidth - Number of contextial regions in the X direction (min 2, max DEF_MAX_WIDTH)
 *   rgnHeight - Number of contextial regions in the Y direction (min 2, max DEF_MAX_HEIGHT)
 *   ranges - Number of greybins for histogram ("dynamic range")
 *   float fCliplimit - Normalized cliplimit (higher values give more contrast)
 * The number of "effective" greylevels in the output image is set by ranges; selecting
 * a small value (eg. 128) speeds up processing and still produce an output image of
 * good quality. The output image will have the same minimum and maximum value as the input
 * image. A clip limit smaller than 1 results in standard (non-contrast limited) AHE.
 */

bool RAWImageToolKit::ApplyCLAHE( unsigned short* pImage,
                                  unsigned int imgWidth, unsigned int imgHeight,
                                  unsigned short Min, unsigned short Max,
                                  unsigned int rgnWidth, unsigned int rgnHeight,
                                  unsigned int ranges, float fCliplimit )
{
    /* counters */
    unsigned int cnt_x;
    unsigned int cnt_y;

    /* size of context. reg. and subimages */
    unsigned int subszW = 0;
    unsigned int subszH = 0;
    unsigned int subImgW = 0;
    unsigned int subImgH = 0;

    /* auxiliary variables interpolation routine */
    unsigned int cnt_xL = 0;
    unsigned int cnt_xR = 0;
    unsigned int cnt_yU = 0;
    unsigned int cnt_yB = 0;

    /* clip limit and region pixel count */
    unsigned long clipLimit = 0;
    unsigned long pixelCnts = 0;

    /* pointer to image */
    unsigned short* pImgPtr = NULL;

    /* lookup table used for scaling of input image */
    unsigned short aLUT[ RAWIMGTK_MAX_D_ARRAY_SZ ] = {0};

    /* pointer to histogram and mappings*/
    unsigned long* pulHist = NULL;
    unsigned long* pMapArray = NULL;

    /* auxiliary pointers interpolation */
    unsigned long* pulLU = NULL;
    unsigned long* pulLB = NULL;
    unsigned long* pulRU = NULL;
    unsigned long* pulRB = NULL;

    if ( rgnWidth > CLAHE_MAX_REG_W )
		rgnWidth = CLAHE_MAX_REG_W;

    if ( rgnHeight > CLAHE_MAX_REG_H )
		rgnHeight = CLAHE_MAX_REG_H;

    if ( imgWidth % rgnWidth )
		return false;       /// x-resolution no multiple of rgnWidth

    if ( imgHeight & rgnHeight )
		return false;       /// y-resolution no multiple of rgnHeight

    if ( Max >= RAWIMGTK_MAX_D_ARRAY_SZ )
		return false;       /// maximum too large

    if ( Min >= Max )
		return false;       /// minimum equal or larger than maximum

    if ( rgnWidth < 2 )
        rgnWidth = 2;

    if ( rgnHeight < 2 )
		rgnHeight = 2;

    if ( fCliplimit == 1.0f )
		return true;	    /// is OK, immediately returns original image.

    if ( fCliplimit < 1.0f )
        fCliplimit = 1.0f;

    if ( ranges == 0 )
		ranges = 128;	    /// default value when not specified

    pMapArray = new unsigned long[ rgnWidth * rgnHeight * ranges ];

    if ( pMapArray == NULL )
		return false;

    subszW    = imgWidth/rgnWidth; subszH = imgHeight/rgnHeight;  /* Actual size of contextual regions */
    pixelCnts = (unsigned long)subszW * (unsigned long)subszH;

    if(fCliplimit > 0.0)
	{
		/* Calculate actual cliplimit	 */
		clipLimit = (unsigned long) (fCliplimit * (subszW * subszH) / ranges);
        clipLimit = (clipLimit < 1UL) ? 1UL : clipLimit;
    }
    else
	{
		/* Large value, do not clip (AHE) */
		clipLimit = 1UL<<14;
	}

    CLAHE_MakeLut(aLUT, Min, Max, ranges);	  /* Make lookup table for mapping of greyvalues */

    /* Calculate greylevel mappings for each contextual region */
    for ( cnt_y=0, pImgPtr=pImage; cnt_y<rgnHeight; cnt_y++ )
	{
		for ( cnt_x=0; cnt_x<rgnWidth; cnt_x++, pImgPtr+=subszW )
		{
			pulHist = &pMapArray[ranges * (cnt_y * rgnWidth + cnt_x)];

			CLAHE_MakeHistogram( pImgPtr,
                                 imgWidth, subszW, subszH,
                                 pulHist, ranges, aLUT );
			CLAHE_ClipHistogram( pulHist, ranges, clipLimit);
			CLAHE_MapHistogram( pulHist, Min, Max, ranges, pixelCnts);
		}

		pImgPtr += (subszH - 1) * imgWidth;		  /* skip lines, set pointer */
    }

    /* Interpolate greylevel mappings to get CLAHE image */
    for (pImgPtr = pImage, cnt_y = 0; cnt_y <= rgnHeight; cnt_y++)
	{
		if (cnt_y == 0)
		{					  /* special case: top row */
			subImgH = subszH >> 1;  cnt_yU = 0; cnt_yB = 0;
		}
		else
		{
			if (cnt_y == rgnHeight)
			{				  /* special case: bottom row */
				subImgH = subszH >> 1;	cnt_yU = rgnHeight-1;	 cnt_yB = cnt_yU;
			}
			else
			{					  /* default values */
				subImgH = subszH; cnt_yU = cnt_y - 1; cnt_yB = cnt_yU + 1;
			}
		}
		for (cnt_x = 0; cnt_x <= rgnWidth; cnt_x++)
		{
			if (cnt_x == 0)
			{				  /* special case: left column */
				subImgW = subszW >> 1; cnt_xL = 0; cnt_xR = 0;
			}
			else
			{
				if (cnt_x == rgnWidth)
				{			  /* special case: right column */
					subImgW = subszW >> 1;  cnt_xL = rgnWidth - 1; cnt_xR = cnt_xL;
				}
				else
				{					  /* default values */
					subImgW = subszW; cnt_xL = cnt_x - 1; cnt_xR = cnt_xL + 1;
				}
			}

			pulLU = &pMapArray[ranges * (cnt_yU * rgnWidth + cnt_xL)];
			pulRU = &pMapArray[ranges * (cnt_yU * rgnWidth + cnt_xR)];
			pulLB = &pMapArray[ranges * (cnt_yB * rgnWidth + cnt_xL)];
			pulRB = &pMapArray[ranges * (cnt_yB * rgnWidth + cnt_xR)];

			CLAHE_Interpolate( pImgPtr,
                               imgWidth, pulLU, pulRU, pulLB, pulRB,
                               subImgW, subImgH, aLUT );

			pImgPtr += subImgW;			  /* set pointer on next matrix */
		}

		pImgPtr += (subImgH - 1) * imgWidth;
    }

    delete[] pMapArray;

	return true;
}

