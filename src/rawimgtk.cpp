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

#define RAWImageToolKitSwapUS( _a_, _b_ )   float t=_a_; _a_=_b_; _b_=t;

////////////////////////////////////////////////////////////////////////////////

bool RAWImageToolKit::FlipHorizontal( float* ptr, unsigned w, unsigned h )
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
                RAWImageToolKitSwapUS( ptr[ w * ( h - 1 - cnth ) + cntw ],
                                       ptr[ w * cnth + cntw ] );
            }
        }

        return true;
    }

    return false;
}

bool RAWImageToolKit::FlipVertical( float* ptr, unsigned w, unsigned h )
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

bool RAWImageToolKit::Rotate90( float* ptr, unsigned* w, unsigned* h )
{
    unsigned cur_w = *w;
    unsigned cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        unsigned src_x = 0;
        unsigned src_y = 0;
        unsigned new_w = cur_h;
        unsigned new_h = cur_w;

        float* tempbuff = new float[ new_w * new_h ];

        if ( tempbuff == NULL )
            return false;

        memset( tempbuff, 0, new_w * new_h * sizeof( float ) );

        unsigned cntw = 0;
        unsigned cnth = 0;


        #pragma omp parallel for private(cnth)
        for( cntw=new_w-1; cntw>0; cntw-- )
        {
            for( cnth=0; cnth<new_h; cnth++ )
            {
                unsigned pos1 = new_w * cnth + cntw;
                unsigned pos2 = cur_w * ( new_w - cntw - 1 ) + cnth;

				tempbuff[ pos1 ] = ptr[ pos2 ];
            }
        }

        memcpy( ptr, tempbuff,  new_w * new_h * sizeof( float ) );

        delete[] tempbuff;

        *w = new_w;
        *h = new_h;

        return true;
    }

    return false;
}

bool RAWImageToolKit::Rotate180( float* ptr, unsigned* w, unsigned* h )
{
    unsigned cur_w = *w;
    unsigned cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        unsigned cntw = 0;
        unsigned cnth = 0;
		unsigned imgmax = *w * *h;
		unsigned cntmax = imgmax / 2;

        #pragma omp parallel for
        for( unsigned cnt=0; cnt<cntmax; cnt++ )
        {
			RAWImageToolKitSwapUS( ptr[ cnt ], ptr[ imgmax - cnt ] );
        }

        return true;
    }

    return false;
}

bool RAWImageToolKit::Rotate270( float* ptr, unsigned* w, unsigned* h )
{
    unsigned cur_w = *w;
    unsigned cur_h = *h;

    if ( ( cur_w > 0 ) && ( cur_h > 0 ) )
    {
        unsigned src_x = 0;
        unsigned src_y = 0;
        unsigned new_w = cur_h;
        unsigned new_h = cur_w;

        float* tempbuff = new float[ new_w * new_h ];

        if ( tempbuff == NULL )
            return false;

        memset( tempbuff, 0,  new_w * new_h * sizeof( float ) );

        unsigned cntw = 0;
        unsigned cnth = 0;

        #pragma omp parallel for private(cnth)
        for( cntw=0; cntw<new_w; cntw++ )
        {
            for( cnth=new_h-1; cnth>0; cnth-- )
            {
                unsigned pos1 = new_w * cnth + cntw;
                unsigned pos2 = cur_w * cntw + new_h - cnth;

                tempbuff[ pos1 ] = ptr[ pos2 ];
            }
        }

        memcpy( ptr, tempbuff,  new_w * new_h * sizeof( float ) );

        delete[] tempbuff;

        *w = new_w;
        *h = new_h;

        return true;
    }

    return false;
}

float ritk_min4(float a, float b, float c, float d)
{
   float mn = a;
   if(mn > b) mn = b;
   if(mn > c) mn = c;
   if(mn > d) mn = d;
   return mn;
}

float ritk_max4(float a, float b, float c, float d)
{
   float mx = a;
   if(mx < b) mx = b;
   if(mx < c) mx = c;
   if(mx < d) mx = d;
   return mx;
}

bool RAWImageToolKit::RotateFree( float* ptr, unsigned* w, unsigned* h, float** newptr, float degree, float background )
{
    if ( ( ptr == NULL ) || ( w == 0 ) || ( h == 0 ) )
        return false;

    float fdeg = ( degree - ( degree / 360.0f ) ) / 360.0f * 100.0f;

    int img_w = *w;
    int img_h = *h;

	float CtX = ( (float) img_w ) / 2.0f;
	float CtY = ( (float) img_h ) / 2.0f;

	float cA = (float)cos( fdeg );
	float sA = (float)sin( fdeg );

	float x1 = CtX + (-CtX) * cA - (-CtY) * sA;
	float x2 = CtX + (img_w - CtX) * cA - (-CtY) * sA;
	float x3 = CtX + (img_w - CtX) * cA - (img_h - CtY) * sA;
	float x4 = CtX + (-CtX) * cA - (img_h - CtY) * sA;

	float y1 = CtY + (-CtY) * cA + (-CtX) * sA;
	float y2 = CtY + (img_h - CtY) * cA + (-CtX) * sA;
	float y3 = CtY + (img_h - CtY) * cA + (img_w - CtX) * sA;
	float y4 = CtY + (-CtY) * cA + (img_w - CtX) * sA;

	int OfX = ((int)floor(ritk_min4(x1, x2, x3, x4)));
	int OfY = ((int)floor(ritk_min4(y1, y2, y3, y4)));

	int dstW = ((int)ceil(ritk_max4(x1, x2, x3, x4))) - OfX;
	int dstH = ((int)ceil(ritk_max4(y1, y2, y3, y4))) - OfY;

    // Now new image !
    float* obuff = new float[ dstW * dstH ];

    if ( obuff == NULL )
        return false;

    #pragma omp parellel for
    for( int cnt=0; cnt<(dstW * dstH); cnt++ )
    {
        obuff[ cnt ] = background;
    }

    int stepY = 0;
    int stepX = 0;

    // pointer to destination.
    float* dst = obuff;

	#pragma omp parellel for private( stepX )
	for ( stepY = 0; stepY<dstH; stepY++ )
	{
		for ( stepX = 0; stepX<dstW; stepX++ )
		{
            float CtX2 = CtX - OfX;
            float CtY2 = CtY - OfY;

            float orgX = ( cA*(stepX-CtX2) + sA*(stepY-CtY2)) + CtX;
            float orgY = (-sA*(stepX-CtX2) + cA*(stepY-CtY2)) + CtY;

			int iorgX  = (int) orgX;
			int iorgY  = (int) orgY;

            float diffX = (orgX - iorgX);
            float diffY = (orgY - iorgY);

            if ((orgX >= 0) && (orgY >= 0) && (orgX < img_w-1) && (orgY < img_h-1))
            {
                float* pd = &obuff[ stepY * dstW + stepX ];
                float* ps = &ptr[ iorgX + iorgY * img_w ];

                // Doing interpolated pixel calculation .
                float pv[4] = {0.0f};

                pv[0] = ps[ 0 ];
                pv[1] = ps[ 1 ]; /// Right pixel
                pv[2] = ps[ img_w ]; /// Below pixel
                pv[3] = ps[ img_w + 1 ]; /// Right below pixel.

                float \
                pf  = pv[0] * ( 1.0f - diffX ) * ( 1.0f - diffY ) +
                      pv[1] *          diffX   * ( 1.0f - diffY ) +
                      pv[2] * ( 1.0f - diffX ) *          diffY +
                      pv[3] *          diffX   *          diffY;

                *pd = (float)pf;
            }
		}

		dst += dstW;
	}

	*newptr = obuff;
	*w      = dstW;
	*h      = dstH;

	return true;
}

bool RAWImageToolKit::Crop( float* ptr, unsigned w, unsigned h, float** newptr, unsigned cx, unsigned cy, unsigned cw, unsigned ch )
{
    if ( ( ptr == NULL ) || ( w == 0 ) || ( h == 0 ) )
        return false;

    if ( ( cx == 0 ) || ( cx > w ) || ( cw > w ) ||
         ( cy == 0 ) || ( cy > h ) || ( ch > h ) )
    return false;

    unsigned rsx = cx;
    unsigned rsy = cy;
    unsigned rw  = cw;
    unsigned rh  = ch;

    if ( w < ( rw + rsx ) )
    {
        rw = w - rsx;
    }

    if ( h < ( rh + rsy ) )
    {
        rh = h - rsy;
    }

    float* obuff = new float[ rw * rh ];

    if ( obuff != NULL )
    {
        unsigned srcw = w;
        unsigned srch = h;
        unsigned dsth = srch - rsy;
        unsigned cnty;
        unsigned cntx;

        #pragma omp parellel for private( cntx )
        for( cnty=rsy; cnty<srch; cnty++ )
        {
            for( cntx=rsx; cntx<srcw; cntx ++ )
            {
                obuff[ (cnty - rsy) * rw + ( cntx - rsx) ] = ptr[ cnty * srcw + cntx ];
            }
        }

        *newptr = obuff;

        return true;
    }

    return false;
}

bool RAWImageToolKit::CropCenter( float* ptr, unsigned w, unsigned h, float** newptr, unsigned cw, unsigned ch )
{
    if ( ( ptr == NULL ) || ( w == 0 ) || ( h == 0 ) ||
         ( cw == 0 ) || ( cw > w ) || ( ch == 0 ) || ( ch > h ) )
        return false;

    unsigned c_l = ( w - cw ) / 2;
    unsigned c_t = ( h - ch ) / 2;

    return Crop( ptr, w, h, newptr, c_l, c_t, cw, ch );
}

bool RAWImageToolKit::AdjustGamma( float* ptr, unsigned arraysz, double gamma )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    float LUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};
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
        LUT[ cnt ] = (float)floor( newlvl + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT );
}

bool RAWImageToolKit::AdjustBrightness( float* ptr, unsigned arraysz, double perc )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    float LUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};
    const \
    double bscaled  = ( 100.0 + perc ) / 100.0;
    double bval     = 0.0;

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<RAWIMGTK_MAX_D_ARRAY_SZ; cnt++ )
    {
        bval = (double)cnt * bscaled;
        bval = MAX( 0.0, MIN( bval, RAWIMGTK_MAX_F_VAL ) );
        LUT[ cnt ] = (float)floor( bval + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT );
}

bool RAWImageToolKit::AdjustContrast( float* ptr, unsigned arraysz, double perc )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) )
        return false;

    float LUT[RAWIMGTK_MAX_D_ARRAY_SZ] = {0};
    const \
    double bscaled  = ( 100.0 + perc ) / 100.0;
    double bval     = 0.0;

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<RAWIMGTK_MAX_D_ARRAY_SZ; cnt++ )
    {
        bval = (double)RAWIMGTK_HALF_F_VAL + ( (double)cnt - RAWIMGTK_HALF_F_VAL ) * bscaled;
        bval = MAX( 0.0, MIN( bval, RAWIMGTK_MAX_F_VAL ) );
        LUT[ cnt ] = (float)floor( bval + 0.5 );
    }

    return AdjustCurve( ptr, arraysz, LUT );
}

bool RAWImageToolKit::AdjustCurve( float* ptr, unsigned arraysz, float* LUT )
{
    if ( ( ptr == NULL ) || ( arraysz == 0 ) || ( LUT == NULL ) )
        return false;

    #pragma omp parallel for
    for( unsigned cnt=0; cnt<arraysz; cnt++ )
    {
        //ptr[ cnt ] = LUT[ ptr[ cnt ] ];
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
void CLAHE_MakeHistogram ( float* pImage,
                           unsigned int imgWidth,
                           unsigned int rgnszW, unsigned int rgnszH,
                           unsigned long* pHisto,
                           unsigned int greyLvl, float* pLUT )
{
    float* pImgPtr = NULL;
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
			//pHisto[ pLUT[ (float)*pImage++ ] ]++;
		}

		pImgPtr += imgWidth;
		pImage = &pImgPtr[-(long long)rgnszW];
    }
}

/* CLAHE_MapHistogram() :
 * This function calculates the equalized lookup table (mapping) by
 * cumulating the input histogram. Note: lookup table is rescaled in range [Min..Max].
 */
void CLAHE_MapHistogram ( unsigned long* pHisto,
                          float Min, float Max,
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
void CLAHE_MakeLut( float * pLUT, float Min, float Max, unsigned int ranges )
{
    const float BinSize = (float) (1 + (Max - Min) / ranges);

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
void CLAHE_Interpolate( float * pImage,
                        int imgWidth, unsigned skipWidth,
                        unsigned long * pMapLU, unsigned long * pMapRU, unsigned long * pMapLB,  unsigned long * pMapRB,
                        unsigned int subszW, unsigned int subszH, float * pLUT)
{
    const unsigned int incSz = imgWidth-subszW; /* Pointer increment after processing row */
    float GreyValue;
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
				/*
				GreyValue = pLUT[*pImage];

				*pImage++ = (float )( ( invcoefH *
												 ( invcoefW*pMapLU[GreyValue] +  coefW * pMapRU[GreyValue] ) +
						                         coefH * (invcoefW * pMapLB[GreyValue] +
												 coefW * pMapRB[GreyValue]) )
                                               / normFactor );

                pImage += skipWidth;
                */
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
			  coefH++, invcoefH--, pImage+=incSz )
		{
			 for ( coefW = 0, invcoefW = subszW;
			       coefW < subszW;
			       coefW++, invcoefW-- )
			{
				/* get histogram bin value */
				/*
				GreyValue = pLUT[*pImage];

				*pImage++ = (float)( ( invcoefH *
				                                ( invcoefW * pMapLU[GreyValue] + coefW * pMapRU[GreyValue] ) +
												coefH * (invcoefW * pMapLB[GreyValue] +
												coefW * pMapRB[GreyValue]) )
                                              >> shifts );
                pImage += skipWidth;
                */
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

bool RAWImageToolKit::ApplyCLAHE( float* pImage,
                                  unsigned int imgWidth, unsigned int imgHeight,
                                  float Min, float Max,
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
    unsigned int skipW = 0;

    /* auxiliary variables interpolation routine */
    unsigned int cnt_xL = 0;
    unsigned int cnt_xR = 0;
    unsigned int cnt_yU = 0;
    unsigned int cnt_yB = 0;

    /* clip limit and region pixel count */
    unsigned long clipLimit = 0;
    unsigned long pixelCnts = 0;

    /* pointer to image */
    float* pImgPtr = NULL;

    /* lookup table used for scaling of input image */
    float aLUT[ RAWIMGTK_MAX_D_ARRAY_SZ ] = {0};

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

    /*
    if ( imgWidth % rgnWidth )
		return false;       /// x-resolution no multiple of rgnWidth

    if ( imgHeight % rgnHeight )
		return false;       /// y-resolution no multiple of rgnHeight
    */

    if ( Min >= Max )
		return false;       /// minimum equal or larger than maximum

    if ( rgnWidth < 2 )
        rgnWidth = 2;

    if ( rgnHeight < 2 )
		rgnHeight = 2;

    // Make it skip when image size not divided done by zero.
    skipW = imgWidth % rgnWidth;
    if ( skipW > 0 )
    {
        skipW++;
    }

    if ( fCliplimit == 1.0f )
		return true;	    /// is OK, immediately returns original image.

    if ( fCliplimit < 1.0f )
        return false;

    if ( fCliplimit < 2.0f )
        fCliplimit = 2.1f;

    if ( ranges == 0 )
		ranges = 65536;	    /// default value when not specified

    pMapArray = new unsigned long[ rgnWidth * rgnHeight * ( ranges + 1 ) ];

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

		/* skip lines, set pointer */
		pImgPtr += (subszH - 1) * imgWidth + skipW;
    }

    /* Interpolate greylevel mappings to get CLAHE image */
    for (pImgPtr = pImage, cnt_y = 0; cnt_y <= rgnHeight; cnt_y++)
	{
		if (cnt_y == 0)
		{
		    /* special case: top row */
			subImgH = subszH >> 1;
			cnt_yU = 0;
			cnt_yB = 0;
		}
		else
		{
			if (cnt_y == rgnHeight)
			{
			    /* special case: bottom row */
				subImgH = subszH >> 1;
				cnt_yU = rgnHeight-1;
				cnt_yB = cnt_yU;
			}
			else
			{
			    /* default values */
				subImgH = subszH;
				cnt_yU = cnt_y - 1;
				cnt_yB = cnt_yU + 1;
			}
		}
		for (cnt_x = 0; cnt_x <= rgnWidth; cnt_x++)
		{
			if (cnt_x == 0)
			{
			    /* special case: left column */
				subImgW = subszW >> 1;
				cnt_xL = 0;
				cnt_xR = 0;
			}
			else
			{
				if (cnt_x == rgnWidth)
				{
				    /* special case: right column */
					subImgW = subszW >> 1;
					cnt_xL = rgnWidth - 1;
					cnt_xR = cnt_xL;
				}
				else
				{
				    /* default values */
					subImgW = subszW;
					cnt_xL = cnt_x - 1;
					cnt_xR = cnt_xL + 1;
				}
			}

			pulLU = &pMapArray[ranges * (cnt_yU * rgnWidth + cnt_xL)];
			pulRU = &pMapArray[ranges * (cnt_yU * rgnWidth + cnt_xR)];
			pulLB = &pMapArray[ranges * (cnt_yB * rgnWidth + cnt_xL)];
			pulRB = &pMapArray[ranges * (cnt_yB * rgnWidth + cnt_xR)];

			CLAHE_Interpolate( pImgPtr,
                               imgWidth, 0,
                               pulLU, pulRU, pulLB, pulRB,
                               subImgW, subImgH, aLUT );

			/* set pointer on next matrix */
			pImgPtr += subImgW;
		}

		if ( skipW > 0 )
        {
            pImgPtr += (subImgH - 1) * imgWidth + skipW - 1;
        }
        else
        {
            pImgPtr += (subImgH - 1) * imgWidth;
        }

    }

    delete[] pMapArray;

	return true;
}

