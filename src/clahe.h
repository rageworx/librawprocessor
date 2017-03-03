#ifndef CLAHE_NATIVE_FEATURE
#ifndef __CLAHE_H__
#define __CLAHE_H__

/*
 * ANSI C code from the article
 * "Contrast Limited Adaptive Histogram Equalization"
 * by Karel Zuiderveld, karel@cv.ruu.nl
 * in "Graphics Gems IV", Academic Press, 1994
 *
 *
 *  These functions implement Contrast Limited Adaptive Histogram Equalization.
 *  The main routine (CLAHE) expects an input image that is stored contiguously in
 *  memory;  the CLAHE output image overwrites the original input image and has the
 *  same minimum and maximum values (which must be provided by the user).
 *  This implementation assumes that the X- and Y image resolutions are an integer
 *  multiple of the X- and Y sizes of the contextual regions. A check on various other
 *  error conditions is performed.
 *
 *  #define the symbol BYTE_IMAGE to make this implementation suitable for
 *  8-bit images. The maximum number of contextual regions can be redefined
 *  by changing uiMAX_REG_X and/or uiMAX_REG_Y; the use of more than 256
 *  contextual regions is not recommended.
 *
 *  The code is ANSI-C and is also C++ compliant.
 *
 *  Author: Karel Zuiderveld, Computer Vision Research Group,
 *	     Utrecht, The Netherlands (karel@cv.ruu.nl)
 */

#define     uiNR_OF_GREY 0xFFFF

typedef unsigned short kz_pixel_t;	 /* for 12 bit-per-pixel images (default) */

#ifdef __cplusplus
    extern "C" {
#endif // __cplusplus

/******** Prototype of CLAHE function. Put this in a separate include file. *****/
int CLAHE(kz_pixel_t* pImage, unsigned int uiXRes, unsigned int uiYRes, kz_pixel_t Min,
	      kz_pixel_t Max, unsigned int uiNrX, unsigned int uiNrY,
	      unsigned int uiNrBins, float fCliplimit);

#ifdef __cplusplus
    };
#endif // __cplusplus

#endif
#endif /// of CLAHE_NATIVE_FEATURE
