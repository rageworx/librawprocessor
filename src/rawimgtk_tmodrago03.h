#ifndef __RAWIMGTK_TMODRAGO03_H__
#define __RAWIMGTK_TMODRAGO03_H__

/*******************************************************************************
** " Alogrithm of Logarithmic mapping operator "
** -----------------------------------------------------------------------------
**  Source code re-programmed from FreeImage 3 Library,
**                                           Tone mapping operator (Drago, 2003)
**                       FreeImage3 lib. Design and implementation by
**                                           - Herv?Drolon (drolon@infonie.fr)
**   Referenced:
**      [1] F. Drago, K. Myszkowski, T. Annen, and N. Chiba,
**   Adaptive Logarithmic Mapping for Displaying High Contrast Scenes,
**                                                            Eurographics 2003.
**   Re-programmed for librawprocessor:
**      Raphael Kim ( rageworx@gmail.com, rage.kim@gmail.com )
*******************************************************************************/

namespace RAWImageToolKit
{
    /***
    * parameters :
    *  - gamma    : Gamma correction (gamma > 0).
    *               1 means no correction, 2.2 in the original paper.
    *  - exposure : Exposure parameter (0 means no correction, 0 in the original paper)
    ***/
    bool tmoDrago03( unsigned short *src, unsigned srcsz, float compressf,
                     double gamma, double exposure);
};

#endif /// of __RAWIMGTK_TMODRAGO03_H__
