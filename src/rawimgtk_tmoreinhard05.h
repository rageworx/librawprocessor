#ifndef __RAWIMGTK_TMOREINHARD05_H__
#define __RAWIMGTK_TMOREINHARD05_H__

/*******************************************************************************
** " Global and/or local tone mapping operator "
** -----------------------------------------------------------------------------
**  Source code re-programmed from FreeImage 3 Library,
**                                       Tone mapping operator (Reinhard 2005 )
**                       FreeImage3 lib. Design and implementation by
**                           - Herv?Drolon (drolon@infonie.fr)
**                           - Mihail Naydenov (mnaydenov@users.sourceforge.net)
**   Referenced:
**       [1] Erik Reinhard and Kate Devlin,
**                'Dynamic Range Reduction Inspired by Photographer Physiology',
**          IEEE Transactions on Visualization and Computer Graphics, 11(1),
**                                                                 Jan/Feb 2005.
**       [2] Erik Reinhard,
**                    'Parameter estimation for photographic tone reproduction',
**                    Journal of Graphics Tools, vol. 7, no. 1, pp. 45?1, 2003.
**
**   Re-programmed for librawprocessor:
**      Raphael Kim ( rageworx@gmail.com, rage.kim@gmail.com )
*******************************************************************************/

namespace RAWImageToolKit
{
    /***
    * parameters :
    *  - intensity Overall intensity in range [-8:8] : default to 0
    *  - contrast Contrast in range [0.3:1) : default to 0
    *  - adaptation Adaptation in range [0:1] : default to 1
    *  - color_correction Color correction in range [0:1] : default to 0
    ***/
    bool tmoReinhard2005( unsigned short* src, unsigned srcsz, float normalf,
                          float intensity, float contrast, float adaptation, float colorcorrection );
};

#endif /// of __RAWIMGTK_TMOREINHARD05_H__
