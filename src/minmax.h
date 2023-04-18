#ifndef __MINMAX_H__
#define __MINMAX_H__

#ifndef MIN
    #define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef MAX
    #define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#define _MIN_F_     ( -1e6F )
#define _MAX_F_     ( +1e6F )

#endif /// of __MINMAX_H__
