#ifndef __RAWFILTERS_H__
#define __RAWFILTERS_H__

#include "rawprocessor.h"

namespace RAWImageFilterKit
{
    #define RAWFILTER_PRESET_BLUR       "blur"
    #define RAWFILTER_PRESET_BLURMORE   "blurmore"
    #define RAWFILTER_PRESET_EDGE       "edgeonly"
    #define RAWFILTER_PRESET_EDGE2      "edgemore"
    #define RAWFILTER_PRESET_EDGE3      "edgestrong"
    #define RAWFILTER_PRESET_SHARPEN    "sharpen"
    #define RAWFILTER_PRESET_SHARPEN2   "sharpenhard"
    #define RAWFILTER_PRESET_SHARPENSB  "sharpensubtle"
    #define RAWFILTER_PRESET_UNSHARPENM "unsharpenmask"

    RAWProcessor::FilterConfig* GetPresetFilter( const char* name );
    void RemoveFilter( RAWProcessor::FilterConfig* fp );
}

#endif /// of __RAWFILTERS_H__
