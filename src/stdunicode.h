#ifndef __STDUNICODE_H__
#define __STDUNICODE_H__
////////////////////////////////////////////////////////////////////////////////
//
//  stdunicode
//  ==========================================================================
//  (C)Copyright 2013, Raphael Kim (rage.kim@gmail.com)
//
////////////////////////////////////////////////////////////////////////////////

#ifdef RAWPROCESSOR_USE_LOCALTCHAR
#include "tchar.h"
#endif /// of RAWPROCESSOR_USE_LOCALTCHAR
#include <wchar.h>

char*    convertW2M(const wchar_t* src);
wchar_t* convertM2W(const char* src);

#endif /// of __STDUNICODE_H__
