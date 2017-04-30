#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "stdunicode.h"

char* convertW2M(const wchar_t* src)
{
    static char* convM = NULL;

    if ( src == NULL )
    {
        if ( convM != NULL )
        {
            delete[] convM;
            convM = NULL;
        }
    }

    int sz  = wcslen( src );

    if ( convM != NULL )
    {
        delete[] convM;
        convM = NULL;
    }

    convM = new char[ sz + 1 ];
    memset( convM, 0, sz + 1 );

    wcstombs( convM, src, sz );

    return convM;
}

wchar_t* convertM2W(const char* src)
{
    static wchar_t* convW = NULL;

    if ( src == NULL )
    {
        if ( convW != NULL )
        {
            delete[] convW;
            convW = NULL;
        }
    }

    int sz  = strlen( src );

    if ( convW != NULL )
    {
        delete[] convW;
        convW = NULL;
    }

    convW = new wchar_t[ sz + 1 ];
    memset( convW, 0, ( sz + 1 ) * sizeof(wchar_t) );

    mbstowcs( convW, src, sz * sizeof(wchar_t) );

    return convW;
}
