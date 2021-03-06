#ifdef USE_OMP
#include <omp.h>
#endif // USE_OMP

#include "rawscale.h"
#include "minmax.h"

/// Clamp function
template <class T> T CLAMP(const T &value, const T &min_value, const T &max_value)
{
	return ((value < min_value) ? min_value : (value > max_value) ? max_value : value);
}

RawScaleWeightsTable::RawScaleWeightsTable( RAWGenericFilter* pFilter, unsigned uDstSize,
                                            unsigned uSrcSize )
 : _WeightTable( NULL ),
   _WindowSize( 0 ),
   _LineLength( uDstSize )
{
    if ( pFilter != NULL )
    {

        unsigned    u;
        double      dWidth;
        double      dFScale         = 1.0;
        const \
        double      dFilterWidth    = pFilter->GetWidth();
        const \
        double      dScale = double(uDstSize) / double(uSrcSize);

        if( dScale < 1.0 )
        {
            dWidth  = dFilterWidth / dScale;
            dFScale = dScale;
        }
        else
        {
            dWidth= dFilterWidth;
        }

        _WindowSize = 2 * (int)ceil(dWidth) + 1;

        _WeightTable = new Contribution[ _LineLength + 1 ];

        if ( _WeightTable != NULL)
        {
            for( u=0; u<_LineLength; u++ )
            {
                _WeightTable[ u ].Weights = new double[ _WindowSize + 1 ];
            }

            const double dOffset = ( 0.5 / dScale ) - 0.5;

            for( u=0; u<_LineLength; u++ )
            {
                const double dCenter = (double)u / dScale + dOffset;

                int iLeft  = MAX( 0, (int)floor (dCenter - dWidth) );
                int iRight = MIN( (int)ceil (dCenter + dWidth), int(uSrcSize) - 1 );

                if( ( iRight - iLeft + 1 ) > int(_WindowSize) )
                {
                    if( iLeft < ( int(uSrcSize) - 1 / 2 ) )
                    {
                        iLeft++;
                    }
                    else
                    {
                        iRight--;
                    }
                }

                _WeightTable[ u ].Left  = iLeft;
                _WeightTable[ u ].Right = iRight;

                int iSrc = 0;
                double dTotalWeight = 0;

                for( iSrc=iLeft; iSrc<=iRight; iSrc++ )
                {
                    const double weight = dFScale *
                                          pFilter->Filter( dFScale * (dCenter - (double)iSrc) );

                    if ( _WeightTable[ u ].Weights != NULL )
                    {
                        _WeightTable[ u ].Weights[ iSrc - iLeft ] = weight;
                    }
                    dTotalWeight += weight;
                }

                if( ( dTotalWeight > 0 ) && ( dTotalWeight != 1 ) )
                {
                    for( iSrc = iLeft; iSrc <= iRight; iSrc++ )
                    {
                        if ( _WeightTable[ u ].Weights != NULL )
                        {
                            _WeightTable[ u ].Weights[ iSrc-iLeft ] /= dTotalWeight;
                        }
                    }

                    iSrc = iRight - iLeft;

                    if ( _WeightTable[ u ].Weights != NULL )
                    {
                        while( _WeightTable[ u ].Weights[ iSrc ] == 0 )
                        {
                            _WeightTable[ u ].Right--;
                            iSrc--;

                            if( _WeightTable[ u ].Right == _WeightTable[ u ].Left )
                                break;
                        }
                    }
                }
            }
        } /// of if ( _WeightTable != NULL)
    } /// of if ( pFilter != NULL )
}

RawScaleWeightsTable::~RawScaleWeightsTable()
{
    if ( _WeightTable != NULL )
    {
        for( unsigned u = 0; u < _LineLength; u++ )
        {
            if ( _WeightTable[u].Weights != NULL )
            {
                delete[] _WeightTable[u].Weights;
            }
        }

        delete[] _WeightTable;
    }
}

double RawScaleWeightsTable::getWeight( unsigned dst_pos, unsigned src_pos )
{
    if ( dst_pos < _LineLength )
    {
        //int sz = _WeightTable[dst_pos].Right - _WeightTable[dst_pos].Left;
        if ( src_pos < _WindowSize )
        {
            return _WeightTable[dst_pos].Weights[src_pos];
        }
    }

    return 0.0;
}

unsigned RawScaleWeightsTable::getLeftBoundary( unsigned dst_pos )
{
    return _WeightTable[dst_pos].Left;
}

unsigned RawScaleWeightsTable::getRightBoundary( unsigned dst_pos )
{
    return _WeightTable[dst_pos].Right;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

RAWResizeEngine::RAWResizeEngine( RAWGenericFilter* filter )
 : _pFilter( filter )
{
}

unsigned RAWResizeEngine::scale( const unsigned short* src, unsigned src_width, unsigned src_height,
                                 unsigned dst_width, unsigned dst_height, unsigned short** dst )
{
	if ( src == NULL)
        return 0;

    if ( ( src_width == 0 ) || ( src_width == 0 ) || ( dst_width == 0 ) || ( dst_height == 0 ) )
        return 0;

    unsigned imgsz = dst_width * dst_height;

    if ( *dst != NULL )
    {
        delete[] *dst;
    }

    *dst = new unsigned short[ imgsz ];

    if ( *dst == NULL )
    {
        return 0;
    }

	if ( ( src_width == dst_width ) && ( src_height == dst_height ) )
	{
        if ( *dst != NULL )
        {
            unsigned cpsz = src_width * src_height * sizeof( unsigned short );
            memcpy( *dst, src, cpsz );
            return cpsz;
        }
	}

	if ( dst_width <= src_width )
	{
        unsigned short* tmp_buff = NULL;

        if ( src_width != dst_width )
        {
            if ( src_height != dst_height )
            {
                tmp_buff = new unsigned short[ dst_width * src_height ];

                if ( tmp_buff == NULL )
                {
                    delete[] *dst;
                    return 0;
                }
            }
            else
            {
                tmp_buff = *dst;
            }

            horizontalFilter( src, src_height, src_width,
                             0, 0, tmp_buff, dst_width );

        }
        else
        {
            tmp_buff = (unsigned short*)src;
        }

        if ( src_height != dst_height )
        {
            verticalFilter( tmp_buff, dst_width, src_height, 0, 0,
                            *dst, dst_width, dst_height );
        }

		if ( ( tmp_buff != src ) && ( tmp_buff != *dst ) )
		{
		    delete[] tmp_buff;
		    tmp_buff = NULL;
		}

	}
	else    /// == ( dst_width > src->w() )
	{
		unsigned short* tmp_buff = NULL;

		if ( src_height != dst_height )
		{
			if ( src_width != dst_width )
			{
			    tmp_buff = new unsigned short[ src_width * dst_height ];
				if ( tmp_buff == NULL )
				{
					delete[] dst;
					return 0;
				}
			}
			else
            {
				tmp_buff = *dst;
			}

			verticalFilter( src, src_width, src_height,
                            0, 0, tmp_buff, dst_width, dst_height );

		}
		else
		{
			tmp_buff = (unsigned short*)src;
		}

		if ( src_width != dst_width )
		{
			horizontalFilter( tmp_buff, dst_height, src_width,
                              0, 0, *dst, dst_width );
		}

		if ( ( tmp_buff != src ) && ( tmp_buff != *dst ) )
		{
			delete[] tmp_buff;
			tmp_buff = NULL;
		}
	}

	if ( dst != NULL )
	{
        return imgsz * sizeof( unsigned short );
	}

	return 0;
}

void RAWResizeEngine::horizontalFilter( const unsigned short* src, const unsigned height, const unsigned src_width,
                                        const unsigned src_offset_x, const unsigned src_offset_y, unsigned short* dst, const unsigned dst_width )
{
	// allocate and calculate the contributions
	RawScaleWeightsTable weightsTable( _pFilter, dst_width, src_width );

	unsigned y = 0;
	unsigned x = 0;
	unsigned i = 0;

	#pragma omp parallel for private(x,i)
    for ( y = 0; y < height; y++)
    {
        const
        unsigned short* src_bits = &src[ ( ( y + src_offset_y ) * src_width ) + src_offset_x  ];
        unsigned short* dst_bits = &dst[ y * dst_width ];

        // scale each row
        for ( x = 0; x < dst_width; x++)
        {
            // loop through row
            const unsigned          iLeft  = weightsTable.getLeftBoundary(x);			/// retrieve left boundary
            const unsigned          iLimit = weightsTable.getRightBoundary(x) - iLeft;	/// retrieve right boundary
            const unsigned short*   pixel  = src_bits + iLeft;
            double                  gray   = 0.0;

            // for(i = iLeft to iRight)
            for ( i = 0; i <= iLimit; i++)
            {
                // scan between boundaries
                // accumulate weighted effect of each neighboring pixel
                const double    weight = weightsTable.getWeight(x, i);
                double          dpixel = *pixel;

                gray += ( weight * dpixel );
                pixel++;
            }

            *dst_bits = (unsigned short)CLAMP<long long>((long long)(gray + 0.5f), 0, 0xFFFFFFFF);
            dst_bits++;
        }
    }
}

/// Performs vertical image filtering
void RAWResizeEngine::verticalFilter( const unsigned short* src, unsigned width, unsigned src_height, unsigned src_offset_x, unsigned src_offset_y,
                                      unsigned short* dst, const unsigned dst_width, unsigned dst_height)
{
	// allocate and calculate the contributions
	RawScaleWeightsTable weightsTable( _pFilter, dst_height, src_height );

    //unsigned dst_pitch = dst_width * src_bpp;
    unsigned        dst_pitch = width;
    unsigned short* dst_base  = dst;
    unsigned        src_pitch = width;

	unsigned y = 0;
	unsigned x = 0;
	unsigned i = 0;

    #pragma omp parallel for private(y,i)
    for ( x = 0; x < width; x++)
    {
        // work on column x in dst
        const unsigned  index = x;
        unsigned short* dst_bits = dst_base + index;

        // scale each column
        for ( y = 0; y < dst_height; y++)
        {
            const
            unsigned short* src_base  = &src[ ( src_offset_y * width ) +
                                              ( src_offset_y * src_pitch + src_offset_x ) ];
            // loop through column
            const unsigned          iLeft       = weightsTable.getLeftBoundary(y);				// retrieve left boundary
            const unsigned          iLimit      = weightsTable.getRightBoundary(y) - iLeft;	// retrieve right boundary
            const unsigned short*   src_bits    = src_base + ( iLeft * src_pitch + index );
            double                  gray        = 0.0;

            for ( i = 0; i <= iLimit; i++)
            {
                // scan between boundaries
                // accumulate weighted effect of each neighboring pixel
                const double    weight = weightsTable.getWeight(y, i);
                double          dpixel = *src_bits;

                gray += weight * dpixel;
                src_bits += src_pitch;
            }

            // clamp and place result in destination pixel
            *dst_bits = (unsigned short)CLAMP<long long>((long long) (gray + 0.5f), 0, 0xFFFFFFFF);
            dst_bits += dst_pitch;
        }
    }
}
