# librawprocessor #

A library for 16bit Gray scaled RAW image processor in GCC w/ FLTK 1.3.x, or any different system.

This library is designed for GCC(MinGW) with FLTK 1.3.x
Reading 12 to 14bit RAW image and make its threshold cut-off image into 8bit another RAW image.

### This library supports these functions. ###
 1. RAW read by image size of height.
 1. Get weight informations.
 1. Threshold cut-off (for lose useless image datas, or make it easy to processing)
 1. Export to 8bit image (vector array)
 1. Export to 16bit image (vector array)

### And still these functions are not embodied. ###
 1. Some Matrix to rotate, reverse image.
 1. ~~Export to a file.~~ 
    Emboded to other component, Fl_ImageViewer (for other project)
 1. ~~Colored RAW image.~~ 
    This library designed for medical raws.

### Usage ###
 1. Include and link (or compiled with in code) librawprocessor.
 1. Then Just declare what you want to use 'new' class.
~~~~
	RAWProcessor* rawproc = new RawProcessor( "test.raw", 1566 );
    if ( rawproc != NULL )
    {
    	if ( rawproc->isLoaded() == true )
    	{
        	// some works ...
            
            rawproc->Unload();
            delete rawproc;
            rawproc = NULL;
        }
    }
~~~~
 1. Or
~~~~
	RAWProcessor* rawproc = new RawProcessor();
    if ( rawproc != NULL )
    {
    	rawproc->Load( 'test.raw', RAWProcessor::LMATRIX_NONE, 1566 );
    	if ( rawproc->isLoaded() == true )
    	{
        	// some works ...

			rawproc->Unload();
            delete rawproc;
            rawproc = NULL;
        }
	}
~~~~
 1. Or from memory..
~~~~
	RAWProcessor* rawproc = new RawProcessor();
    if ( rawproc != NULL )
    {
    	rawproc->LoadFromMemory( refrawmem, refrawmemsz, RAWProcessor::LMATRIX_NONE, 1566 );
    	if ( rawproc->isLoaded() == true )
    	{
        	// some works ...

            delete rawproc;
            rawproc = NULL;
		}
	}
~~~~
 1. Extract raw image to thresholded (as 0 to 1204) 8 bit pixels with not reversed levels.
~~~~
	WeightAnalysisReport wreport;
    vector<unsigned char> downscaled;

    wreport.timestamp = GetTickCount(); /// it doesn't need but doing for debug.
    wreport.threshold_wide_min = 0;
    wreport.threshold_wide_max = 1204;

    rawproc->Get8bitThresholdedImage( wreport, &downscaled, false );
~~~~
