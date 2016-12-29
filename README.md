# librawprocessor #
* A library for 16bit Gray scaled RAW image processor in GCC for any different system.
* This library is designed for GCC(MinGW) with FLTK 1.3.3 or 1.3.4.
* Reading 12 to 16bit RAW image and make its threshold cut-off image into 8bit another RAW image.

### This library supports these functions. ###
 1. RAW read by image size of height, width automatically calculated and padded with Zero.
 1. Get weight informations.
 1. Threshold cut-off (for lose useless image datas, or make it easy to processing)
 1. Export to 8bit pixel array(vector array)
 1. Export to 16bit pixel array(vector array)
 1. Reverse all pixels
 1. Flip H/V and Rotate in 90,180,270 degrees.

### And still these functions are not embodied. ###
 1. Nothing in now.

### Reserved to be included ###
* OpenMP

### Latest update ###
* Version marked as 0.9.26.90
* Now supporting Pixel filterings with :
  1. ApplyFilter();
  2. CloneWithFilter();
* Basically provides these filters:
  1. blur
  2. blur more
  3. edge
  4. edge more
  5. sharpen
  6. sharpen more
  7. sharpen subtle
  8. unsharpen mask (in development)
* User customizable filter availed, see WIKI.
* Median filter included:
  1. ApplyMedianFilter();

### Previous updates ###
* Enhanced for these methods:
  1. Get16bitThresholdedImage();
  1. Get8bitThresholdedImage();
* Changed some integer type to unsigned.
* Added function for Linear, Rectangle, Polygon regional histogram data.
  1. GetLinearPixels();
  1. GetRectPixels();
  1. GetPolygonPixels();
* Can getting information about average, variance, deviation from this method:
  1. GetAnalysisFromPixels();
* Fixed a bug Rescaled result applied again transform.
* Changed some method names:
  1. getWidth() -> Width()
  1. getHeight() -> Height()
  1. rescale() -> Rescale()
* Added method:
  1. SaveToFile()
* Reverse(Invert) all pixels with bit-range ( 9 to 16 )
* LoadingMatrix changed to Transform (unsigned int) definitions.
* Supporting these transforms :
  1. Flip horizontal.
  1. Flip vertical.
  1. Rotate clock wide 90,180,270
  1. Rotate reverse clock wide 90,180,270
* Included RAW image filtered scaling.
* Rescale method supports these: BiLinear/BiCubic/BSpline/Lanczos3
* Source refer to Free Image Project (http://freeimage.sourceforge.net/)
* Removed class methods not handles TCHAR casting.
* Compatible with 64bit gcc compilers.
* Some array struct enhanced for be optimized with AVX SIMD.</br>
  It may easily enhanced by compiler optimized option as like:
  ````
  -mavx (or) -march=core-i7avx
  ````

### Usage & Help ###
 1. [See wiki page](https://github.com/rageworx/librawprocessor/wiki)