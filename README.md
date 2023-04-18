# librawprocessor (float32 version)#
* A library for processing gray scaled medical RAW images.
* Supports multi-platforms: Windows, Linux and MacOS
* Reading 8 to 32bit any type of gray scaled RAW image and make its downscale window size into 8bit another RAW image, or reverse.

### This library supports these functions. ###
 1. RAW read by image size of height, width automatically calculated and padded with Zero.
 1. Get weight informations to calculating window size.
 1. Windowing - threshold cut-off (for lose useless image datas, or make it easy to processing)
 1. CLAHE( Contrast Limited Adaptive Histogram Equalization )
 1. Export to 8bit pixel array(vector array)
 1. Export to 16bit pixel array(vector array)
 1. Reverse all pixels ( black to white, or white to black )
 1. Flip H/V and Rotate in 90,180,270 degrees.
 1. 2 different ways for tone mapping (known as H.D.R image)

### Makefile rule ###
 1. Refer to HOW2MAKE.

### And still these functions are not embodied. ###
 1. Nothing in now.

### Reserved to be included ###
* Nothing now, bug fix if found.

### Latest update ###
* Version marked as 0.9.35.150.
* Included byte swapping

### Previous updates ###
* Fixed some bugs in AdjustingXXX methods.
* Enhanced CLAHE region divider, no failure at running, but result may affected.
* Fixed getting version information.
* Added RotateFree, RotateFree in cropped for original image size.
* Changed method name Reverse() to Invert()
* Changed method name ReverseAuto() to InvertAuto()
* Fixed a bug of Rotate180()
* Added Low Frequency Filter (may used for removing artifacts)
* Added Edge Enhancement Filter
* Added Anistropic Filter
* Optimized CLAHE algorithm works.
* Included prototype of CLAHE, not optimized.
* Fixed wrong range threshold method.
* Fixed some issues
  1. Tone mapping divide & multiply fixed different range.
* Supporting 2 different methods of Hight Dynamic Range tone mapping with this methods:
  1. bool AdjustToneMapping();
* Now supporting openmp, compiler may need support this.
* Added ReverseAuto();
* Fixed a bug in resize to same image.
* Fixed a bug of AdjustGamma();
* Supports Gamma, Brightness, Contrast with these methods:
  1. bool AdjustGamma( float gamma );
  2. bool AdjustBrightness( float percent );
  3. bool AdjustContrast( float percent );
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

### Usage & Help ###
 1. [See wiki page](https://github.com/rageworx/librawprocessor/wiki)
