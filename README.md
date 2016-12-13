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

### And still these functions are not embodied. ###
 1. Some Matrix to rotate, reverse image.
 1. ~~Export to a file.~~ 
    Emboded to other component, Fl_ImageViewer (for other project)

### Latest update ###
* Included RAW image filtered scaling.
* Source refer to Free Image Project (http://freeimage.sourceforge.net/)

### Previous updates ###
* Removed class methods not handles TCHAR casting.
* Compatible with 64bit gcc compilers.
* Some array struct enhanced for be optimized with AVX SIMD.</br>
  It may easily enhanced by compiler optimized option as like:
  ````
  -mavx (or) -march=core-i7avx
  ````

### Reserved to be included ###
* Brightness control.
* Filtered resize method, supports BiLinear/BiCubic/Lanczos3

### Usage & Help ###
 1. [See wiki page](https://github.com/rageworx/librawprocessor/wiki)