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
 1. [See wiki page](https://github.com/rageworx/librawprocessor/wiki)