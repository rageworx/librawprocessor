librawprocessor
===============

A library for 16bit Gray scaled RAW image processor in GCC w/ FLTK 1.3.x.

This library is designed for GCC(MinGW) with FLTK 1.3.x
Reading 12 to 14bit RAW image and make its threshold cut-off image into 8bit another RAW image.

This library supports these functions.
 1. RAW read by image size of height.
 2. Get weight informations.
 3. Threshold cut-off (for lose useless image datas, or make it easy to processing)
 2. Export to 8bit image (vector array)
 3. Export to 16bit image (vector array)

And still these functions are not embodied.
 1. Export to a file.
 2. Colored RAW image.
