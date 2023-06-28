# Image-Abstraction

## ABOUT

* IPOL website: TO_COMPLETE
* Version number: TO_COMPLETE
* Date: June 2023
* Author    : Noura Faraj <noura.faraj@umontpellier.fr>, Lucía Bouza <lucia.bouza-hegeurte@u-paris.fr> , Julie Delon <julie.delon@u-paris.fr>
* Copyright : (C) 2022 IPOL Image Processing On Line http://www.ipol.im/
* Licence   : GPL v3+, see GPLv3.txt

## OVERVIEW

This source code provides a C++ implementation of the framework for the structured abstraction of images.

### Contents

* main.cpp: Main Script to parse arguments and run the task. 
* Segmentation.h and cpp:  Class to implement segmentation on input image.
* TreeOfShapes.h and TreeOfShapes.cpp:  Class to implement abstraction process.
* src/tree_of_shapes.h and src/tree_of_shapes.cpp:  Class to implement initialization of shapes and default parameters.
* Image-Abstraction.pro: file used for compilation with qmake. 
* src/segment: Library used by the solution, implemented by Pedro Felzenszwalb.
* src/mw3, src/flst: Library used by the solution, implemented by other authors. [link](http://dev.ipol.im/git/nil/megawave.git/)
* src/kdtree: Library used by the solution, implemented by David M. Mount and Sunil Arya [link](https://www.cs.umd.edu/~mount/ANN/)

## USER GUIDE

The code as is runs in C++ with the following dependencies:

### Dependencies

* build-essentials
* qt4-qmake
* qt4-dev-tools


### Usage

#### 1. Compilation

```
cd ~/Image-Abstraction
<Path_to_qmake>/qmake -makefile .
make 
```

#### 2. Run

Executable: `./image_abstraction`

Parameters:
- 1: file name of image. 
- 2: Task to perform: Abstraction: 0; Watercolor: 1; Shaking: 2; Shape smoothing: 3; Style transfer: 4;
- 3: Synthesis model to use: orignal shapes: 0; ellipse: 1; rectangle: 2; circle 3;  dictionary: 4 (just for Style Transfer);
- 4: Segmentation of input image: No: 0; Yes: 1;
- 5: Alternative model: If we use mask to mix models, this parameter is use to choose the second model: orignal shapes: 0; ellipse: 1; rectangle: 2; circle 3;  dictionary: 4 (just for Style Transfer);
- 6 - 14: Advance options for shape abstraction task
- 15 - 23: Advanced options for watercolor task
- 24 - 32: Advanced options for shaking task
- 33 - 42: Advanced options for filtering task
- 42 - 55: Advanced options for style transfer task
- 56: Style image. In case we select the task Style Transfer (4), we need to indicate the file name of the Style image. Absolute or relative path can be used.
- 57: Mask: In case we will use mixed models or segmentation, we need to indicate the file name of the mask. 

#### 3. Examples

Apply Abstraction with ellipse to image input_0.png and default parameters.
`./image_abstraction input_0.png 0 1 false 0 false 1 0.01 100 3 0.5 0 1 0 false 0 0 100 3 0.6 0 0 0 false 0 0 100 3 0.6 0 0 0 false 0 0.1 100 3 0.8 0 0 0 false 0 0 100 3 0 0 1 0.2 2 1 0 0 0 input_1.png mask_0.png`

Apply Style Transfer to image input_0.png, with Style image input_1.png and default parameters. 
`./image_abstraction input_0.png 4 0 false 0 false 1 0.01 100 3 0.5 0 1 0 false 0 0 100 3 0.6 0 0 0 false 0 0 100 3 0.6 0 0 0 false 0 0.1 100 3 0.8 0 0 0 false 0 0 100 3 0 0 1 0.2 2 1 0 0 0 input_1.png mask_0.png`

Apply Watercolor to image bordeauxResize.jpg, leaving original shapes, and apply segmentation.
`./image_abstraction input_0.png 1 0 true 0 false 1 0.01 100 3 0.5 0 1 0 false 0 0 100 3 0.6 0 0 0 false 0 0 100 3 0.6 0 0 0 false 0 0.1 100 3 0.8 0 0 0 false 0 0 100 3 0 0 1 0.2 2 1 0 0 0 input_1.png mask_0.png`


## ABOUT THIS FILE

Copyright 2022 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.

