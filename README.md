# Image-Abstraction

## ABOUT

* IPOL website: TO_COMPLETE
* Version number: TO_COMPLETE
* Date: October 2022
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
* src/segment, src/mw3, src/kdtree and src/flst: Modules used by the solutions, implemented by other authors. 
* Image-Abstraction.pro: file used for compilation with qmake. 

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
- 1: file name of image. Can be an absolute path like **/mnt/data/lbouza/Image-Abstraction-Modif/bordeauxResize.jpg** or relative path like **bordeauxResize.jpg**
- 2: Task to perform: Abstraction: 0; Watercolor: 1; Shaking: 2; Shape smoothing: 3; Style transfer: 4;
- 3: Synthesis model to use: orignal shapes: 0; ellipse: 1; rectangle: 2; circle 3;  dictionary: 4 (use just for Style Transfer), random: 5 (not use);
- 4: Segmentation of input image: No: 0; Yes: 1;
- 5: Style image. In case we select the task Style Transfer (4), we need to indicate the file name of the Style image. Absolute or relative path can be used.

#### 3. Examples

Apply Abstraction with ellipse to image bordeauxResize.jpg, but not use Segmentation:
`./image_abstraction bordeauxResize.jpg 0 1 0`

Apply Style Transfer to image bordeauxResize.jpg, with Style image VanGogh.jpg
`./image_abstraction bordeauxResize.jpg 4 4 0 VanGogh.jpg`

Apply Watercolor to image bordeauxResize.jpg, leaving original shapes, and apply segmentation.
`./image_abstraction bordeauxResize.jpg 1 0 1`


## ABOUT THIS FILE

Copyright 2022 IPOL Image Processing On Line http://www.ipol.im/

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.

