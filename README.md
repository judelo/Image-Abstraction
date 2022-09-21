# Image-Abstraction
The code depend only on Qt4.

## Compilation

cd ~/Image-Abstraction
/usr/lib/x86_64-linux-gnu/qt4/bin/qmake -makefile .
make 

## Run

Executable: ./image_abstraction
Parameters:
- 1: file name of image. Can be an absolute path like **/mnt/data/lbouza/Image-Abstraction-Modif/bordeauxResize.jpg** or relative path like **bordeauxResize.jpg**
- 2: Task to perform: Abstraction: 0; Watercolor: 1; Shaking: 2; Shape smoothing: 3; Style transfer: 4;
- 3: Synthesis model to use: orignal shapes: 0; ellipse: 1; rectangle: 2; circle 3;  dictionary: 4 (use just for Style Transfer), random: 5 (not use);
- 4: Segmentation of input image: No: 0; Yes: 1;
- 5: Style image. In case we select the task Style Transfer (4), we need to indicate the file name of the Style image. Absolute or relative path can be used.

## Examples

Apply Abstraction with ellipse to image bordeauxResize.jpg, but not use Segmentation:
./image_abstraction bordeauxResize.jpg 0 1 0

Apply Style Transfer to image bordeauxResize.jpg, with Style image VanGogh.jpg
./image_abstraction bordeauxResize.jpg 4 4 0 VanGogh.jpg

Apply Watercolor to image bordeauxResize.jpg, leaving original shapes, and apply segmentation.
./image_abstraction bordeauxResize.jpg 1 0 1

