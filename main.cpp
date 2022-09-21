/*
Main function to parse the parameters and run the task. 

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include "Segmentation.h"
//#include "AbstractionProcess.h"
#include "TreeOfShapes.h"
#include <iostream>
#include <QImage>


Cfimage cfimages_from_qimage( const QImage &input_image  ){

    int nx = input_image.width(),ny = input_image.height();
    Cfimage out = mw_change_cfimage(NULL,ny,nx);

    float * red = (float *)malloc(nx*ny*sizeof(float));
    float * green = (float *)malloc(nx*ny*sizeof(float));
    float * blue = (float *)malloc(nx*ny*sizeof(float));

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            QColor color = input_image.pixel( x, y );
            int comp = y * nx + x;
            red[comp]= color.red();
            green[comp]= color.green();
            blue[comp]= color.blue();
        }
    }
    free(out->red);
    free(out->green);
    free(out->blue);
    out->red = red;
    out->green = green;
    out->blue = blue;

    return out;
}

int main(int argc, char *argv[])
{
    
    // Read input parameters
    char * fileName = argv[1];    // something like: "/mnt/data/lbouza/Image-Abstraction-Modif/bordeauxResize.jpg"
    char * mode_char = argv[2];   //Task Abstraction: 0; Watercolor:1; Shaking: 2; Shape smoothing:3; Style transfer:4;
    char * model_char = argv[3]; //synthesis model: orignal shape: m=0; ellipse: m=1; rectangle: m=2; circle m=3,  dictionary m=4 (use just for styletransfer), random: m=5 (not use);
    char * seg_char = argv[4];   //Segmentation of input image: No 0, Yes 1;
    char * style_fileName = argv[5];   // something like: "/mnt/data/lbouza/Image-Abstraction-Modif/VanGogh.jpg"
    
    int mode = atoi(mode_char);
    int model = atoi(model_char);
    int seg = atoi(seg_char);
    
    // Load Image
    QImage image(fileName);
    
    // Update image if segmentation is selected. 
    // _segParameters.c = 500;  _segParameters.min_size = 50; _segParameters.sigma = 0.5;
    if (seg==1){
       Segmentation * segmentation = new Segmentation (image);
       image = segmentation->segment( 0.5, 500, 50);
       image.save("Segment.png");
       image = segmentation->getResult();
       image.save("After_Segment.png");
    };
    
    //AbstractionProcess * imageAbstractionProcess = new AbstractionProcess(image);    
    TreeOfShapes * TOS = new TreeOfShapes(cfimages_from_qimage(image));

    QImage resulting_image;
    bool tree_recomputed;

    // Load parameters depending on the task
    TOSParameters _TOSParameters = getAbstractionTOSParameters();
    
    if( mode==0 ){
        _TOSParameters = getAbstractionTOSParameters();
        std::cout << "Abstraction " << std::endl;
    } else if( mode==1 ){
        std::cout << "Watercolor " << std::endl;
        _TOSParameters =  getWaterColorTOSParameters();
    } else if( mode==2 ){
        std::cout << "Shaking " << std::endl;
        _TOSParameters =  getShapeShakingTOSParameters();
    } else if( mode==3 ){
        std::cout << "Shape smoothing " << std::endl;
        _TOSParameters =  getShapeSmoothingTOSParameters();
    } else if( mode==4 ){
        std::cout << "Style transfer " << std::endl;
        _TOSParameters =  getStyleTransferTOSParameters();
        DictionaryParameters _dictionaryParameters = getDefaultDictionaryParameters();

        // load dictionary of Style image
        QImage image_dict(style_fileName);
        TreeOfShapes * dictionary = new TreeOfShapes(cfimages_from_qimage(image_dict));
        dictionary->compute_tree( getDefaultTOSParameters(), true);
        //resulting_image = imageAbstractionProcess->render(_TOSParameters, tree_recomputed, _dictionaryParameters, dictionary);
        resulting_image = TOS->render(_TOSParameters, tree_recomputed, _dictionaryParameters, dictionary);
    };
     

    if (mode!=4){
       // Select model
       _TOSParameters.model = model;  
       //resulting_image = imageAbstractionProcess->render(_TOSParameters, tree_recomputed);
       resulting_image = TOS->render(_TOSParameters, tree_recomputed);
    };
    
    resulting_image.save("result.png");

}
