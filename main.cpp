/*
Main function to parse the parameters and run the task. 

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include "Segmentation.h"
#include "TreeOfShapes.h"
#include <iostream>
#include <QImage>

// Auxiliary Function to get cfimage from a QImage
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
    char * file_name = argv[1];          // Something like: "/mnt/data/lbouza/Image-Abstraction-Modif/bordeauxResize.jpg"
    char * mode_char = argv[2];         // Task Abstraction: 0; Watercolor:1; Shaking: 2; Shape smoothing:3; Style transfer:4;
    char * model_char = argv[3];        // Synthesis model: orignal shape: m=0; ellipse: m=1; rectangle: m=2; circle m=3,  dictionary m=4 (use just for styletransfer), random: m=5 (not use);
    char * seg_char = argv[4];          // Segmentation of input image: No 0, Yes 1;
    char * style_file_name = argv[5];    // Something like: "/mnt/data/lbouza/Image-Abstraction-Modif/VanGogh.jpg"
    
    int mode = atoi(mode_char);
    int model = atoi(model_char);
    int seg = atoi(seg_char);
    
    // Load Image
    QImage image(file_name);
    
    // Update image if segmentation is selected. 
    if (seg==1){
       Segmentation * segmentation = new Segmentation (image);
       image = segmentation->segment( 0.5, 500, 50); // segParameters.c = 500;  segParameters.min_size = 50; segParameters.sigma = 0.5;
       image.save("Segment.png");
       image = segmentation->getResult();
       image.save("After_Segment.png");
    };
       
    TreeOfShapes * TOS = new TreeOfShapes(cfimages_from_qimage(image));

    QImage resulting_image;
    bool tree_recomputed;

    // Load parameters depending on the task
    TOSParameters TOSParameters = getAbstractionTOSParameters();
    
    if( mode==0 ){
        TOSParameters = getAbstractionTOSParameters();
        std::cout << "Abstraction " << std::endl;
    } else if( mode==1 ){
        std::cout << "Watercolor " << std::endl;
        TOSParameters =  getWaterColorTOSParameters();
    } else if( mode==2 ){
        std::cout << "Shaking " << std::endl;
        TOSParameters =  getShapeShakingTOSParameters();
    } else if( mode==3 ){
        std::cout << "Shape smoothing " << std::endl;
        TOSParameters =  getShapeSmoothingTOSParameters();
    } else if( mode==4 ){
        std::cout << "Style transfer " << std::endl;
        TOSParameters =  getStyleTransferTOSParameters();
        DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters();

        // Load dictionary of Style image
        QImage image_dict(style_file_name);
        TreeOfShapes * dictionary = new TreeOfShapes(cfimages_from_qimage(image_dict));
        dictionary->compute_tree( getDefaultTOSParameters(), true);
        resulting_image = TOS->render(TOSParameters, tree_recomputed,  dictionary, dictionaryParameters);
    };
     
    if (mode!=4){
       // Select model
       TOSParameters.model = model;  
       resulting_image = TOS->render(TOSParameters, tree_recomputed);
    };
    
    resulting_image.save("result.png");
}
