/*
Main function to parse the parameters and run the task. 

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include "TreeOfShapes.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <QImage>
#include <QColor>

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
    char * file_name = argv[1];
    int mode = atoi(argv[2]);
    int model = atoi(argv[3]);
    int alternative_model = atoi(argv[4]);

    std::stringstream ss(argv[5]);
    bool options_shapeabstraction;
    ss >> std::boolalpha >> options_shapeabstraction;
    int color_sketch_shapeabstraction=atoi(argv[6]);
    int minarea_shapeabstraction=atof(argv[7]);
    int maxarea_shapeabstraction=atof(argv[8]);
    int scaleratio_shapeabstraction=atoi(argv[9]);
    float threshold_shapeabstraction=atof(argv[10]);
    float eps_shapeabstraction=atof(argv[11]);
    int renderOrder_shapeabstraction=atoi(argv[12]);
    float alpha_shapeabstraction=atof(argv[13]);

    std::stringstream ss1(argv[14]);
    bool options_watercolor;
    ss1 >> std::boolalpha >> options_watercolor;
    int color_sketch_watercolor=atoi(argv[15]);
    int minarea_watercolor=atof(argv[16]);
    int maxarea_watercolor=atof(argv[17]);
    int scaleratio_watercolor=atoi(argv[18]);
    float threshold_watercolor=atof(argv[19]);
    float eps_watercolor=atof(argv[20]);
    int renderOrder_watercolor=atoi(argv[21]);
    float alpha_watercolor=atof(argv[22]);

    std::stringstream ss2(argv[23]);
    bool options_shaking;
    ss2 >> std::boolalpha >> options_shaking;
    int color_sketch_shaking=atoi(argv[24]);
    int minarea_shaking=atof(argv[25]);
    int maxarea_shaking=atof(argv[26]);
    int scaleratio_shaking=atoi(argv[27]);
    float threshold_shaking=atof(argv[28]);
    float eps_shaking=atof(argv[29]);
    int renderOrder_shaking=atoi(argv[30]);
    float alpha_shaking=atof(argv[31]);

    std::stringstream ss3(argv[32]);
    bool options_shapefiltering;
    ss3 >> std::boolalpha >> options_shapefiltering;
    int color_sketch_shapefiltering=atoi(argv[33]);
    int minarea_shapefiltering=atof(argv[34]);
    int maxarea_shapefiltering=atof(argv[35]);
    int scaleratio_shapefiltering=atoi(argv[36]);
    float threshold_shapefiltering=atof(argv[37]);
    float eps_shapefiltering=atof(argv[38]);
    int renderOrder_shapefiltering=atoi(argv[39]);
    float alpha_shapefiltering=atof(argv[40]);

    std::stringstream ss4(argv[41]);
    bool options_styletransfer;
    ss4 >> std::boolalpha >> options_styletransfer;
    int color_sketch_styletransfer=atoi(argv[42]);
    int minarea_styletransfer=atof(argv[43]);
    int maxarea_styletransfer=atof(argv[44]);
    int scaleratio_styletransfer=atoi(argv[45]);
    float threshold_styletransfer=atof(argv[46]);
    float eps_styletransfer=atof(argv[47]);
    int renderOrder_styletransfer=atoi(argv[48]);
    float alpha_styletransfer=atof(argv[49]);
    int modelDictionary = atoi(argv[50]); 
    int mcolor = atoi(argv[51]);
    int equal = atoi(argv[52]);
    float kappaDict = atof(argv[53]);

    char * dictionary_file_name = argv[54];
    char * mask_file_name = argv[55];  
    
    // Load Image
    QImage image(file_name);

    // Load Mask
    QImage image_mask(mask_file_name);
       
    TreeOfShapes * TOS = new TreeOfShapes(cfimages_from_qimage(image));

    QImage resulting_image;

    // Load parameters depending on the task
    TOSParameters TOSParameters = getAbstractionTOSParameters();
    
    if( mode==0 ){
        TOSParameters = getAbstractionTOSParameters();
        std::cout << "Shape abstraction " << std::endl;
        TOSParameters.color_sketch = color_sketch_shapeabstraction;
        if (options_shapeabstraction){
            TOSParameters.order = renderOrder_shapeabstraction;
            TOSParameters.alpha = alpha_shapeabstraction; 
            TOSParameters.minarea = minarea_shapeabstraction;
            TOSParameters.maxarea = maxarea_shapeabstraction;
            TOSParameters.eps = eps_shapeabstraction;
            TOSParameters.threshold = threshold_shapeabstraction;
            TOSParameters.ns = scaleratio_shapeabstraction;
        }
    } else if( mode==1 ){
        std::cout << "Watercolor " << std::endl;
        TOSParameters =  getWaterColorTOSParameters();
        TOSParameters.color_sketch = color_sketch_watercolor; 
        if (options_watercolor){
            TOSParameters.order = renderOrder_watercolor;
            TOSParameters.alpha = alpha_watercolor;
            TOSParameters.minarea = minarea_watercolor;
            TOSParameters.maxarea = maxarea_watercolor;
            TOSParameters.eps = eps_watercolor;
            TOSParameters.threshold = threshold_watercolor;
            TOSParameters.ns = scaleratio_watercolor;
        }
    } else if( mode==2 ){
        std::cout << "Shaking " << std::endl;
        TOSParameters =  getShapeShakingTOSParameters();
        TOSParameters.color_sketch = color_sketch_shaking; 
        if (options_shaking){
            TOSParameters.order = renderOrder_shaking;
            TOSParameters.alpha = alpha_shaking;
            TOSParameters.minarea = minarea_shaking;
            TOSParameters.maxarea = maxarea_shaking;
            TOSParameters.eps = eps_shaking;
            TOSParameters.threshold = threshold_shaking;
            TOSParameters.ns = scaleratio_shaking;
        }
    } else if( mode==3 ){
        std::cout << "Shape filtering " << std::endl;
        TOSParameters =  getShapeSmoothingTOSParameters();
        TOSParameters.color_sketch = color_sketch_shapefiltering; 
        if (options_shapefiltering){
            TOSParameters.order = renderOrder_shapefiltering;
            TOSParameters.alpha = alpha_shapefiltering;
            TOSParameters.minarea = minarea_shapefiltering;
            TOSParameters.maxarea = maxarea_shapefiltering;
            TOSParameters.eps = eps_shapefiltering;
            TOSParameters.threshold = threshold_shapefiltering;
            TOSParameters.ns = scaleratio_shapefiltering;
        }
    } else if( mode==4 ){
        std::cout << "Style transfer " << std::endl;
        TOSParameters =  getStyleTransferTOSParameters();
        model= 4;
        TOSParameters.color_sketch = color_sketch_styletransfer; 
        if (options_styletransfer){
            TOSParameters.order = renderOrder_styletransfer;
            TOSParameters.alpha = alpha_styletransfer;
            TOSParameters.minarea = minarea_styletransfer;
            TOSParameters.maxarea = maxarea_styletransfer;
            TOSParameters.eps = eps_styletransfer;
            TOSParameters.threshold = threshold_styletransfer;
            TOSParameters.ns = scaleratio_styletransfer;
        }
    };

    TOSParameters.model = model; 
    
    if (model == 4) { 
        // Load dictionary parameters
        DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters();

        if (options_styletransfer){
            dictionaryParameters.randS = modelDictionary;
            dictionaryParameters.mcolor = mcolor;
            dictionaryParameters.equal = equal; 
            dictionaryParameters.kappaDict = kappaDict; 
        };

        // Load dictionary of dictionary image
        QImage image_dict(dictionary_file_name);

        if (image_dict.isNull()){
           std::cout << "A Style image it is necessary" << std::endl; 
           std::ofstream demo_failure;
           demo_failure.open ("demo_failure.txt");
           demo_failure << "A Style image it is necessary.\n";
           demo_failure.close();
           return 0;
        };

        TreeOfShapes * dictionary = new TreeOfShapes(cfimages_from_qimage(image_dict));
        dictionary->compute_tree( getDefaultTOSParameters(), true);

        // Run abstraction
        resulting_image = TOS->render(TOSParameters, image_mask, alternative_model, dictionary, dictionaryParameters);
    } else {
        // Run abstraction
        resulting_image = TOS->render(TOSParameters, image_mask, alternative_model);
    };

    resulting_image.save("result.png");   
}
