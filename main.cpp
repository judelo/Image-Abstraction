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
    char * file_name = argv[1];               // Something like: "/mnt/data/lbouza/Image-Abstraction-Modif/bordeauxResize.jpg"
    char * mode_char = argv[2];               // Task Abstraction: 0; Watercolor:1; Shaking: 2; Shape smoothing:3; Style transfer:4;
    char * model_char = argv[3];              // Synthesis model: orignal shape: m=0; ellipse: m=1; rectangle: m=2; circle m=3,  dictionary m=4, random: m=5 (not use);
    std::stringstream ss(argv[4]);
    char * color_sketch_char = argv[5];       // Keep meaningful boundaries: No 0, Yes 1;
    char * renderOrder_char = argv[6];        //rendering order of the shapes: top->down: o=0 ; large->small: o=1; random: o=2"
    char * alpha_char = argv[7];              // Transparency (between 0 and 1)
    char * modelDictionary_char = argv[8];    // Selection model: randS=0, randomly select shapes; randS=1, select shapes according to elongation, compactness and scale; randS=2, select shapes according to elongation, compactness, scale and color",
    char * mcolor_char = argv[9];             // Select de source of color
    char * equal_char = argv[10];              // Scaling shape with equal aspect ratio or not
    char * kappaDict_char = argv[11];          // Compactness parameter of the attribute filtering on the transferred image
    char * mpixel_char = argv[12];             // minimal area (in pixel) for FLST
    char * maxarea_char = argv[13];             // maximal area (in pixel) for FLST
    char * ScaleRatio_char = argv[14];          // "scale ratio order for color filtering",
    char * Threshold_char = argv[15];           // "threshold for color filtering",
    char * eps_char = argv[16];                 // "-log10(max number of false alarms)",
    char * dictionary_file_name = argv[17];    // Something like: "/mnt/data/lbouza/Image-Abstraction-Modif/VanGogh.jpg"
    char * mask_file_name = argv[18];    

    int mode = atoi(mode_char);
    int model = atoi(model_char);
    int renderOrder = atoi(renderOrder_char);
    float alpha = atof(alpha_char);
    int color_sketch = atoi(color_sketch_char);
    int modelDictionary = atoi(modelDictionary_char); 
    int mcolor = atoi(mcolor_char);
    int equal = atoi(equal_char);
    float kappaDict = atof(kappaDict_char);
    int mpixel = atoi(mpixel_char);
    int maxarea = atoi(maxarea_char);
    int ScaleRatio = atoi(ScaleRatio_char);
    int Threshold = atoi(Threshold_char);
    int eps = atoi(eps_char);
    bool advanceOptions;
    ss >> std::boolalpha >> advanceOptions;
    
    // Load Image
    QImage image(file_name);

    // Load Mask
    QImage image_mask(mask_file_name);
       
    TreeOfShapes * TOS = new TreeOfShapes(cfimages_from_qimage(image));

    QImage resulting_image;
    QImage background;
    bool tree_recomputed = false;

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
        if (model!=4){
           std::cout << "For Style Transfer, model has to be dictionary." << std::endl;
           std::ofstream demo_failure("demo_failure.txt");
           demo_failure << "For Style Transfer, model has to be dictionary.\n";
           demo_failure.close();
           return 0;
        };
    };

    TOSParameters.model = model; 

    if (advanceOptions){
       TOSParameters.order = renderOrder;
       TOSParameters.alpha = alpha;
       TOSParameters.color_sketch = color_sketch; 
       TOSParameters.mpixel = mpixel;
       TOSParameters.maxarea = maxarea;
       TOSParameters.eps = eps;
       TOSParameters.threshold = Threshold;
       TOSParameters.ns = ScaleRatio;
    };
    
    if (model == 4) { 
        std::cout << "Model Dictionary" << std::endl; 

        if (mode!=4){
           std::cout << "Use Dictionary only for Style transfer." << std::endl;
           std::ofstream demo_failure("demo_failure.txt");
           demo_failure << "Use Dictionary only for Style transfer.\n";
           demo_failure.close();
           return 0;
        };

        // Load dictionary parameters
        DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters();

        if (advanceOptions){
            dictionaryParameters.randS = modelDictionary;
            dictionaryParameters.mcolor = mcolor;
            dictionaryParameters.equal = equal; 
            dictionaryParameters.kappaDict = kappaDict; 
            std::cout << "mcolor " << mcolor << std::endl;
        };

        // Load dictionary of dictionary image
        QImage image_dict(dictionary_file_name);

        if (image_dict.isNull()){
           std::cout << "An image for dictionary it is necessary" << std::endl; 
           std::ofstream demo_failure;
           demo_failure.open ("demo_failure.txt");
           demo_failure << "An image for dictionary it is necessary.\n";
           demo_failure.close();
           return 0;
           return 0;
        };

        TreeOfShapes * dictionary = new TreeOfShapes(cfimages_from_qimage(image_dict));
        dictionary->compute_tree( getDefaultTOSParameters(), true);

        // Run abstraction
        resulting_image = TOS->render(TOSParameters, tree_recomputed,  image_mask, background, dictionary, dictionaryParameters);
    } else {
        // Run abstraction
        background = TOS->renderOrigShapesBackground(TOSParameters, tree_recomputed, image_mask);
        resulting_image = TOS->render(TOSParameters, tree_recomputed, image_mask, background);
    };

    resulting_image.save("result.png");
    
}
