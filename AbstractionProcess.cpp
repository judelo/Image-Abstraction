/*
 Module to implement abstraction process.

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include "AbstractionProcess.h"
#include "tree_of_shapes.h"
#include <iostream>
#include <sys/time.h>
#include <QColor>

extern void flst();

AbstractionProcess::AbstractionProcess( std::string fileNameIn ){

    // read input image
    _imgin = cfimageread(fileNameIn.c_str());

    // set default input options
    _image_loaded = true;
    _tosParameters = getDefaultTOSParameters();
    _tree_computed = false;
}

Cfimage AbstractionProcess::cfimages_from_qimage( const QImage &input_image  ){

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

AbstractionProcess::AbstractionProcess( const QImage &imageIn ){

    // read input image
    _imgin = cfimages_from_qimage(imageIn);

    // set default input options
    _image_loaded = true;
    _tosParameters = getDefaultTOSParameters();
    _tree_computed = false;

    _treeOfShapes = new TreeOfShapes( cfimages_from_qimage(imageIn) );
}

void AbstractionProcess::init(Cfimage inputImg, Shapes &pTree){

    Shape pShape;
    Fimage imgIntensity;

    if  ( ((imgIntensity = mw_new_fimage()) == NULL) ||
          (mw_alloc_fimage(imgIntensity,inputImg->nrow,inputImg->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    if((pTree = mw_new_shapes()) == NULL)
        mwerror(FATAL, 1,
                "fgrain --> Not enough memory to allocate the tree of shapes");
    if  ( ((_NormOfDu = mw_new_fimage()) == NULL) ||
          (mw_alloc_fimage(_NormOfDu,inputImg->nrow,inputImg->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");
    /*=======================================*/
    /*==    Compute Intensity imag        ===*/
    /*=======================================*/
    for( int i= 0; i< inputImg->ncol; i++)
        for( int j= 0; j< inputImg->nrow; j++)
        {
            //imgIntensity->gray[j*imgin->ncol + i] = (imgin->blue[j*imgin->ncol + i] * 255.0 );
            imgIntensity->gray[j*inputImg->ncol + i] = (int)(inputImg->blue[j*inputImg->ncol + i]
                    + inputImg->red[j*inputImg->ncol + i]
                    + inputImg->green[j*inputImg->ncol + i])/3;
        }
    float fzero = 0.; int nsize = 3;
    fderiv(imgIntensity,NULL,NULL,NULL,NULL,NULL,NULL,_NormOfDu,NULL,&fzero,&nsize);

    /*=======================================*/
    /*=== Compute FLST on Intensity image ===*/
    /*=======================================*/
    int minArea = 5;
    flst(&minArea, imgIntensity, pTree);

    /*=======================================*/
    /*==    Initialization                ===*/
    /*=======================================*/
    shapeInitialize(pTree);

    /*=======================================*/
    /*==   Assign color to each shape     ===*/
    /*=======================================*/
    // assign color to each shape
    Shape * ppShapeOfPixel = pTree->smallest_shape;

    for(int i= 0; i< pTree->ncol*pTree->nrow; i++){
        pShape = ppShapeOfPixel[i];

        ((Info*)(pShape->data))->r += inputImg->red[i];
        ((Info*)(pShape->data))->g += inputImg->green[i];
        ((Info*)(pShape->data))->b += inputImg->blue [i];
        ((Info*)(pShape->data))->n_cc++;
    }

    for(int i= 0; i< pTree->nb_shapes; i++){
        pShape = pTree->the_shapes + i;

        ((Info*)(pShape->data))->r /= (float)((Info*)(pShape->data))->n_cc;
        ((Info*)(pShape->data))->g /= (float)((Info*)(pShape->data))->n_cc;
        ((Info*)(pShape->data))->b /= (float)((Info*)(pShape->data))->n_cc;
        //printf("%f, %f\n", ((Info*)(pShape->data))->b, pShape->value/255.0);
    }

    mw_delete_fimage(imgIntensity);
}

AbstractionProcess::~AbstractionProcess(){
    if(_tree_computed){
        mw_delete_shapes(_pTree);
    }
}

Cfimage AbstractionProcess::cfimageread(const char* name){

    QImage image(name);
    return cfimages_from_qimage(image);
}

/*
QImage AbstractionProcess::render(TOSParameters tosParameters, bool &tree_recomputed, DictionaryParameters dictionaryParameters, TreeOfShapes * dictionnary){

    Ccimage imgsyn =NULL;
    if( tosParameters.model == 4 )
        imgsyn = _treeOfShapes->render( tosParameters, tree_recomputed, dictionnary, dictionaryParameters );
    else
        imgsyn = _treeOfShapes->render( tosParameters, tree_recomputed );
    QImage result_image( QSize(imgsyn->ncol, imgsyn->nrow), QImage::Format_RGB32 );

    _tosParameters = tosParameters;

    for( int j= 0; j< imgsyn->nrow; j++)
        for( int i= 0; i< imgsyn->ncol; i++)
        {
            int comp = j*imgsyn->ncol + i;

            QColor color (imgsyn->red[comp], imgsyn->green[comp], imgsyn->blue[comp]);
            result_image.setPixel(i, j , qRgb(color.red(), color.green(), color.blue()));
        }

    if( imgsyn != NULL )
        mw_delete_ccimage(imgsyn);
    return result_image;
}
*/


