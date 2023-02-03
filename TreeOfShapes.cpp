/*
 Class to implement abstraction process.

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#include "TreeOfShapes.h"
#include "tree_of_shapes.h"
#include <iostream>
#include <sys/time.h>
#include <QColor>
#include <queue>
#include <algorithm>
#include "flst_boundary.c"

extern void flst();
#define PI 3.1415926
#define _MIN(x, y) ( (x)<(y) ? (x) : (y) )
#define _MAX(x, y) ( (x)>(y) ? (x) : (y) )

int TreeOfShapes::_tree_count = 0;

TreeOfShapes::TreeOfShapes( Cfimage imageIn ){
    
    // Read input image
    _imgin = imageIn;

    // Set default input options   
    _tosParameters = getDefaultTOSParameters();
    _dictionaryParameters = getDefaultDictionaryParameters();
    _tree_computed = false;
    _tree_recomputed = false;
    _tree_id = _tree_count++;
    _large_to_small_index = NULL;
    _large_to_small_index_computed = false;
    _texture_image_loaded = false;
    _use_kdtree = false;
    _ArrayPixelsMask = (Point_plane) malloc(imageIn->ncol * imageIn->nrow * sizeof(struct point_plane));
    _len_ArrayPixelsMask = 0;
}


void TreeOfShapes::init(Cfimage inputImg, Shapes &pTree){

    std::cout << "TreeOfShapes::init(Cfimage inputImg, Shapes &pTree)" << std::endl;
    Shape pShape;
    Fimage imgIntensity;

    if  ( ((imgIntensity = mw_new_fimage()) == NULL) || (mw_alloc_fimage(imgIntensity,inputImg->nrow,inputImg->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");
    if((pTree = mw_new_shapes()) == NULL)
        mwerror(FATAL, 1, "fgrain --> Not enough memory to allocate the tree of shapes");
    if  ( ((_NormOfDu = mw_new_fimage()) == NULL) || (mw_alloc_fimage(_NormOfDu,inputImg->nrow,inputImg->ncol) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    // Compute Intensity image
    for( int i= 0; i< inputImg->ncol; i++)
        for( int j= 0; j< inputImg->nrow; j++){
            imgIntensity->gray[j*inputImg->ncol + i] = (int)(inputImg->blue[j*inputImg->ncol + i] 
                    + inputImg->red[j*inputImg->ncol + i] + inputImg->green[j*inputImg->ncol + i])/3;
        }

    float fzero = 0.; 
    int nsize = 3;
    fderiv(imgIntensity,NULL,NULL,NULL,NULL,NULL,NULL,_NormOfDu,NULL,&fzero,&nsize);
     
    // Compute FLST on Intensity image   
    int minArea = 5;
    flst(&minArea, imgIntensity, pTree);

    // Initialization                 
    shapeInitialize(pTree);  
     
    // Assign color to each shape
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
    }

    mw_delete_fimage(imgIntensity);
     
    // Compute shape attribute
    std::cout << "Compute shape attributes" << std::endl;
    compute_shape_attribute();
}


// Delete tree of shapes
TreeOfShapes::~TreeOfShapes(){
    if( _tree_computed ){
        mw_delete_shapes(_pTree);
        mw_delete_fimage(_NormOfDu);
        mw_delete_cfimage(_imgin);
        if( _large_to_small_index_computed )
            mw_delete_fsignal(_large_to_small_index);
        for (std::map<int, Fsignal>::iterator it = _dictionary_selections.begin(); it !=  _dictionary_selections.end(); ++it)
            mw_delete_fsignal( it->second );
        if( _texture_image_loaded )
            mw_delete_cfimage(_texture_image);
    }
}


// Compute the orientation elongation and engienvalues of shapes
void TreeOfShapes::shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k, float *pX0, float *pY0){

    float size, a11, a20, a02, x0, y0, sumx, sumy, lambda1, lambda2;
    int i;

    size = (float)pShape->area;
    sumx = 0;
    sumy = 0;

    for(i = 0; i< size; i++){
        sumx += (double)((pShape->pixels+i)->x);
        sumy += (double)((pShape->pixels+i)->y);
    }
    x0 = sumx /(float)size;
    y0 = sumy /(float)size;

    a11 = 0;
    a20 = 0;
    a02 = 0;
    for(i = 0; i< size; i++){
        a11 += (float) ((pShape->pixels+i)->x - x0)*((pShape->pixels+i)->y - y0);
        a20 += (float) ((pShape->pixels+i)->x - x0)*((pShape->pixels+i)->x - x0);
        a02 += (float) ((pShape->pixels+i)->y - y0)*((pShape->pixels+i)->y - y0);
    }
    a11 = a11 /(float)pow(size,1);
    a20 = a20 /(float)pow(size,1)+ 1.0/12.0;
    a02 = a02 /(float)pow(size,1)+ 1.0/12.0;

    lambda1 =0.5*( a02 + a20 + sqrt((a20-a02)*(a20-a02) + 4*a11*a11));
    lambda2 =0.5*( a02 + a20 - sqrt((a20-a02)*(a20-a02) + 4*a11*a11));
    
    // Compute the orientation and engienvalues of shapes
    *out_ori = (0.5*atan2(2*a11,(a20-a02)));
    *out_e = lambda1;
    *out_k = lambda2;

    *pX0 = x0;
    *pY0 = y0;
}


// Index the mn-order parent of the pShape
Shape TreeOfShapes::m_order_parent(Shape pShape, int mn, bool dict){
    Shape pShapeTemp;
    pShapeTemp = pShape;

    for(int t=0; t<mn; t++){
        if(pShapeTemp->parent == NULL)
            break; 
        if(dict)
            for(pShapeTemp=pShape->parent; ((Info*)(pShapeTemp->data))->show != 1 && pShapeTemp !=NULL; pShapeTemp=pShapeTemp->parent){;}
        else 
            pShapeTemp = pShapeTemp->parent;
    }

    return pShapeTemp;
}


// Index the tree by the Breadth-first order
void TreeOfShapes::top2bottom_index_tree(Fsignal t2b_index){
    int queArr[_pTree->nb_shapes];
    int i, qInd, *head, *rear;
    Shape pShape, pShapeTemp;

    mw_change_fsignal(t2b_index, _pTree->nb_shapes);

    for(i=0; i< _pTree->nb_shapes; i++) {
        pShape = _pTree->the_shapes + i;
        ((Info*)(pShape->data))->index = i;
        queArr[i] = -i;
    }

    head = queArr;
    rear = queArr;
    queArr[0] = 0;
    rear++;
    t2b_index->values[0] = 0;
    i=0; 
    qInd =1;

    while(head != rear){
        t2b_index->values[i++] = (float) *head;
        pShape = _pTree->the_shapes + *head;
        head++;

        for(pShapeTemp = pShape->child; pShapeTemp != NULL; pShapeTemp = pShapeTemp->next_sibling){
            queArr[qInd++] = (((Info*)(pShapeTemp->data))->index);
            rear++;
        }
    }
}


// Compute boundingbox of each shape on the tree                                   
void TreeOfShapes::shape_boundingbox(Shape pShape){

    float xmin, xmax, ymin, ymax, theta, x0temp, y0temp, x, y, xr, yr;
    int size, ncol, nrow;

    ncol = _pTree->ncol;
    nrow = _pTree->nrow;
    size = pShape->area;
    theta = ((Info*)(pShape->data))->oren;
    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;
    xmin = ncol;  xmax = -ncol;
    ymin = nrow;  ymax = -nrow;

    for(int i = 0; i< size; i++){
        x = (float)((pShape->pixels+i)->x - x0temp);
        y = (float)((pShape->pixels+i)->y - y0temp);

        // Scaling_rotation_dict2(&x, &y, &theta, &xr, &yr);
        xr = (int)( x*cos(theta) + y*sin(theta));
        yr = (int)( y*cos(theta) - x*sin(theta));

        xmin = _MIN(xr, xmin);
        ymin = _MIN(yr, ymin);
        xmax = _MAX(xr, xmax);
        ymax = _MAX(yr, ymax);
    }

    ((Info*)(pShape->data))->boundingbox[0] = (int) xmin;
    ((Info*)(pShape->data))->boundingbox[1] = (int) xmax;
    ((Info*)(pShape->data))->boundingbox[2] = (int) ymin;
    ((Info*)(pShape->data))->boundingbox[3] = (int) ymax;
}


// Compute boundingbox of each shape on the tree  
void TreeOfShapes::tree_boundingbox(){
    Shape pShape;

    for(int i = _pTree->nb_shapes-1; i>=0; i--){
        pShape = _pTree->the_shapes + i;
        shape_boundingbox(pShape);
    }
}


// Compute shape attributes for the tree of shapes: orientation, engienvalues of shapes, color, contrast.  
void TreeOfShapes::compute_shape_attribute(){
    float oren, lamb1, lamb2, x0, y0, length;
    Shape pShape;

    _average_r = 0.;
    _average_g = 0.;
    _average_b = 0.;
    _maxArea = 0;

    for(int i = _pTree->nb_shapes-1; i>=0; i--){
        pShape = _pTree->the_shapes + i;

        shape_orilam(pShape, &oren, &lamb1, &lamb2, &x0, &y0);

        ((Info*)(pShape->data))->lambda1 = lamb1;
        ((Info*)(pShape->data))->lambda2 = lamb2;
        ((Info*)(pShape->data))->oren = oren;
        ((Info*)(pShape->data))->x0 = x0;
        ((Info*)(pShape->data))->y0 = y0;

        pShape->removed = 0;

        _average_r += ((Info*)(pShape->data))->r ;
        _average_g += ((Info*)(pShape->data))->g ;
        _average_b += ((Info*)(pShape->data))->b ;

        if(i != 0){
            Flist pBoundary = NULL;
            pBoundary = mw_change_flist(pBoundary, 4*pShape->area+1, 0, 2);
            flst_boundary(_pTree, pShape, pBoundary);
            ((Info*)(pShape->data))->contrast = fabs(min_contrast(pBoundary,&length,_NormOfDu));
            _maxArea = std::max( pShape->area, _maxArea );
        }
    }

    _average_r /= _pTree->nb_shapes;
    _average_g /= _pTree->nb_shapes;
    _average_b /=  _pTree->nb_shapes;
}


// Compute Median filter and Gaussian Blur
void TreeOfShapes::MedianFilterAndGaussianBlur(float left, float right, float top, float bottom, 
                                               Cimage imgShapeLabelSyn,Fimage imgShapeBlurSyn,
                                               Fsignal gaussKernel, int *median){

    int x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain;

    // Median Filter  
    MedSize = (int)((*median)/2.0);
    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++){
            numMedain = 0;
            for(iKer = - MedSize; iKer <= MedSize; iKer++)
                for(jKer = - MedSize; jKer <= MedSize; jKer++){
                    xKer = x + iKer;
                    yKer = y + jKer;

                    if(xKer<0 || xKer>= imgShapeLabelSyn->ncol || yKer<0 || yKer>= imgShapeLabelSyn->nrow )
                        continue;

                    imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] += (float)(imgShapeLabelSyn->gray[yKer*imgShapeLabelSyn->ncol + xKer]);
                    numMedain++;
                }
            if( imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] < ((float) numMedain)/2.0 )
                imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 0.0;
            else
                imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 1.0;
        }

    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++){
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = (int) imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];
            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] = 0.0;
        }

    // Add Gaussian Blur  
    KerSize = (int) ( sqrt( (double) gaussKernel->size) /2.0 );
    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++){
            for(iKer = -KerSize; iKer <= KerSize; iKer++)
                for(jKer = -KerSize; jKer <= KerSize; jKer++){
                    xKer = x + iKer;
                    yKer = y + jKer;

                    if(xKer<0 || xKer>= imgShapeLabelSyn->ncol || yKer<0 || yKer>= imgShapeLabelSyn->nrow )
                        continue;

                    imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] += 
                    gaussKernel->values[(iKer + KerSize)*KerSize + (jKer + KerSize)]*
                            (float)(imgShapeLabelSyn->gray[yKer*imgShapeLabelSyn->ncol + xKer]);
                }
       }
} 

// Synthesis by Shape Shaking for Ellipses, Rectangles or Circular shapes
// Before the Shaking, smooth the shape with a gaussian kernel or a median filter  
void TreeOfShapes::synshape(int model, Shape pShape,
                                  Ccimage imgsyn,
                                  Cimage imgShapeLabelSyn,
                                  Fimage imgShapeBlurSyn,
                                  Fsignal gaussKernel,
                                  int *median,
                                  float *alpha){

    int xi, yi, x, y;
    float ALPHA, BETA, a, b, x0temp, y0temp, top, right, left, bottom, phi, xi_e, yi_e;
    float xShift, yShift, theta, tR, tG, tB, TR, TG, TB, tr, tg, tb;
    bool condition;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    x0temp += xShift;
    y0temp += yShift;
    phi = ((Info*)(pShape->data))->oren;

    // Compute limits of shape before shaking
    if (model == 1){ //Ellipse
        a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
        b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
        left   = _MAX(0, x0temp - a);
        right  = _MIN(_pTree->ncol -1, x0temp + a);
        top    = _MAX(0, y0temp - a);
        bottom = _MIN(_pTree->nrow - 1, y0temp + a);
    } else if (model == 2){ //Rectangle
        a = (sqrt(3.0 * ((Info*)(pShape->data))->lambda1));
        b = (sqrt(3.0 * ((Info*)(pShape->data))->lambda2));
        float bLimit = sqrt(2.0)*a;
        left   = _MAX(0, x0temp - bLimit);
        right  = _MIN(_pTree->ncol -1, x0temp + bLimit);
        top    = _MAX(0, y0temp - bLimit);
        bottom = _MIN(_pTree->nrow - 1, y0temp + bLimit);
    }else if (model == 3){ //Circle
        a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
        b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
        left   = _MAX(0, x0temp - b);
        right  = _MIN(_pTree->ncol -1, x0temp + b);
        top    = _MAX(0, y0temp - a);
        bottom = _MIN(_pTree->nrow - 1, y0temp + b);
    }

    TR  = ((Info*)(pShape->data))->r;
    TG  = ((Info*)(pShape->data))->g;
    TB  = ((Info*)(pShape->data))->b;

    for( xi= ceil(left); xi<= right; xi++)
        for( yi= ceil(top); yi<= bottom; yi++){
            xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
            yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);

            if (model == 1) //Ellipse
                condition = ( xi_e*xi_e/(a*a) + yi_e*yi_e/(b*b) <= 1 );
            else if (model == 2){ //Rectangle
                condition = ( xi_e >= -a && xi_e <= +a && yi_e >= -b && yi_e <= +b );
            }else if (model == 3){ //Circle
                condition = ( xi_e*xi_e/(b*b) + yi_e*yi_e/(b*b) <= 1 );
            };

            if( condition ){
                if ((xi<0 || xi>= imgsyn->ncol || yi<0 || yi>= imgsyn->nrow ))
                    continue;

                imgShapeLabelSyn->gray[yi*imgShapeLabelSyn->ncol + xi] = 1;

                left   = _MIN(xi, left);
                top    = _MIN(yi, top);
                right  = _MAX(xi, right);
                bottom = _MAX(yi, bottom);
            }
        }

    MedianFilterAndGaussianBlur(left, right, top, bottom, imgShapeLabelSyn,imgShapeBlurSyn,gaussKernel, median);
    
    // Synthesis  
    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++){
            if(imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] == 0)
                continue;

            BETA = imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            tr = ((float) imgsyn->red[y*imgsyn->ncol + x])*(1-BETA)   + BETA*TR;
            tg = ((float) imgsyn->green[y*imgsyn->ncol + x])*(1-BETA) + BETA*TG;
            tb = ((float) imgsyn->blue[y*imgsyn->ncol + x])*(1-BETA)  + BETA*TB;

            tR = ((float) imgsyn->red[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tr;
            imgsyn->red[y*imgsyn->ncol + x]   = (int) tR;
            tG = ((float) imgsyn->green[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tg;
            imgsyn->green[y*imgsyn->ncol + x] = (int) tG; 
            tB = ((float) imgsyn->blue[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tb;
            imgsyn->blue[y*imgsyn->ncol + x]  = (int) tB;

            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x]  = 0.0;
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 0;
        }
}


// Synthesis by Shape Shaking for Original, Ellipses, Rectangles or Circular shapes
void TreeOfShapes::synshape(int model, Shape pShape,Ccimage imgsyn, float *alpha){
    
    int xi, yi, i, aux;
    float a, b, x0temp, y0temp, top, right, left, bottom, ALPHA, shLambda, shTR, shTG, shTB;
    float phi, x, y, xr, yr, xi_e, yi_e, xShift, yShift, theta, tR, tG, tB;
    bool condition;
    shLambda = 0.3;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    x0temp += xShift;
    y0temp += yShift;
    phi = ((Info*)(pShape->data))->oren;

    // Compute limits of shape before shaking
    if (model == 0){ //Original. We compute the actual point in for loop. 
        left   = 0;
        right  = _pTree->ncol -1;
        top    = 0;
        bottom = _pTree->nrow - 1;
        x0temp = (((Info*)(pShape->data))->x0);
        y0temp = (((Info*)(pShape->data))->y0);
    } else if (model == 1){ //Ellipse
        a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
        b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
        left   = _MAX(0, x0temp - a);
        right  = _MIN(_pTree->ncol -1, x0temp + a);
        top    = _MAX(0, y0temp - a);
        bottom = _MIN(_pTree->nrow - 1, y0temp + a);
    } else if (model == 2){ //Rectangle
        a = (sqrt(3.0 * ((Info*)(pShape->data))->lambda1));
        b = (sqrt(3.0 * ((Info*)(pShape->data))->lambda2));
        float bLimit = sqrt(2.0)*a;
        left   = _MAX(0, x0temp - bLimit);
        right  = _MIN(_pTree->ncol -1, x0temp + bLimit);
        top    = _MAX(0, y0temp - bLimit);
        bottom = _MIN(_pTree->nrow - 1, y0temp + bLimit);
    }else if (model == 3){ //Circle
        a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
        b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
        left   = _MAX(0, x0temp - b);
        right  = _MIN(_pTree->ncol -1, x0temp + b);
        top    = _MAX(0, y0temp - a);
        bottom = _MIN(_pTree->nrow - 1, y0temp + b);
    }

    shTR  = shLambda*((Info*)(pShape->data))->r;
    shTG  = shLambda*((Info*)(pShape->data))->g;
    shTB  = shLambda*((Info*)(pShape->data))->b;

    // Synthesis  
    i = 0;
    for( xi= ceil(left); xi<= right; xi++){
        for( yi= ceil(top); yi<= bottom; yi++){

            if (model == 0){ //Original shapes
                x = (float)((pShape->pixels+i)->x);
                y = (float)((pShape->pixels+i)->y);

                xr = (x - x0temp)*cos(theta) + (y - y0temp)*sin(theta);
                yr = (y - y0temp)*cos(theta) - (x - x0temp)*sin(theta);
                // Point of shape: (xi_e, yi_e)
                xi_e = floor(xShift + x0temp + xr);
                yi_e = floor(yShift + y0temp + yr);

                condition = xi_e>= 0 && xi_e< _pTree->ncol && yi_e>= 0 && yi_e< _pTree->nrow; 
                i++;
                if (i == pShape->area)
                    break; 
                
            } else {
                xi_e = ((float)xi - x0temp)*cos(phi+theta) + ((float)yi - y0temp)*sin(phi+theta);
                yi_e = ((float)yi - y0temp)*cos(phi+theta) - ((float)xi - x0temp)*sin(phi+theta);
                
                if (model == 1) //Ellipse
                    condition = ( xi_e*xi_e/(a*a) + yi_e*yi_e/(b*b) <= 1 );
                else if (model == 2) //Rectangle
                    condition = ( xi_e >= -a && xi_e <= +a && yi_e >= -b && yi_e <= +b );
                else if (model == 3) //Circle
                    condition = ( xi_e*xi_e/(b*b) + yi_e*yi_e/(b*b) <= 1 );
            }

            if(condition){
                if (model ==0)
                   aux = yi_e * _pTree->ncol + xi_e;
                else
                   aux = yi*_pTree->ncol + xi;

                tR = ((float) imgsyn->red[aux])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->r;
                imgsyn->red[aux] = (int) tR; 
                tG = ((float) imgsyn->green[aux])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->g;
                imgsyn->green[aux] = (int) tG; 
                tB = ((float) imgsyn->blue[aux])*ALPHA + (1-ALPHA)*((Info*)(pShape->data))->b;
                imgsyn->blue[aux] = (int) tB;
            }
        }
        if ((i == pShape->area) && (model == 0))
            break; 
    }
}

// Synthesis of shape using model dictionary and blur
void TreeOfShapes::synShapeDict(Shape pShapeDict, Shape pShape,
                                Ccimage imgsyn,Cfimage imgDict, 
                                Cfimage imgShapeColorSyn,Cimage imgShapeLabelDict, Cimage imgShapeLabelSyn,Fimage imgShapeBlurSyn,
                                Fsignal gaussKernel,int *median,
                                float *alpha,
                                int *equal, int *mcolor){

    int i, x, y, iKer, jKer, KerSize, MedSize, xKer, yKer, numMedain, index;
    float xi, yi, xr, yr, xt, yt, x0temp, y0temp, x0tempDict, y0tempDict;
    float tR, tG, tB, ALPHA, tr, tg, tb, BETA, TR, TG, TB;
    float xBound_l, xBound_r, yBound_t, yBound_b;
    float SCALEx, SCALEy, theta, thetaDict;
    float top, right, left, bottom;
    float a, b, bLimit, tempx;
    ALPHA = *alpha;

    xBound_l = (int) ((Info*)(pShapeDict->data))->boundingbox[0];
    xBound_r = (int) ((Info*)(pShapeDict->data))->boundingbox[1];
    yBound_t = (int) ((Info*)(pShapeDict->data))->boundingbox[2];
    yBound_b = (int) ((Info*)(pShapeDict->data))->boundingbox[3];

    // Begin Label pShapeDict 
    for(i=0; i < pShapeDict->area; i++){
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabelDict->gray[y*imgShapeLabelDict->ncol + x] = 1;
    }

    if(*equal == 1){
        SCALEx = sqrt( (double) pShape->area / (double) pShapeDict->area );
        SCALEy = SCALEx;
    } else{
        SCALEx =  sqrt( ((Info*)(pShape->data))->lambda1) /  sqrt(((Info*)(pShapeDict->data))->lambda1 );
        SCALEy =  sqrt( ((Info*)(pShape->data))->lambda2) /  sqrt(((Info*)(pShapeDict->data))->lambda2 );
    }

    a = SCALEx * (xBound_r - xBound_l);
    b = SCALEy * (yBound_b - yBound_t);

    bLimit = (1.0*sqrt(a*a + b*b)/2.0);

    x0temp = ((Info*)(pShape->data))->x0;
    y0temp = ((Info*)(pShape->data))->y0;
    x0tempDict =  ((Info*)(pShapeDict->data))->x0;
    y0tempDict =  ((Info*)(pShapeDict->data))->y0;

    theta  = ((Info*)(pShape->data))->oren;
    thetaDict = ((Info*)(pShapeDict->data))->oren;

    left   = _MAX(0, x0temp - bLimit);
    right  = _MIN(imgsyn->ncol - 1, x0temp + bLimit );
    top    = _MAX(0, y0temp - bLimit);
    bottom = _MIN(imgsyn->nrow - 1, y0temp + bLimit);

    // Transformation  
    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++){
            index = y*imgShapeLabelSyn->ncol + x;

            xt = (x - x0temp);
            yt = (y - y0temp);

            // Rotate to the major axis for scaling
            xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
            yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

            // Inverse-transformation
            xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
            yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);

            if(xr <0 || xr >=imgShapeLabelDict->ncol || yr <0 || yr >=imgShapeLabelDict->nrow)
                continue;
            
            if( imgShapeLabelDict->gray[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)] == 1 ){
                imgShapeLabelSyn->gray[index]  = 1;
                imgShapeColorSyn->red[index]   =  imgDict->red[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
                imgShapeColorSyn->green[index] =  imgDict->green[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
                imgShapeColorSyn->blue[index]  =  imgDict->blue[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
            }
            if(imgShapeLabelSyn->gray[index] == 0){     
                imgShapeColorSyn->red[index]   =  imgDict->red[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
                imgShapeColorSyn->green[index] =  imgDict->green[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
                imgShapeColorSyn->blue[index]  = imgDict->blue[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
            }
        }

    MedianFilterAndGaussianBlur(left, right, top, bottom, imgShapeLabelSyn,imgShapeBlurSyn,gaussKernel, median);

    std::cout <<"b"<< std::endl;

    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++){
            index = y*imgShapeLabelSyn->ncol + x;

            if(imgShapeLabelSyn->gray[index] == 0 ){
                xt = (x - x0temp);
                yt = (y - y0temp);

                // Rotate to the major axis for scaling
                xi = ( (1/SCALEx)* (xt*cos(theta) + yt*sin(theta)));
                yi = ( (1/SCALEy)* (yt*cos(theta) - xt*sin(theta)));

                // Inverse-transformation
                xr = ( (xi*cos(thetaDict) - yi*sin(thetaDict)) + x0tempDict);
                yr = ( (yi*cos(thetaDict) + xi*sin(thetaDict)) + y0tempDict);

                if(xr <0 || xr >=imgShapeLabelDict->ncol || yr <0 || yr >=imgShapeLabelDict->nrow)
                   continue;

                imgShapeColorSyn->red[index]   =  imgDict->red[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
                imgShapeColorSyn->green[index] =  imgDict->green[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
                imgShapeColorSyn->blue[index]  =  imgDict->blue[(int)(yr)*imgShapeLabelDict->ncol + (int)(xr)];
            }
        }

    std::cout <<"c"<< std::endl;

    if(*mcolor == 2){
        TR = ((Info*)(pShapeDict->data))->r;
        TG = ((Info*)(pShapeDict->data))->g;
        TB = ((Info*)(pShapeDict->data))->b;
    }

    for(x = ceil(left); x <= right; x++)
        for(y = ceil(top); y <= bottom; y++){
            if(imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] == 0)
                continue;

            BETA = imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            if(*mcolor == 1){
                TR  = imgShapeColorSyn->red[y*imgShapeLabelSyn->ncol + x];
                TG  = imgShapeColorSyn->green[y*imgShapeLabelSyn->ncol + x];
                TB  = imgShapeColorSyn->blue[y*imgShapeLabelSyn->ncol + x];
            }

            tr = ((float) imgsyn->red[y*imgsyn->ncol + x])*(1-BETA)   + BETA*TR;
            tg = ((float) imgsyn->green[y*imgsyn->ncol + x])*(1-BETA) + BETA*TG;
            tb = ((float) imgsyn->blue[y*imgsyn->ncol + x])*(1-BETA)  + BETA*TB;

            tR = ((float) imgsyn->red[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tr;
            imgsyn->red[y*imgsyn->ncol + x]   = (int) tR; 
            tG = ((float) imgsyn->green[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tg;
            imgsyn->green[y*imgsyn->ncol + x] = (int) tG;  
            tB = ((float) imgsyn->blue[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tb;
            imgsyn->blue[y*imgsyn->ncol + x]  = (int) tB;  
            
            // Reset indicator images
            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x]   = 0.0;
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 0;
        }
    
    // Reset indicator image
    for(i=0; i < pShapeDict->area; i++){
        x = (pShapeDict->pixels+i)->x;
        y = (pShapeDict->pixels+i)->y;
        imgShapeLabelDict->gray[y*imgShapeLabelDict->ncol + x] = 0;
    }
}

// Synthesis by Shape Shaking witn original shapes
// Before the Shaking, smooth the shape with a gaussian kernel or a median filter  
void TreeOfShapes::synshapeOriginal(Shape pShape,
                                    Ccimage imgsyn,
                                    Cimage imgShapeLabelSyn,
                                    Fimage imgShapeBlurSyn,
                                    Fsignal gaussKernel,
                                    int *median,
                                    float *alpha){

    int i, xi, yi, x, y;
    float xr, yr, x0temp, y0temp, ALPHA, BETA;
    float xShift, yShift, theta, tR, tG, tB, tr, tg, tb, TR, TG, TB;
    float top, right, left, bottom;

    ALPHA = *alpha;
    x0temp = (((Info*)(pShape->data))->x0);
    y0temp = (((Info*)(pShape->data))->y0);

    xShift = (((Info*)(pShape->data))->xShift);
    yShift = (((Info*)(pShape->data))->yShift);
    theta  = (((Info*)(pShape->data))->rotation);

    TR  = ((Info*)(pShape->data))->r;
    TG  = ((Info*)(pShape->data))->g;
    TB  = ((Info*)(pShape->data))->b;

    right = 0.0;  bottom= 0.0;
    left = (float)(imgShapeLabelSyn->ncol-1);
    top  = (float)(imgShapeLabelSyn->nrow-1);
    for(i=0; i< pShape->area; i++){
        x = (pShape->pixels+i)->x;
        y = (pShape->pixels+i)->y;

        xr = (x - x0temp)*cos(theta) + (y - y0temp)*sin(theta);
        yr = (y - y0temp)*cos(theta) - (x - x0temp)*sin(theta);

        xi = floor(xShift + x0temp + xr);
        yi = floor(yShift + y0temp + yr);

        if(xi<0 || xi>= imgShapeLabelSyn->ncol || yi<0 || yi>= imgShapeLabelSyn->nrow )
            continue;

        imgShapeLabelSyn->gray[yi*imgShapeLabelSyn->ncol + xi] = 1;

        left   = _MIN(xi, left);
        top    = _MIN(yi, top);
        right  = _MAX(xi, right);
        bottom = _MAX(yi, bottom);
    }

    MedianFilterAndGaussianBlur(left, right, top, bottom, imgShapeLabelSyn,imgShapeBlurSyn,gaussKernel, median);

    // Synthesis  
    for(x = left; x <= right; x++)
        for(y = top; y <= bottom; y++){
            if(imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x] == 0)
                continue;

            BETA = imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x];

            tr = ((float) imgsyn->red[y*imgsyn->ncol + x])*(1-BETA)   + BETA*TR;
            tg = ((float) imgsyn->green[y*imgsyn->ncol + x])*(1-BETA) + BETA*TG;
            tb = ((float) imgsyn->blue[y*imgsyn->ncol + x])*(1-BETA)  + BETA*TB;

            tR = ((float) imgsyn->red[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tr;
            imgsyn->red[y*imgsyn->ncol + x]   = (int) tR; 

            tG = ((float) imgsyn->green[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tg;
            imgsyn->green[y*imgsyn->ncol + x] = (int) tG;  

            tB = ((float) imgsyn->blue[y*imgsyn->ncol + x])*ALPHA + (1-ALPHA)*tb;
            imgsyn->blue[y*imgsyn->ncol + x]  = (int) tB;  

            imgShapeBlurSyn->gray[y*imgShapeBlurSyn->ncol + x]  = 0.0;
            imgShapeLabelSyn->gray[y*imgShapeLabelSyn->ncol + x] = 0;
        }
}


// Sort shapes of the tree
void TreeOfShapes::sortShapes(Fsignal t2b_index){

    std::priority_queue < std::pair<int, int>, std::deque< std::pair<int, int> > , std::less<std::pair<int, int> > > AreaShapeIDQueue;
    Shape pShape;

    for(int i=0; i< _pTree->nb_shapes; i++) {
        pShape = _pTree->the_shapes + i;
        AreaShapeIDQueue.push( std::make_pair(pShape->area, i));
    }

    int i = 0;
    while(!AreaShapeIDQueue.empty()){
        std::pair<int, int> area_id_pair = AreaShapeIDQueue.top();
        AreaShapeIDQueue.pop();
        t2b_index->values[i++] = area_id_pair.second;
    }
}


// Generate a 1D Gaussian kernel. Used for generate a 2D Gaussian kernel
Fsignal TreeOfShapes::sgauss(float *std, Fsignal out, int *size){  
    int i,n;
    double sum, v;

    n = *size;
    v = (0.5*(double)(n-1)) / (double)(*std);
    v = 0.5*v*v/log(10.);

    out = mw_change_fsignal(out, n);
    if (!out) mwerror(FATAL,1,"Not enough memory.");

    out->shift = -0.5*(float)(n-1);

    if (n==1){
        out->values[0]=1.0;
    }
    else{
        // Store Gaussian signal  
        for(i=(n+1)/2; i--; ){
            v = ((double)i+(double)out->shift)/(double)(*std);
            out->values[i] = out->values[n-1-i] = (float)exp(-0.5*v*v);
        }
        // Normalize to get unit mass  
        for (sum=0.0,i=n; i--; )
            sum += (double)out->values[i];
        for (i=n; i--; )
            out->values[i] /= (float)sum;
    }

    return(out);
}


// Generate a 2D Gaussian kernel
Fsignal TreeOfShapes::Sgauss(float *std, Fsignal out, int *size){
    Fsignal sgaussX, sgaussY;
    int i, j, n;
    float sum= 0.0;
    n = *size;

    if  ( ((sgaussX = mw_new_fsignal()) == NULL) || (mw_alloc_fsignal(sgaussX, n) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");
    if  ( ((sgaussY = mw_new_fsignal()) == NULL) || (mw_alloc_fsignal(sgaussY, n) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    sgaussX = sgauss(std, sgaussX, size);
    sgaussY = sgauss(std, sgaussY, size);
    out = mw_change_fsignal(out, n*n);

    for (i = 0; i<n; i++)
        for (j = i; j<n; j++){
            out->values[j*n + i] = sgaussX->values[i] * sgaussY->values[j];
            out->values[i*n + j] = out->values[j*n + i];
        }
    for (i = 0; i<n*n; i++)
         sum += out->values[i];
    for (i = 0; i<n*n; i++)
        out->values[i] /= sum;

    mw_delete_fsignal(sgaussX);
    mw_delete_fsignal(sgaussY);
    return(out);
}

// Recursively compute area of shapes when holes are removed 
void TreeOfShapes::get_shapes_truearea(Shape s, Shape root, int *truearea){
    int index = s-root;
    truearea[index] = s->area;
    Shape t = mw_get_first_child_shape(s);
    while(t){
        get_shapes_truearea(t,root,truearea);
        truearea[index] -= t->area;
        t = mw_get_next_sibling_shape(t);
    }
}


void TreeOfShapes::filter_shapes( Cfimage out, char *local, float *eps){

    std::cout <<"TreeOfShapes::filter_shapes()::begin"<< std::endl;
    Fimage Fv;
    float lgeo,step,hstep,tcleanup,farea, std, *gray,*red,*green,*blue;
    int prec,visit,nrow,ncol,i,j,indexshape,*truearea;
    char all;
    Shapes tree;
    Shape s,t;

    tree = mw_new_shapes();
    Fv = mw_change_fimage(NULL,_imgin->nrow,_imgin->ncol);

    if(!Fv ) mwerror(FATAL,1,"Not enough memory.\n");

    lgeo = 10.; step = 1.; prec = 2; hstep = 0.01; std=0.5;
    tcleanup = 1.; all=(char) 1; visit = 100;

    nrow = _imgin->nrow;
    ncol = _imgin->ncol;
    for(i=0;i<ncol*nrow;i++)
        Fv->gray[i] = (_imgin->red[i]+_imgin->green[i]+_imgin->blue[i])/3.;
    
    // Function from MegaWave module. Extract (local or not) meaningful boundaries from FV
    ll_boundaries2(Fv,eps,NULL,&step,&prec,&std,&hstep,NULL,&visit,NULL,NULL,tree);

    // Compute recursively integral of gray level saturation cos and sin of hue 
    red = (float*)calloc(tree->nb_shapes,sizeof(float));
    green = (float*)calloc(tree->nb_shapes,sizeof(float));
    blue = (float*)calloc(tree->nb_shapes,sizeof(float));
    gray = (float*)calloc(tree->nb_shapes,sizeof(float));
    if(!(gray && red  && green && blue))
        mwerror(FATAL,1,"Not enough memory.\n");

    // Integrate grey level saturation and hue
    for(i=0;i<ncol;i++)
        for(j=0;j<nrow;j++){
            s = mw_get_smallest_shape(tree,i,j);
            indexshape = s-tree->the_shapes;
            gray[indexshape] += Fv->gray[i+j*ncol];
            red[indexshape] += _imgin->red[i+j*ncol];
            green[indexshape] += _imgin->green[i+j*ncol];
            blue[indexshape] += _imgin->blue[i+j*ncol];
        }

    // Recursively compute area of meaningful shapes when holes are removed 
    truearea = (int*)malloc(sizeof(int)*tree->nb_shapes);
    if(!truearea)     
       mwerror(FATAL,1,"Not enough memory.\n");
    get_shapes_truearea(tree->the_shapes,tree->the_shapes,truearea);

    for(i=0;i<tree->nb_shapes;i++){
        s = tree->the_shapes+i;
        farea = (float) (s->area);
        t = mw_get_first_child_shape(s);
        while(t){
            farea -= (float)t->area;
            t = mw_get_next_sibling_shape(t);
        }

        gray[i] /= (float)truearea[i];
        red[i] /= (float)truearea[i];
        green[i] /= (float)truearea[i];
        blue[i] /= (float)truearea[i];
    }

    // Replace gray saturation and hue by average
    out = mw_change_cfimage(out,_imgin->nrow,_imgin->ncol);
    for(i=0;i<ncol;i++)
        for(j=0;j<nrow;j++){
            s = mw_get_smallest_shape(tree,i,j);
            indexshape = s-tree->the_shapes;
            out->red[i+j*ncol]= red[indexshape];
            out->green[i+j*ncol]= green[indexshape];
            out->blue[i+j*ncol]= blue[indexshape];
        }

    free(gray);
    free(red);
    free(green);
    free(blue);
    free(truearea);
}


// Select shapes to remove of the image 
void TreeOfShapes::filter_image(int *ns,float *threshold,float *minarea,float *maxarea, int totalSize, float *k){
    // Declare variables here
    float thre, CONTR, mpixel, maxpixel;
    Shape pShape, pShapeTemp;
    thre = *threshold;

    // Define size in pixels for filtering shapes
    mpixel = totalSize * (*minarea) / 100;
    maxpixel = totalSize * (*maxarea) / 100;

    // Filtering the image
    for(int i = 0; i<=_pTree->nb_shapes-1; i++){
        pShape = _pTree->the_shapes + i;

        if(pShape->parent == NULL)
            continue;

        pShapeTemp = m_order_parent(pShape, *ns);
        ((Info*)(pShape->data))->attribute[0] = ((float) pShape->area)/((float) pShapeTemp->area);

        // Update attribute 0 of m-order-parent shape. 
        if( ((Info*)(pShapeTemp->data))->attribute[0] <= ((float) pShape->area)/((float) pShapeTemp->area))
            ((Info*)(pShapeTemp->data))->attribute[0] = ((float) pShape->area)/((float) pShapeTemp->area);
        
        // compute distance of color between shape and parent shape (Contrast)
        CONTR = sqrt(pow((((Info*)(pShape->data))->r - ((Info*)(pShape->parent->data))->r), 2.0) +
                     pow((((Info*)(pShape->data))->g - ((Info*)(pShape->parent->data))->g), 2.0) +
                     pow((((Info*)(pShape->data))->b - ((Info*)(pShape->parent->data))->b), 2.0));

        // Conditions to remove shape:
        // Area smaller than mpixel or larger than maxpixel or
        // (area/area_m-order_parent)*(Contrast between shape and parent) is less than threshold or. Why??
        // (area/area_grand_parent) < compactness parameter 
        if(pShape->area <= mpixel || maxpixel < pShape->area || (((Info*)(pShape->data))->attribute[0])*CONTR<= thre || ((((Info*)(pShape->data))->attribute[0])< (*k)))
            pShape->removed = 1;
        else
            pShape->removed = 0;
        
        // First shape it's not removed
        if(i ==0)
            pShape->removed = 0;
    };
}

// Selects random number smaller than (*M). 
int TreeOfShapes::random_number(int *M){
    int i, size, select_i;
    float temp = ((float)rand())/RAND_MAX;
    Fsignal pb=0;
    size = *M;
    select_i = size-1;

    pb = mw_change_fsignal(pb,size);
    mw_clear_fsignal(pb,0.0);

    // Sampling by uniform distribution
    for(i=0; i< size; i++)
        pb->values[i] = 1/size;

    for(i=1; i< size; i++){
        pb->values[i] += pb->values[i-1];
        if( temp <= pb->values[i]){
            select_i = i;
            break;
        }
    }

    mw_delete_fsignal(pb);
    return select_i;
}


// Function to Compute kd_tree in case there are more than 1000 shapes in the tree.
void TreeOfShapes::computeKdTree(float average_r, float average_g, float average_b ){
    _use_kdtree = false;
    float lambda1, lambda2, elongDict, kappaDict, scaDict;

    // Compute kd_tree in case there are more than 1000 shapes in the tree. 
    if( _pTree->nb_shapes > 1000 ){
        std::cout <<"Building KD tree" << std::endl;
        _annTree = BasicANNkdTree(6);
        Shape pShape;

        std::vector<std::vector<float>> values;
        // Compute values for all shapes to push into the kd_tree. 
        for(int i=1; i < _pTree->nb_shapes; i++){
            pShape = _pTree->the_shapes + i;
            lambda1 = ((Info*)(pShape->data))->lambda1;
            lambda2 = ((Info*)(pShape->data))->lambda2;
            elongDict = lambda2 / lambda1;
            kappaDict = ((float) pShape->area)/(sqrt(lambda2*lambda1)*4*PI);
            scaDict = ((float) pShape->area);

            std::vector<float> value;
            value.push_back(elongDict);
            value.push_back(kappaDict);
            value.push_back(scaDict);
            value.push_back(((Info*)(pShape->data))->r/average_r);
            value.push_back(((Info*)(pShape->data))->g/average_g);
            value.push_back(((Info*)(pShape->data))->b/average_b);
            values.push_back(value);
        }

        _annTree.build(values);
        _use_kdtree = true;
    }
}

// Select Shape according to the distance of its attributes
Shape TreeOfShapes::selectShapeDict(Shape pShape,
                                    float *paDict,
                                    int *randS,
                                    int &index,
                                    float average_r,
                                    float average_g,
                                    float average_b){

    Shape pShapeDict, pShapeTemp;
    float sca, scaDict, pa, elong, kappa, elongDict, kappaDict, Dist, minDist, lambda1, lambda2;
    int i, temp;

    lambda1 = ((Info*)(pShape->data))->lambda1;
    lambda2 = ((Info*)(pShape->data))->lambda2;
    elong = lambda2 / lambda1;
    kappa = ((float) pShape->area)/(sqrt(lambda2*lambda1)*4*PI);
    sca = ((float) pShape->area);
    index = 1; 
    minDist = 10000.0;
    
    if(*randS == 0){ // randS=0, randomly select shapes; 
        temp = _pTree->nb_shapes -1;
        index = random_number(&temp) + 1;
    } else{
        // randS=1, select shapes according to elongation, compactness and scale;  
        // randS=2, select shapes according to elongation, compactness, scale and color
        if (!_use_kdtree){
            for(i= 1; i<_pTree->nb_shapes; i++){

                // Attribute filtering. Index the 3th parent of the shape to compute compactness.  
                pShapeDict = _pTree->the_shapes + i;
                pShapeTemp =  m_order_parent(pShapeDict, 3, true);
                if((((float) pShapeDict->area)/((float) pShapeTemp->area)) < *paDict)
                    continue;

                lambda1 = ((Info*)(pShapeDict->data))->lambda1;
                lambda2 = ((Info*)(pShapeDict->data))->lambda2;
                elongDict = lambda2 / lambda1;
                kappaDict = ((float) pShapeDict->area)/(sqrt(lambda2*lambda1)*4*PI);
                scaDict = ((float) pShapeDict->area);
                
                // compute distance according to elongation, compactness and scale 
                Dist = pow((elong - elongDict), 2.0) + pow((kappa - kappaDict), 2.0) + pow((1 - _MIN(sca/scaDict, scaDict/sca)), 2.0);
                
                // if randS==2, Add information of color to the distance
                if(*randS == 2){
                    Dist +=  ( pow((1 - _MIN(((Info*)(pShape->data))->r/((Info*)(pShapeDict->data))->r,
                                             ((Info*)(pShapeDict->data))->r/((Info*)(pShape->data))->r)), 2.0) +
                               pow((1 - _MIN(((Info*)(pShape->data))->g/((Info*)(pShapeDict->data))->g,
                                             ((Info*)(pShapeDict->data))->g/((Info*)(pShape->data))->g)), 2.0) +
                               pow((1 - _MIN(((Info*)(pShape->data))->b/((Info*)(pShapeDict->data))->b,
                                             ((Info*)(pShapeDict->data))->b/((Info*)(pShape->data))->b)), 2.0) )/3.0;
                }
                
                // Update closest shape
                if(minDist > Dist){
                    minDist = Dist;
                    index = i;
                }
            }
        } else{
            // Get the k=10 nearest neighbors using kd_tree. 
            std::vector<float> value;
            value.push_back(elong);
            value.push_back(kappa);
            value.push_back(sca);
            value.push_back(((Info*)(pShape->data))->r/average_r);
            value.push_back(((Info*)(pShape->data))->g/average_g);
            value.push_back(((Info*)(pShape->data))->b/average_b);

            int k = std::min(10, _pTree->nb_shapes);
            ANNidxArray neighbors = new ANNidx[ k ];
            ANNdistArray neighbors_sqr_dists= new ANNdist[ k ];
            _annTree.knearest(value, k, neighbors, neighbors_sqr_dists);

            for(i= 1; i<k; i++){

                // Attribute filtering. Index the 3th parent of the shape to compute compactness  
                pShapeDict = _pTree->the_shapes + (int)neighbors[i];
                pShapeTemp =  m_order_parent(pShapeDict, 3, true);
                if((((float) pShapeDict->area)/((float) pShapeTemp->area)) < *paDict)
                    continue;

                lambda1 = ((Info*)(pShapeDict->data))->lambda1;
                lambda2 = ((Info*)(pShapeDict->data))->lambda2;
                elongDict = lambda2 / lambda1;
                kappaDict = ((float) pShapeDict->area)/(sqrt(lambda2*lambda1)*4*PI);
                scaDict = ((float) pShapeDict->area);

                // Compute distance according to elongation, compactness and scale 
                Dist = pow((elong - elongDict), 2.0) + pow((kappa - kappaDict), 2.0) + pow((1 - _MIN(sca/scaDict, scaDict/sca)), 2.0);

                // randS=2. Add information of color to the distance
                if(*randS == 2){
                    Dist +=  ( pow((1 - _MIN(((Info*)(pShape->data))->r/((Info*)(pShapeDict->data))->r,
                                             ((Info*)(pShapeDict->data))->r/((Info*)(pShape->data))->r)), 2.0) +
                               pow((1 - _MIN(((Info*)(pShape->data))->g/((Info*)(pShapeDict->data))->g,
                                             ((Info*)(pShapeDict->data))->g/((Info*)(pShape->data))->g)), 2.0) +
                               pow((1 - _MIN(((Info*)(pShape->data))->b/((Info*)(pShapeDict->data))->b,
                                             ((Info*)(pShapeDict->data))->b/((Info*)(pShape->data))->b)), 2.0) )/3.0;
                }
                
                // Update closest shape 
                if(minDist > Dist){
                    minDist = Dist;
                    index = (int)neighbors[i];
                }
            }
        }
    }

    pShapeDict = _pTree->the_shapes + index;
    return pShapeDict;
}


Shape TreeOfShapes::getShape(int index){
    if( index < _pTree->nb_shapes )
        return _pTree->the_shapes + index;
    return NULL;
}


void TreeOfShapes::compute_tree( TOSParameters tosParameters, bool dictionary ){

    std::cout <<"Compute_tree started"<< std::endl;

    if( !_tree_computed || _tosParameters.color_sketch != tosParameters.color_sketch ||
            ( _tosParameters.color_sketch == 1 & _tosParameters.eps != tosParameters.eps ) ){

        if(tosParameters.color_sketch == 1){
            char local=NULL ; // "local boundaries, default NULL"
            Cfimage out = mw_change_cfimage(NULL,_imgin->nrow,_imgin->ncol);
            filter_shapes(out,&local,&tosParameters.eps);
            init(out, _pTree);
            mw_delete_cfimage(out);
        }else if( _tosParameters.color_sketch == 1 & tosParameters.color_sketch == 0 || !_tree_computed)
            init(_imgin, _pTree);

        _tree_computed = true;
        _tree_recomputed = true;
        _use_kdtree = false;

        for (std::map<int, Fsignal>::iterator it = _dictionary_selections.begin(); it !=  _dictionary_selections.end(); ++it)
            mw_delete_fsignal( it->second );
        
        _dictionary_selections.clear();
    }
    
    // Compute shape attribute if dictionary
    if( dictionary ){
        std::cout << "Compute shape attributes if dictionary" << std::endl;
        compute_shape_attribute();
        tree_boundingbox();
    }
}

// Compute list of pixels of mask (mask select shapes to transform by alternative shapes)  
void TreeOfShapes::compute_list_pixels_mask(QImage image_mask){
    QColor color_ij;
    Point_plane pCurrentPoint;
    _len_ArrayPixelsMask = 0;

    for( int i= 0; i< image_mask.width() ; i++)
        for( int j= 0; j< image_mask.height(); j++){
            color_ij =image_mask.pixel( i, j );       
            if (!(color_ij.red() == 0 &&  color_ij.blue() == 0 && color_ij.green() == 0)){
                pCurrentPoint = &_ArrayPixelsMask[_len_ArrayPixelsMask];
                pCurrentPoint->x = i;
                pCurrentPoint->y = j;
                _len_ArrayPixelsMask = _len_ArrayPixelsMask +1;   
            };
        };
}

// Function to add a random shift to each shape of the tree. Additionaly add also a random rotation to each shape (smode =1)
void TreeOfShapes::shift_shapes(float *shift, float *theta, int mode){
    int i,j;
    Shape pShape;
    float SHIFT, THETA, tempx, tempy, tempt, tempShiftx, tempShifty, a, b, phi, CONTRAST, length;
    SHIFT = *shift;
    THETA = *theta;

    srand( time(NULL) );
    for(i = _pTree->nb_shapes-1; i>= 0; i--){
        pShape = _pTree->the_shapes + i;

        if(pShape == NULL )
            continue;

        if(i==0){
            ((Info*)(pShape->data))->xShift = 0;
            ((Info*)(pShape->data))->yShift = 0;
            ((Info*)(pShape->data))->rotation = 0;
            continue;
        }

        if( (rand()%10) >5 )
            tempx = -((float)rand())/RAND_MAX;
        else
            tempx = ((float)rand())/RAND_MAX;

        if( (rand()%10) >5 )
            tempy = -((float)rand())/RAND_MAX;
        else
            tempy = ((float)rand())/RAND_MAX;

        if (mode == 0){ // shaking type: uniform shaking: smode=0
            ((Info*)(pShape->data))->xShift = tempx*SHIFT;
            ((Info*)(pShape->data))->yShift = tempy*SHIFT;
        } else{ // shaking type: dominant shaking: smode=1
            if( (rand()%10) >5 )
                tempt = -((float)rand())/RAND_MAX;
            else
                tempt = ((float)rand())/RAND_MAX;

            a = 2.0 * sqrt(((Info*)(pShape->data))->lambda1);
            b = 2.0 * sqrt(((Info*)(pShape->data))->lambda2);
            phi = ((Info*)(pShape->data))->oren;

            tempShiftx = tempx*SHIFT;
            tempShifty = tempy*SHIFT*pow((b/a),2);

            ((Info*)(pShape->data))->xShift = tempShiftx*cos(phi) + tempShifty*sin(phi);
            ((Info*)(pShape->data))->yShift = tempShifty*cos(phi) - tempShiftx*sin(phi);

            ((Info*)(pShape->data))->rotation = tempt*THETA*PI*(b/a);

            Flist pBoundary = NULL;
            pBoundary = mw_change_flist(pBoundary, 4*pShape->area+1, 0, 2);
            flst_boundary(_pTree, pShape, pBoundary);
            CONTRAST = min_contrast(pBoundary,&length,_NormOfDu);

            ((Info*)(pShape->data))->xShift *= 1/pow(CONTRAST, 0.5);
            ((Info*)(pShape->data))->yShift *= 1/pow(CONTRAST, 0.5);
            ((Info*)(pShape->data))->rotation *= 1/pow(CONTRAST, 0.5);
        }
    }
}


QImage TreeOfShapes::render(TOSParameters tosParameters, bool segmentWithMask, int alternative_model, TreeOfShapes *tosDictionary, DictionaryParameters dictionaryParameters ){
    
    std::cout <<"TreeOfShapes::Abstraction started"<< std::endl;

    //Declare variables
    int i,j, modelToUse, shape_id, totalSize;
    Shape pShape, pShapeTemp, pShapeDict;  
    Cimage imgShapeLabel, imgShapeLabelDict;
    Cfimage imgShapeColorSyn, imgDict;
    Fimage imgShapeBlur;
    Fsignal t2b_index, gaussKernel, dictionary_correspondance;

    // Define synthesis image
    Ccimage imgsyn = mw_change_ccimage(imgsyn, _imgin->nrow, _imgin->ncol);

    //Step 1: Decomposition. 
    compute_tree(tosParameters, false);

    // Image filtering    
    std::cout << "Image filtering" << std::endl;
    totalSize = _imgin->nrow * _imgin->ncol;
    filter_image(&tosParameters.ns,&tosParameters.threshold, &tosParameters.minarea, &tosParameters.maxarea, totalSize, &tosParameters.kappa);

    // Select the rendering order 
    std::cout << "Rendering order " << tosParameters.order <<std::endl;

    if  ( ((t2b_index = mw_new_fsignal()) == NULL) ||(mw_alloc_fsignal(t2b_index,_pTree->nb_shapes) == NULL) )
        mwerror(FATAL,1,"Not enough memory.\n");

    if(tosParameters.order == 0)
        top2bottom_index_tree(t2b_index);
    else{ // tosParameters.order == 1     
        if ( ((_large_to_small_index = mw_new_fsignal()) == NULL) || (mw_alloc_fsignal(_large_to_small_index,_pTree->nb_shapes) == NULL) )
            mwerror(FATAL,1,"Not enough memory.\n");
        sortShapes(_large_to_small_index);
        _large_to_small_index_computed = true;     
        mw_copy_fsignal_values(_large_to_small_index, t2b_index);
    } 
  
    // Special resources needed if blur is selected
    if (tosParameters.blur == 1){
        // Define auxiliar images for abstraction
        if  ( ((imgShapeLabel = mw_new_cimage()) == NULL) || (mw_alloc_cimage(imgShapeLabel, _imgin->nrow, _imgin->ncol) == NULL) )
            mwerror(FATAL,1,"Not enough memory.\n");
        if  ( ((imgShapeBlur = mw_new_fimage()) == NULL) || (mw_alloc_fimage(imgShapeBlur, _imgin->nrow, _imgin->ncol) == NULL) )
            mwerror(FATAL,1,"Not enough memory.\n");

        imgShapeLabel = mw_change_cimage(imgShapeLabel, _imgin->nrow, _imgin->ncol);
        imgShapeBlur  = mw_change_fimage(imgShapeBlur, _imgin->nrow, _imgin->ncol);
        mw_clear_cimage(imgShapeLabel,0);
        mw_clear_fimage(imgShapeBlur,0.0);

        // Compute a Gaussian kernel
        std::cout << "Compute a Gaussian kernel" << std::endl;     
        if  ( ((gaussKernel = mw_new_fsignal()) == NULL) || (mw_alloc_fsignal(gaussKernel, tosParameters.kerSize*tosParameters.kerSize) == NULL) )
            mwerror(FATAL,1,"Not enough memory.\n");
        gaussKernel = Sgauss(&tosParameters.kerStd, gaussKernel, &tosParameters.kerSize);

        // Sepcial resources needed if dictionary is selected
        if( tosParameters.model == 4 ){
            // Define auxiliar images for abstraction
            imgDict = tosDictionary->getCfImage();
            if  ( ((imgShapeColorSyn = mw_new_cfimage()) == NULL) || (mw_alloc_cfimage(imgShapeColorSyn, _imgin->nrow, _imgin->ncol) == NULL) )
                mwerror(FATAL,1,"Not enough memory.\n");
            if  ( ((imgShapeLabelDict = mw_new_cimage()) == NULL) || (mw_alloc_cimage(imgShapeLabelDict, imgDict->nrow, imgDict->ncol) == NULL) )
                mwerror(FATAL,1,"Not enough memory.\n");

            imgShapeColorSyn = mw_change_cfimage(imgShapeColorSyn, _imgin->nrow, _imgin->ncol);
            imgShapeLabelDict = mw_change_cimage(imgShapeLabelDict, imgDict->nrow, imgDict->ncol);
            mw_clear_cimage(imgShapeLabelDict,0);

            // Compute kd-tree to perform efficient search to compute matching shapes between shapes of images - shapes of dictionary.  
            if  ( ((dictionary_correspondance = mw_new_fsignal()) == NULL) ||(mw_alloc_fsignal(dictionary_correspondance,_pTree->nb_shapes) == NULL) )
                mwerror(FATAL,1,"Not enough memory.\n");
            if  ( ((_dictionary_selections[ tosDictionary->getTreeId() ] = mw_new_fsignal()) == NULL) || (mw_alloc_fsignal(_dictionary_selections[ tosDictionary->getTreeId() ],_pTree->nb_shapes) == NULL) )
                mwerror(FATAL,1,"Not enough memory.\n");
            mw_clear_fsignal(_dictionary_selections[ tosDictionary->getTreeId() ],-1.0);
            tosDictionary->computeKdTree(_average_r, _average_g, _average_b);
        }
    }

    // Add a random shift to each shape 
    std::cout << "Image Shaking" << std::endl;   
    shift_shapes(&tosParameters.shift, &tosParameters.theta, tosParameters.smodel);
     
    // Iterate in shapes
    std::cout << "Iterate in shapes: " <<  _pTree->nb_shapes << " shapes." << std::endl;

    for(i=0; i < _pTree->nb_shapes; i++) {
        pShape = _pTree->the_shapes + (int)t2b_index->values[i];
         
        if((int)t2b_index->values[i] == 0 ) { // Background shape: Rectangle and average color.
           
            //before: if (tosParameters.model == 4 && (dictionaryParameters.mcolor == 1 || dictionaryParameters.mcolor ==2)){
            if (tosParameters.model == 4 && (dictionaryParameters.color_background==0)){
                // Take color from dictionary for background
                pShapeDict = tosDictionary->getShape(0);
                ((Info*)(pShape->data))->r = ((Info*)(pShapeDict->data))->r;
                ((Info*)(pShape->data))->g = ((Info*)(pShapeDict->data))->g;
                ((Info*)(pShape->data))->b = ((Info*)(pShapeDict->data))->b;
            }
            float ALPHA = 0.0;
            synshape(2, pShape, imgsyn, &ALPHA);              
        } 
        else if(pShape->removed != 1){
                // Verify if some point of the mask touch the shape. 
                modelToUse = tosParameters.model;
                if (!segmentWithMask)
                    for (j=0; j<_len_ArrayPixelsMask; j++)
                        if (point_in_shape((&_ArrayPixelsMask[j])->x, (&_ArrayPixelsMask[j])->y, pShape, _pTree)){
                            modelToUse = alternative_model;
                            break;
                        };

                // Modification of shape according to model
                if (modelToUse < 4){ // Rendering Model: Original, Rectangle, Ellipse or Circular
                    if(tosParameters.blur == 0)
                       synshape(modelToUse, pShape, imgsyn, &tosParameters.alpha);
                    else if (modelToUse == 0)
                        synshapeOriginal(pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha);
                    else 
                        synshape(modelToUse, pShape, imgsyn, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha);
                }
                else{ // modelToUse ==4 -> Rendering Model: Dictionary
                    pShapeDict = tosDictionary->selectShapeDict(pShape, &dictionaryParameters.kappaDict, &dictionaryParameters.randS, shape_id, _average_r, _average_g, _average_b);
                    dictionary_correspondance->values[(int)t2b_index->values[i]] = shape_id;
                    synShapeDict( pShapeDict, pShape, imgsyn, imgDict, imgShapeColorSyn, imgShapeLabelDict, imgShapeLabel, imgShapeBlur, gaussKernel, &tosParameters.median, &tosParameters.alpha, &dictionaryParameters.equal, &dictionaryParameters.mcolor);
                    //std::cout << "3 " << std::endl; 
                }
        }
    }

    // Compute Resuting image
    std::cout << "Compute Resuting image" << std::endl; 
    QImage result_image( QSize(imgsyn->ncol, imgsyn->nrow), QImage::Format_RGB32 );
    for( int j= 0; j< imgsyn->nrow; j++)
        for( int i= 0; i< imgsyn->ncol; i++){
            int comp = j*imgsyn->ncol + i;
            QColor color (imgsyn->red[comp], imgsyn->green[comp], imgsyn->blue[comp]);
            result_image.setPixel(i, j , qRgb(color.red(), color.green(), color.blue()));
        }
 
    // Delete image and signals
    std::cout << "Delete auxiliar images and signals" << std::endl;
    
    mw_delete_fsignal(t2b_index);
    if (imgsyn != NULL)
        mw_delete_ccimage(imgsyn);
    if (tosParameters.blur == 1){
        mw_delete_cimage(imgShapeLabel);
        mw_delete_fimage(imgShapeBlur);
        mw_delete_fsignal(gaussKernel);
        if (tosParameters.model == 4){
            mw_delete_cfimage(imgShapeColorSyn);
            mw_delete_cimage(imgShapeLabelDict);
        }
    }
 
    return result_image;
}