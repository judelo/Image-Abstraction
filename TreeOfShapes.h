/*
 Module to implement abstraction process.

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef TREEOFSHAPES_H
#define TREEOFSHAPES_H

#include <QImage>

#include "mw3.h"
#include "mw3-modules.h"
#include "stdio.h"
#include "stdlib.h"
#include "tree_of_shapes.h"
#include "jmkdtree.h"
#include <cfloat>
#include <map>

class TreeOfShapes
{
public:

    static int _tree_count;
    TreeOfShapes( Cfimage imageIn );
    ~TreeOfShapes();
    QImage render(TOSParameters tosParameters, bool &tree_recomputed, QImage image_mask, int alternative_model, TreeOfShapes *tosDictionary=NULL, DictionaryParameters dictionaryParameters=getDefaultDictionaryParameters() );
    void compute_tree( TOSParameters tosParameters, bool dictionary=false );
    void computeKdTree(float average_r, float average_g, float average_b );
    Cfimage getCfImage(){ if( _texture_image_loaded ) return _texture_image; else return _imgin; }
    Shape selectShapeDict(Shape pShape, float *paDict, int *randS, int &index, float average_r, float average_g, float average_b);
    Shape getShape(int index);
    int getTreeId(){ return _tree_id; }
    int getMaxArea(){ return _maxArea; }
protected:
    bool _tree_computed;
    bool _texture_image_loaded;
    bool _tree_recomputed;
    bool _large_to_small_index_computed;
    bool _use_kdtree;
    BasicANNkdTree _annTree;
    int _tree_id;
    int _maxArea;
    Point_plane  _ArrayPixelsMask;
    int _len_ArrayPixelsMask;

    Cfimage _imgin;
    Cfimage _texture_image;
    Shapes _pTree;
    Fimage _NormOfDu;
    Fsignal _large_to_small_index;
    std::map<int, Fsignal> _dictionary_selections;

    TOSParameters _tosParameters;
    DictionaryParameters _dictionaryParameters;

    float _average_r;
    float _average_g;
    float _average_b;

    void init(Cfimage inputImg, Shapes &pTree);
    void sortShapes(Fsignal t2b_index);
    void fgrain_side(int MinArea, float *in, int nx, int ny, float *out, int sideflag);
    Shape m_order_parent(Shape pShape, int *mn, bool dict = false);
    void Order(Fsignal t2b_index, int *p, int *q);
    void shape_orilam(Shape pShape,float *out_ori, float *out_e, float *out_k);
    void shape_orilam(Shape pShape, float *out_ori, float *out_e, float *out_k, float *pX0, float *pY0);
    void compute_shape_attribute();
    float min_contrast(Shape pShape);
    void synshapeCircle(Shape pShape,
                        Ccimage imgsyn,
                        float *alpha,
                        int *relief,
                        float *reliefOrentation, float *reliefHeight);
    void synshapeCircle(Shape pShape,
                        Ccimage imgsyn,
                        Cimage imgShapeLabelSyn,
                        Fimage imgShapeBlurSyn,
                        Fsignal gaussKernel,
                        int *median,
                        float *alpha,
                        int *relief,
                        float *reliefOrentation, float *reliefHeight);
    void synshapeEllipse(Shape pShape,
                         Ccimage imgsyn,
                         float *alpha,
                         int *relief,
                         float *reliefOrentation, float *reliefHeight);
    void synshapeEllipse(Shape pShape,
                         Ccimage imgsyn,
                         Cimage imgShapeLabelSyn,
                         Fimage imgShapeBlurSyn,
                         Fsignal gaussKernel,
                         int *median,
                         float *alpha,
                         int *relief,
                         float *reliefOrentation, float *reliefHeight);
    void synshapeRect(Shape pShape,
                      Ccimage imgsyn,
                      float *alpha,
                      int *relief,
                      float *reliefOrentation, float *reliefHeight);
    void synshapeRect(Shape pShape,
                      Ccimage imgsyn,
                      Cimage imgShapeLabelSyn,
                      Fimage imgShapeBlurSyn,
                      Fsignal gaussKernel,
                      int *median,
                      float *alpha,
                      int *relief,
                      float *reliefOrentation, float *reliefHeight);
    void synshapeOriginal(Shape pShape,
                          Ccimage imgsyn,
                          Cimage imgShapeLabelSyn,
                          Fimage imgShapeBlurSyn,
                          Fsignal gaussKernel,
                          int *median,
                          float *alpha,
                          int *relief,
                          float *reliefOrentation, float *reliefHeight);
    void synshapeOriginal( Shape pShape,
                           Ccimage imgsyn,
                           float *alpha,
                           int *relief,
                           float *reliefOrentation, float *reliefHeight);

    void top2bottom_index_tree(Fsignal t2b_index);

    Fsignal sgauss(float *std, Fsignal out, int *size);
    Fsignal Sgauss(float *std, Fsignal out, int *size);

    void compute_shape_attribute(int *ns);
    void filter_image(int *ns,
                      float *alpha,
                      int *mpixel,
                      int *maxpixel);
    void filter_shapes( Cfimage out,
                        char *local,
                        float *eps);
    void get_shapes_truearea(Shape s, Shape root,
                             int *truearea);

    int random_number(int *M);

    void shape_boundingbox(Shape pShape);
    void tree_boundingbox();

    void synShapeDict(Shape pShapeDict, Shape pShape,
                      Ccimage imgsyn,
                      Cfimage imgDict, Cfimage imgShapeColorSyn,
                      Cimage imgShapeLabel, Cimage imgShapeLabelSyn,
                      Fimage imgShapeBlurSyn,
                      Fsignal gaussKernel,
                      int *median,
                      float *alpha,
                      int *equal, int *mcolor, int *relief,
                      float *reliefOrentation, float *reliefHeight);
    void synShapeDict(Shape pShapeDict, Shape pShape,
                      Ccimage imgsyn,
                      Cfimage imgDict, Cfimage imgShapeColorSyn,
                      Cimage imgShapeLabel, Cimage imgShapeLabelSyn,
                      float *alpha,
                      int *equal, int *mcolor, int *relief,
                      float *reliefOrentation, float *reliefHeight);
};

#endif // TREEOFSHAPES_H
