/*
 Module to implement abstraction process.

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef ABSTRACTIONGRAINPROCESS_H
#define ABSTRACTIONGRAINPROCESS_H

#include <QImage>

#include "mw3.h"
#include "mw3-modules.h"
#include "stdio.h"
#include "stdlib.h"
#include "tree_of_shapes.h"
#include <cfloat>
#include "TreeOfShapes.h"

enum Abstraction_Mode { COLOR_SKETCH=0, FILTER_COLOR=1, SHOW_TREE=2, SYNTEXTURE_COLOR=3, SYNTEXTURE_COLOR_WA=4,
                        SYNTEXTURE_COLOR_DICT=5, SYNTEXTURE_COLOR_DICT2=6, SYNTEXTURE_COLOR_V2=7,
                        SYNTEXTURE_COLOR_V3=8, SYNTEXTURE_COLOR_V4=9,
                        SYNTEXTURE_COLOR_DICT2_OUTPUT=10, SYNTEXTURE_COLOR_DICT3=11, SYNTEXTURE_COLOR_MULT=12,
                        SYNTEXTURE_COLOR_TT=13};

class AbstractionProcess
{
public:
    AbstractionProcess(){ _image_loaded = false; }
    AbstractionProcess( std::string fileNameIn );
    AbstractionProcess( const QImage &imageIn );
    ~AbstractionProcess();

    //QImage render(TOSParameters tosParameters, bool &tree_recomputed, DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters(), TreeOfShapes * dictionnary=NULL);

protected:
    bool _tree_computed;
    Cfimage _imgin;
    bool _image_loaded;

    TreeOfShapes *_treeOfShapes;

    Shapes _pTree;

    Fimage   _NormOfDu;

    TOSParameters _tosParameters;

    void init(Cfimage inputImg, Shapes &pTree);

    Cfimage cfimageread(const char* name);

    Cfimage cfimages_from_qimage( const QImage &input_image  );

};

#endif // ABSTRACTIONGRAINPROCESS_H
