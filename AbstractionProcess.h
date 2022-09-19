/*--------------------------- MegaWave2 Module -----------------------------*/
/* mwcommand
  name = {fgrain_side};
  version = {"1.2"};
  author = {"Pascal Monasse, Frederic Guichard, G. Facciolo"};
  function = {"Grain filter of an image"};
  usage = {
    'a': [min_area=20]-> pMinArea   "Min area of grains we keep",
    image_in -> pFloatImageInput    "Input fimage",
    image_out <- pFloatImageOutput  "Output fimage"
    };
*/
/*----------------------------------------------------------------------
 v1.2 (04/2007): simplified header (LM)
 v1.3 (2013): portable version (GF)
 v1.4 (2014): fgrain_side: removes only upper (bright) or lower (dark) courves
----------------------------------------------------------------------*/

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

    QImage render(TOSParameters tosParameters, bool &tree_recomputed, DictionaryParameters dictionaryParameters = getDefaultDictionaryParameters(), TreeOfShapes * dictionnary=NULL);

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
