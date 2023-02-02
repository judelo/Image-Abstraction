#!/usr/bin/env bash
set -e

task=${1}
model=${2}
segmentWithMask=${3}
alternative_model=${4}
options_shapeabstraction=${5}
color_sketch_shapeabstraction=${6}
minarea_shapeabstraction=${7}
maxarea_shapeabstraction=${8}
scaleratio_shapeabstraction=${9}
threshold_shapeabstraction=${10}
eps_shapeabstraction=${11}
renderOrder_shapeabstraction=${12}
alpha_shapeabstraction=${13}
options_watercolor=${14}
color_sketch_watercolor=${15}
minarea_watercolor=${16}
maxarea_watercolor=${17}
scaleratio_watercolor=${18}
threshold_watercolor=${19}
eps_watercolor=${20}
renderOrder_watercolor=${21}
alpha_watercolor=${22}
options_shaking=${23}
color_sketch_shaking=${24}
minarea_shaking=${25}
maxarea_shaking=${26}
scaleratio_shaking=${27}
threshold_shaking=${28}
eps_shaking=${29}
renderOrder_shaking=${30}
alpha_shaking=${31}
options_shapefiltering=${32}
color_sketch_shapefiltering=${33}
minarea_shapefiltering=${34}
maxarea_shapefiltering=${35}
scaleratio_shapefiltering=${36}
threshold_shapefiltering=${37}
eps_shapefiltering=${38}
renderOrder_shapefiltering=${39}
alpha_shapefiltering=${40}
options_styletransfer=${41}
color_sketch_styletransfer=${42}
minarea_styletransfer=${43}
maxarea_styletransfer=${44}
scaleratio_styletransfer=${45}
threshold_styletransfer=${46}
eps_styletransfer=${47}
renderOrder_styletransfer=${48}
alpha_styletransfer=${49}
modelDictionary=${50}
mcolor=${51}
color_background=${52}
equal=${53}
kappaDict=${54}

echo $bin/image_abstraction input_0.png $task $model $segmentWithMask $alternative_model $options_shapeabstraction $color_sketch_shapeabstraction $minarea_shapeabstraction $maxarea_shapeabstraction $scaleratio_shapeabstraction $threshold_shapeabstraction $eps_shapeabstraction $renderOrder_shapeabstraction $alpha_shapeabstraction $options_watercolor $color_sketch_watercolor $minarea_watercolor $maxarea_watercolor $scaleratio_watercolor $threshold_watercolor $eps_watercolor $renderOrder_watercolor $alpha_watercolor $options_shaking $color_sketch_shaking $minarea_shaking $maxarea_shaking $scaleratio_shaking $threshold_shaking $eps_shaking $renderOrder_shaking $alpha_shaking $options_shapefiltering $color_sketch_shapefiltering $minarea_shapefiltering $maxarea_shapefiltering $scaleratio_shapefiltering $threshold_shapefiltering $eps_shapefiltering $renderOrder_shapefiltering $alpha_shapefiltering $options_styletransfer $color_sketch_styletransfer $minarea_styletransfer $maxarea_styletransfer $scaleratio_styletransfer $threshold_styletransfer $eps_styletransfer $renderOrder_styletransfer $alpha_styletransfer $modelDictionary $mcolor $color_background $equal $kappaDict input_1.png mask_0.png

$bin/image_abstraction input_0.png $task $model $segmentWithMask $alternative_model $options_shapeabstraction $color_sketch_shapeabstraction $minarea_shapeabstraction $maxarea_shapeabstraction $scaleratio_shapeabstraction $threshold_shapeabstraction $eps_shapeabstraction $renderOrder_shapeabstraction $alpha_shapeabstraction $options_watercolor $color_sketch_watercolor $minarea_watercolor $maxarea_watercolor $scaleratio_watercolor $threshold_watercolor $eps_watercolor $renderOrder_watercolor $alpha_watercolor $options_shaking $color_sketch_shaking $minarea_shaking $maxarea_shaking $scaleratio_shaking $threshold_shaking $eps_shaking $renderOrder_shaking $alpha_shaking $options_shapefiltering $color_sketch_shapefiltering $minarea_shapefiltering $maxarea_shapefiltering $scaleratio_shapefiltering $threshold_shapefiltering $eps_shapefiltering $renderOrder_shapefiltering $alpha_shapefiltering $options_styletransfer $color_sketch_styletransfer $minarea_styletransfer $maxarea_styletransfer $scaleratio_styletransfer $threshold_styletransfer $eps_styletransfer $renderOrder_styletransfer $alpha_styletransfer $modelDictionary $mcolor $color_background $equal $kappaDict input_1.png mask_0.png