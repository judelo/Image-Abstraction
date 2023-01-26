#!/usr/bin/env bash
set -e

task=${1}
model=${2}
alternative_model=${3}
options_shapeabstraction=${4}
color_sketch_shapeabstraction=${5}
minarea_shapeabstraction=${6}
maxarea_shapeabstraction=${7}
scaleratio_shapeabstraction=${8}
threshold_shapeabstraction=${9}
eps_shapeabstraction=${10}
renderOrder_shapeabstraction=${11}
alpha_shapeabstraction=${12}
options_watercolor=${13}
color_sketch_watercolor=${14}
minarea_watercolor=${15}
maxarea_watercolor=${16}
scaleratio_watercolor=${17}
threshold_watercolor=${18}
eps_watercolor=${19}
renderOrder_watercolor=${20}
alpha_watercolor=${21}
options_shaking=${22}
color_sketch_shaking=${23}
minarea_shaking=${24}
maxarea_shaking=${25}
scaleratio_shaking=${26}
threshold_shaking=${27}
eps_shaking=${28}
renderOrder_shaking=${29}
alpha_shaking=${30}
options_shapefiltering=${31}
color_sketch_shapefiltering=${32}
minarea_shapefiltering=${33}
maxarea_shapefiltering=${34}
scaleratio_shapefiltering=${35}
threshold_shapefiltering=${36}
eps_shapefiltering=${37}
renderOrder_shapefiltering=${38}
alpha_shapefiltering=${39}
options_styletransfer=${40}
color_sketch_styletransfer=${41}
minarea_styletransfer=${42}
maxarea_styletransfer=${43}
scaleratio_styletransfer=${44}
threshold_styletransfer=${45}
eps_styletransfer=${46}
renderOrder_styletransfer=${47}
alpha_styletransfer=${48}
modelDictionary=${49}
mcolor=${50}
color_background=${51}
equal=${52}
kappaDict=${53}

echo $bin/image_abstraction input_0.png $task $model $alternative_model $options_shapeabstraction $color_sketch_shapeabstraction $minarea_shapeabstraction $maxarea_shapeabstraction $scaleratio_shapeabstraction $threshold_shapeabstraction $eps_shapeabstraction $renderOrder_shapeabstraction $alpha_shapeabstraction $options_watercolor $color_sketch_watercolor $minarea_watercolor $maxarea_watercolor $scaleratio_watercolor $threshold_watercolor $eps_watercolor $renderOrder_watercolor $alpha_watercolor $options_shaking $color_sketch_shaking $minarea_shaking $maxarea_shaking $scaleratio_shaking $threshold_shaking $eps_shaking $renderOrder_shaking $alpha_shaking $options_shapefiltering $color_sketch_shapefiltering $minarea_shapefiltering $maxarea_shapefiltering $scaleratio_shapefiltering $threshold_shapefiltering $eps_shapefiltering $renderOrder_shapefiltering $alpha_shapefiltering $options_styletransfer $color_sketch_styletransfer $minarea_styletransfer $maxarea_styletransfer $scaleratio_styletransfer $threshold_styletransfer $eps_styletransfer $renderOrder_styletransfer $alpha_styletransfer $modelDictionary $mcolor $color_background $equal $kappaDict input_1.png mask_0.png

$bin/image_abstraction input_0.png $task $model $alternative_model $options_shapeabstraction $color_sketch_shapeabstraction $minarea_shapeabstraction $maxarea_shapeabstraction $scaleratio_shapeabstraction $threshold_shapeabstraction $eps_shapeabstraction $renderOrder_shapeabstraction $alpha_shapeabstraction $options_watercolor $color_sketch_watercolor $minarea_watercolor $maxarea_watercolor $scaleratio_watercolor $threshold_watercolor $eps_watercolor $renderOrder_watercolor $alpha_watercolor $options_shaking $color_sketch_shaking $minarea_shaking $maxarea_shaking $scaleratio_shaking $threshold_shaking $eps_shaking $renderOrder_shaking $alpha_shaking $options_shapefiltering $color_sketch_shapefiltering $minarea_shapefiltering $maxarea_shapefiltering $scaleratio_shapefiltering $threshold_shapefiltering $eps_shapefiltering $renderOrder_shapefiltering $alpha_shapefiltering $options_styletransfer $color_sketch_styletransfer $minarea_styletransfer $maxarea_styletransfer $scaleratio_styletransfer $threshold_styletransfer $eps_styletransfer $renderOrder_styletransfer $alpha_styletransfer $modelDictionary $mcolor $color_background $equal $kappaDict input_1.png mask_0.png