#!/usr/bin/env bash
set -e

task=${1}
model=${2}
options=${3}
segmentation=${4}
color_sketch=${5}
renderOrder=${6}
alpha=${7}
modelDictionary=${8}
mcolor=${9}
equal=${10}
kappaDict=${11}
minSize=${12}
mpixel=${13}

echo $bin/image_abstraction input_0.png $task $model $options $segmentation $color_sketch $renderOrder $alpha $modelDictionary $mcolor $equal $kappaDict $minSize $mpixel input_1.png

$bin/image_abstraction input_0.png $task $model $options $segmentation $color_sketch $renderOrder $alpha $modelDictionary $mcolor $equal $kappaDict $minSize $mpixel input_1.png
