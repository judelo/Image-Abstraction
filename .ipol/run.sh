#!/usr/bin/env bash
set -e

task=${1}
model=${2}
options=${3}
color_sketch=${4}
renderOrder=${5}
alpha=${6}
modelDictionary=${7}
mcolor=${8}
equal=${9}
kappaDict=${10}
minSize=${11}
mpixel=${12}
scaleratio=${13}
threshold=${14}
eps=${15}

echo $bin/image_abstraction input_0.png $task $model $options $color_sketch $renderOrder $alpha $modelDictionary $mcolor $equal $kappaDict $mpixel $maxarea $scaleratio $threshold $eps input_1.png

$bin/image_abstraction input_0.png $task $model $options  $color_sketch $renderOrder $alpha $modelDictionary $mcolor $equal $kappaDict $mpixel $maxarea $scaleratio $threshold $eps input_1.png
