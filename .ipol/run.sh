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
maxarea=${13}

echo $bin/image_abstraction input_0.png $task $model $options $color_sketch $renderOrder $alpha $modelDictionary $mcolor $equal $kappaDict $mpixel $maxarea input_1.png

$bin/image_abstraction input_0.png $task $model $options  $color_sketch $renderOrder $alpha $modelDictionary $mcolor $equal $kappaDict $mpixel $maxarea input_1.png
