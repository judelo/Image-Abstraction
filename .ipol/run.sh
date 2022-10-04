#!/usr/bin/env bash
set -e

task=${1}
model=${2}
options=${3}
segmentation=${4}
color_sketch=${5}
renderOrder=${6}
alpha=${7}

echo $bin/image_abstraction input_0.png $task $model $options $segmentation $color_sketch $renderOrder $alpha input_1.png

$bin/image_abstraction input_0.png $task $model $options $segmentation $color_sketch $renderOrder $alpha input_1.png
