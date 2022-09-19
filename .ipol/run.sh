#!/usr/bin/env bash
set -e

task=${1}
model=${2}
segmentation=${3}

echo ./image_abstraction input_0.png $task $model $segmentation input_1.png

cp /workdir/bin/image_abstraction /workdir/exec/image_abstraction

./image_abstraction input_0.png $task $model $segmentation input_1.png
