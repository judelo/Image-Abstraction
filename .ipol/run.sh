#!/usr/bin/env bash
set -e

task=${1}
model=${2}
segmentation=${3}

echo ./image_abstraction input_0.png $task $model $segmentation input_1.png

echo ls /workdir/bin/
ls /workdir/bin/

mv /workdir/bin/image_abstraction /workdir/exec/image_abstraction

echo ls /workdir/exec
ls /workdir/exec

echo ls pwd
ls $pwd

pwd

./image_abstraction input_0.png $task $model $segmentation input_1.png
