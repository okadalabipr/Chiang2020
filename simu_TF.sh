#!/bin/bash

pwd=`pwd`

while IFS=$'\t' read -r -a tempTerm
do
    ./bin/simu_TF ${pwd}/input_TF/${tempTerm[0]}_param ${pwd}/output_TF/${tempTerm[0]}_output.txt
done < para_dep/TFlist

