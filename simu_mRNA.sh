#!/bin/bash

pwd=`pwd`

synr=("A" "B" "C")
tf=("S" "E" "Q" "I")
mrna=("M10" "M15" "M30" "M60" "M120" "M180")

for i in `seq 0 1 2`; do
  for j in `seq 0 1 5`; do
    for m in `seq 0 1 3`; do
      for n in `seq ${m} 1 3`; do
	fileName=${synr[$i]}_${mrna[$j]}_${tf[$m]}${tf[$n]}
	./bin/simu_mRNA_${synr[$i]} ${pwd}/input_mRNA/${fileName}_param ${pwd}/output_mRNA/${fileName}_output.txt
      done
    done
  done
done

