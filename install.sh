#!/bin/bash

pwd=`pwd`

cd cvode/src
make

cd ../../model_TF
make
make clean

cd ../model_mRNA
cp simulation_A.C simulation.C
make
make clean
mv ../bin/simu_mRNA ../bin/simu_mRNA_A

cp simulation_B.C simulation.C
make
make clean
mv ../bin/simu_mRNA ../bin/simu_mRNA_B

cp simulation_C.C simulation.C
make
make clean
mv ../bin/simu_mRNA ../bin/simu_mRNA_C

rm simulation.C

cd ../cvode/src
make clean

cd ../..
rm lib/libsundials_cvode.a

if [ ! -d "input_TF" ]; then
    mkdir input_TF
else
   rm -rf input_TF/*
fi

if [ ! -d "input_mRNA" ]; then
    mkdir input_mRNA
else
   rm -rf input_mRNA/*
fi

if [ ! -d "output_TF" ]; then
    mkdir output_TF
else
   rm -rf output_TF/*
fi

if [ ! -d "output_mRNA" ]; then
    mkdir output_mRNA
else
   rm -rf output_mRNA/*
fi

perl prod_input_TF.pl
perl prod_input_mRNA.pl

chmod 0755 simu_TF.sh
chmod 0755 simu_mRNA.sh

