#!/bin/bash

input=griffin_core.in
mesh=core.e
xslib=../OpenMC/isoxml.xml
ncpu=6

mkdir tmprun
cp $input ./tmprun/
ln -s $mesh ./tmprun/mesh.e
ln -s $xslib ./tmprun/isoxml
cd tmprun
mpirun -n $ncpu griffin-opt -i $input | tee griffin.out
