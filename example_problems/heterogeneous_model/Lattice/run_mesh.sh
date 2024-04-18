#!/bin/bash


input=lattice.i
ncpu=6

mkdir tmp_$$
cd tmp_$$
cp ../$input .
mpirun -n $ncpu griffin-opt -i $input --mesh-only | tee griffin.out

