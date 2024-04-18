#!/bin/bash
#PBS -N griffin
#PBS -P neams
#PBS -l select=1:ncpus=40:mem=170gb
#PBS -l walltime=1:00:00
#PBS -j oe
#PBS -o output_griffin
#PBS -e error_griffin

module unload griffin
module load use.moose moose-apps griffin

input=prac.i
ncpu=40

cd $PBS_O_WORKDIR
mkdir tmp_$$
cd tmp_$$
cp $PBS_O_WORKDIR/$input .
mpirun -n $ncpu griffin-opt -i $input --mesh-only | tee griffin.out

