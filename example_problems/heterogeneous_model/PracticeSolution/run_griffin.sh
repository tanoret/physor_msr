#!/bin/bash
#PBS -N griffin
#PBS -P neams
#PBS -l select=2:ncpus=48:mem=170gb
#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o output_griffin
#PBS -e error_griffin

module load use.moose moose-apps griffin

input=griffin_prac.in
mesh=lattice.e
xslib=../OpenMC/isoxml.xml
ncpu=96

cd $PBS_O_WORKDIR
mkdir tmprun
cp $PBS_O_WORKDIR/$input ./tmprun/
ln -s $mesh ./tmprun/mesh.e
ln -s $xslib ./tmprun/isoxml
cd tmprun
mpirun -n $ncpu griffin-opt -i $input | tee griffin.out
