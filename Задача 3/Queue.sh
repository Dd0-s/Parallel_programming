#!/bin/bash

#PBS -l walltime=00:10:00,nodes=1:ppn=4
#PBS -N Bogdanov
#PBS -q batch

cd $PBS_O_WORKDIR
./a.out
