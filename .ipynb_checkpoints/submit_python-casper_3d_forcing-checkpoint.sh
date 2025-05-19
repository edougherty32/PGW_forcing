#!/bin/bash

#PBS -N get_delta_3d
#PBS -A P48500028
#PBS -l select=1:ncpus=2
#PBS -l walltime=24:00:00
#PBS -q casper
#PBS -j oe                    

source /glade/u/home/doughert/miniconda3/bin/activate pangeo3

ulimit -s unlimited

python lens2_pgw_forcing_3D.py

