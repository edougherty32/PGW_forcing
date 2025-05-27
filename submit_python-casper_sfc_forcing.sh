#!/bin/bash

#PBS -N get_delta_sfc
#PBS -A 
#PBS -l select=1:ncpus=2:mem=200GB
#PBS -l walltime=24:00:00
#PBS -q casper
#PBS -j oe                    

source /glade/u/home/doughert/miniconda3/bin/activate pangeo3

ulimit -s unlimited

python lens2_pgw_forcing_sfc.py

