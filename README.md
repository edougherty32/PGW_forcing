This code is for generating PGW forcing for the mean of LENS2 ensembles. Here are the following steps:
1. Edit lens2_pgw_forcing_3D.py 
  This script takes the monthly mean difference between a historical and future period from the ensemble mean of LENS2 for 3D variables
  Change historical period, future period, and month of difference, and bounds for domain
  Exports netcdfs of all variables
2. > qsub submit_python-casper_3d_forcing.sh
  This runs lens2_pgw_forcing_3D.py
3. Edit lens2_pgw_forcing_sfc.py 
  Same as above, but for 2D variables
4. > qsub submit_python-casper_sfc_forcing.sh
  This runs lens2_pgw_forcing_sfc.py
5. Edit + run ERA5_perturb_sfc_new.ipynb
  This script perturbs the ERA5 files with the PGW forcing and outputs modified ERA5 file to run MPAS
horizontally reproject LENS2 to ERA5 
6. Edit ERA5_perturb_3D_new.py
  This script perturbs the ERA5 files with the PGW forcing and outputs modified ERA5 file to run MPAS
  Vertically interpolate LENS2 to ERA5 and horizontally reproject LENS2 to ERA5 
> qsub submit_python-casper_ptb_3D.sh
