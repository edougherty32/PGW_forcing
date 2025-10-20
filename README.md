## This set of scripts generate PGW forcing from the ensemble mean of LENS2 that is added to ERA5 reanalysis data. It is assumed that you have already downloaded ERA5 data on pressure levels. The final output are the perturbed ERA5 files output in the format needed to run WRF/MPAS. 

### Please run the scripts in the following order:
**1. Edit + run lens2_pgw_forcing_3D.ipynb**
- This script takes the monthly mean difference between a historical and future period from the ensemble mean of LENS2 for 3D variables
- Change historical period, future period, and month of difference, and bounds for domain
- Exports netcdfs of all variables

**2. Edit + run lens2_pgw_forcing_sfc.ipynb**
- Same as above, but for 2D variables

**3. Edit + run ERA5_perturb_sfc_new.ipynb**
- This script perturbs the ERA5 files with the PGW forcing and outputs modified ERA5 file to run MPAS
horizontally reproject LENS2 to ERA5 

**4. Edit ERA5_perturb_3D_new.py**
- This script perturbs the ERA5 files with the PGW forcing and outputs modified ERA5 file to run MPAS
- Vertically interpolate LENS2 to ERA5 and horizontally reproject LENS2 to ERA5 

**5. > qsub submit_python-casper_ptb_3D.sh**
- Need to submit this via job because it much more memory intensive

## Please cite Dougherty et al. (2023) if you use this code in your project
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022JD037537







