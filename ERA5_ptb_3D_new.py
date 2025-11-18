### Script for adding CESM-LENS perturbation to ERA5 3D files 
### date created: 13 May 2025
### author: doughert@ucar.edu

import math
import numpy as np
import pandas as pd
import matplotlib as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# import netCDF4 as nc
# from netCDF4 import Dataset, num2date
#from datetime import datetime, date, timedelta
import glob
import xarray as xr

### change year, month, day as well as export path
yr_mo = '201709'
days = ['15', '16', '17','18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29','30',]
export_dir = '/'
input_dir = '/'

#### var list from CESM delta outputs - DON'T CHANGE BELOW HERE
var_list = ['T', 'U',  'V', 'Z3', 'RELHUM', ]

#### list ERA5 variables that will be perturbed
era5_pl = '/glade/campaign/collections/rda/data/ds633.0/e5.oper.an.pl/'
var_era = ['128_130_t.ll025sc', '128_131_u.ll025uv','128_132_v.ll025uv', '128_129_z.ll025sc', '128_157_r.ll025sc', ]
#era5_name = ['SSTK', 'SKT', 'SP', 'MSL', 'CI', 'SD']
param_ids = ["130.128", "131.128", "132.128", "129.128", "157.128",]

### make definitions
#### interpolate CESM deltas to ERA5 levels

import wrf

def vert_interp(var_3d):
    # this defines shape of cesm variable with vertical 
    # dimension matching era5 (levs = 37)
    var_interp = np.zeros((len(era5_lvls), len(var_3d[0]), len(var_3d[0][0]))) 
    
    for x in range(len(var_3d[0])):
        for y in range(len(var_3d[0][0])):
            var_interp[:,x,y] = wrf.interp1d(var_3d[:,x, y], cesm_levs, era5_lvls, missing=np.nan)

    return(var_interp)

#### replace original ERA5 variable with perturbed variable
def replace_era5_var(orig_era5, var, ptb_era5):
    era5_var_list = []

    for i in range(len(ptb_era5)):
        era5_rp = orig_era5[i].assign({var:ptb_era5[i]})
        era5_var_list.append(era5_rp)
           
    return(era5_var_list)

#### reformat files in the way MPAS expects for input
def fmt_files(var, pid):
    from cdo import Cdo
    cdo = Cdo(cdo='/glade/u/home/doughert/miniconda3/envs/pangeo3/bin/cdo')
    
    ### select first variable
    substr = '.'+yr_mo
    add_str3 = 'var1'

    for file in sorted(glob.glob(export_dir+'era5.pl.'+var+'.*.nc')):
        stridx = file.index(substr)
        newfname = str(file[:stridx+9]) + add_str3 + str(file[stridx+9:]) 
        newfname = newfname[:-2]+'grb'
        print(newfname)
        cdo.selparam("-1", input=file, output=newfname, options='-f grb')

    ## change z-axis/levels to Pa
    zaxis_file = 'plev_zaxis'
    substr = 'var1'
    add_str3 = 'zaxis'

    for file in sorted(glob.glob(export_dir+'era5.pl.'+var+'.*var1*.grb')):
        stridx = file.index(substr)
        newfname = str(file[:stridx]) + add_str3 + str(file[stridx+3:]) 
        #print(newfname)
        cdo.setzaxis(export_dir+zaxis_file, input=file, output=newfname)
     
    ## change parameter id code 
    substr = 'zaxis'
    add_str3 = 'pid'
    
    for file in sorted(glob.glob(export_dir+'era5.pl.'+var+'.*zaxis1*.grb')):
        stridx = file.index(substr)
        newfname = str(file[:stridx]) + add_str3 + str(file[stridx+6:]) 
        #print(newfname)
        cdo.setparam(pid, input=file, output=newfname, )

    ## change dtype to P16
    substr = 'pid'
    add_str3 = 'final'
    
    for file in sorted(glob.glob(export_dir+'era5.pl.'+var+'.*pid*.grb')):
        stridx = file.index(substr)
        newfname = str(file[:stridx]) + add_str3 + str(file[stridx+3:]) 
        #print(newfname)
        cdo.copy(input=file, output=newfname, options="-b P16")

### this routine opens delta from LENS2, opens ERA5 for same variable, interpolates ERA5 and CESM to match
### then adds delta to ERA5, and export file with perturbed variables

import pyresample

for i in range(len(var_list)):
    print(var_list[i])
    ## open deltas
    delta_var = xr.open_dataset(input_dir+'LENS2-Sept_2070-2100_1991-2021_delta_'+var_list[i]+'.nc')[var_list[i]]
    
    ## convert geopotential height to geopotential
    if var_list[i]=='Z3':
        delta_var = delta_var*9.8
    cesm_lat = delta_var.lat
    cesm_lon = delta_var.lon
    cesm_levs = delta_var.lev
    
    ## open ERA5
    era5_var_days = []
    for d in days:
        file = xr.open_dataset(era5_pl+yr_mo+'/e5.oper.an.pl.'+var_era[i]+'.'+yr_mo+d+'00_'+yr_mo+d+'23.nc')
        era5_var_days.append(file)

    era5_lat = era5_var_days[0].latitude
    era5_lon = era5_var_days[0].longitude
    era5_lvls = era5_var_days[0].level.values
    era5_lon2d, era5_lat2d  = np.meshgrid(era5_lon.values, era5_lat.values)

    ## vertically interpolate
    delta_var_vert_interp = vert_interp(delta_var)
    # replace the last level with the surface data 
    delta_var_vert_interp[36] = delta_var[31]
    
    ## horizontally reproject data
    cesm_lon2d, cesm_lat2d = np.meshgrid(cesm_lon.values, cesm_lat.values)
    # change CESM lons to be from -180 to 180 instead of 0 to 360
    cesm_lon2d_new = pyresample.utils.check_and_wrap(cesm_lon2d, cesm_lat2d)[0]
    # create presample object of the CESM data
    orig_def = pyresample.geometry.SwathDefinition(cesm_lon2d_new, cesm_lat2d)
    # change ERA-5 lons to be from -180 to 180 instead of 0 to 360
    era5_lon2d_new = pyresample.utils.check_and_wrap(era5_lon2d, era5_lat2d)[0]
    # presample object of ERA-5 data
    targ_def = pyresample.geometry.SwathDefinition(era5_lon2d_new, era5_lat2d)

    re_cesm_delta_var = []
    for z in range(len(delta_var_vert_interp)):
        re_cesm_delta_var_3d_ar = pyresample.kd_tree.resample_nearest(orig_def, delta_var_vert_interp[z], targ_def, radius_of_influence=110000, fill_value=np.nan)
        re_cesm_delta_var.append(re_cesm_delta_var_3d_ar)   
        
    ## add delta to ERA5 for each day of interest
    if var_list[i]=='T' or var_list[i]=='U' or var_list[i]=='V':
        era5_times = len(era5_var_days[0][var_list[i]])
        era5_var_ptb = [era5_var_days[d][var_list[i]]+([np.nan_to_num(re_cesm_delta_var, nan=0.0)]*era5_times) for d in range(len(era5_var_days))]
    elif var_list[i]=='Z3':
        era5_times = len(era5_var_days[0]['Z'])
        era5_var_ptb = [era5_var_days[d]['Z']+([np.nan_to_num(re_cesm_delta_var, nan=0.0)]*era5_times) for d in range(len(era5_var_days))]
    elif var_list[i]=='RELHUM':
        era5_times = len(era5_var_days[0]['R'])
        era5_var_ptb = [era5_var_days[d]['R']+([np.nan_to_num(re_cesm_delta_var, nan=0.0)]*era5_times) for d in range(len(era5_var_days))]

    ## replace ERA5 with perturbed variable
    if var_list[i]=='T' or var_list[i]=='U' or var_list[i]=='V' or var_list[i]=='Z3':
        era5_var_change = replace_era5_var(era5_var_days, var_list[i], era5_var_ptb)
    else:
        era5_var_change = replace_era5_var(era5_var_days, 'R', era5_var_ptb)

    ## # export modified ERA5
    for d in range(len(era5_var_change)):
        era5_var_change[d].to_netcdf(export_dir+'era5.pl.'+var_list[i]+'.'+yr_mo+days[d]+'.nc')
        
    fmt_files(var_list[i], param_ids[i])
