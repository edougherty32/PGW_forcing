### Script for adding CESM-LENS2 perturbation to ERA5 data
### date created: 27 Nov. 2023
### author: Erin Dougherty (doughert@ucar.edu) 


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


### define time slices - CHANGE month and 30-year slices (optional)
def is_hist(year, month): 
    return (year >= 1991) & (year <= 2021) & (month ==9)

def is_future(year, month):
    return (year >= 2070) & (year <= 2100) & (month ==9)

### set bounds for domain - CHANGE based on area covered by MPAS
llat = -25
ulat = 60
llon = -74
rlon = 70

### change where you will output the data
outdir = '/glade/scratch/doughert/CESM/LENS2/'

######## DON'T CHANGE BELOW HERE ########

### set path to monthly data
cesm_dir= '/glade/campaign/cgd/cesm/CESM2-LE/timeseries/atm/proc/tseries/month_1/'
cesm_land= '/glade/campaign/cgd/cesm/CESM2-LE/timeseries/lnd/proc/tseries/month_1/'

### 2D variables
cesm_dir_sst = cesm_dir+'SST/'
cesm_dir_ts = cesm_dir+'TS/'
cesm_dir_tsoil = cesm_land+'TSOI/'
cesm_dir_ps = cesm_dir+'PS/'
cesm_dir_psl = cesm_dir+'PSL/'
cesm_dir_phi = cesm_dir+'PHIS/'
cesm_dir_ice = cesm_dir+'ICEFRAC/'
cesm_dir_snow = cesm_dir+'SNOWHLND/'


### open historical and future variables
hist_dates = ['199001', '200001', '201001', '201501']
future_dates = ['206501', '207501', '208501', '209501']

def cesm_hist(var_dir, varname):
    hist_var = []

    for dir in glob.glob(var_dir+'b.e21.*'):
        for c, item in enumerate(hist_dates): 
            for name in glob.glob(dir):
                if item in name:
                    file = xr.open_mfdataset(name)[varname]
                    #change longitude from 0-360, to -180 to 180
                    file['_longitude_adjusted'] = xr.where(file['lon'] > 180, file['lon']-360, file['lon'])
                    file = (file.swap_dims({'lon': '_longitude_adjusted'}).sel(**{'_longitude_adjusted': sorted(file._longitude_adjusted)}).drop('lon'))
                    file = file.rename({'_longitude_adjusted': 'lon'})
                    # sub select 
                    file_sub = file.sel(time=is_hist(file['time.year'], file['time.month']), lat=slice(llat, ulat), lon=slice(llon, rlon))
                    hist_var.append(file_sub)

    return(hist_var)

def cesm_future(var_dir, varname):
    future_var = []

    for dir in glob.glob(var_dir+'b.e21.*'):
        for c, item in enumerate(future_dates): 
            for name in glob.glob(dir):
                if item in name:
                    file = xr.open_mfdataset(name)[varname]
                    #change longitude from 0-360, to -180 to 180
                    file['_longitude_adjusted'] = xr.where(file['lon'] > 180, file['lon']-360, file['lon'])
                    file = (file.swap_dims({'lon': '_longitude_adjusted'}).sel(**{'_longitude_adjusted': sorted(file._longitude_adjusted)}).drop('lon'))
                    file = file.rename({'_longitude_adjusted': 'lon'})
                    # sub select 
                    file_sub = file.sel(time=is_future(file['time.year'], file['time.month']), lat=slice(llat, ulat), lon=slice(llon, rlon))
                    future_var.append(file_sub)

    return(future_var)


import dask
# this avoids creating large chunks
dask.config.set({"array.slicing.split_large_chunks": True})

hist_sst = cesm_hist(cesm_dir_sst, 'SST')
hist_ts = cesm_hist(cesm_dir_ts, 'TS')
hist_tsoil = cesm_hist(cesm_dir_tsoil, 'TSOI')
hist_ps = cesm_hist(cesm_dir_ps, 'PS')
hist_psl = cesm_hist(cesm_dir_psl, 'PSL')
hist_phi = cesm_hist(cesm_dir_phi, 'PHIS')
hist_ice = cesm_hist(cesm_dir_ice, 'ICEFRAC')
hist_snow = cesm_hist(cesm_dir_snow, 'SNOWHLND')

future_sst = cesm_future(cesm_dir_sst, 'SST')
future_ts = cesm_future(cesm_dir_ts, 'TS')
future_tsoil = cesm_future(cesm_dir_tsoil, 'TSOI')
future_ps = cesm_future(cesm_dir_ps, 'PS')
future_psl = cesm_future(cesm_dir_psl, 'PSL')
future_phi = cesm_future(cesm_dir_phi, 'PHIS')
future_ice = cesm_future(cesm_dir_ice, 'ICEFRAC')
future_snow = cesm_future(cesm_dir_snow, 'SNOWHLND')

### concat all members together/take mean
hist_sst_all = xr.concat(hist_sst, dim='member_time')
hist_sst_mean = hist_sst_all.mean(dim=['member_time', 'time'])
hist_ts_all = xr.concat(hist_ts, dim='member_time')
hist_ts_mean = hist_ts_all.mean(dim=['member_time', 'time'])
hist_tsoil_all = xr.concat(hist_tsoil, dim='member_time')
hist_tsoil_mean = hist_tsoil_all.mean(dim=['member_time', 'time'])
hist_ps_all = xr.concat(hist_ps, dim='member_time')
hist_ps_mean = hist_ps_all.mean(dim=['member_time', 'time'])
hist_psl_all = xr.concat(hist_psl, dim='member_time')
hist_psl_mean = hist_psl_all.mean(dim=['member_time', 'time'])
hist_phi_all = xr.concat(hist_phi, dim='member_time')
hist_phi_mean = hist_phi_all.mean(dim=['member_time', 'time'])
hist_ice_all = xr.concat(hist_ice, dim='member_time')
hist_ice_mean = hist_ice_all.mean(dim=['member_time', 'time'])
hist_snow_all = xr.concat(hist_snow, dim='member_time')
hist_snow_mean = hist_snow_all.mean(dim=['member_time', 'time'])

future_sst_all = xr.concat(future_sst, dim='member_time')
future_sst_mean = future_sst_all.mean(dim=['member_time', 'time'])
future_ts_all = xr.concat(future_ts, dim='member_time')
future_ts_mean = future_ts_all.mean(dim=['member_time', 'time'])
future_tsoil_all = xr.concat(future_tsoil, dim='member_time')
future_tsoil_mean = future_tsoil_all.mean(dim=['member_time', 'time'])
future_ps_all = xr.concat(future_ps, dim='member_time')
future_ps_mean = future_ps_all.mean(dim=['member_time', 'time'])
future_psl_all = xr.concat(future_psl, dim='member_time')
future_psl_mean = future_psl_all.mean(dim=['member_time', 'time'])
future_phi_all = xr.concat(future_phi, dim='member_time')
future_phi_mean = future_phi_all.mean(dim=['member_time', 'time'])
future_ice_all = xr.concat(future_ice, dim='member_time')
future_ice_mean = future_ice_all.mean(dim=['member_time', 'time'])
future_snow_all = xr.concat(future_snow, dim='member_time')
future_snow_mean = future_snow_all.mean(dim=['member_time', 'time'])

### delete variables to free up some memory
del hist_sst_all
del hist_ts_all
del hist_tsoil_all
del hist_phi_all
del hist_ice_all
del hist_snow_all

del future_sst_all
del future_ts_all
del future_tsoil_all
del future_phi_all
del future_ice_all
del future_snow_all

### future - historical difference
delta_sst_mean = future_sst_mean - hist_sst_mean
delta_ts_mean = future_ts_mean - hist_ts_mean
delta_tsoil_mean = future_tsoil_mean - hist_tsoil_mean
delta_ps_mean = future_ps_mean - hist_ps_mean
delta_psl_mean = future_psl_mean - hist_psl_mean
delta_phi_mean = future_phi_mean - hist_phi_mean
delta_ice_mean = future_ice_mean - hist_ice_mean
delta_snow_mean = future_snow_mean - hist_snow_mean

### load variables into memory
delta_sst_mean_ar = delta_sst_mean.compute()
delta_ts_mean_ar = delta_ts_mean.compute()
delta_tsoil_mean_ar = delta_tsoil_mean.compute()
delta_ps_mean_ar = delta_ps_mean.compute()
delta_psl_mean_ar = delta_psl_mean.compute()
delta_phi_mean_ar = delta_phi_mean.compute()
delta_ice_mean_ar = delta_ice_mean.compute()
delta_snow_mean_ar = delta_snow_mean.compute()

### export to NETCDFs
delta_sst_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_sst.nc', compute=False)
delta_ts_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_ts.nc', compute=False)
delta_tsoil_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_tsoil.nc', compute=False)
delta_ps_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_ps.nc', compute=False)
delta_psl_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_psl.nc', compute=False)
delta_phi_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_phi.nc', compute=False)
delta_ice_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_ice.nc', compute=False)
delta_snow_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_snow.nc', compute=False)

