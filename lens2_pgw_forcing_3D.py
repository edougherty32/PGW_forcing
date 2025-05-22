#!/usr/bin/env python
# coding: utf-8

# ### Script for calculating CESM-LENS2 climate forcing to force MPAS simulations
# ### date created: 10 Nov 2023
# ### author: Erin Dougherty (doughert@ucar.edu)

# In[1]:


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


# ### set path to monthly data and export dir

# In[2]:


cesm_dir= '/glade/campaign/cgd/cesm/CESM2-LE/timeseries/atm/proc/tseries/month_1/'
outdir = ''


# ### 3D variables

# In[3]:


cesm_dir_t = cesm_dir+'T/'
cesm_dir_u = cesm_dir+'U/'
cesm_dir_v = cesm_dir+'V/'
cesm_dir_rh = cesm_dir+'RELHUM/'
cesm_dir_q = cesm_dir+'Q/'
cesm_dir_z = cesm_dir+'Z3/'


# ### define time slices- CHANGE month of forcing (changing 30-year slice is optional)

# In[4]:


def is_hist(year, month):
    return (year >= 1991) & (year <= 2021) & (month ==9)


# In[5]:


def is_future(year, month):
    return (year >= 2070) & (year <= 2100) & (month ==9)


# ### set bounds for domain-CHANGE if you want this over different location

# In[6]:


llat = -25
ulat = 60
llon = -74
rlon = 70


# ### open historical and future 3D variables

# In[8]:


hist_dates = ['199001', '200001', '201001', '201501']
future_dates = ['206501', '207501', '208501', '209501']


# In[9]:


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


# In[10]:


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


# In[ ]:


print('opening historical and future variables')


# In[11]:


import dask
# this avoids creating large chunks
dask.config.set({"array.slicing.split_large_chunks": True})

hist_t = cesm_hist(cesm_dir_t, 'T')
hist_u = cesm_hist(cesm_dir_u, 'U')
hist_v = cesm_hist(cesm_dir_v, 'V')
hist_q = cesm_hist(cesm_dir_q, 'Q')
hist_rh = cesm_hist(cesm_dir_rh, 'RELHUM')
hist_z = cesm_hist(cesm_dir_z, 'Z3')


# In[15]:


import dask
# this avoids creating large chunks
dask.config.set({"array.slicing.split_large_chunks": True})

future_t = cesm_future(cesm_dir_t, 'T')
future_u = cesm_future(cesm_dir_u, 'U')
future_v = cesm_future(cesm_dir_v, 'V')
future_q = cesm_future(cesm_dir_q, 'Q')
future_rh = cesm_future(cesm_dir_rh, 'RELHUM')
future_z = cesm_future(cesm_dir_z, 'Z3')


# In[16]:


print('concatenating across times')


# ### concat members/times together and take mean

# In[17]:


hist_t_all = xr.concat(hist_t, dim='member_time')
hist_t_mean = hist_t_all.mean(dim=['member_time', 'time'])

hist_u_all = xr.concat(hist_u, dim='member_time')
hist_u_mean = hist_u_all.mean(dim=['member_time', 'time'])

hist_v_all = xr.concat(hist_v, dim='member_time')
hist_v_mean = hist_v_all.mean(dim=['member_time', 'time'])

hist_q_all = xr.concat(hist_q, dim='member_time')
hist_q_mean = hist_q_all.mean(dim=['member_time', 'time'])

hist_rh_all = xr.concat(hist_rh, dim='member_time')
hist_rh_mean = hist_rh_all.mean(dim=['member_time', 'time'])

hist_z_all = xr.concat(hist_z, dim='member_time')
hist_z_mean = hist_z_all.mean(dim=['member_time', 'time'])


# In[22]:


future_t_all = xr.concat(future_t, dim='member_time')
future_t_mean = future_t_all.mean(dim=['member_time', 'time'])

future_u_all = xr.concat(future_u, dim='member_time')
future_u_mean = future_u_all.mean(dim=['member_time', 'time'])

future_v_all = xr.concat(future_v, dim='member_time')
future_v_mean = future_v_all.mean(dim=['member_time', 'time'])

future_q_all = xr.concat(future_q, dim='member_time')
future_q_mean = future_q_all.mean(dim=['member_time', 'time'])

future_rh_all = xr.concat(future_rh, dim='member_time')
future_rh_mean = future_rh_all.mean(dim=['member_time', 'time'])

future_z_all = xr.concat(future_z, dim='member_time')
future_z_mean = future_z_all.mean(dim=['member_time', 'time'])


# In[ ]:


print('taking difference')


# ### future - historical difference 

# In[28]:


delta_t_mean = future_t_mean - hist_t_mean
delta_u_mean = future_u_mean - hist_u_mean
delta_v_mean = future_v_mean - hist_v_mean
delta_z_mean = future_z_mean - hist_z_mean
delta_q_mean = future_q_mean - hist_q_mean
delta_rh_mean = future_rh_mean - hist_rh_mean


# ### load variables into memory

# In[ ]:


delta_t_mean_ar = delta_t_mean.compute()
delta_u_mean_ar = delta_u_mean.compute()
delta_v_mean_ar = delta_v_mean.compute()
delta_z_mean_ar = delta_z_mean.compute()
delta_q_mean_ar = delta_q_mean.compute()
delta_rh_mean_ar = delta_rh_mean.compute()


# In[33]:


print('variables computed')


# ### export to netcdfs

# In[ ]:


delta_t_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_T.nc', compute=False)
delta_u_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_U.nc', compute=False)
delta_v_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_V.nc', compute=False)
delta_z_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_Z.nc', compute=False)
delta_q_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_Q.nc', compute=False)
delta_rh_mean_ar.to_netcdf(outdir+'LENS2-Sept_2070-2100_1991-2021_delta_RH.nc', compute=False)


# In[ ]:


print('netcdfs exported')

