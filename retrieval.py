import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors

deg = unichr(176)

fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/SABER/NIGHTLY/atox_athy_night_YY2014_V5.3_fixedfnight_SV2.nc', 'r', format='NETCDF4')
o3 = fname.variables['o3conc'][:]
lats = fname.variables['lat'][:]
pressure = fname.variables['pressure'][:]
alt = fname.variables['alt'][:]
year = fname.variables['year'][:]
day = fname.variables['day'][:]
orbit = fname.variables['orbit'][:]
event = fname.variables['event'][:]
fname.close

o3_bin = np.zeros([16,36])
alt_bin = np.zeros([16,36])

def make_retrieval_arrays(tracer, set_year, set_day):
    tracer_bin = np.zeros([16,36])
    for j in range(0,16):
        for k in range(0,36):            
            k_min = (k*5) - 90
            k_max = (k*5) - 85      
            time_lat_bin = np.append(np.where((lats > k_min) & (lats < k_max) & (year == set_year) & (day == set_day)), True)        
            tracer_scans = tracer[time_lat_bin, j]
            tracer_scans_mean = np.mean(tracer_scans)
            tracer_bin[j,k] = tracer_scans_mean
    return tracer_bin      

def pressure_plot_2d(tracer_bin, name, units):
    if name == 'ozone':
        diffs = [1.e+7, 1.75e+7, 2.5e+7, 3.25e+7, 4.e+7, 4.75e+7, 5.5e+7, 6.25e+7, 7.e+7, 7.75e+7, 8.5e+7, 9.25e+7, 1.e+8, 1.75e+8, 2.5e+8, 3.25e+8, 4.e+8, 4.75e+8, 5.5e+8, 6.25e+8, 7.e+8, 7.75e+8, 8.5e+8, 9.25e+8, 1.e+9]
    x, y = np.meshgrid(lats_bins, pressure)
    plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.xticks(np.arange(-90, 120, 30), fontsize=12)
    plt.xlabel('Latitude [%s]' %deg, fontsize=12)
    plt.yscale('log')
    plt.yticks((1.e-4, 4.e-4, 7.e-4, 1.e-3, 4.e-3, 7.e-3, 1.e-2), fontsize=12)
    plt.gca().invert_yaxis()
    plt.ylabel('Pressure [hPa]', fontsize=12)
    plt.title('%s, DOY=%s' %(set_year, set_day), fontsize=14)
    return

def altitude_plot_2d(tracer_bin, name, units):
    if name == 'ozone':
        diffs = [1.e+7, 1.75e+7, 2.5e+7, 3.25e+7, 4.e+7, 4.75e+7, 5.5e+7, 6.25e+7, 7.e+7, 7.75e+7, 8.5e+7, 9.25e+7, 1.e+8, 1.75e+8, 2.5e+8, 3.25e+8, 4.e+8, 4.75e+8, 5.5e+8, 6.25e+8, 7.e+8, 7.75e+8, 8.5e+8, 9.25e+8, 1.e+9]
    x, y = np.meshgrid(lats_bins, pressure)
    y = alt_bin
    plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.xticks(np.arange(-90, 120, 30), fontsize=12)
    plt.xlabel('Latitude [%s]' %deg, fontsize=12)
    plt.yticks([80, 85, 90, 95, 100], fontsize=12)
    plt.ylabel('Altitude [km]', fontsize=12)
    plt.title('%s, DOY=%s' %(set_year, set_day), fontsize=14)  

set_year = 2014
set_day = 6

o3_bin = make_retrieval_arrays(o3, set_year, set_day)
alt_bin = make_retrieval_arrays(alt, set_year, set_day)
lats_bins = np.arange(-87.5, 92.5, 5)

#pressure_plot_2d(o3_bin, 'ozone', '$\mathregular{cm^{-3}}$')
altitude_plot_2d(o3_bin, 'ozone', '$\mathregular{cm^{-3}}$')

plt.show()