import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import math
from matplotlib import colors

deg = unichr(176)
set_year = 2014
set_day = 6

fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/SABER/NIGHTLY/atox_athy_night_YY%s_V5.3_fixedfnight_SV2.nc' %set_year, 'r', format='NETCDF4')
o3 = fname.variables['o3conc'][:]
o = fname.variables['qatox'][:]*(1.e+6)
h = fname.variables['qathy'][:]*(1.e+6)
T = fname.variables['ktemp'][:]
lats = fname.variables['lat'][:]
pressure = fname.variables['pressure'][:]
alt = fname.variables['alt'][:]
year = fname.variables['year'][:]
day = fname.variables['day'][:]
orbit = fname.variables['orbit'][:]
event = fname.variables['event'][:]
fname.close
lats_bin_mid_saber = np.arange(-87.5, 92.5, 5)

o3_bin = np.zeros([16,36])
o_bin = np.zeros([16,36])
h_bin = np.zeros([16,36])
T_bin = np.zeros([16,36])
alt_bin = np.zeros([16,36])

def make_retrieval_arrays(tracer, set_year, set_day):
    tracer_bin = np.zeros([16,36])
    for j in range(0,16):
        for k in range(0,36):            
            k_min = (k*5) - 90
            k_max = (k*5) - 85      
            time_lat_bin = np.append(np.where((lats > k_min) & (lats < k_max) & (year == set_year) & (day == set_day)), True)        
            tracer_scans = tracer[time_lat_bin, j]
            tracer_scans_good = np.array([])
            for i in range(0, len(tracer_scans)):
                if tracer_scans[i] > 0:
                    tracer_scans_good = np.append(tracer_scans_good, tracer_scans[i])
            if len(tracer_scans_good) < 2:
                tracer_scans_good = 0
            tracer_scans_mean = np.mean(tracer_scans_good)
            tracer_bin[j,k] = tracer_scans_mean                               
    return tracer_bin      

def calc_cos_factor(tracer_bin, lowlat, highlat):
    tracer_weighted = np.zeros(16)    
    for j in range (0,16):    
        sig_cos_x = 0
        sig_cos = 0
        for k in range (lowlat, highlat):
            sig_cos_x = sig_cos_x + (math.cos(math.radians(lats_bin_mid_saber[k])) * tracer_bin[j][k])
            sig_cos = sig_cos + math.cos(math.radians(lats_bin_mid_saber[k]))         
            if  k == (highlat - 1):
                tracer_weighted[j] = sig_cos_x / sig_cos
    return tracer_weighted

def pressure_plot_2d(tracer_bin, name, units):
    x, y = np.meshgrid(lats_bin, pressure)
    if name == 'ozone':
        diffs = [1.e+7, 1.75e+7, 2.5e+7, 3.25e+7, 4.e+7, 4.75e+7, 5.5e+7, 6.25e+7, 7.e+7, 7.75e+7, 8.5e+7, 9.25e+7, 1.e+8, 1.75e+8, 2.5e+8, 3.25e+8, 4.e+8, 4.75e+8, 5.5e+8, 6.25e+8, 7.e+8, 7.75e+8, 8.5e+8, 9.25e+8, 1.e+9]
        plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    else:
        plt.contourf(x[:,:], y[:,:], tracer_bin, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
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
    x, y = np.meshgrid(lats_bin, pressure)
    y = alt_bin
    if name == 'ozone':
        diffs = [1.e+7, 1.75e+7, 2.5e+7, 3.25e+7, 4.e+7, 4.75e+7, 5.5e+7, 6.25e+7, 7.e+7, 7.75e+7, 8.5e+7, 9.25e+7, 1.e+8, 1.75e+8, 2.5e+8, 3.25e+8, 4.e+8, 4.75e+8, 5.5e+8, 6.25e+8, 7.e+8, 7.75e+8, 8.5e+8, 9.25e+8, 1.e+9]
        plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    if name == 'atomic_oxygen':
        diffs = [1.e+0, 2.5e+0, 4.e+0, 5.5e+0, 7.e+0, 8.5e+0, 1.e+1, 2.5e+1, 4.e+1, 5.5e+1, 7.e+1, 8.5e+1, 1.e+2, 2.5e+2, 4.e+2, 5.5e+2, 7.e+2, 8.5e+2, 1.e+3, 2.5e+3, 4.e+3, 5.5e+3, 7.e+3, 8.5e+3, 1.e+4, 2.5e+4, 4.e+4, 5.5e+4, 7.e+4, 8.5e+4, 1.e+5]
        plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    if name == 'atomic_hydrogen':
        diffs = [1.e-2, 2.5e-2, 4.e-2, 5.5e-2, 7.e-2, 8.5e-2, 1.e-1, 2.5e-1, 4.e-1, 5.5e-1, 7.e-1, 8.5e-1, 1.e+0, 2.5e+0, 4.e+0, 5.5e+0, 7.e+0, 8.5e+0, 1.e+1]
        plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    if name == 'temperature':
        diffs = np.arange(150,225,5)
        plt.contourf(x[:,:], y[:,:], tracer_bin, diffs, cmap=plt.get_cmap('viridis'))
    else:
        plt.contourf(x[:,:], y[:,:], tracer_bin, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.xticks(np.arange(-90, 120, 30), fontsize=12)
    plt.xlabel('Latitude [%s]' %deg, fontsize=12)
    plt.ylim(77,100)
    plt.yticks([80, 85, 90, 95, 100], fontsize=12)
    plt.ylabel('Altitude [km]', fontsize=12)
    plt.title('%s, DOY=%s' %(set_year, set_day), fontsize=14)  

def altitude_plot_1d(name, tracer_weighted, alt_weighted, lowlat, highlat, units, plot_no):
    lowlat_no = int((lowlat * 5) - 90)
    highlat_no = int((highlat * 5) - 90)
    if plot_no > 5:
        plot_no = plot_no - 6
    plt.subplot(gs1[plot_no])
    plt.title('%s%s to %s%s' %(lowlat_no, deg, highlat_no, deg), fontsize=14)
    x = tracer_weighted[::-1]
    y = alt_weighted[::-1]
    plt.plot(x, y, color='k', label='saber')
    plt.ylim(77,100)
    if plot_no == 0:
        plt.ylabel('Altitude [km]', fontsize=12)
        plt.tick_params(labelbottom='off')
    if plot_no == 1:
        plt.tick_params(labelleft='off')
        plt.tick_params(labelbottom='off')
    if plot_no == 2:
        plt.tick_params(labelleft='off')
        plt.tick_params(labelbottom='off')
    if plot_no == 3:
        plt.xlabel('%s [%s]' %(name, units), fontsize=12)
        plt.ylabel('Altitude [km]', fontsize=12)
    if plot_no == 4:
        plt.xlabel('%s [%s]' %(name, units), fontsize=12)
        plt.tick_params(labelleft='off')
    if plot_no == 5:
        plt.xlabel('%s [%s]' %(name, units), fontsize=12)
        plt.tick_params(labelleft='off')
    if name == 'ozone':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlim(0,5.e+8)
    if name == 'atomic_oxygen':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlim(0,50000)
    if name == 'atomic_hydrogen':
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.xlim(0,5)
    if name == 'temperature':
        plt.xlim(150,220)
    if plot_no == 2:
        plt.legend(loc=1)
    return

def setup_altitude_plot_1d(tracer_bin, name, units):
    plt.suptitle('%s, DOY=%s, night' %(set_year, set_day), fontsize=16)
    for i in range(0,6):    
        lowlat = i * 6
        highlat = (i * 6) + 6
        tracer_weighted = calc_cos_factor(tracer_bin, lowlat, highlat)
        alt_weighted = calc_cos_factor(alt_bin, lowlat, highlat)
        altitude_plot_1d(name, tracer_weighted, alt_weighted, lowlat, highlat, units, i)

o3_bin = make_retrieval_arrays(o3, set_year, set_day)
o_bin = make_retrieval_arrays(o, set_year, set_day)
h_bin = make_retrieval_arrays(h, set_year, set_day)
T_bin = make_retrieval_arrays(T, set_year, set_day)
alt_bin = make_retrieval_arrays(alt, set_year, set_day)
lats_bin = np.arange(-87.5, 92.5, 5)

# 1D plot
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.1, hspace=0.1)
#setup_altitude_plot_1d(o3_bin, 'ozone', '$\mathregular{cm^{-3}}$')
setup_altitude_plot_1d(o_bin, 'atomic_oxygen', 'ppmv')
#setup_altitude_plot_1d(h_bin, 'atomic_hydrogen', 'ppmv')
#setup_altitude_plot_1d(T_bin, 'temperature', 'K')

# 2D plot
#pressure_plot_2d(o3_bin, 'ozone', '$\mathregular{cm^{-3}}$')
#altitude_plot_2d(o3_bin, 'ozone', '$\mathregular{cm^{-3}}$')
#altitude_plot_2d(o_bin, 'atomic_oxygen', 'ppmv')
#altitude_plot_2d(h_bin, 'atomic_hydrogen', 'ppmv')
#altitude_plot_2d(T_bin, 'temperature', 'K')

plt.show()