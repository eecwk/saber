import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import math
from matplotlib import colors

deg = unichr(176)
set_year = 2014
set_month = 'April'

fname = netCDF4.Dataset('/nfs/a265/earfw/CHRIS/SABER/SABER_Temp_O3_%s%s_v2.0.nc' %(set_month, set_year), 'r', format='NETCDF4')
o3_96 = fname.variables['O3_96'][:]
#o3_127 = fname.variables['O3_127'][:]
#T = fname.variables['ktemp'][:]
lats = fname.variables['tplatitude'][:]
#pressure = fname.variables['pressure'][:]
alt = fname.variables['tpaltitude'][:]
date = fname.variables['date'][:]
time = fname.variables['time'][:]
fname.close
levels = 374

o3_96_bin = np.zeros([levels,36])
o3_127_bin = np.zeros([levels,36])
T_bin = np.zeros([levels,36])
lats_bin = np.zeros([levels,36])
alt_bin = np.zeros([levels,36])
time_bin = np.zeros([levels,36])

def make_retrieval_arrays(tracer):
    tracer_bin = np.zeros([levels,36])
    for j in range(0,levels):
        for k in range(0,36):            
            k_min = (k*5) - 90
            k_max = (k*5) - 85      
            time_lat_bin = np.append(np.where((lats[:,j] > k_min) & (lats[:,j] < k_max)), True)    
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
    tracer_weighted = np.zeros(levels)    
    for j in range (0,levels):    
        sig_cos_x = 0
        sig_cos = 0
        for k in range (lowlat, highlat):
            sig_cos_x = sig_cos_x + (math.cos(math.radians(lats[k,j])) * tracer_bin[j][k])
            sig_cos = sig_cos + math.cos(math.radians(lats[k,j]))         
            if  k == (highlat - 1):
                tracer_weighted[j] = sig_cos_x / sig_cos
    return tracer_weighted

def altitude_plot_2d(tracer_bin, name, units):
    x = lats_bin
    y = alt_bin
    plt.contourf(x[:,:], y[:,:], tracer_bin, norm=colors.LogNorm(), cmap=plt.get_cmap('viridis'))
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('%s [%s]' %(name, units), fontsize=12)
    plt.xticks(np.arange(-90, 120, 30), fontsize=12)
    plt.xlabel('Latitude [%s]' %deg, fontsize=12)
    plt.ylim(60,160)
    plt.ylabel('Altitude [km]', fontsize=12)
    plt.title('%s %s' %(set_month, set_year), fontsize=14)  

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
    plt.ylim(60,160)
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
    if name == 'ozone_96':
        plt.xlim(0,10)
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
    plt.suptitle('%s %s' %(set_month, set_year), fontsize=14)
    for i in range(0,6):    
        lowlat = i * 6
        highlat = (i * 6) + 6
        tracer_weighted = calc_cos_factor(tracer_bin, lowlat, highlat)
        alt_weighted = calc_cos_factor(alt_bin, lowlat, highlat)
        altitude_plot_1d(name, tracer_weighted, alt_weighted, lowlat, highlat, units, i)

o3_96_bin = make_retrieval_arrays(o3_96)
#T_bin = make_retrieval_arrays(T)
alt_bin = make_retrieval_arrays(alt)
lats_bin = make_retrieval_arrays(lats)
time_bin = make_retrieval_arrays(time)

# 1D plot
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11,8))
gs1 = gridspec.GridSpec(2, 3)
gs1.update(wspace=0.1, hspace=0.1)
setup_altitude_plot_1d(o3_96_bin*1.e+6, 'ozone_96', 'ppmv')

# 2D plot
#altitude_plot_2d(o3_96_bin*1.e+6, 'ozone_96', 'ppmv')

plt.show()