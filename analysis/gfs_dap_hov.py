import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.mlab as mlab
import os
from pydap.client import open_url
from datetime import datetime, timedelta

prefix_path = os.path.dirname(os.path.abspath(__file__)) + '/'
NPZPATH = prefix_path + 'npz/'

def compute_theta(temp, pres):
    ''' Compute potential temperature '''
    return temp * (100000. / pres) ** 0.286

def get_cmap():
    ''' Setup a custom colortable. '''
    cdict = {'red':     ((0.00, 240/255., 220/255.),
                         (0.25,  40/255.,  20/255.),
                         (0.50, 225/255., 255/255.),
                         (0.75, 150/255., 150/255.),
                         (1.00, 255/255., 255/255.)),

             'green':   ((0.00, 240/255., 220/255.),
                         (0.25,   0/255.,  50/255.),
                         (0.50, 255/255., 255/255.),
                         (0.75,   0/255.,  35/255.),
                         (1.00, 225/255., 240/255.)),

             'blue':    ((0.00, 255/255., 255/255.),
                         (0.25, 160/255., 150/255.),
                         (0.50, 255/255., 170/255.),
                         (0.75,   0/255.,  35/255.),
                         (1.00, 225/255., 240/255.))}
    return mpl.colors.LinearSegmentedColormap('cool2warm', cdict, 256)

def uniquify(list):
    unique = []
    seen = {}
    for item in list:
        if item in seen: continue
        seen[item] = 1
        unique.append(item)
    return unique


def parse_date(dt):
    ''' Parse a string and return a datetime object '''
    return datetime.strptime(dt, '%Y%m%d%H')

def get_dataset(url):
    
    try:
        data = Dataset(url)
        lats = data.variables['lat'][:]
        lons = data.variables['lon'][:]
        var = data.variables[variable_to_plot][:]
    
        u_lat_idx = np.where(lats == lat_upper)
        l_lat_idx = np.where(lats == lat_lower)
        
        print u_lat_idx[0][0]
        print l_lat_idx[0][0]
        
        bounded_vars = []
        
        for lat_idx in np.arange(u_lat_idx[0][0], l_lat_idx[0][0]+1, 1):
            bounded_vars.append(var[0][lat_idx])
    except:    
        return np.nan
    
    return bounded_vars

def get_data_mean_lat(array, lats):
        
    u_lat_idx = np.where(lats == lat_upper)
    l_lat_idx = np.where(lats == lat_lower)
    print array.shape
    print l_lat_idx[0].shape, u_lat_idx[0].shape
    bounded_vars = array[0][l_lat_idx[0][0]:u_lat_idx[0][0],:]
    print "BOUNDED_VARS:",bounded_vars
    print "BOUNDED_VARS:",bounded_vars.shape
    mean_lat = np.ma.mean(bounded_vars, axis=0)
    print mean_lat.shape
    print mean_lat.max()
    print mean_lat
    return np.ma.mean(bounded_vars, axis=0)

#def get_url(dt):
###    ''' Create the URL '''
#    date = dt.strftime('%Y%m%d')
#    hour = dt.strftime('%H')
#    url = "http://motherlode.ucar.edu:8080/thredds/dodsC/grib/NCEP/GFS/Global_onedeg/files/GFS_Global_onedeg_" + date + "_" + hour + "00.grib2"
    #url = 'http://motherlode.ucar.edu:8080/thredds/dodsC/fmrc/NCEP/GFS/Global_onedeg/runs/NCEP-GFS-Global_onedeg_RUN_' + date[0:4] + '-' + date[4:6] + '-' + date[6:8] + 'T' + hour + ':00:00Z'
#    return url, date, hour

def get_url(dt):
    ''' Create the URL '''
    date = dt.strftime('%Y%m%d')
    hour = dt.strftime('%H')
    url = r'http://nomads.ncep.noaa.gov:9090'
    url += '/dods/gfs_0p50/gfs%s/gfs_0p50_%sz' % (date, hour)
    return url, date, hour

variable_to_plot = 'v-component_of_wind_tropopause'


#30 N to 60 N
MS2KTS = 1.94384449
lat_upper = 60
lat_lower = 30
LONLATS = prefix_path + "gfsHD.lonlats.npz"
MISSING = 9.99900026e+20
CLON = 255.

dt = datetime.utcnow()
dt = dt.replace(hour=0)
BDATE = datetime.strftime(dt - timedelta(seconds=24*60*60*12), '%Y%m%d%H')
EDATE = datetime.strftime(dt - timedelta(seconds=60*60*12), '%Y%m%d%H')
print BDATE, EDATE
#BDATE = '2012060100'

#This part of the code loops through all of the dates
bdt = parse_date(BDATE)
edt = parse_date(EDATE)

delta = int((edt - bdt + timedelta(hours=6)).total_seconds())

plottable_data = []
plottable_date = []
plottable_hour = []

initial = True

# Read in lons and lats from lcoal file if
try:
    lonlats = np.load(LONLATS)
    lons = lonlats['lons']
    lats = lonlats['lats']
except:
    lons = np.arange(0, 360.5, .5)
    lats = np.arange(-90, 90.5, .5)
    lons, lats = np.meshgrid(lons, lats)
    np.savez_compressed(LONLATS, lons=lons, lats=lats)

tmplons = np.arange(0, 360, 0.5)

for deltat in range(-21600, delta, 21600):
    if initial:
        deltat += 21600
    dt = bdt + timedelta(seconds=deltat)
    url, date, hour = get_url(dt)
    print date, hour
    npzfile = os.path.join(NPZPATH, 'dt.%s%s.npz' % (date, hour))
    
    if os.path.isfile(npzfile):
        dobj = np.load(npzfile)
        p = dobj['p']
        t = dobj['t']
        u = dobj['u']
        v = dobj['v']
        thta = dobj['thta']
        # If local file doesn't exist, try to download it from NOMADS server
    else:
        print url
        dobj = open_url(url)
        pres1p5pv = dobj.pres2pv[0,:,:].array[:]
        tmp1p5pv = dobj.tmp2pv[0,:,:].array[:]
        ugrd1p5pv = dobj.ugrd2pv[0,:,:].array[:]
        vgrd1p5pv = dobj.vgrd2pv[0,:,:].array[:]
        
        u = ugrd1p5pv
        v = vgrd1p5pv
        t = tmp1p5pv
        p = pres1p5pv
        #print t, p
        thta = compute_theta(t, p)
        #print thta
        np.savez_compressed(npzfile, u=u, v=v, t=t, p=p, thta=thta)
    print p
    p = np.ma.masked_where(p > 9e10, p)
    t = np.ma.masked_where(t > 9e10, t)
    print p,t
    u = np.ma.masked_where(np.abs(u) > 1E10, u)
    v = np.ma.masked_where(np.abs(v) > 1E10, v)
    thta = np.ma.masked_where(thta > 1E10, thta)
    plottable_data.append(get_data_mean_lat(v, lats))
    print v.shape
    print plottable_data[0].shape
    plottable_hour.append(hour)
    plottable_date.append(datetime.strftime(dt, '%Y-%m-%d'))

fig = plt.figure(figsize=(12,9))
#ax = fig.subplot(111)
print plottable_data
print np.asarray(plottable_data).shape
cnt = plt.contourf(np.array(plottable_data),levels=np.arange(-60,62,2),cmap=get_cmap(),extend='both')

indexes = []
unique = uniquify(plottable_date)
for day in unique:
    indexes.append(plottable_date.index(day))

lon_locations = [0, 30, 60, 120, 180, 210, 240, 270, 300, 330, 360]
lon_pts = []
lon_labels = ['0$^{\circ}$', '30$^{\circ}$E', '60$^{\circ}$E', '90$^{\circ}$E', '120$^{\circ}$E', '150$^{\circ}$E', '180$^{\circ}$' ,'150$^{\circ}$W', '120$^{\circ}$W', '90$^{\circ}$W', '60$^{\circ}$W', '30$^{\circ}$W']
print "BLAH"
for l in lon_locations:
    print lons
    idx = np.where(lons == float(l))[0]
    print idx

    lon_pts.append(lons[0,idx])
	#$rint idx
lon_pts = np.asarray(lon_pts)
#print lon_pts.shape
print "STOP"
print lon_pts, lon_labels
lon_pts = np.arange(0, 721, 60)
plt.xticks(lon_pts, lon_labels)
plt.yticks(indexes, unique)
plt.gca().invert_yaxis()
plt.grid(True)
plt.axvline(x=230*2, color='k', linestyle='-', linewidth=3)
plt.axvline(x=300*2, color='k', linestyle='-', linewidth=3)
plt.xlabel("Longitude")
plt.colorbar(extend='both')
plt.title("GFS 1 Degree Tropopause V-wind (m/s) \nAveraged between " + str(lat_upper) + " and " + str(lat_lower) + "$^{\circ}$N",fontsize=18)
plt.tight_layout()
OUTDIR = '/data/soundings/http/blumberg/'

plt.savefig(OUTDIR + 'dt_hov.png', bbox_inches='tight')


