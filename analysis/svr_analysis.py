import datetime as dt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytz
import glob
import os
import scipy.ndimage
from mpl_toolkits.basemap import Basemap
from pydap.client import open_url
from siphon.tds import TDSCatalog
from netCDF4 import Dataset

# Adapted from the example at:
# https://github.com/lesserwhirls/siphon/tree/master/examples/notebooks/radar

def radar_colormap():
    nws_reflectivity_colors = [
    #"#646464", # ND
    "#000000", # ND
    "#ccffff", # -30
    "#cc99cc", # -25
    "#996699", # -20
    "#663366", # -15
    "#cccc99", # -10
    "#999966", # -5
    "#646464", # 0
    "#04e9e7", # 5
    "#019ff4", # 10
    "#0300f4", # 15
    "#02fd02", # 20
    "#01c501", # 25
    "#008e00", # 30
    "#fdf802", # 35
    "#e5bc00", # 40
    "#fd9500", # 45
    "#fd0000", # 50
    "#d40000", # 55
    "#bc0000", # 60
    "#f800fd", # 65
    "#9854c6", # 70
    "#fdfdfd" # 75
    ]

    cmap = mpl.colors.ListedColormap(nws_reflectivity_colors)

    return cmap
#path_to_imgs = './'

max_no = 95

path_to_imgs = '/data/soundings/http/blumberg/nexrad_rap/'
#cur_imgs = np.sort(glob.glob(path_to_imgs + '*.png'))
#if 'radar_' + str(max_no) + '.png' in cur_imgs:
#    # remove the last image
#    os.system('/usr/bin/rm ' + path_to_imgs + 'radar_' + str(max_no) + '.png')

# update the names of all the other images
#for img in cur_imgs:
#    no = int(img.split('_')[1].split('.')[0])
#    old_path = path_to_imgs + 'radar_' + str(no) + '.png'
#    new_path = path_to_imgs + 'radar_' + str(no + 1) + '.png'
#    print old_path
#    print new_path
#    os.system('/usr/bin/mv ' + old_path + ' ' + new_path)

today = dt.datetime.utcnow() 
today = today.replace(minute=today.minute - today.minute % 15)
end = today - dt.timedelta(days=1)
#delta = dt.timedelta(minutes=15)

## GET RAP DATA
today_rap = today.replace(minute=0) - dt.timedelta(hours=1)
# Make all of today's graphics
rap_url = 'http://atm.ucar.edu//thredds/catalog/grib/NCEP/RAP/CONUS_13km/RR_CONUS_13km_{}.grib2/catalog.xml'.format(today_rap.strftime('%Y%m%d_%H%M'))
rap_cat = TDSCatalog(rap_url)
key = rap_cat.datasets.keys()[0]
latestDs = rap_cat.datasets[key]

dataset = open_url(latestDs.accessUrls['OPENDAP'])

stride = 7
def getData(dataset, var_name, index, stride=7):
    data = dataset[var_name][1,index,::stride,::stride]
    grid = Dataset('/data/soundings/blumberg/programs/wx_maps/utils/13km_latlon.nc')
    grid_lat = grid.variables['lat'][::stride,::stride].T
    grid_lon = grid.variables['lon'][::stride,::stride].T
    elevation = grid.variables['Geopotential_height_surface'][::stride,::stride].T
    return data[var_name][:].squeeze(), grid_lat, grid_lon, elevation

CAPE, cape_grid_lat, cape_grid_lon, cape_grid_elev = getData(dataset, 'Convective_available_potential_energy_pressure_difference_layer',0, stride=1)
CAPE = scipy.ndimage.gaussian_filter(CAPE.T, 2)
sfc_u,grid_lat, grid_lon, grid_elev = getData(dataset , 'u-component_of_wind_height_above_ground', 0)
sfc_v = getData(dataset , 'v-component_of_wind_height_above_ground', 0)[0]

# Compute the 0-6 km SHEAR
u_6km = np.empty(sfc_u.shape)
v_6km = np.empty(sfc_u.shape)
all_levels = np.arange(len(dataset['isobaric'][:]))
height_msl = getData(dataset, 'Geopotential_height_isobaric', all_levels)[0]
v_comp = getData(dataset, 'v-component_of_wind_isobaric', all_levels)[0]
u_comp = getData(dataset, 'u-component_of_wind_isobaric', all_levels)[0]
height_agl = height_msl - grid_elev
for idx in np.ndenumerate(u_6km):
    v = v_comp[:,idx[0][0], idx[0][1]]
    u = u_comp[:,idx[0][0], idx[0][1]]
    h = height_agl[:,idx[0][0], idx[0][1]]
    v_6km[idx[0]] = np.interp(6000, h[::-1], v[::-1])
    u_6km[idx[0]] = np.interp(6000, h[::-1], u[::-1])
u_shear = 1.94384 * (u_6km - sfc_u)
v_shear = 1.94384 * (v_6km - sfc_v)
mag = np.sqrt(np.power(u_shear,2) + np.power(v_shear,2))
u_shear = np.ma.masked_where(mag < 30, u_shear)
v_shear = np.ma.masked_where(mag < 30, v_shear)

print "SHEAR:",u_shear.shape, v_shear.shape
print "CAPE:", CAPE.shape
print "GRID:", grid_lat.shape, grid_lon.shape

## GET NEXRAD MOSAIC DATA
url = "http://atm.ucar.edu/thredds/catalog/nexrad/composite/gini/n0r/1km/{}/catalog.xml".format(today.strftime("%Y%m%d"))
cat = TDSCatalog(url)
names = cat.datasets.keys()
names.sort()
name_fmt = 'Level3_Composite_n0r_1km_%Y%m%d_%H%M.gini'
idx = np.where(np.asarray(names) == today.strftime(name_fmt))

latest = names[idx[0]]
latestDs = cat.datasets[latest]

dataset = open_url(latestDs.accessUrls['OPENDAP'])

# get basic info from the file regarding projection attributes
# see the following for more info on projections as described here:
#   http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218
#http://www.wmo.int/pages/prog/www/WDM/Guides/Guide-binary-2.html [see LAMBERT CONFORMAL SECANT OR TANGENT CONE GRIDS]
#   http://www.unidata.ucar.edu/mailing_lists/archives/netcdf-java/2006/msg00006.html [starndard parallels in CDM]
proj_attributes = dataset['LambertConformal'].attributes
rsphere = proj_attributes['earth_radius']

# lat_0	: center of desired map domain (in degrees) [Basemap]
# CDM : 'latitude_of_projection_origin'
# NCO : where Dx and Dy are defined
# :
lat_0 = proj_attributes['latitude_of_projection_origin']

# lon_0	: center of desired map domain (in degrees) [Basemap]
# CDM : 'longitude_of_central_meridian'
# NCO : Lov
#
# Lov - The orientation of the grid; i.e., the east longitude value
# of the meridian which is parallel to the y-axis (or columns of the
#  grid) along which latitude increases as the y-coordinate increases.
# (Note: The orientation longitude may, or may not, appear within a
# particular grid.)
#
lon_0 = proj_attributes['longitude_of_central_meridian'] # Lov

# lat_1, lat_2 : 1st and second parallels [Basemap]
# CDM : 'standard_parallel' - this attr contains both 1st and 2nd
# NCO : ??? Not sure where this shows up on the nco page
lat_1 = proj_attributes['standard_parallel']

# hardcoded from catalog metadata - will add metadata into catalog class
# to grab this programatically
llcrnrlon = -120.02283
llcrnrlat = 23.0132
urcrnrlon = -63.5105
urcrnrlat = 47.30635

xstride = 1
ystride = 1

# download x and y coords and convert them from
# km to m, offset for use in basemap
x = dataset['x'][:] * 1000.
y = dataset['y'][:] * 1000;

width = x.max() - x.min()
height = y.max() - y.min()

data = dataset["Reflectivity"][0,0::ystride,0::xstride]
x = x[0::xstride]
y = y[0::ystride]
data = np.squeeze(data)
data = np.ma.masked_where(data < 10, data)

## PLOT ALL THE DATA

time = dataset["time"][:][0] / 1000.
dt_obj = dt.datetime.fromtimestamp(time, pytz.utc)
title = "1 km NEXRAD Mosaic for " + dt.datetime.strftime(dt_obj , '%Y/%m/%d %H:%M UTC')
title += '\n' + today_rap.strftime('%H UTC 13km RAP Analysis: 0-6 km BWD and MUCAPE')

fig = plt.figure(figsize=(27,30))

ax = fig.add_subplot(1,1,1)

m = Basemap(projection='lcc', lat_0 = lat_0, lon_0 = lon_0, lat_1 = lat_1,
      llcrnrlon = llcrnrlon, llcrnrlat = llcrnrlat,
      urcrnrlat = urcrnrlat, urcrnrlon = urcrnrlon,
      area_thresh = 1000., rsphere = rsphere, resolution='l')
#grid_lon, grid_lat = m.makegrid(len(x),len(y))
#print grid_lon, grid_lat

x = (m.llcrnrx - x.min()) + x
y = (m.llcrnry - y.min()) + y

cmap = radar_colormap()
cmap.set_bad('k')
norm = mpl.colors.Normalize(vmin=-35, vmax=80)
ax.text(0.5, .97, title, transform=plt.gca().transAxes,fontsize=16, horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.9, boxstyle='round'))                                                     
#cax = m.pcolormesh(x, y, data, cmap=cmap, norm=norm)
cax = m.imshow(data, extent = (x.min(), x.max(), y.min(), y.max()), cmap=cmap, norm=norm, origin="upper")
m.drawcoastlines(color='#FFFFFF')
m.drawstates(color='#FFFFFF')
m.drawcountries(color='#FFFFFF')
m.drawcounties(color='#FFFFFF', linewidth=0.09)
cbar = m.colorbar(cax)

x_grid, y_grid = m(cape_grid_lon, cape_grid_lat)
c = m.contour(x_grid, y_grid, CAPE, np.arange(500,6500,500), cmap='spring_r')
plt.clabel(c, fmt='%4.0f')
stride = 10
#m.barbs(x_grid[::stride, ::stride], y_grid[::stride, ::stride], u_shear[::stride,::stride], v_shear[::stride,::stride], mag[::stride,::stride], cmap='RdPu') 
m.barbs(grid_lon, grid_lat, u_shear.T, v_shear.T, mag.T, clim=[30,70], latlon=True, cmap='RdPu') 

plt.tight_layout()
# Make the newest graphic.
out_filename = path_to_imgs + 'nexrad_rap_composite_' + today.strftime("%Y%m%d_%H%M.png")
plt.savefig(out_filename, bbox_inches='tight')

# Delete the oldest graphic.
old_filename = path_to_imgs + 'nexrad_rap_composite_' + end.strftime("%Y%m%d_%H%M.png")
os.system('/usr/bin/rm ' + old_filename)

print "Making looper"
path_to_looper = '/data/soundings/blumberg/programs/wx_maps/looper/make_looper.py'
os.system('/data/soundings/anaconda/bin/python ' + path_to_looper + ' ' + path_to_imgs + ' nexrad_rap')

