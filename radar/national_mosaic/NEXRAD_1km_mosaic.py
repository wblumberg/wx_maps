import datetime as dt
import matplotlib.pyplot as plt
plt.use('Agg')
import numpy as np
import matplotlib as mpl
import pytz
import glob

from mpl_toolkits.basemap import Basemap
from pydap.client import open_url
from siphon.tds import TDSCatalog

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
path_to_imgs = './'
'''
max_no = 95

path_to_imgs = '/data/soundings/blumberg/nexrad1km/imgs/'
cur_imgs = np.sort(glob.glob(path_to_imgs + '*.png'))
if 'radar_95.png' in cur_imgs:
    # remove the last image
    os.system('/usr/bin/rm ' + path_to_imgs + 'radar_95.png')

# update the names of all the other images
for img in cur_imgs:
    no = int(img.split('_')[1].split('.')[0])
    os.system('/usr/bin/mv ' + path_to_imgs + 'radar_' + str(no + 1) + '.png')
'''

today = dt.datetime.utcnow() 
today = today.replace(minute=today.minute - today.minute % 15)
#end = today - dt.timedelta(days=1)
#delta = dt.timedelta(minutes=15)


# Make all of today's graphics

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
data = np.ma.masked_where(data == -30, data)


time = dataset["time"][:][0] / 1000.
dt_obj = dt.datetime.fromtimestamp(time, pytz.utc)
title = "1 km NEXRAD Mosaic for " + dt.datetime.strftime(dt_obj , '%Y/%m/%d %H:%M UTC')

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
ax.text(0.5, .95, title, transform=plt.gca().transAxes,fontsize=16, horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.9, boxstyle='round'))                                                     
#cax = m.pcolormesh(x, y, data, cmap=cmap, norm=norm)
cax = m.imshow(data, extent = (x.min(), x.max(), y.min(), y.max()), cmap=cmap, norm=norm, origin="upper")
m.drawcoastlines(color='#FFFFFF')
m.drawstates(color='#FFFFFF')
m.drawcountries(color='#FFFFFF')
m.drawcounties(color='#FFFFFF', linewidth=0.09)
cbar = m.colorbar(cax)
plt.tight_layout()
plt.savefig(path_to_imgs + 'radar_0.png')



