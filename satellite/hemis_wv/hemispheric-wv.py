import matplotlib as mpl
mpl.use('Agg')
import Nio
import matplotlib.pyplot as plt
import sys
import numpy as np
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as colors
#from colormaps import khcm

def from_ascii(filename, name):
    palette = open(filename)
    lines = palette.readlines()
    carray = np.zeros([len(lines), 3])
    for num, line in enumerate(lines):
        carray[num, :] = [float(val) for val in line.strip().split()]
    carray = carray
    cmap = mpl.colors.ListedColormap(carray/255, name=name)
    mpl.cm.register_cmap(name=name, cmap=cmap)
 
def grayify_cmap(cmap):
    """Return a grayscale version of the colormap"""
    cmap = plt.cm.get_cmap(cmap)
    colors = cmap(np.arange(cmap.N))

    # convert RGBA to perceived greyscale luminance
    # cf. http://alienryderflex.com/hsp.html
    RGB_weight = [0.299, 0.587, 0.114]
    luminance = np.sqrt(np.dot(colors[:, :3] ** 2, RGB_weight))
    colors[:, :3] = luminance[:, np.newaxis]

    return cmap.from_list(cmap.name + "_grayscale", colors, cmap.N)


from_ascii("kh_cm.txt", 'kh_cm')
hemisxy = np.load('nhemispheric_mosaic.npz')
m = Basemap(projection='npstere', lat_0=60.,lat_ts=30.,  boundinglat=0, lon_0=-105., area_thresh=1000., k_0=.933)
cmap = grayify_cmap(plt.get_cmap('afmhot'))
cmap.set_bad('black', 1)
 
# Get the current time
now = datetime.utcnow()
i = 0
no_imgs = 48
sat_type = sys.argv[1]
now = now - timedelta(hours = now.hour % 3.)
last_sat_img_string = datetime.strftime(now, '%Y%m%d_%H00')
 
while i < no_imgs:
    data_link = 'http://motherlode.ucar.edu/thredds/dodsC/satellite/' + sat_type + '/NHEM-MULTICOMP_1km/current/NHEM-MULTICOMP_1km_' + sat_type + '_' + last_sat_img_string + '.gini'
    
    try:
        #data = Nio.open_file(data_link)
        data = open_url(data_link)
    except Exception,e:
        print e
        now = now - timedelta(hours = 3.)
        last_sat_img_string = datetime.strftime(now, '%Y%m%d_%H00')       
        continue
    print data
    print "Drawing for: ", last_sat_img_string
    print data["IR_WV"][:]
    data_length = a.variables['IR_WV'].shape[1]
    wv = np.vstack((a.variables['IR_WV'][0,0:data_length -1,:], a.variables['IR_WV'][0,data_length,:]))
    my_dpi = 96
    
    plt.figure(figsize=(1000/my_dpi, 900/my_dpi), dpi=my_dpi)
    plt.gca().set_axis_bgcolor('black')

    lons = hemisxy['lon']
    lats = hemisxy['lat']
    x,y = m(lons, lats)

    wv = np.ma.masked_where(wv>=-15, wv)

    plt.pcolormesh(x,y, wv,cmap=plt.get_cmap('kh_cm'), vmin=-130, vmax=-10)
    m.drawcoastlines(color='#3399FF')
    m.drawstates(color='#3399FF')
    m.drawcountries(color='#3399FF')
    
    plt.text(0.5, .98, 'Hemispheric Water Vapor (Imager 6.7/6.5 micron IR WV) Mosaic for ' + datetime.strftime(now, '%Y/%m/%d %H00 UTC'), transform=plt.gca().transAxes,horizontalalignment='center', bbox=dict(facecolor='white', alpha=0.7))
    plt.text(0,0, "Made by Greg Blumberg (wblumberg@ou.edu)", fontsize=10, transform=plt.gca().transAxes, color='yellow')
    plt.tight_layout()
    plt.savefig('wvhemis_' + str(i) + '.png', bbox_inches='tight', pad_inches=0)
    #plt.show()
    plt.close()
    i += 1

    now = now - timedelta(hours = 3.)
    last_sat_img_string = datetime.strftime(now, '%Y%m%d_%H00')
