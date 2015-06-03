#!/usr/bin/env python
# encoding: utf-8
from __future__ import division
import os
import datetime
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic
from mpl_toolkits.axes_grid1 import ImageGrid
from pydap.client import open_url


# CONSTANT & VARIABLES THAT CAN BE SET
MS2KTS = 1.94384449
MISSING = 9.999E20
CLON = 255.
WINDBARBS = False
PVU = "1.5"
dt = datetime.utcnow()
dt = dt.replace(hour=0)
BDATE = datetime.strftime(dt - timedelta(seconds=24*60*60*12), '%Y%m%d%H')
EDATE = datetime.strftime(dt - timedelta(seconds=60*60*12), '%Y%m%d%H')
IMGPATH = r'./imgs'
IMGPATH2 = r'./imgs/winds'
NPZPATH = r'./npz'
LONLATS = r'gfsHD.lonlats.npz'


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


def draw_map_background(m):
    ''' Setup the map background '''
    m.drawcoastlines(linewidth=2, color='#444444', zorder=6)
    m.drawcountries(linewidth=1, color='#444444', zorder=5)
    m.drawstates(linewidth=0.66, color='#444444', zorder=4)


def compute_theta(temp, pres):
    ''' Compute potential temperature '''
    return temp * (100000. / pres) ** 0.286


def parse_date(dt):
    ''' Parse a string and return a datetime object '''
    return datetime.datetime.strptime(dt, '%Y%m%d%H')


def get_url(dt):
    ''' Create the URL '''
    date = dt.strftime('%Y%m%d')
    hour = dt.strftime('%H')
    url = r'http://nomads.ncep.noaa.gov:9090'
    url += '/dods/gfs_hd/gfs_hd%s/gfs_hd_%sz' % (date, hour)
    return url, date, hour


def add_validation_time(ax, dt, lon, lat, pvu):
    pos = ax.get_position().bounds
    l, b, w, h = pos
    ax2 = plt.axes([l+0.025, b+0.025, 0.2, 0.1], xticks=[], yticks=[],
        frameon=False)
    date = dt.strftime('%d %b %Y @ %H%M UTC')
    ax2.text(0.05, 0.5,
        'DT: %s PVU\n' \
        'Latitude: %s$^{\circ}$N\n' \
        'Longitude: %s$^{\circ}$E\n' \
        r'GFS $0.5^{\circ} \! \times \, 0.5^{\circ}$ Resolution' \
        '\nAnalysis: %s' \
        % (pvu, lat, lon, date), ha='left', va='center', fontsize=9,
        bbox=dict(boxstyle='round, pad=1, rounding_size=0.25', fc="white",
        ec="k", lw=1))
    plt.axes(ax)




if __name__ == '__main__':

    # Get the colormap and set up the figure and map instrance
    cmap = get_cmap()
    fig = plt.figure(figsize=(9, 9))
    #m = Basemap(projection='npstere', boundinglat=12.5, lon_0=CLON,
    #    resolution='l', area_thresh=5000., lat_ts=0)
    m = Basemap(projection='npstere', boundinglat=12.5, lon_0=CLON,
        resolution='l', area_thresh=5000., lat_ts=0, round=True)
    print "Making map."
    # make sure directories exist; if not, create them
    if not os.path.isdir(IMGPATH):
        os.makedirs(IMGPATH)
    if not os.path.isdir(IMGPATH2):
        os.makedirs(IMGPATH2)
    if not os.path.isdir(NPZPATH):
        os.makedirs(NPZPATH)

    # Read in lons and lats from lcoal file if
    try:
        lonlats = np.load(LONLATS)
        lons = lonlats['lons']
        lats = lonlats['lats']
    except:
        lons = np.arange(0, 360, 0.5)
        lats = np.arange(-90, 90.001, 0.5)
        lons, lats = np.meshgrid(lons, lats)
        np.savez_compressed(LONLATS, lons=lons, lats=lats)
    
    print "Loaded lats, lons."
    tmplons = np.arange(0, 360, 0.5)
    x, y = m(lons, lats)

    # Set the initial and end dates
    bdt = parse_date(BDATE)
    edt = parse_date(EDATE)
    delta = int((edt - bdt + datetime.timedelta(hours=6)).total_seconds())
    print "Begin DT:", bdt
    print "End DT:", edt
    print delta
    print edt - bdt
#    delta = datetime.timedelta(hours=6)
#    print np.arange(bdt, edt, delta)
    import os
    print os.path.dirname(os.path.abspath(__file__))
    initial = True
    for deltat in range(-21600, delta, 21600):
        if initial:
            deltat += 21600
        dt = bdt + datetime.timedelta(seconds=deltat)
        
        url, date, hour = get_url(dt)
        print "plotting for:",dt
        npzfile = os.path.join(NPZPATH, 'dt.%s%s.npz' % (date, hour))

        # Use local file if it has been downloaded
        if os.path.isfile(npzfile):
            dobj = np.load(npzfile)
            p = dobj['p']
            t = dobj['t']
            u = dobj['u']
            v = dobj['v']
           
            # Fix the grids with the add cyclic command
            if len(p) == 3:
                p, lons2 = addcyclic(p[0], tmplons)
            else:
                p, lons2 = addcyclic(p, tmplons)

            if len(t.shape) == 3:
                t, lons2 = addcyclic(t[0], tmplons)
            else:
                t, lons2 = addcyclic(t, tmplons)

            if len(u.shape) == 3:
                u, lons2 = addcyclic(u[0], tmplons)
            else:
                u, lons2 = addcyclic(u, tmplons)

            if len(v.shape) == 3:
                v, lons2 = addcyclic(v[0], tmplons)
            else:
                v, lons2 = addcyclic(v, tmplons)

            thta = compute_theta(t, p)
         
        # If local file doesn't exist, try to download it from NOMADS server
        else:
            dobj = open_url(url)
            pres1p5pv = dobj.pres1p5pv[0,:,:].array[:]
            tmp1p5pv = dobj.tmp1p5pv[0,:,:].array[:]
            ugrd1p5pv = dobj.ugrd1p5pv[0,:,:].array[:] * MS2KTS
            vgrd1p5pv = dobj.vgrd1p5pv[0,:,:].array[:] * MS2KTS

            print pres1p5pv.shape, tmp1p5pv.shape, ugrd1p5pv.shape, vgrd1p5pv.shape
            if len(pres1p5pv.shape) == 3:
                p, lons2 = addcyclic(pres1p5pv[0], tmplons)
            else:
                p, lons2 = addcyclic(pres1p5p, tmplons)

            if len(tmp1p5pv.shape) == 3:
                t, lons2 = addcyclic(tmp1p5pv[0], tmplons)
            else:
                t, lons2 = addcyclic(tmp1p5pv, tmplons)

            if len(ugrd1p5pv.shape) == 3:
                u, lons2 = addcyclic(ugrd1p5pv[0], tmplons)
            else:
                u, lons2 = addcyclic(ugrd1p5pv, tmplons)

            if len(vgrd1p5pv.shape) == 3:
                v, lons2 = addcyclic(vgrd1p5pv[0], tmplons)
            else:
                v, lons2 = addcyclic(vgrd1p5pv, tmplons)

            thta = compute_theta(t, p)
            np.savez_compressed(npzfile, u=u, v=v, t=t, p=p, thta=thta)

        # Rotate the winds to make sure they jive with the plotting grid
        print v.shape, u.shape, lons.shape, lats.shape
        #u, v = m.rotate_vector(u[0], v[0], lons, lats)

        # Mask missing data so they aren't plotted
        p = ma.masked_where(p==MISSING, p)
        t = ma.masked_where(t==MISSING, t)
        #u = ma.masked_where(np.abs(u) > 1E10, u)
        #v = ma.masked_where(np.abs(v) > 1E10, v)
        thta = ma.masked_where(thta > 1E10, thta)

        # Clear the figure so we don't duplicate images
        plt.clf()

        # Setup the axes grids
        #grid = ImageGrid(fig, 111, nrows_ncols = (1, 1), direction="row",
        #       axes_pad = 0.5, add_all=True, label_mode = "1",
        #       share_all = True, cbar_location="right", cbar_mode="each",
        #       cbar_size="3%", cbar_pad="1%")

        #ax1 = grid[0]
        draw_map_background(m)


        # Setup the color fill bounds and the skip frequency of wind barbs
        maxa = 100; offset = 330
        iskip = 15; jskip = 20

        # Color filled contours
        cflevs = np.arange(-maxa, maxa+1e-10, 4) + offset

        # What contours will be labled on the colorbar
        cnlevs = np.arange(-maxa, maxa+1e-10, 8) + offset
        # Create the color filled contour
        cf = m.contourf(x, y, thta[0], cmap=cmap, levels=cflevs, extend='both', zorder=1)

        # Create the colorbar
        cbar = plt.colorbar(cf, extend='both')
        cbar.set_ticks(cnlevs)
        clabs = ['%i K' % f for f in cnlevs]
        cbar.ax.set_yticklabels(clabs, size=10)

        # Setup the title
        date2 = dt.strftime('%d %b %Y @ %H%M UTC')
        plt.title('Dynamic Tropopause (%s PVU)\n' \
            'Potential Temperature (K)\n' \
            '%s' % (PVU, date2), fontsize=18)

        # Hack to fix colorbar not plotting on first time through
        if initial:
            initial = False
            continue

        # Save the image
        imgfile = os.path.join(IMGPATH, 'dt_%03i_%s_%sUTC.png' % (CLON,
            date, hour))
        #plt.set_frame_on(False)
        plt.savefig(imgfile, dpi=100, bbox_inches='tight')

        imgfile2 = os.path.join(IMGPATH2, 'dt_%03i_%s_%sUTC.png' % (CLON,
            date, hour))

        WINDBARBS=False
        # Create the image with wind barbs (if chosen)
        if WINDBARBS:
            ax1.set_title('Dynamic Tropopause (%s PVU)\n' \
                'Potential Temperature (K) & Wind Barbs (kts)\n' \
                '%s' % (PVU, date2), fontsize=18)

            m.barbs(x[::iskip, ::jskip], y[::iskip, ::jskip],
                u[::iskip, ::jskip], v[::iskip, ::jskip], pivot='middle',
                length=5.5, zorder=10, ax=ax1)
            plt.savefig(imgfile2, dpi=100, bbox_inches='tight')
