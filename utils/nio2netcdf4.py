#!/usr/bin/env python
import Nio
from netCDF4 import Dataset
import sys
import numpy as NP

# convert a file read by Nio file to netCDF (NETCDF4_CLASSIC format
# with zlib compression) using PyNIO and netCDF4.

if len(sys.argv) < 2:
    print """
        usage:  python grib2nc4.py <nio file name>
        a netCDF file (in NETCDF4_CLASSIC format, with zlib compression)
        called <nio file name>.nc will be created
        
        requires PyNIO to read nio files.
        """
    sys.exit(1)

niofile = sys.argv[1]
f = Nio.open_file(niofile)
new_filename = niofile.split('/')[-1]
nc = Dataset(new_filename + '.nc','w')

for name,att in f.__dict__.items():
    setattr(nc,name,att)

print 'dimensions (name, length):'
for dimname,dim in f.dimensions.items():
    nc.createDimension(dimname,dim)
    print dimname, dim

print 'variables (name, type, dimensions, units):'
for varname,var in f.variables.items():
    niodata = NP.array(var[:])
    print varname,var.typecode(),var.dimensions
    varo = nc.createVariable(varname,var.typecode(),var.dimensions)
    for name, att in var.__dict__.items():
        setattr(varo,name,att)
    varo[:] = niodata[:]

f.close()
nc.close()

'http://motherlode.ucar.edu/thredds/dodsC/satellite/WV/NHEM-MULTICOMP_1km/current/NHEM-MULTICOMP_1km_WV_20121024_2100.gini'
