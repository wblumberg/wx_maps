from netCDF4 import Dataset
import os
import sys
import numpy as np

file_path = sys.argv[1]
d = Dataset(file_path)
print d.variables.keys()

dim = {}
dim['lat'] = np.array(d.variables['lat'][:])
dim['lon'] = np.array(d.variables['lon'][:])
dim['elevation'] = np.array(d.variables['Geopotential_height_surface'][:])
print dim
d.close()
np.savez('40kmRAPgrid.npz', **dim)


