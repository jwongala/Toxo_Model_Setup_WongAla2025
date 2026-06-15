#!/usr/bin/env python3
# coding: utf-8

# # Particle Tracking Script using OceanParcels to model transport and connectivity of T. gondii 
# 
# This is the version of the script includes:
# - base version of the script with particles bouncing perpendicular to the coast line
# - nesting of the Kaneohe and Oahu South Shore ROMS
# 
# by Jennifer Wong-Ala
# Original script put together by Gabi Mukai with aditional code from  Johanna Wren
# 
# Date last updated: 06/28/2023

## Load in libraries
# Includes some extra stuff from the tutorial that I didn't use


import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
# import feather
from scipy import interpolate

from parcels import FieldSet, ParticleSet, JITParticle, ScipyParticle, DiffusionUniformKh, Variable, Field, GeographicPolar, Geographic, ErrorCode, plotTrajectoriesFile, NestedField, AdvectionRK4
from datetime import timedelta as delta

from operator import attrgetter

# get_ipython().run_line_magic('matplotlib', 'inline')
#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#from matplotlib.colors import ListedColormap
#from matplotlib.lines import Line2D
#from copy import copy
#import cmocean


## Initial conditions (start & end dates, simulation depth, etc...)


file_name="toots_test_06282023_oocysts_model_2019.zarr"
path2 = '/home/pi/wongalaj/Ciannelli_Lab/wongalaj/TOOTS/toots_test_06282023_oocysts_model_2019.nc' # give path to where .nc file is
feather_path = 'toots_test_06282023_oocysts_model_2019.feather' # assign name of feather file
# nc_file = "toots_test_06282023_oocysts_model_2019.nc"

startDate = '2019-01-01' # YYYY-MM-DD 
endDate = '2019-12-31' # YYYY-MM-DD

run_days = 365  # day 2018: 240, 2019: 365, 2020: 366, 2021: 365

# 2018 runs=  2018-05-05 - 2018-12-31 
# 2019 - 2021 runs = 20xx-01-01 - 20xx-12-31

simDepthMHI = 1

simDepthKANE = 1

simDepthOSS = 1 

kh = 10   # This is the eddy diffusivity in m2/s

pld = 120   # in days

# ## Create MHI ROMS fieldset (online)


## MHI BASE

file2 = 'https://pae-paha.pacioos.hawaii.edu/erddap/griddap/roms_hiig' # MHI ROMS (BASE)

ds2 = xr.open_dataset(file2)  # this puts the opendap data into a xarray dataset

# ds2.load() # https://github.com/pydata/xarray/issues/593 # work around to prevent IndexError

myDat2 = ds2.isel(**{'depth': simDepthMHI}).sel(**{'time': slice(startDate,endDate)}) # subset based on time and depth layer

variables = {'U': 'u', 'V': 'v'}
dimensions = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}

U=Field.from_xarray(myDat2.u, 'U', dimensions, allow_time_extrapolation=True)
V=Field.from_xarray(myDat2.v, 'V', dimensions, allow_time_extrapolation=True)


## 3. Generate fieldset
fieldset = FieldSet(U, V)


# ## Displacement

# ### Code to make landmask (MHI)

# In[ ]:


def make_landmask1(fielddata):
    """Returns landmask where land = 1 and ocean = 0
    fielddata is a netcdf file.
    """
    datafile = Dataset(fielddata)

    landmask = datafile.variables['u'][0, 0]
    landmask = np.ma.masked_invalid(landmask)
    landmask = landmask.mask.astype('int')

    return landmask


# ### Import velocity field for MHI ROMS

# In[ ]:


file_path2 = '/nfs7/CEOAS/Ciannelli_Lab/wongalaj/TOOTS/roms_hiig_2018_05_05_685a_c38d_b3fe.nc'


# ### Make landmask for MHI ROMS
# "land = 1" and "ocean = 0"


landmask_mhi = make_landmask1(file_path2)

# ## Detect the coast
# We can detect the edges between land and ocean nodes by computing the Laplacian with the 4 nearest neighbors [i+1,j], [i-1,j], [i,j+1] and [i,j-1]:
# 
# 
# ∇2landmask=∂xxlandmask+∂yylandmask,
# 
# and filtering the positive and negative values. This gives us the location of coast nodes (ocean nodes next to land) and shore nodes (land nodes next to the ocean).
# 
# Additionally, we can find the nodes that border the coast/shore diagonally by considering the 8 nearest neighbors, including [i+1,j+1], [i-1,j+1], [i-1,j+1] and [i-1,j-1].



def get_coastal_nodes(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes directly
    next to land. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')

    return coastal

def get_shore_nodes(landmask):
    """Function that detects the shore nodes, i.e. the land nodes directly
    next to the ocean. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the shore nodes, the shore nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore



def get_coastal_nodes_diagonal(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes where 
    one of the 8 nearest nodes is land. Computes the Laplacian of landmask
    and the Laplacian of the 45 degree rotated landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    coastal = np.ma.masked_array(landmask, mask_lap > 0)
    coastal = coastal.mask.astype('int')
    
    return coastal
    
def get_shore_nodes_diagonal(landmask):
    """Function that detects the shore nodes, i.e. the land nodes where 
    one of the 8 nearest nodes is ocean. Computes the Laplacian of landmask 
    and the Laplacian of the 45 degree rotated landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the shore nodes, the shore nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap += np.roll(landmask, (-1,1), axis=(0,1)) + np.roll(landmask, (1, 1), axis=(0,1))
    mask_lap += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (1, -1), axis=(0,1))
    mask_lap -= 8*landmask
    shore = np.ma.masked_array(landmask, mask_lap < 0)
    shore = shore.mask.astype('int')

    return shore



coastal_mhi = get_coastal_nodes_diagonal(landmask_mhi)
shore_mhi = get_shore_nodes_diagonal(landmask_mhi)


# #### Check if coastal and shore nodes was created



# coastal_mhi
# shore_mhi


# ## Assigning coastal velocities
# For the displacement kernel we define a velocity field that pushes the particles back to the ocean. This velocity is a vector normal to the shore.
# 
# For the shore nodes directly next to the ocean, we can take the simple derivative of landmask and project the result to the shore array, this will capture the orientation of the velocity vectors.
# 
# For the shore nodes that only have a diagonal component, we need to take into account the diagonal nodes also and project the vectors only onto the inside corners that border the ocean diagonally.
# 
# Then to make the vectors unitary, we normalize them by their magnitude.


def create_displacement_field(landmask, double_cell=False):
    """Function that creates a displacement field 1 m/s away from the shore.

    - landmask: the land mask dUilt using `make_landmask`.
    - double_cell: Boolean for determining if you want a double cell.
      Default set to False.

    Output: two 2D arrays, one for each camponent of the velocity.
    """
    shore = get_shore_nodes(landmask)
    shore_d = get_shore_nodes_diagonal(landmask) # bordering ocean directly and diagonally
    shore_c = shore_d - shore                    # corner nodes that only border ocean diagonally
    
    Ly = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0) # Simple derivative
    Lx = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
    
    Ly_c = np.roll(landmask, -1, axis=0) - np.roll(landmask, 1, axis=0)
    Ly_c += np.roll(landmask, (-1,-1), axis=(0,1)) + np.roll(landmask, (-1,1), axis=(0,1)) # Include y-component of diagonal neighbours
    Ly_c += - np.roll(landmask, (1,-1), axis=(0,1)) - np.roll(landmask, (1,1), axis=(0,1))
    
    Lx_c = np.roll(landmask, -1, axis=1) - np.roll(landmask, 1, axis=1)
    Lx_c += np.roll(landmask, (-1,-1), axis=(1,0)) + np.roll(landmask, (-1,1), axis=(1,0)) # Include x-component of diagonal neighbours
    Lx_c += - np.roll(landmask, (1,-1), axis=(1,0)) - np.roll(landmask, (1,1), axis=(1,0))
    
    v_x = -Lx*(shore)
    v_y = -Ly*(shore)
    
    v_x_c = -Lx_c*(shore_c)
    v_y_c = -Ly_c*(shore_c)
    
    v_x = v_x + v_x_c
    v_y = v_y + v_y_c

    magnitude = np.sqrt(v_y**2 + v_x**2)
    # the coastal nodes between land create a problem. Magnitude there is zero
    # I force it to be 1 to avoid problems when normalizing.
    ny, nx = np.where(magnitude == 0)
    magnitude[ny, nx] = 1

    v_x = v_x/magnitude
    v_y = v_y/magnitude

    return v_x, v_y



# create displacement field for MHI ROMS
v_x_c, v_y_c = create_displacement_field(landmask_mhi)


# ## Calculate the distance to the shore
# In this tutorial, we will only displace particles that are within some distance (smaller than the grid size) to the shore.
# 
# For this we map the distance of the coastal nodes to the shore: Coastal nodes directly neighboring the shore are 1dx away. Diagonal neighbors are 2‾√dx away. The particles can then sample this field and will only be displaced when closer than a threshold value. This gives a crude estimate of the distance.



def distance_to_shore(landmask, dx=1):
    """Function that computes the distance to the shore. It is based in the
    the `get_coastal_nodes` algorithm.

    - landmask: the land mask dUilt using `make_landmask` function.
    - dx: the grid cell dimension. This is a crude approxsimation of the real
    distance (be careful).

    Output: 2D array containing the distances from shore.
    """
    ci = get_coastal_nodes(landmask) # direct neighbours
    dist = ci*dx                     # 1 dx away
    
    ci_d = get_coastal_nodes_diagonal(landmask) # diagonal neighbours
    dist_d = (ci_d - ci)*np.sqrt(2*dx**2)       # sqrt(2) dx away
        
    return dist+dist_d


# create land mask for MHI ROMS
d_2_s_c = distance_to_shore(landmask_mhi)


# ### Add displacement for u and v


u_displacement = v_x_c
v_displacement = v_y_c


fieldset.add_field(Field('dispU', data=u_displacement,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))
fieldset.add_field(Field('dispV', data=v_displacement,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))

fieldset.dispU.units = GeographicPolar()
fieldset.dispV.units = Geographic()


# ### Add landmask and distance to shore


fieldset.add_field(Field('landmask', landmask_mhi,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))

fieldset.add_field(Field('distance2shore', d_2_s_c,
                         lon=fieldset.U.grid.lon, lat=fieldset.U.grid.lat,
                         mesh='spherical'))


# ### Add eddy diffusivity

# Add even diffusivity to the fieldset
fieldset.add_constant_field('Kh_zonal', kh, mesh='spherical')
fieldset.add_constant_field('Kh_meridional', kh, mesh='spherical')   


# Add PLD to fieldset
fieldset.add_constant('pld', (pld*86400))


# fieldset.distance2shore.show(domain={'N':21.3, 'S':21, 'E':-156.5, 'W':-157.5})


# ## Particle and Kernels
# The distance to shore, used to flag whether a particle must be displaced, is stored in a particle Variable d2s. To visualize the displacement, the zonal and meridional displacements are stored in the variables dU and dV.
# 
# To write the displacement vector to the output before displacing the particle, the set_displacement kernel is invoked after the advection kernel. Then only in the next timestep are particles displaced by displace, before resuming the advection.

# ### Define DisplacementParticle


class DisplacementParticle(JITParticle):
    dU = Variable('dU')
    dV = Variable('dV')
    d2s = Variable('d2s', initial=1e3)
    age = Variable('age', dtype=np.float32, initial=0.)
    releaseSite = Variable('releaseSite', dtype=np.int32)
    IslandReleaseSite = Variable('IslandReleaseSite', dtype=np.int32)


# ### Define set_displacement


def set_displacement(particle, fieldset, time):
    particle.d2s = fieldset.distance2shore[time, particle.depth,
                               particle.lat, particle.lon]
    if  particle.d2s < 0.5:
        dispUab = fieldset.dispU[time, particle.depth, particle.lat,
                               particle.lon]
        dispVab = fieldset.dispV[time, particle.depth, particle.lat,
                               particle.lon]
        particle.dU = dispUab
        particle.dV = dispVab
    else:
        particle.dU = 0.
        particle.dV = 0.
    


# ### Define displace


def displace(particle, fieldset, time):    
    if  particle.d2s < 0.5:
        particle.lon += particle.dU*particle.dt
        particle.lat += particle.dV*particle.dt
        


# ## MERGE Particle 

# Use this to have particles interact with one another
# class MergeParticle(ScipyInteractionParticle):
#     nearest_neighbor = Variable('nearest_neighbor', dtype=np.int64, to_write=False)
#     mass = Variable('mass', initial=1, dtype=np.float32)
#     dU = Variable('dU')
#     dV = Variable('dV')
#     d2s = Variable('d2s', initial=1e3)
#     age = Variable('age', dtype=np.float32, initial=0.)
#     releaseSite = Variable('releaseSite', dtype=np.int32)
#     IslandReleaseSite = Variable('IslandReleaseSite', dtype=np.int32)


## Delete Particle
# For when it goes past the boundary


def DeleteParticle(particle, fieldset, time):
    print('deleted particle')
    particle.delete()


# Age Particle
# I want the model to record age at each timestep so it can be used for analyses over time

# def Age(particle, fieldset, time):
#    particle.age += particle.dt


def Ageing(particle, fieldset, time):
    particle.age += particle.dt
    if particle.age >= fieldset.pld:
        particle.delete()
    
    
### MHI ROMS Boundaries


# add min and max lon and lat for the ROMS
fieldset.add_constant('min_lon', -163.83070)

fieldset.add_constant('max_lon', -152.51930)

fieldset.add_constant('min_lat', 17.01843)

fieldset.add_constant('max_lat', 23.98238)


## Simulation

### Number of particle released per location

# constant release number
nrepeat = 500 # how many times do you want locations to repeat 500


# Read in the seeding location file
source_loc = pd.read_csv("all_release_locs_oocysts_final.csv") # read in file I make 


### Release particles only in a specific area


# kane_source = source_loc[source_loc["ID"].isin([30, 41, 42, 43, 44, 45, 29, 46, 28, 25, 40, 39, 38, 37, 36, 26])]
# print(kane_source)

# Kaneohe = 30, 41, 42, 43, 44, 45, 29, 46, 28
# Oahu SS = 25,40,39,38,37,36,26


### Hydrological model particle forcing 
 
# Create for loop to release number of particles based on oocysts hydrological model output 

toots_norm = source_loc.oocysts_normalized #assign column of oocysts_normalized to toots-norm
# print(toots_norm)

par_release_prop = nrepeat * toots_norm # multiple nrepeat to the toots_norm data to get the number of particles that need to be added in addition to the 500 minimum

# print(par_release_prop)

release_par_final = par_release_prop + nrepeat # add the additional particles to nrepeat (500) so all locations release a minimum of 500 particles (max of 1000)
# print(release_par_final[21:45])
# print (np.finfo(release_par_final).max) # look at min of toots_norm
# print (np.finfo(par_release_prop).max) # look at min of toots_norm
# print (np.finfo(toots_norm).min) # look at min of toots_norm0


# create empty np.array to append data too
lon_fin = np.array([]) 
lat_fin = np.array([])
site_fin = np.array([])
islandsite_fin = np.array([])

# type(lon1)
# print(lon_fin)

for i in range(0,117):
#    print (i) # what i-th is loop at 
    
     tmp = release_par_final[i] # get out number of particles to be released at that specific location using [i]

     lon1 = np.repeat(source_loc.lon[i], tmp) # subset out the variables I need to be repeated in the pset and replicate it by number of particles needed
     lat1 = np.repeat(source_loc.lat[i], tmp)
     site1 = np.repeat(source_loc.ID[i], tmp)
     islandsite1 = np.repeat(source_loc.island_release[i], tmp)
    
     lon_fin = np.append(lon_fin, lon1) # bind all of the initial conditions data to one final vector for the model psete input
     lat_fin = np.append(lat_fin, lat1)
     site_fin = np.append(site_fin, site1)
     islandsite_fin = np.append(islandsite_fin, islandsite1)
    

# print(lat_fin.shape) # look at dimension (shape) of  numpy array
# type(lon1) # look at what is class of object
# print(lon1) # test to make sure it is working


# Release location from the file read in above
# load in lon and lat from file
habilon = lon_fin # np.repeat(source_loc.lon, nrepeat) # np.repeat(infile[1],nrepeat) #add repeat if applicable 
habilat = lat_fin # np.repeat(source_loc.lat, nrepeat) # np.repeat(infile[0],nrepeat)
habisite = site_fin # np.repeat(source_loc.ID, nrepeat)
islandsite = islandsite_fin # np.repeat(source_loc.island_release, nrepeat)


# ## Constant release of particles 
# Assign pset variables 

# habilon = np.repeat(source_loc.lon, nrepeat) # repeat lon (and other variables) by whatever nrepeat is 
# habilat = np.repeat(source_loc.lat, nrepeat) 
# habisite = np.repeat(source_loc.ID, nrepeat)
# islandsite = np.repeat(source_loc.island_release, nrepeat)


## Define the pset and associated variables


# Time interval between particle release (in seconds)

release_days = 15  # how often to release particles (original == 15 days)

release_int = release_days*86400 # 15 days converted to seconds

# Start date for release (if you want it different from the first day of the currents in the fielset)
#start_date = datetime(2000, 1, 16)


## Define the pset
pset = ParticleSet.from_list(fieldset=fieldset, #parameter found in ParticleSet
                             pclass=DisplacementParticle, #parameter found in ParticleSet
                             lon=habilon, #parameter found in ParticleSet
                             lat=habilat, #parameter found in ParticleSet
                             releaseSite=habisite, 
                             IslandReleaseSite=islandsite,
                             repeatdt=release_int) #parameter found in ParticleSet
 


## Assign kernels, label output file, define run_days


kernels = pset.Kernel(displace) + pset.Kernel(AdvectionRK4) + pset.Kernel(DiffusionUniformKh) + pset.Kernel(set_displacement) + pset.Kernel(Ageing)

output_file = pset.ParticleFile(name=file_name, outputdt=delta(hours=4))


## Execute model


# don't print depth
pset.set_variable_write_status('depth', False)

pset.execute(kernels, 
             runtime=delta(days=run_days), 
             dt=delta(minutes=15),
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
             output_file=output_file)



# now stop the repeated release
pset.repeatdt = None

# now continue running for the remaining length of the PLD
pset.execute(kernels,
            runtime=delta(days=pld+1),
            dt=delta(minutes=15), 
            recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
            output_file=output_file)

output_file.close()


## Convert .nc to .feather file 

import pyarrow.feather as feather

# reading zarr file to ncdf 

dat_from_zarr = xr.open_zarr(store = file_name)
dat_from_zarr.to_netcdf(path2)

ds = xr.open_dataset(path2) # load in the .nc file

df = ds.to_dataframe() # convert ds to a dataframe called df

feather.write_feather(df,feather_path) # write df as a feather file


