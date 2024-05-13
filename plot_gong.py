plot_gong.py:
GONG PFSS extrapolation 
=======================
Calculating PFSS solution for a GONG synoptic magnetic field map. 
"""
############################################################################### 
# First, import required modules
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.coordinates import SkyCoord 
from sunpy.data import manager
import pfsspy
import pfsspy.tracing as tracing
from pfsspy.sample_data import get_gong_map 
#@manager.require('gong1',
#        'https://gong2.nso.edu/oQR/zqs/202009/mrzqs200901/'
#        'mrzqs200901t1304c2234_022.fits.gz') 
@manager.require('gong_map',
          'https://gong2.nso.edu/oQR/zqs/201207/mrzqs120723/' 
          'mrzqs120723t0554c2126_241.fits.gz', 
          'b85e18b3fee0d7565930c01a222c3962b62f6b91ffacb7d3e357beee89cdf142')
def get_my_gong_map(name, uri, file_name): 
return manager.get(name)
############################################################################### 
# Load a GONG magnetic field map
manager.skip_hash_check()
#manager.require('gong_map',
#         'https://gong2.nso.edu/oQR/zqs/202009/mrzqs200901/'
#         'mrzqs200901t1304c2234_022.fits.gz',
#         'aad927d8f617f32b72255b862c4910f13640fc7ca13edf982'
#         '88cd0735a2db6a0')
gong_fname = get_my_gong_map('gong_map',
          'https://gong2.nso.edu/oQr/zqs/201207/mrzqs120723/',
          'mrzqs120723t0554c2126_241.fits.gz') 
print(gong_fname)
gong_map = sunpy.map.Map(gong_fname)
############################################################################### 
# The PFSS solution is calculated on a regular 3D grid in (phi, s, rho), where
# rho = ln(r), and r is the standard spherical radial coordinate. We need to
# define the number of rho grid points, and the source surface radius.
nrho = 35 
rss = 2.5

 ############################################################################### 
# From the boundary condition, number of radial grid points, and source
# surface, we now construct an Input object that stores this information
pfss_in = pfsspy.Input(gong_map, nrho, rss)


def set_axes_lims(ax): 
ax.set_xlim(0, 360) 
ax.set_ylim(0, 180)
############################################################################### 
# Using the Input object, plot the input field
m = pfss_in.map
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot()
plt.colorbar()
ax.set_title('Input field')
set_axes_lims(ax)
fig.savefig('plot1.png') 
############################################################################### 
# Now calculate the PFSS solution
pfss_out = pfsspy.pfss(pfss_in)

 ############################################################################### 
# Using the Output object we can plot the source surface field, and the
# polarity inversion line.
ss_br = pfss_out.source_surface_br
# Create the figure and axes
fig = plt.figure()
ax = plt.subplot(projection=ss_br)


# Plot the source surface map
ss_br.plot()
# Plot the polarity inversion line
ax.plot_coord(pfss_out.source_surface_pils[0])
# Plot formatting
plt.colorbar()
ax.set_title('Source Surface Magnetic Field')
set_axes_lims(ax)
fig.savefig('plot2.png') 
############################################################################### 
# It is also easy to plot the magnetic field at an arbitrary height within
# the PFSS solution.
# Get the radial magnetic field at a given height
ridx = 15
br = pfss_out.bc[0][:, :, ridx]
# Create a sunpy Map object using output WCS
br = sunpy.map.Map(br.T, pfss_out.source_surface_br.wcs)

# Get the radial coordinate
r = np.exp(pfss_out.grid.rc[ridx])
# Create the figure and axes 
fig = plt.figure()
ax = plt.subplot(projection=br)
# Plot the source surface map
br.plot(cmap='RdBu')
# Plot formatting
plt.colorbar()
ax.set_title('$B_{r}$ ' + f'at r={r:.2f}' + '$r_{\\odot}$') set_axes_lims(ax)
fig.savefig('plot3.png')
############################################################################### 
# Finally, using the 3D magnetic field solution we can trace some field lines.
# In this case 64 points equally gridded in theta and phi are chosen and
# traced from the source surface outwards.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


tracer = tracing.FortranTracer()
r = 1.2 * const.R_sun
lat = np.linspace(-np.pi / 2, np.pi / 2, 8, endpoint=False) 
lon = np.linspace(0, 2 * np.pi, 8, endpoint=False)

lat, lon = np.meshgrid(lat, lon, indexing='ij') 
lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad


seeds = SkyCoord(lon, lat, r, frame=pfss_out.coordinate_frame)


field_lines = tracer.trace(seeds, pfss_out)


for field_line in field_lines:
  color = {0: 'black', -1: 'tab:blue', 1: 'tab:red'}.get(field_line.polarity) coords = field_line.coords
  coords.representation_type = 'cartesian'
  ax.plot(coords.x / const.R_sun,
    coords.y / const.R_sun, 
    coords.z / const.R_sun, 
    color=color, linewidth=1)
ax.set_title('PFSS Solution') 
fig.savefig('plot4.png') 
plt.show()


# sphinx_gallery_thumbnail_number = 4
