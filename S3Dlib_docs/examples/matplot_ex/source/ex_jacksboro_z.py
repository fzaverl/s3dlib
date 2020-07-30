
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import s3dlib.surface as s3d

#.. Matplotlib Examples: Datagrid Geometry, 2

# 1. Define function to examine .....................................

with np.load('data/jacksboro_fault_dem.npz') as dem :
    z = dem['elevation']
    nrows, ncols = z.shape
    x = np.linspace(dem['xmin'], dem['xmax'], ncols)
    y = np.linspace(dem['ymin'], dem['ymax'], nrows)
    x, y = np.meshgrid(x, y)

region = np.s_[5:50, 5:50]
x, y, z = x[region], y[region], z[region]

datagrid = np.flip(z,0)

# 2. Setup and map surfaces .........................................
rez=5
ls = s3d.elev_azim_2vector(90,-135)

surface = s3d.PlanarSurface(rez, cmap='gist_earth')
surface.map_geom_from_datagrid( datagrid )
surface.map_cmap_from_op( lambda xyz : xyz[2] )
surface.shade(0.6,direction=ls).hilite(0.5,direction=ls, focus=0.5)
surface.scale_dataframe(x,y,datagrid)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure()
fig.text(0.975,0.975,str(surface), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-84.415,-84.375), ylim=(36.690,36.740), zlim=(350,700) )
ax.xaxis.set_major_locator(LinearLocator(5))
ax.yaxis.set_major_locator(LinearLocator(6))
ax.zaxis.set_major_locator(LinearLocator(8))
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')

ax.add_collection3d(surface)

plt.show()