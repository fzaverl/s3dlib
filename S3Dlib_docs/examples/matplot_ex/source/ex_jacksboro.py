
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Matplotlib Examples: Datagrid Geometry

# 1. Define function to examine .....................................

z=np.load('data/jacksboro_fault_dem.npz')['elevation']
datagrid = np.flip( z[5:50, 5:50], 0 )

# 2. Setup and map surfaces .........................................
rez=5
cmu.rgb_cmap_gradient([0.25,0.15,0], [1,.9,.75], 'cardboard' )
ls = s3d.elev_azim_2vector(90,-135)

surface = s3d.PlanarSurface(rez,cmap='cardboard')
surface.map_geom_from_datagrid( datagrid )
surface.map_cmap_from_normals(direction=ls)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure()
fig.text(0.975,0.975,str(surface), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(0,1) )
ax.xaxis.set_major_locator(LinearLocator(5))
ax.yaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_locator(LinearLocator(6))

ax.add_collection3d(surface)

plt.show()