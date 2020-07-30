import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Datagrid Geometric Mapping

# 1. Define functions to examine ....................................

Z=np.load('data/jacksboro_fault_dem.npz')['elevation']
datagrid = np.flip(Z,0)

# 2. Setup and map surfaces .........................................
rez = 6
cmu.rgb_cmap_gradient( [0.25,0.15,0], [1,.9,.75], 'cardboard' )

plate = s3d.PolarSurface(rez,basetype='squ_c')
plate.map_geom_from_datagrid(datagrid, 0.075)
plate.transform(scale=1.2, rotate=s3d.eulerRot(45,45,180))
plate.map_cmap_from_normals(cmap='cardboard',direction=[-1,1,1])

tube = s3d.CylindricalSurface(rez, basetype='squ_s')
tube.map_geom_from_datagrid(datagrid, 0.15)
tube.transform(rotate=s3d.eulerRot(135,0))
tube.map_cmap_from_normals(cmap='cardboard',direction=[-1,1,1])

# 3. Construct figures, add surfaces, and plot .....................

fig = plt.figure(figsize=plt.figaspect(0.5))
fig.text(0.47,0.975,str(plate), ha='right', va='top', fontsize='smaller')
fig.text(0.9,0.975,str(tube), ha='right', va='top', fontsize='smaller', multialignment='right')

ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')
minmax = (-1,1)
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax2.set(xlim=minmax, ylim=minmax, zlim=minmax)
ticks=[-1,-.5,0,.5,1]
ax1.set_xticks(ticks)
ax1.set_yticks(ticks)
ax1.set_zticks(ticks)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_zticks(ticks)

ax1.add_collection3d(plate)
ax2.add_collection3d(tube)

plt.show()