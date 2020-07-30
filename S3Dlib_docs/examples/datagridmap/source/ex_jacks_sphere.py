import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Spherical Datagrid Mapping 

# 1. Define functions to examine ....................................

Z=np.load('data/jacksboro_fault_dem.npz')['elevation']
datagrid = np.flip(Z,0)

# 2. Setup and map surfaces .........................................
rez=6

fault = s3d.SphericalSurface(rez, cmap='gist_earth')
fault.map_cmap_from_datagrid(datagrid)
fault.shade(direction=[0,1,1])

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1) )
fig.text(0.975,0.975,str(fault), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = fig.add_subplot(111, projection='3d')
minmax = (-0.8,0.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_title('Jacksboro Fault')
ax.set_axis_off()
ax.view_init(0,90)

ax.add_collection3d(fault)

fig.tight_layout()
plt.show()