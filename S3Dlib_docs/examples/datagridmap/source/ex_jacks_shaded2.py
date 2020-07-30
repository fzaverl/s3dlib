import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Geometric and Color Datagrid Mapping, 2

# 1. Define function to examine .....................................

Z=np.load('data/jacksboro_fault_dem.npz')['elevation']
datagrid = np.flip(Z,0)

# 2. Setup and map surfaces .........................................
rez=6

surface = s3d.PlanarSurface(rez, basetype='oct1', cmap='gist_earth')
surface.map_cmap_from_datagrid( datagrid )
surface.map_geom_from_datagrid( datagrid, scale=0.2 ).shade(contrast=1.3)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=(7.5,2.5))
fig.text(0.975,0.975,str(surface), ha='right', va='top',
        fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(0,0.25) )
minc, maxc = surface.bounds['vlim']
cbar=plt.colorbar(surface.cBar_ScalarMappable, ax=ax,
        ticks=np.linspace(minc,maxc,5), shrink=0.6, pad=-.08  )
cbar.set_label('Elevation', rotation=270, labelpad = 15)
ax.set_axis_off()
ax.set_proj_type('ortho')
ax.view_init(elev=70, azim=60)

ax.add_collection3d(surface)
ax.add_collection3d(s3d.PlanarSurface(color='k'))

fig.tight_layout()
plt.show()