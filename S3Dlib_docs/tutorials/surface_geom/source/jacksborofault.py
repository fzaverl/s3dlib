import matplotlib.pyplot as plt
import numpy as np
import s3dlib.surface as s3d

# 1. Define function to examine .....................................

dem = np.load('jacksboro_fault_dem.npz')
z = dem['elevation']
datagrid = z[5:50, 5:50]

# 2. Setup and map surfaces .........................................

surface = s3d.PlanarSurface(5)
surface.map_geom_from_datagrid( datagrid )
surface.shade()

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(0.75))
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(0,1) )
ax.view_init(elev=30, azim=120)

ax.add_collection3d(surface)

plt.show()
