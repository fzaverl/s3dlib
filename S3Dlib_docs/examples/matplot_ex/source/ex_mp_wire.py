
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import axes3d
import s3dlib.surface as s3d

#.. Matplotlib Examples: Datagrid Wireframe plot

# 1. Define function to examine .....................................

X, Y, Z = axes3d.get_test_data()

# 2. Setup and map surfaces .........................................
rez=3

surface = s3d.PlanarSurface(rez, basetype='oct1')
surface.map_geom_from_datagrid( Z )

# 3. Construct figure, add surface, plot ............................

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(0,1) )
ax.xaxis.set_major_locator(LinearLocator(5))
ax.yaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_locator(LinearLocator(5))

ax.add_collection3d(surface.edges)

plt.show()