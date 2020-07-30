import matplotlib.pyplot as plt
import numpy as np
import s3dlib.surface as s3d

# 1. Define function to examine ....................................

def planarfunc(xyz) :
    x,y,z = xyz
    r = np.sqrt( x**2 + y**2)
    Z = np.sin( 6.0*r )/2
    return x,y,Z

# 2. Setup and map surface .........................................

surface = s3d.PlanarSurface(4)
surface.map_geom_from_op( planarfunc )

# 3. Construct figure, add surface, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))

ax.add_collection3d(surface.edges)

plt.show()
