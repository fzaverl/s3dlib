import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

# 1. Define function to examine ....................................

def polarfunc(rtz) :
    r,t,z = rtz
    z = np.sin( 6.0*r )/2
    return r,t,z

# 2. Setup and map surface .........................................

surface = s3d.PolarSurface(4)
surface.map_geom_from_op( polarfunc )
surface.shade()

# 3. Construct figure, add surface, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))

ax.add_collection3d(surface)

plt.show()
