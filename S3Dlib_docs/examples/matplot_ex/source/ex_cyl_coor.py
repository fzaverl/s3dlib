import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Matplotlib Examples: 3D surface in polar coordinates

# 1. Define function to examine ....................................

def polarfunc(rtz) :
    r,t,z = rtz
    R = 1.25*r  # radial direction scaled [0,1] -> [0,1.25]
    Z = ((R**2 - 1)**2)
    return R,t,Z

# 2. Setup and map surface .........................................

surface = s3d.PolarSurface(6)
surface.map_geom_from_op( polarfunc ).shade()
#surface.map_cmap_from_op( lambda rtz : rtz[2] , 'YlGnBu_r').shade(.5)

# 3. Construct figure, add surface, and plot ......................

fig = plt.figure()
ax = plt.axes(projection='3d')
maxmin = (-1.3,1.3)
ax.set(xlim=maxmin, ylim=maxmin, zlim=(0,1))

ax.set_xlabel(r'$\phi_\mathrm{real}$')
ax.set_ylabel(r'$\phi_\mathrm{im}$')
ax.set_zlabel(r'$V(\phi)$')

ax.add_collection3d(surface)

plt.show()