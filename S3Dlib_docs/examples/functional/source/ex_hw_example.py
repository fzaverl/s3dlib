import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Hello World Example

# 1. Define function to examine .....................................

def geo_map(xyz) :
    x,y,z = xyz
    X,Y = 3*x, 3*y
    Z1 = np.exp(-X**2 - Y**2)
    Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
    Z = Z1-Z2
    return x,y,Z

# 2. Setup and map surfaces .........................................
rez = 6

surface = s3d.PlanarSurface(rez)
surface.map_geom_from_op( geo_map )
surface.map_cmap_from_op(lambda xyz : xyz[2])
surface.shade().hilite(.5)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(0.75))
fig.text(0.975,0.975,str(surface), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
maxmin = ( -0.8,0.8 )
ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
plt.colorbar(surface.cBar_ScalarMappable, ax=ax,  shrink=0.6 )
ax.set_axis_off()
ax.view_init(20,-55)

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()