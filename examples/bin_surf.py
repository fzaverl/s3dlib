import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

# 1. Define function to examine .....................................
def twistFunction(rtz,twists=6) :
    r,t,z = rtz
    phi = twists*t/2
    w = 0.33*z 
    R = 1 + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

# 2. Setup and map surfaces .........................................
bcmap = cmu.binary_cmap('silver', 'sandybrown', name='slvr_brwn' )

surface = s3d.CylindricalSurface(5, basetype='squ_s', cmap=bcmap)
surface.map_geom_from_op( twistFunction )

# 3. Construct figures, add surface, plot ...........................
minmax = (-.8,.8)
fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()

ax.view_init(azim=-70)
surface.map_cmap_from_normals(direction=ax)
ax.add_collection3d(surface.shade(ax=ax).hilite(ax=ax))

fig.tight_layout(pad=0)
plt.show()