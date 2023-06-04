import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

# 1. Define function to examine ....................................

def townsend(xyz) :
    x,y,z = xyz
    A = -np.cos( (x-0.1)*y )**2
    B = -np.sin(3*x + y)
    return x, y, (A + x*B)

# 2. Setup and map surface .........................................

plane = s3d.PlanarSurface(6).domain( 2.25, (-2.5, 1.75) )
plane.map_geom_from_op(townsend)
plane.map_cmap_from_op()     # default: z-direction

contours = plane.contourLineSet(20)
contours.map_to_plane( -5 )  # default: xy plane
contours.set_linewidth(1)

# 3. Construct figure & axes, add surface & contours, show .........

fig = plt.figure()
ax = plt.axes(projection='3d', proj_type='ortho')
ax.set(xlabel='X',ylabel='Y',zlabel='Z')

s3d.auto_scale(ax,plane,contours)
ax.add_collection3d(plane.shade(0.5).hilite(.5))
ax.add_collection3d(contours)

fig.tight_layout(pad=0)
plt.show()