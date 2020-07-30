import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Multiple Geometric Maps

# 1. Define functions to examine ....................................

def torusFunc(rtz) :
    r,t,z = rtz
    ratio = .2
    Z = ratio*np.sin(z*np.pi)
    R = r + ratio*np.cos(z*np.pi)
    return R,t,Z

def knot(rtz) :
    r,t,z = rtz
    R = r*( (1+np.cos(5*t))/2 + 0.65*(1+np.cos(np.pi + 5*t))/2 )
    Z = z +  0.25*np.sin(5*t)
    return R,2*t,Z

# 2. Setup and map surfaces .........................................
rez = 5
revhsv = cmu.reversed_cmap('hsv')

torus = s3d.CylindricalSurface(rez)
torus.map_cmap_from_op( lambda xyz : xyz[1] , revhsv)
torus.map_geom_from_op( torusFunc )
torus.map_geom_from_op( knot )
torus.shade().hilite(.5)

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,str(torus), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-.8,.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.set_axis_off()
ax.view_init(azim=0)

ax.add_collection3d(torus)

fig.tight_layout()
plt.show()