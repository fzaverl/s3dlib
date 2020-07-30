import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Polar Coordinates to XYZ

# 1. Define function to examine .....................................

def dinisurf(rtz) :
    r,t,z = rtz
    a, b = 1, 0.2
    T = 2*t
    x = a*np.cos(T)*np.sin(r)
    y = a*np.sin(T)*np.sin(r)
    z = a*(np.cos(r) + np.log(np.tan(r/2))) + b*T
    return x,y,z

# 2. Setup and map surfaces .........................................
rez = 4

surface = s3d.PolarSurface(rez, basetype='hex_c', minrad=0.01)
surface.map_cmap_from_op( lambda rtz: rtz[0] , cmap='inferno' )
surface.map_geom_from_op( dinisurf, returnxyz=True )

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975, "Dini's Surface", ha='right', va='top', fontsize='larger', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-.75,.75), ylim=(-.75,.75), zlim=(-3,1) )
ax.set_axis_off()
ax.view_init(elev=20)
ax.add_collection3d(surface)

fig.tight_layout()
plt.show()