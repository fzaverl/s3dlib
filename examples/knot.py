import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

# 1. Define function to examine .....................................
def knot(rtz) :
    r,t,z = rtz
    rho,zeta,delta = 0.25, 0.3, 0.3
    R = (1-rho)*(1-delta*np.sin(3*t)) + rho*np.cos(z*np.pi) 
    Z = rho*np.sin(z*np.pi) + zeta*np.cos(3*t)
    return R, 2*t, Z

# 2. Setup and map surface  .........................................
surface = s3d.CylindricalSurface(6)
surface.map_cmap_from_op(lambda c: c[1],'hsv')
surface.map_geom_from_op( knot )

# 3. Construct figure & axes, add surface, show .....................
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set( xlim=(-1,1),ylim=(-1,1),zlim=(-1,1) )
ax.add_collection3d(surface.shade().hilite(.75))
plt.show()