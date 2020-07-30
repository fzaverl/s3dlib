import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
from cmap_xtra import Lab_cmap_gradient

#.. Multiple Geometric Maps 2

# 1. Define functions to examine ....................................

def Dodecahedron(rez) :
    v,f = s3d.SphericalSurface.get_dodecahedron()
    surface = s3d.Surface3DCollection(v,f)
    surface.triangulate(rez)
    surface.name = 'dodecahedron'
    return surface

def burst(xyz) :
    r,t,p = s3d.SphericalSurface.coor_convert(xyz, tocart=False)
    R = np.power(r,4)
    Rtp = np.vstack((R,t,p))
    XYZ = s3d.SphericalSurface.coor_convert(Rtp, tocart=True)
    return XYZ

def radialColor(xyz) :
    r,t,p = s3d.SphericalSurface.coor_convert(xyz, tocart=False)
    return r

# 2. Setup and map surfaces .........................................
rez=5
cmap = Lab_cmap_gradient('maroon', 'yellow')

surface = Dodecahedron(rez)
surface.map_geom_from_op(burst)
surface.map_cmap_from_op(radialColor, cmap)
surface.shade().hilite(.7,focus=2)

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975, str(surface), \
    ha='right', va='top', fontsize='small', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-.7,.7), ylim=(-.7,.7), zlim=(-.7,.7))
ax.set_axis_off()
ax.view_init(elev=35, azim=-55)

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()