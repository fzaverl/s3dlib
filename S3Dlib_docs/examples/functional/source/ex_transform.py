import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Coordinate Transform

# 1. Define function to examine .....................................

def dipole(rtp):
    r,t,p = rtp
    R = np.abs(r*np.cos(p))
    return R,t,p

# 2. Setup and map surfaces .........................................
rez = 4
cmap = cmu.rgb_cmap_gradient([0.25,0.15,0,.1], [1,1,1,0.1])

demo = s3d.SphericalSurface(rez)
demo.map_geom_from_op(dipole)
demo.transform(scale=[2,2,1], rotate=s3d.eulerRot(20,30,15))
demo.map_cmap_from_normals(cmap=cmap)

# 3. Construct figures, add dasurfaces, and plot ....................

fig = plt.figure(figsize=plt.figaspect(1))
ax = fig.add_subplot(111, projection='3d')
s3d.standardAxis( ax, negaxis=False )

ax.add_collection3d(demo)
ax.add_collection3d(demo.get_transformAxis())

fig.tight_layout()
plt.show()

