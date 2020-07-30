import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Surface Addition

# 1. Define function to examine .....................................

def torusFunc(rtz) :
    # surface geometry f(V) -> V
    r,t,z = rtz
    ratio = .5
    Z = ratio*np.sin(z*np.pi)
    R = r + ratio*np.cos(z*np.pi)
    return R,t,Z

# 2. Setup and map surfaces .........................................
rez = 5
offset = 0.5
posOff, negOff = [offset,0,0], [-offset,0,0]
cmu.rgb_cmap_gradient( [0.25,0.15,0], [1,.9,.75], 'cardboard' )

torus_Z = s3d.CylindricalSurface(rez,basetype='squ' )
torus_Z.map_geom_from_op(torusFunc)
torus_Z.transform(translate=negOff)

torus_X = s3d.CylindricalSurface(rez,basetype='squ' )
torus_X.map_geom_from_op(torusFunc)
torus_X.transform(translate=posOff, rotate=s3d.eulerRot(0,90))

links = torus_X + torus_Z
links.map_cmap_from_normals(direction=[1,1,1], cmap='cardboard')

# 3. Construct figure, add surfaces, and plot ......................
fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,str(links), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = fig.add_subplot(111, projection='3d')
minmax = (-1.5,1.5)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)

ax.add_collection3d(links)

fig.tight_layout()
plt.show()