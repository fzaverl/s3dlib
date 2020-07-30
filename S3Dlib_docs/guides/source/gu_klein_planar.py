import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Klein Bottle, PlanarSurface construction

# 1. Define function to examine ....................................

def klein_planar(xyz) :
    x,y,z = xyz
    t = -np.pi*(x+1)
    p = np.pi*(y+1)/2
    #.........
    u = p
    v = t
    cU, sU = np.cos(u), np.sin(u)
    cV, sV = np.cos(v), np.sin(v)
    x = -(2/15)*cU* \
        (  ( 3 )*cV + \
           ( -30 + 90*np.power(cU,4) - 60*np.power(cU,6) + 5*cU*cV )*sU \
        )
    y = -(1/15)*sU* \
        (  ( 3 - 3*np.power(cU,2) -48*np.power(cU,4) +48*np.power(cU,6) )*cV + \
           (-60 + ( 5*cU - 5*np.power(cU,3) - 80*np.power(cU,5) + 80*np.power(cU,7) )*cV  )*sU \
        )
    z = (2/15)*( 3 + 5*cU*sU )*sV
    return x,y,z

# 2. Setup and map surface .........................................
rez=6
cmap = cmu.mirrored_cmap('viridis')
cmap = cmu.alpha_cmap(cmap,0.7)

surface = s3d.PlanarSurface(rez,basetype='oct1', linewidth=0 )
surface.map_geom_from_op( klein_planar, returnxyz=True )
surface.map_cmap_from_normals(cmap=cmap, direction=[1,1,1])
surface.transform(s3d.eulerRot(0,-90),translate=[0,0,2])

# 3. Construct figure, add surface plot ............................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975, "Klein Bottle", \
    ha='right', va='top', fontsize='larger', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-1.5,1.5)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.view_init(elev=20, azim=-125)

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()