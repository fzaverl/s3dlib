import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Figure 8 Klein Bottle

# 1. Define function to examine ....................................

def fig8(rtp) :
    r,t,p = rtp
    R=2
    v = 2*p
    Q = ( R + np.cos(t/2)*np.sin(v) - np.sin(t/2)*np.sin(2*v) )
    x = Q*np.cos(t)
    y = Q*np.sin(t)
    z = np.sin(t/2)*np.sin(v) + np.cos(t/2)*np.sin(2*v)
    return x,y,z

# 2. Setup and map surface .........................................
rez=7
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboardMrrd',mirrored=True)

surface = s3d.SphericalSurface(rez,basetype='octa_c')
surface.map_geom_from_op( fig8, returnxyz=True )
surface.map_cmap_from_normals(cmap='cardboardMrrd')

# 3. Construct figure, add surface plot ............................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975, "Figure 8 Immersion of the Klein Bottle", \
    ha='right', va='top', fontsize='larger', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-2,2)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.view_init(elev=35, azim=-60)

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()