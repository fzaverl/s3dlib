import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Shade: Combined Example 2

# 1. Define function to examine .....................................

def randfunc(rtp) :
    r,t,p = rtp
    sigma = 0.005
    R = r + sigma*np.random.rand( len(r) )
    return R,t,p

# 2. Setup and map surfaces .........................................
rez = 5

surface = s3d.SphericalSurface(rez,color='orange')
surface.map_geom_from_op(randfunc).shade(0.2).hilite(0.5)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=(3,3))
info = 'orange-color, shade, hilite'
fig.text(0.975,0.975,info, ha='right', va='top', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-.8,.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.view_init( azim=-55)

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()