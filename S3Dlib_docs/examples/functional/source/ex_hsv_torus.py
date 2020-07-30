import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Cyclic Colormapped (angular)

# 1. Define functions to examine ....................................

def torusFunc(rtz) :
    r,t,z = rtz
    ratio = .45
    Z = ratio*np.sin(z*np.pi)
    R = r + ratio*np.cos(z*np.pi)
    return R,t,Z

# 2. Setup and map surfaces .........................................
rez= 5
cmap = cmu.reversed_cmap(cmu.hue_cmap())

torus = s3d.CylindricalSurface(rez).map_geom_from_op( torusFunc )
torus.map_cmap_from_op( lambda xyz : xyz[1] , cmap ).shade(.5)

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(0.75))
fig.text(0.975,0.975,str(torus), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax.xaxis.set_major_locator(LinearLocator(5))
ax.yaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_locator(LinearLocator(5))
minc = torus.bounds['vlim'][0]
maxc = torus.bounds['vlim'][1]
plt.colorbar(torus.cBar_ScalarMappable, ax=ax, ticks=np.linspace(minc,maxc,5), shrink=0.6 )

ax.add_collection3d(torus)

fig.tight_layout()
plt.show()