import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Matplot: Color Mapping Normals

# 1. Define function to examine .....................................

def pringle(rtz) :
    r,t,z = rtz
    xy = r*r*np.cos(t)*np.sin(t)
    Z = 2*np.sin(-xy )
    return r,t,Z

# 2. Setup and map surfaces .........................................
cboard = cmu.rgb_cmap_gradient([0.25,0.15,0], [1,.9,.75], 'cardboard' )

saddle = s3d.PolarSurface(4)
saddle.map_geom_from_op( pringle )
saddle.map_cmap_from_normals( cmap=cboard, direction=[1,1,1] )

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(0.8))
fig.text(0.975,0.975,str(saddle), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax.xaxis.set_major_locator(LinearLocator(5))
ax.yaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_locator(LinearLocator(5))
ax.set_title('pringle surface')

ax.add_collection3d(saddle)

plt.show()