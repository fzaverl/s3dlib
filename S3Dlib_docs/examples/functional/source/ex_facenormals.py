import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Face Normals Vector Field

# 1. Define functions to examine ....................................

def twistFunction(rtz,twists) :
    r,t,z = rtz
    # Note: sliced surface needed due to discontinuity @ t=0 if twists is odd
    thickness = 0.5
    w = thickness*z 
    phi = 0.5*t*twists
    R = 1 + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

# 2. Setup and map surfaces .........................................

twist = s3d.CylindricalSurface(2, basetype='squ_s', color=[0,.5,.5,.5])
twist.map_geom_from_op( lambda rtz : twistFunction(rtz,2) ).shade()
twist.set_edgecolor([0,0,0,0])
vf = twist.facenormals(scale=0.3,color='saddlebrown',width=0.75)

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(0.8))
info = str(twist) +'\n'+ str(vf)
fig.text(0.975,0.975,info, ha='right', va='top', fontsize='smaller', multialignment='right')
ax1 = fig.add_subplot(111, projection='3d')
ax1.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax1.xaxis.set_major_locator(LinearLocator(5))
ax1.yaxis.set_major_locator(LinearLocator(5))
ax1.zaxis.set_major_locator(LinearLocator(5))

ax1.add_collection3d(twist)
ax1.add_collection3d(vf)

fig.tight_layout()
plt.show()