import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

# 1. Define functions to examine ....................................

def twistFunction(rtz,twists=6) :
    r,t,z = rtz
    phi = twists*t/2
    w = 0.33*z 
    R = 1 + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

# 2 & 3. Construct figure & axes, add surfaces, show ................

fig = plt.figure(figsize=plt.figaspect(0.6))
for i in range(1,7) :
    ax = fig.add_subplot(2,3,i, projection='3d')
    ax.set(xlim=(-0.8,0.8), ylim=(-0.8,0.8), zlim=(-0.8,0.8) )
    ax.set_title('twists: '+str(i))
    ax.set_axis_off()

    twist = s3d.CylindricalSurface(5, basetype='squ_s', color=[1,.9,.75])
    twist.map_geom_from_op( lambda rtz : twistFunction(rtz,i) )
    twist.shade(direction=[1,1,1],ax=ax).hilite(direction=[1,1,1],ax=ax)

    ax.add_collection3d(twist)

fig.tight_layout()
plt.show()
