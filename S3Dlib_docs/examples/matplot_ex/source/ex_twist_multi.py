import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Matplotlib Examples: Cylinderical Coordinates, multiple surfaces

# 1. Define function to examine .....................................

def twistFunction(rtz,twists) :
    r,t,z = rtz
    # Note: sliced surface needed due to discontinuity @ t=0 if twists is odd
    thickness =  0.33
    w = thickness*z
    phi = 0.5*t*twists
    R = 1 + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

# 2 & 3. Setup surfaces and plot ....................................
rez = 5
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboardMrrd',mirrored=True)

fig = plt.figure(figsize=plt.figaspect(0.6))

for i in range(1,7) :
    ax = fig.add_subplot(2,3,i, projection='3d')
    ax.set(xlim=(-0.8,0.8), ylim=(-0.8,0.8), zlim=(-0.8,0.8) )
    twist = s3d.CylindricalSurface(rez, basetype='squ_s')
    twist.map_geom_from_op( lambda rtz : twistFunction(rtz,i) )
    twist.map_cmap_from_normals(cmap='cardboardMrrd', direction=[1,1,1])
    ax.set_title('twists: '+str(i))
    ax.add_collection3d(twist)
    ax.set_axis_off()

fig.tight_layout()
plt.show()