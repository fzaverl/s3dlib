from matplotlib import pyplot as plt
import numpy as np
import s3dlib.surface as s3d

# 1. Define functions to examine ....................................

def randfunc(rtz) :
    r,t,z = rtz
    sigma = 0.005
    R = r + sigma*np.random.rand( len(r) )
    return R,t,z

def twistFunction(rtz,twists) :
    r,t,z = rtz
    thickness = 0.5
    w = thickness*z 
    phi = 0.5*t*twists
    R = r + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

# 2. Setup and map surfaces .........................................

twist = s3d.CylindricalSurface(5, basetype='squ_s')
twist.map_cmap_from_op( lambda rtz : randfunc(rtz)[0], 'Oranges' )
twist.map_geom_from_op( lambda rtz : twistFunction(rtz,2) )
twist.shade(direction=[1,1,1]).hilite(direction=[0,1,1])

'''
# 2. Setup and map surfaces .........................................

twist = s3d.CylindricalSurface(5, basetype='squ_s')
twist.map_geom_from_op(randfunc)
twist.map_cmap_from_op( lambda rtz : rtz[0], 'Oranges' )
twist.map_geom_from_op( lambda rtz : twistFunction(rtz,2) )
twist.shade().hilite(direction=[0,1,1])

'''
# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(0.8))
info = str(twist) 
fig.text(0.975,0.975,info, ha='right', va='top', fontsize='smaller', multialignment='right')
ax = fig.add_subplot(111, projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax.set_axis_off()

ax.add_collection3d(twist)

fig.tight_layout()
plt.show()