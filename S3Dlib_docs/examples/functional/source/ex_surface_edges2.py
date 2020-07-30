import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu  

#.. Wireframe Plot 2

# 1. Define function to examine .....................................

def fun009(xyz) :
    x,y,z = xyz
    a,b,c = 1,1,15
    lim = 0.001  # use L'Hosital's rule as x -> 0
    A = np.where(np.abs(x)<lim, np.ones(len(x)), np.divide( np.sin(a*x), a*x ) )
    B = np.where(np.abs(y)<lim, np.ones(len(y)), np.divide( np.sin(b*y), b*y ) )
    Z  = c*A*B
    return x,y,Z

def fun012(xyz) :
    x,y,z = xyz
    A = 0.9*np.exp( np.sin(2*x)*np.sin(0.2*y))
    B = 0.9*np.exp( np.sin(2*y)*np.sin(0.2*x))
    Z  = A*B
    return x,y,Z

# 2. Setup and map surfaces .........................................
rez=5
cmap = cmu.hsv_cmap_gradient([1.166,1,1],[0.333,1,1])
mcmap = cmu.mirrored_cmap(cmap)
surface012 = s3d.PlanarSurface(rez,basetype='oct1',linewidth=.3)
surface012.transform(scale=10)
surface012.map_geom_from_op(fun012)
surface012.map_cmap_from_op(lambda xyz : np.abs(fun012(xyz)[2]), cmap )
surface012.shade(.5).hilite(.3)
surface012.set_facecolor([0,0,0,0])

surface009 = s3d.PlanarSurface(rez,basetype='oct2',linewidth=.3)
surface009.transform(scale=10)
surface009.map_geom_from_op(fun009)
surface009.map_cmap_from_op(lambda xyz : np.abs(fun009(xyz)[2]), cmap )
surface009.shade(.5).hilite(.3)
surface009.set_facecolor([0,0,0,0])

# 3. Construct figures, add surface, plot ...........................

minmax = (-8,8)
fig = plt.figure(figsize=plt.figaspect(1/2), facecolor='black')
ax1 = fig.add_subplot(121, projection='3d')
ax1.set_title('      fun009',color='white', ha='left')
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax1.set_axis_off()
ax1.set_facecolor('black')
ax1.set_proj_type('ortho')
ax1.view_init(40,120)

ax1.add_collection3d(surface009)

ax2 = fig.add_subplot(122, projection='3d')
ax2.set_title('      fun012',color='white', ha='left')
ax2.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax2.set_axis_off()
ax2.set_facecolor('black')
ax2.set_proj_type('ortho')
ax2.view_init(40,120)

ax2.add_collection3d(surface012)

fig.tight_layout()
plt.show()