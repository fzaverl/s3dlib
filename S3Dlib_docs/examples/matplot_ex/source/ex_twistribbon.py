import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Matplotlib Examples: Cylindrical Coordinates

# 1. Define functions to examine ....................................

def twistFunction(rtz,twists) :
    r,t,z = rtz
    # Note: sliced surface needed due to discontinuity @ t=0 if twists is odd
    thickness = 0.33
    w = thickness*z 
    phi = 0.5*t*twists
    R = 1 + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

def ribbonFunc(rtz) :
    r,t,z = rtz
    min_radius, max_radius = 0.25, 0.95
    d = (max_radius-min_radius)/2
    R = d + min_radius + d*z
    Z = np.cos(R)*np.sin(3*t)
    return R,t,Z

# 2. Setup and map surfaces .........................................
rez = 5
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboardMrrd',mirrored=True)

twist = s3d.CylindricalSurface(rez, basetype='squ_s')
twist.map_geom_from_op( lambda rtz : twistFunction(rtz,1) )
twist.map_cmap_from_normals(cmap='cardboardMrrd', direction=[1,1,1])
#twist.map_cmap_from_op( lambda xyz : xyz[2] , 'Spectral')

ribbon = s3d.CylindricalSurface(rez, basetype='tri')
ribbon.map_geom_from_op( ribbonFunc )
ribbon.map_cmap_from_normals(cmap='cardboard',direction=[1,1,1])
#ribbon.map_cmap_from_op( lambda xyz : xyz[2] , 'CMRmap')

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(0.5/1.4))
fig.text(0.42,0.975,str(twist), ha='right', va='top', fontsize='smaller', multialignment='right')
fig.text(0.845,0.975,str(ribbon), ha='right', va='top', fontsize='smaller', multialignment='right')

ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')
ax1.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax2.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax1.xaxis.set_major_locator(LinearLocator(5))
ax1.yaxis.set_major_locator(LinearLocator(5))
ax1.zaxis.set_major_locator(LinearLocator(5))
ax2.xaxis.set_major_locator(LinearLocator(5))
ax2.yaxis.set_major_locator(LinearLocator(5))
ax2.zaxis.set_major_locator(LinearLocator(5))
ax1.set_title('Twist')
ax2.set_title('Ribbon')
plt.colorbar(twist.cBar_ScalarMappable, ax=ax1, ticks=np.linspace(-1,1,3), shrink=0.6 )
plt.colorbar(ribbon.cBar_ScalarMappable, ax=ax2, ticks=np.linspace(-1,1,3), shrink=0.6 )

ax1.add_collection3d(twist)
ax2.add_collection3d(ribbon)

#fig.tight_layout()
plt.show()