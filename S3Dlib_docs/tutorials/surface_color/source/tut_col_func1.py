import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

# 1. Define function to examine .....................................

def f_func(xyz) :
    x,y,z = xyz
    r = np.sqrt( x**2 + y**2)
    Z = np.sin( 6.0*r )/2
    return x,y,Z

# 2. Setup and map surfaces .........................................

surface_1 = s3d.PlanarSurface(5, cmap='RdYlGn' )
surface_1.map_geom_from_op( f_func )
surface_1.map_cmap_from_op( lambda xyz: f_func(xyz)[2])

surface_2 = s3d.PlanarSurface(5, cmap='RdYlGn' )
surface_2.map_geom_from_op( f_func )
surface_2.map_cmap_from_op( lambda xyz: xyz[2] )

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=(8,4.5))

ax1 = fig.add_subplot(121, projection='3d')
ax1.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax1.set_title('z-axis: f(x,y)\ncolor: f(x,y)')

ax1.add_collection3d(surface_1)
# .........
ax2 = fig.add_subplot(122, projection='3d')
ax2.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax2.set_proj_type('ortho')
ax2.view_init(elev=90, azim=-90)
ax2.set_title('elev=90, azim=-90\nproj_type=ortho\n\ncolor: f(x,y)')

ax2.add_collection3d(surface_2)

fig.tight_layout()
plt.show()