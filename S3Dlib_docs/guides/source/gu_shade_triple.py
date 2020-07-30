import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Shade: Combined Example

# 1. Define function to examine .....................................

def wavefunc(xyz) :
    x,y,z = xyz
    X = 3*x-1
    Y = 3*y-1
    Z = np.cos( X**2 + Y**2 )/5
    return x,y,Z

# 2. Setup and map surfaces .........................................
rez = 7

wave = s3d.PlanarSurface(rez, basetype='oct1')
wave.map_geom_from_op( wavefunc )
wave.map_cmap_from_normals( 'copper' )
wave.shade()
wave.hilite(focus=2)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(1))
info = str(wave) + '\n' + wave.cmap.name + '-normals, shade, hilite'
fig.text(0.975,0.975,info, ha='right', va='top', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-.8, 0.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.set_axis_off()
ax.view_init( azim=20 )

ax.add_collection3d(wave)

fig.tight_layout()
plt.show()