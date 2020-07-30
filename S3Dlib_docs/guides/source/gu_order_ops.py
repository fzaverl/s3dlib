import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. OOPs, Order of Operations

# 1. Define functions to examine ....................................

def torusFunc(rtz) :
    r,t,z = rtz
    ratio = 0.9
    f = 1 + ratio
    Z = ratio*np.sin(z*np.pi)/f
    R = (r + ratio*np.cos(z*np.pi))/f
    return R,t,Z

# 2. Setup and map surfaces .........................................
rez=5
shading, lumin = 0.5, [0,0,1]
rotation = s3d.eulerRot(105,45)
cmap_func = lambda rtz : 1+rtz[2]

torus_1 = s3d.CylindricalSurface(rez)
torus_1.map_geom_from_op( torusFunc )
torus_1.transform( rotate=rotation)
torus_1.map_cmap_from_op( cmap_func )
torus_1.shade(shading,lumin)
title_1 = r'geom$\Rightarrow$trans$\Rightarrow$color'

torus_2 = s3d.CylindricalSurface(rez)
torus_2.map_cmap_from_op( cmap_func )
torus_2.map_geom_from_op( torusFunc )
torus_2.transform(rotate=rotation)
torus_2.shade(shading,lumin)
title_2 = r'color$\Rightarrow$geom$\Rightarrow$trans'

torus_3 = s3d.CylindricalSurface(rez)
torus_3.transform(rotate=rotation)
torus_3.map_geom_from_op( torusFunc )
torus_3.map_cmap_from_op( cmap_func )
torus_3.shade(shading,lumin)
title_3 = r'trans$\Rightarrow$geom$\Rightarrow$color'

torus_4 = s3d.CylindricalSurface(rez)
torus_4.map_geom_from_op( torusFunc )
torus_4.map_cmap_from_op( cmap_func )
torus_4.transform(rotate=rotation)
torus_4.shade(shading,lumin)
title_4 = r'geom$\Rightarrow$color$\Rightarrow$trans'

torus_5 = s3d.CylindricalSurface(rez)
torus_5.map_cmap_from_op( cmap_func )
torus_5.transform(rotate=rotation)
torus_5.map_geom_from_op( torusFunc )
torus_5.shade(shading,lumin)
title_5 = r'color$\Rightarrow$trans$\Rightarrow$geom'

torus_6 = s3d.CylindricalSurface(rez)
torus_6.transform(rotate=rotation)
torus_6.map_cmap_from_op( cmap_func )
torus_6.map_geom_from_op( torusFunc )
torus_6.shade(shading,lumin)
title_6 = r'trans$\Rightarrow$color$\Rightarrow$geom'

# 3. Construct figures, add surfaces, and plot ......................
minmax = (-0.8, 0.8)
title = [ title_1, title_2, title_3, title_4, title_5, title_6 ]
torus = [ torus_1, torus_2, torus_3, torus_4, torus_5, torus_6 ]

fig = plt.figure(figsize=plt.figaspect(0.666))
for i in range(6) :
    ax = fig.add_subplot(2, 3, i+1, projection='3d')
    ax.set_axis_off()
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
    ax.set_title(title[i])
    ax.add_collection3d(torus[i])

plt.tight_layout()
plt.show()
