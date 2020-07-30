import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Viewport Comparison Plot

# 1. Define functions to examine ....................................

def color_surf(surface) :
    vpBlue   = [0.13, 0.45, 0.30, 0.85]
    vpOrange = [0.95, 0.30, 0.10, 0.55]
    vpGreen  = [0.00, 0.65, 0.06, 1.00]
    surface.map_color_from_image('data/C0b_color.png',viewport=vpBlue)
    surface.map_color_from_image('data/C1b_color.png',viewport=vpOrange)
    surface.map_color_from_image('data/C2b_color.png',viewport=vpGreen)
    surface.shade(0.5,direction=[0,1,1])
    return

# 2. Setup and map surfaces .........................................
rez=6
surf = [None]*4

surf[0] = s3d.PlanarSurface(rez+1, basetype='oct1', facecolor='beige')
surf[1] = s3d.PolarSurface(rez+1,basetype='hex', facecolor='beige')
surf[2] = s3d.CylindricalSurface(rez,facecolor='beige')
surf[3] = s3d.SphericalSurface(rez,facecolor='beige')

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=(6,6) )
minmax = (-0.9,0.9)
for i in range(len(surf)) :
    surface = surf[i]
    ax = fig.add_subplot(2,2,i+1, projection='3d')
    s3d.standardAxis(ax, length=1.5, offset=1,negaxis=False)
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    if i < 2 : ax.view_init(90,-90)
    else :     ax.view_init(15,30)
    color_surf(surface)
    ax.add_collection3d(surface)

fig.tight_layout()
plt.show()