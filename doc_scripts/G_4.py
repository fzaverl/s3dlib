import copy
import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Guides: Shading, Highlighting and Color Mapped Normals -
#   Shading: Depth and Contrast

rez = 4
minmax = (-.7,.7)

#. Figure 1 - Depth .................................................
depthLevel = [0, 0.25, 0.5, 0.75, 1 ]

surface_ref = s3d.SphericalSurface(rez,color='peru')
surface_ref.map_color_from_image('data/earth.png')

width,height = 6.75, 1.6
fig = plt.figure(figsize=(width,height))
for i in range( len(depthLevel) ) :
    surface = copy.copy(surface_ref)
    surface.shade( depthLevel[i], direction=[1,0.2,1])
    ax = fig.add_subplot( 151+i , projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.set_title(str(depthLevel[i]) , fontsize='large')
    ax.set_axis_off()

    ax.add_collection3d(surface)

fig.tight_layout()

#. Figure 2 - Constrast .............................................
contrastLevel = [1,1.25,1.67,2,  1, 0.67, 0.5, 0.25 ]

width,height = 5.25, 3.0
fig = plt.figure(figsize=(width,height))
for i in range( len(contrastLevel) ) :
    surface = s3d.SphericalSurface(rez,color='peru')
    surface.shade( direction=[1,0.2,1], contrast=contrastLevel[i])
    ax = fig.add_subplot( 241+i , projection='3d')
    ax.set_title(str(contrastLevel[i]) , fontsize='large')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.set_axis_off()

    ax.add_collection3d(surface)

fig.tight_layout()

#....................................................................
plt.show()