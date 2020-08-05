import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Guides: Shading, Highlighting and Color Mapped Normals

rez=4
color = 'peru'
width,height,scale = 6.75, 1.84, 1.15
minmax = (-1,1)

#. Figure 1 .........................................................
illum = [ (1,1,1), (0,1,1), (0,0,1), (-1,1,1) ]

width,height = scale*width, scale*height
fig = plt.figure(figsize=(width,height))
for i in range(len(illum)) :
    
    surface = s3d.SphericalSurface(rez,facecolor=color)
    surface.shade(direction=illum[i],contrast=0.7)
    surface.hilite(direction=illum[i],focus=2)

    ax = fig.add_subplot(141+i, projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    s3d.setupAxis(ax, length=1.7, width=2, color='black', offset=1.0, negaxis=True)
    ax.text(0,1.5,2,str(illum[i]), color = 'black', horizontalalignment='center', verticalalignment='center')
    ax.set_axis_off()
    ax.view_init(30, 45)

    ax.add_collection3d(surface)

fig.tight_layout()

#. Figure 2 .........................................................
illum = (1,0,1)

surface = s3d.SphericalSurface(rez,facecolor=color)
surface.shade(contrast=0.7,direction=illum)
surface.hilite(focus=2,direction=illum)

width = height
width,height = scale*width, scale*height
fig = plt.figure(figsize=(width,height))

ax = plt.axes(projection='3d')
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
s3d.setupAxis(ax, length=1.7, width=2, color='black', offset=1.0, negaxis=True)
ax.text(1.5,0,2,str(illum), color = 'black', horizontalalignment='center', verticalalignment='center')
ax.set_axis_off()
ax.add_collection3d(surface)

#....................................................................
plt.show()