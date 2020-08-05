import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Guides: Orientation -
#   Coordinate Views, Illumination Source

isRelative = True  # two different figures generated with bool

color, absdir = 'lightseagreen', [0,1,1]
if isRelative : 
    color = 'yellowgreen'
    absdir = [1,1,1]

rez=3
viewCoor = [ [30,-60], [30,-10], [30,30], [40,140] ]

#. Figure 1 & 2 - ...................................................
fig = plt.figure(figsize=(7.762,2.116))
minmax = (-1,1)
for i in range(len(viewCoor)) :
    ax = fig.add_subplot(141+i, projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    axcolor = 'black'
    if isRelative and i==2 : axcolor = 'firebrick'
    s3d.standardAxis(ax, length=1.7, width=2, color=axcolor, offset=1.0, negaxis=True)
    ax.text(0,0,-2,str(viewCoor[i]), color = 'black', horizontalalignment='center', verticalalignment='center')
    ax.set_axis_off()
    elev, azim = viewCoor[i]
    ax.view_init(elev,azim)
    illum = absdir
    if isRelative : illum = s3d.rtv(absdir,elev, azim)

    s = s3d.SphericalSurface(rez,facecolor=color)
    s.shade(direction=illum,contrast=0.7).hilite(direction=illum,focus=2)
    
    ax.add_collection3d(s)

fig.tight_layout()
plt.show()