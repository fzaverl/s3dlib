import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import s3dlib.surface as s3d

#.. Guides: Orientation - Object Rotations
#   surface object rotations

fourColor_cmap = ListedColormap(['firebrick','forestgreen','yellow','mediumblue'])
rez = 6
illum = [0,1,1]

# 1. Define functions to examine ....................................

def deflate(rtp) :
    r,t,p = rtp
    scale = 0.2
    Rz = np.cos(p)
    Rxys = (1-scale)*np.sin(p) + scale*np.cos(4*t)
    R = np.sqrt( Rz**2 + Rxys**2)
    return R,t,p

# Setup and mapsurfaces .............................................

surface = s3d.SphericalSurface(rez,basetype='octa')
surface.map_cmap_from_op(lambda rtp : rtp[1], fourColor_cmap)
surface.map_geom_from_op(deflate)

minmax, width = (-1,1), 5
# Figure 1, 2 unrotated surfaces.....................................

fig = plt.figure(figsize=(width, width/2))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax2.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax1.set_title('default view' +'\n', fontsize='small')
ax2.set_title('standardAxis' +'\n', fontsize='small')
surf1 = copy.copy(surface).shade(.3).hilite(.7)
surf2 = copy.copy(surface).shade(.3,illum).hilite(.7,illum)
s3d.setupAxis( ax1, offset=1.0 )
ax1.set_axis_off()
s3d.standardAxis( ax2, offset=1.0 )

ax1.add_collection3d(surf1)
ax2.add_collection3d(surf2)

fig.tight_layout()

# Figure 2, 4 rotated surfaces.......................................
theta, phi, psi_v = 60, 45, 50

cnt = [ [True,0], [False,0], [True, psi_v], [False, psi_v]  ]
fig = plt.figure(figsize=(width,width))

for i in range(4) :
    surf = copy.copy(surface)
    isX, psi = cnt[i]
    surf.transform(s3d.eulerRot(theta,phi,psi,useXconv=isX))
    surf.shade(.3,illum).hilite(.7,illum)
    title = 'Y-convension\n'
    if isX : title = 'X-convension\n'
    title2 = r'$\theta$ = {}, $\phi$ = {}, $\psi$ = {}'.format(theta,phi,psi)
    ax = fig.add_subplot(221+i, projection='3d')
    ax.set_title(title+title2, fontsize='small')
    s3d.standardAxis( ax, offset=1.0 )

    ax.add_collection3d(surf)

fig.tight_layout()

#....................................................................
plt.show()