import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu  

# 2. Setup and map surfaces .........................................
rez = 1
C,F,E = 'C2', 'gold', 'C3'

Fa = colors.to_rgba('gold',0.2)
Ea = colors.to_rgba('C3', 0.00)

surf_Fa_s =    s3d.SphericalSurface( rez, color=Fa ).shade()
surf_FaEa_a =  s3d.SphericalSurface( rez, facecolor=Fa, edgecolor=Ea ).shade()
surf_FaW0_a =  s3d.SphericalSurface( rez, color=Fa, linewidth=0 ).shade()

#................
surf_Fs =  s3d.SphericalSurface( rez, color=F ).shade()
surf_Fs.set_facecolor([0,0,0,0])
surf_cmap =  s3d.SphericalSurface( rez, color=F )
surf_cmap.map_cmap_from_normals('hsv')
surf_cmap.set_facecolor([0,0,0,0])

# 3. Construct figure, add surfaces, and plot .....................
minmax = (-.7,.7)

title = [ 'shade', 'hsv normals' ]
surf =  [ surf_Fs, surf_cmap ]

fig = plt.figure(figsize=(4,2))

for i in range(2) :
    ax = fig.add_subplot(121+i, projection='3d')
    ax.set_axis_off()
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
    ax.set_title(title[i])
    ax.add_collection3d(surf[i])

fig.tight_layout()

#....................
title = [ 'F.alpha',  'FE.alpha', 'F.alpha w=0']
surf =  [ surf_Fa_s, surf_FaEa_a, surf_FaW0_a ]

fig = plt.figure(figsize=(5,5/2.8))

for i in range(3) :
    ax = fig.add_subplot(131+i, projection='3d')
    ax.set_axis_off()
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
    ax.set_title(title[i])
    ax.add_collection3d(surf[i])

fig.tight_layout()

plt.show()