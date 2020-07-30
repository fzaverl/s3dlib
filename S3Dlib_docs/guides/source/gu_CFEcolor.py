import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
 
# 2. Setup and map surfaces .........................................
rez = 1
C,F,E = 'C2', 'gold', 'C3'

surf_C =    s3d.SphericalSurface( rez, color=C )
surf_F =    s3d.SphericalSurface( rez, facecolor=F )
surf_E =    s3d.SphericalSurface( rez, edgecolor=E )
surf_CF =   s3d.SphericalSurface( rez, color=C, facecolor=F )
surf_CE =   s3d.SphericalSurface( rez, color=C, edgecolor=E )
surf_FE =   s3d.SphericalSurface( rez, facecolor=F, edgecolor=E )
surf_0 =    s3d.SphericalSurface( rez )
surf_CFE =  s3d.SphericalSurface( rez, color=C, facecolor=F, edgecolor=E )
surf_CFEs = s3d.SphericalSurface( rez, color=C, facecolor=F, edgecolor=E ).shade()

# 3. Construct figure, add surfaces, and plot .....................

title = [    'C',    'F',    'E',    'CF',    'CE',    'FE', 'default',    'CFE',  'CFE.shade']
surf =  [ surf_C, surf_F, surf_E, surf_CF, surf_CE, surf_FE,    surf_0, surf_CFE,   surf_CFEs ]

fig = plt.figure(figsize=(5,5))
minmax = (-.7,.7)

for i in range(9) :
    ax = fig.add_subplot(331+i, projection='3d')
    ax.set_axis_off()
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
    ax.set_title(title[i])
    ax.add_collection3d(surf[i])

fig.tight_layout()
plt.show()