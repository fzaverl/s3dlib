import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Base Face Variations

# 2. Setup and map surface .........................................
rez = 2

surface = s3d.SphericalSurface(rez)
surface.map_cmap_from_op( lambda rtz: rtz[0] , cmap='magma' ).shade(.2)

# 3. Construct figure, add surface plot ............................

fig = plt.figure(figsize=plt.figaspect(0.8))
ax = plt.axes(projection='3d')
minmax = (-0.8,0.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
minc, maxc = surface.bounds['vlim']
plt.colorbar(surface.cBar_ScalarMappable, ax=ax, ticks=np.linspace(minc,maxc,4), shrink=0.6 )
ax.set_axis_off()

ax.add_collection3d(surface)

plt.tight_layout()
plt.show()