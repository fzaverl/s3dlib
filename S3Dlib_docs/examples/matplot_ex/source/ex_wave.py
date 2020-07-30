import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Matplotlib Examples: Shading

# 1. Define function to examine .....................................

def wavefunc(xyz) :
    x,y,z = xyz
    r = np.sqrt( x**2 + y**2)
    Z = np.sin( 6.0*r )/2
    return x,y,Z

# 2. Setup and map surfaces .........................................
rez = 5
fc = [1,.88,.72]
cmap = cmu.rgb_cmap_gradient( 'black', fc )

wave = s3d.PlanarSurface(rez, facecolor=fc , cmap=cmap )
wave.map_geom_from_op( wavefunc ).shade(direction=[1,1,1])
#wave.map_cmap_from_op( lambda xyz : xyz[2] , 'coolwarm')

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(0.75))
fig.text(0.975,0.975,str(wave), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax.xaxis.set_major_locator(LinearLocator(5))
ax.yaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_locator(LinearLocator(5))
plt.colorbar(wave.cBar_ScalarMappable, ax=ax,  shrink=0.6 )

ax.add_collection3d(wave)

plt.show()