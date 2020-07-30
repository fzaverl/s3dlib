import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Percentile Visualization, 2

# 1. Define data to examine .........................................

data = np.loadtxt(open("data/example_stats.csv", "rb"), delimiter=",", skiprows=1)
xa,ya,za = np.transpose(data)
N = len(xa)

# 2. Setup and map surfaces .........................................
colormap=cmu.hsv_cmap_gradient( [0.333,1,.65], [0,1,1], smooth=1.6 )
grey = [.5,.5,.5]
prct = 0.5

surface = s3d.SphericalSurface(3, color=grey)
surface.set_surface_alpha(.05).shade()
disArr_a,t = surface.svd(data,prct)
info = str(N) +', '+"{:.0%}".format(prct) + ', ' + '{:04.2f}'.format( t[0] )

maxdis = max(disArr_a)
colors_a = []
for val in disArr_a :
    if val < 1 : colors_a.append(grey)
    else : colors_a.append(colormap( (val-1)/(maxdis-1) ) )

# 3. Construct figures, add dasurfaces, and plot ....................
minmax,ticks  = (0,20), [0,5,10,15,20]
fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_xlabel('Define Function')
ax.set_ylabel('Setup Surface')
ax.set_zlabel('Construct Figure')
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_zticks(ticks)
ax.set_title(info, horizontalalignment='left')
ax.scatter(xa,ya,za, c=colors_a, edgecolor='k')

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()
