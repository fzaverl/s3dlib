import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import s3dlib.surface as s3d

#.. Principal components analysis (PCA)

# 1. Define data to examine .........................................

np.random.seed(0)
N=3000
y = np.random.normal(scale=0.5, size=N)
x = np.random.normal(scale=0.5, size=N)
z = np.random.normal(scale=0.1, size=N)
a = x + y
b = 2 * y
c = a - b + z
norm = np.sqrt(a.var() + b.var())
a /= norm
b /= norm
data = np.transpose([ a,b,c ])

# 2. Setup and map surfaces .........................................

ellipsoid = s3d.SphericalSurface(3, color='darkgoldenrod', linewidth=0  )
plate     = s3d.PlanarSurface(3, color='darkgoldenrod', linewidth=0  )

# 3. Construct figures, add surfaces, and plot ......................
surfaces = [ plate, ellipsoid ]
elevazim = [ (-75,-80), (45,15) ]

fig = plt.figure(figsize=plt.figaspect(.5))
minmax=(-3,3)
for i in range(2) :
    surface = surfaces[i]
    ea = elevazim[i]
    # setup surfaces .......
    disArr_a,t = surface.svd(data)
    tAxis = surface.get_transformAxis([2.3,2.5,-20])
    surface.transform(scale=2).set_surface_alpha(.2)
    # .....................
    ax = fig.add_subplot(121+i, projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.scatter(a,b,c, c='k', marker='.', s=1)
    ax.view_init(ea[0],ea[1])
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.axes.zaxis.set_ticklabels([])
    ax.set_xlabel('X', labelpad=-10)
    ax.set_ylabel('Y', labelpad=-10)
    ax.set_zlabel('Z', labelpad=-10)
    ax.add_collection(surface)
    ax.add_collection3d(tAxis)

fig.tight_layout()

plt.show()

