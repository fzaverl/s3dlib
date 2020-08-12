import numpy as np
from scipy import special as sp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import s3dlib.surface as s3d

#.. S3Dlib logo construction, figure & background

#. Figure 1 =======================================================

# 1. Define function to examine ...................................

def showquad(xyz) :
    x,y,z = xyz
    return x*y*z > 0

def sphHar(rtp) :
    r, theta, phi = rtp
    m, l = 2,3
    r = sp.sph_harm(m, l, theta, phi).imag
    return r, theta, phi

def sphHar_absR(rtp) :
    r, theta, phi = sphHar(rtp)
    return np.abs(r), theta, phi

# 2. Setup surface ................................................
rez = 6
elev, azim, reldir = 15, -75, [1,1,1]
direction = s3d.rtv(reldir,elev,azim)

slm = [ 0.769, 1.000, 0.000 ]
sor = [ 1.000, 0.404, 0.000 ]
scy = [ 0.000, 1.000, 0.769 ]
sgr = [ 0.027, 0.851, 0.000 ]
cmap = ListedColormap([sor,scy,sgr,slm])

surface = s3d.SphericalSurface(rez, basetype='octa', cmap=cmap)
surface.clip(showquad,True)
surface.map_cmap_from_op(lambda rtp : rtp[1])
surface.map_geom_from_op(sphHar_absR)
surface.shade(.2).hilite(.7, direction=direction)

# 3. Construct figure, add surface, plot ..........................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
maxR = surface.bounds['xlim'][1]
minmax=(-maxR,maxR)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.view_init(elev, azim)

ax.add_collection3d(surface)

fig.tight_layout()

#. Figure 2 =======================================================

# 2. Setup and map surface ........................................

surface = s3d.PolarSurface(4, basetype='hex_c' )
edges = surface.edges
edges.set_color('k')
edges.set_linewidth(.5)

# 3. Construct figure, add surface, and plot ......................

fig = plt.figure(figsize=(8,8))
ax = plt.axes(projection='3d')
minmax=(-.6,.6)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.view_init(90,-90)
ax.set_axis_off()

ax.add_collection3d(edges)

fig.tight_layout()

# =================================================================

plt.show()