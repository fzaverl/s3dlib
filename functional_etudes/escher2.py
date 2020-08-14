import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. #2 - influenced by M.C.Escher - Knots
#   https://mcescher.com/gallery/mathematical/

# 1. Define functions to examine ....................................

wdth = 0.9
twst, twstOff = 1.0 , 0.25
elev, azim = 90, -60
illum = s3d.rtv( [1,-1,1] ,elev,azim)

def flatten(rtz) :
    r,t,z = rtz
    R = 1-z*wdth/2
    return R,t,np.zeros(len(z))

def twistFunction(rtz, twists=twst, toff=twstOff) :
    r,t,z = rtz
    offset = toff*np.pi
    x0 = 1-r
    y0 = z
    r0, t0, temp = s3d.PolarSurface.coor_convert([x0,y0,np.zeros(len(z))])
    t0 = t0 - t*twists + offset
    x1, y1, temp = s3d.PolarSurface.coor_convert([r0,t0,np.zeros(len(z))],True)
    R = 1 - x1
    Z = y1
    return R,t,Z

def Trefoil(rtz) :
    r,t,z = rtz
    rw = 1-wdth/2
    X = rw*(np.sin(t)+2*np.sin(2*t))
    Y = rw*(np.cos(t)-2*np.cos(2*t))
    R0,T,Z = s3d.PolarSurface.coor_convert([X,Y,z])
    R = R0 + r - rw
    Z = z - np.sin(3*t)
    return R,T,Z

# 2. Setup and map surfaces .........................................
rez = 7
aa = colors.rgb_to_hsv( [ 0.659, 0.502, 0.400 ] )
ab = colors.rgb_to_hsv( [ 0.886, 0.706, 0.576 ] )
cmap = cmu.hsv_cmap_gradient(aa,ab)

surface1 = s3d.CylindricalSurface(rez,basetype='tri_s')
surface1.map_geom_from_op(flatten)
surface1.map_geom_from_op(twistFunction)
surface1.map_geom_from_op( Trefoil )

surface2 = s3d.CylindricalSurface(rez,basetype='tri_s')
surface2.transform(scale=[1,1,wdth/2])
surface2.map_geom_from_op(twistFunction)
surface2.map_geom_from_op( Trefoil )

cross = surface1+surface2

cross.map_cmap_from_normals(cmap,direction=illum)
cross.shade(.2,direction=illum).hilite(.8,direction=illum)

# 3. Construct figure, add surfaces, and plot ......................
fig = plt.figure(figsize=plt.figaspect(1))
info, infocolor = 'S3Dlib.org',  [0.898,0.843,0.800] 
text = fig.text(0.05, 0.07, info, color=infocolor, fontsize=45, fontweight='bold'  )
ax = plt.axes(projection='3d')
minmax = (-1.7,1.7)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.set_axis_off()
ax.set_proj_type('ortho')
ax.set_facecolor( [ 0.933, 0.902, 0.859 ] )
ax.view_init(elev,azim)

ax.add_collection3d(cross)

fig.tight_layout()
plt.show()
