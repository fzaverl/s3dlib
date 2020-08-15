import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. #3 - influenced by M.C.Escher - Knots
#   https://mcescher.com/gallery/mathematical/

# 1. Define functions to examine ....................................

wdth = 1
elev, azim = 90, 0
illum = s3d.rtv([1,-1,1],elev,azim)

def Torus(rez,width=wdth ) :
    def fold(rtz,width) :
        r,t,z = rtz
        Z = width*np.sin(z*np.pi)/2
        R = 1 + width*np.cos(z*np.pi)/2
        return R,t,Z
    surface = s3d.CylindricalSurface(rez,basetype='tri_s')
    surface.map_geom_from_op( lambda rtz : fold(rtz,width)    )
    return surface

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
rez = 6
aa = colors.rgb_to_hsv( [ 0.659, 0.502, 0.400 ] )
ab = colors.rgb_to_hsv( [ 0.886, 0.706, 0.576 ] )
cmap = cmu.hsv_cmap_gradient(aa,ab)

surface1 = Torus(rez)
surface1.map_geom_from_op( Trefoil )
surface1.map_cmap_from_normals(cmap,direction=illum)
surface1.set_surface_alpha(.3)

surface2 = Torus(rez,wdth/3)
surface2.map_geom_from_op( Trefoil )
surface2.map_cmap_from_normals(cmap,direction=illum)

ring = surface1+surface2
ring.set_linewidth(0)
ring.shade(.2,direction=illum).hilite(.8,direction=illum)

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

ax.add_collection3d(ring)

fig.tight_layout()
plt.show()
