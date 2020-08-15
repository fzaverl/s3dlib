import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. #1 - influenced by M.C.Escher - Knots
#   https://mcescher.com/gallery/mathematical/

# 1. Define functions to examine ....................................

wdth = 0.75
twst, twstOff = 0.75, 0.25  # MC Escher controls
elev, azim = 90, -30
illum = s3d.rtv([1,-1,1],elev,azim)

def SquareRing(rez, width=wdth) :
    # .....................................................
    def fold(rtz,width,height) :
        r,t,z = rtz
        zeros = np.zeros(len(z))
        width_ar = np.full(len(z),width)
        # fold the cylinder into 4 parts..
        alpha = -2*width*z+width
        alpha = np.where( z <= 0.5, zeros ,     alpha )
        alpha = np.where( z <= 0.0, 2*width*z , alpha )
        alpha = np.where( z <= -.5, -width_ar , alpha )
        beta = height
        beta = np.where( z <= 0.5, 2*height*z,         beta)
        beta = np.where( z <= 0.0, zeros,              beta)
        beta = np.where( z <= -.5, -2*height*z-height, beta)
        R = r + alpha
        R = R + width/2
        Z = beta - height/2
        return R,t,Z
    # .....................................................
    surface = s3d.CylindricalSurface(rez,basetype='tri_s')   
    surface.map_geom_from_op( lambda rtz : fold(rtz,width,width) )
    surface.name = 'ring'
    return surface

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
ba = colors.rgb_to_hsv( [ 0.482, 0.333, 0.267 ] )
bb = colors.rgb_to_hsv( [ 0.855, 0.584, 0.427 ] )
cmap = cmu.hsv_cmap_gradient(ba,bb)

ring = SquareRing(rez)
ring.map_geom_from_op(twistFunction)
ring.map_geom_from_op( Trefoil )
ring.map_cmap_from_normals(cmap,direction=illum)
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
