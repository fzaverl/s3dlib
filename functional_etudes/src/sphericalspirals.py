import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. #1 - influenced by M.C.Escher - Sphere Spirals
#   https://mcescher.com/gallery/mathematical/

# 1. Define function to examine .....................................
ttlcolor =  [0.502, 0.200, 0.278, 0.1] 
fgbgcolor = [0.957, 0.925, 0.882]
outercolor = [0.875, 0.682, 0.388]
innercolor = [0.847, 0.525, 0.471]

elev, azim = 45,30
light_direction = [1,-1,1 ]
illum = s3d.rtv(light_direction,elev,azim)
negillum = -1.0*np.array(illum)

def cliptwo(rtp) :
    r,t,p = rtp
    val = np.full(len(r),False)
    val1 = np.logical_and(t>0.0,t<np.pi/2)
    val2 = np.logical_and(t>np.pi,t<3*np.pi/2)
    val = np.logical_or(val1,val2)
    return val

def opensphere(rtp, rtt = 0.0) :
    r,t,p = rtp
    T = t/2.0 + rtt   
    return r,T,p

def twist(rtp, twist) :
    r,t,p = rtp
    f = np.cos(t) + 1
    b =  2*p/np.pi -1
    g = ( 1 + np.sign(b)*np.power(  np.abs(b), 6 ) )/2
    #T = t - 2*twist*p  .... twist rate is not llinear, so...
    T = t - 2*twists*(p+g*2*np.pi)
    return r,T,p

# 2. Setup and map surfaces .........................................
rez = 7
twists = 1

surface1 = s3d.SphericalSurface(rez, basetype='octa', color=outercolor)
surface1.clip( cliptwo )
surface1.map_geom_from_op(lambda rtp : opensphere(rtp))
surface1.map_geom_from_op(lambda rtp : twist(rtp,twists))

surface2 = s3d.SphericalSurface(rez, basetype='octa', color=outercolor)
surface2.clip( cliptwo )
surface2.map_geom_from_op(lambda rtp : opensphere(rtp, np.pi))
surface2.map_geom_from_op(lambda rtp : twist(rtp,twists))

surfaceA = surface1 + surface2
surfaceA.shade(direction=illum).hilite(.5,direction=illum)

surface1 = s3d.SphericalSurface(rez, basetype='octa', color=innercolor)
surface1.clip( cliptwo )
surface1.map_geom_from_op(lambda rtp : opensphere(rtp))
surface1.map_geom_from_op(lambda rtp : twist(rtp,twists))

surface2 = s3d.SphericalSurface(rez, basetype='octa', color=innercolor)
surface2.clip( cliptwo )
surface2.map_geom_from_op(lambda rtp : opensphere(rtp, np.pi))
surface2.map_geom_from_op(lambda rtp : twist(rtp,twists))

surfaceB = surface1 + surface2
surfaceB.transform(scale=0.95)
surfaceB.shade(direction=negillum).hilite(.5,direction=negillum)

surface = surfaceA + surfaceB

grid = s3d.SphericalSurface(4, basetype='icosa', color='k')
grid.map_cmap_from_normals('binary',illum)
grid.shade(direction=negillum).hilite(direction=negillum)
grid.set_surface_alpha(0.05) 
grid.transform(scale=1.01)

surface = surface + grid

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(1), facecolor=fgbgcolor)
text = fig.text(0.97, 0.87, 'S3Dlib.org', color=ttlcolor, ha='right',
    rotation=90, fontsize=55, fontweight='bold'  )
ax = plt.axes(projection='3d')
maxmin = ( -0.7,0.7 )
ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
ax.set_axis_off()
ax.view_init(elev, azim)
ax.set_facecolor(fgbgcolor)

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()