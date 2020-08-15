import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. #4 - annimated, influenced by M.C.Escher - Knots
#   https://mcescher.com/gallery/mathematical/

# 1. Define functions to examine ....................................

wdth, subW = 0.4, 0.4

elev, azim = 90, -30
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

def shift(rtz, os=0) :
    r,t,z = rtz
    T = np.mod(t+os,np.full(len(t),2*np.pi))
    return T

# 2. Setup and map surfaces .........................................
rez = 6
white2cyan = cmu.hsv_cmap_gradient( [0.5,0,1], [.5,1,1] )
cyan2clear = cmu.hsv_cmap_gradient( [.5,1,1,1], [.5,1,1,0] )
clear = cmu.hsv_cmap_gradient( [.5,1,1,0], [.5,1,1,0] )
cmap1 = cmu.stitch_cmap(clear,cyan2clear, clear,bndry=[0.005,0.2] )
infRev = cmu.reversed_cmap('inferno')
cmap2 = cmu.stitch_cmap(infRev, clear,bndry=[0.5] )

offset = 0

surface1 = Torus(rez)
surf1 = copy.copy(surface1)
surface1.map_cmap_from_op( lambda rtz : shift(rtz,offset), cmap1)
surface1.map_geom_from_op( Trefoil )

surface2 = Torus(rez, subW*wdth)
surf2 = copy.copy(surface2)
surface2.map_cmap_from_op( lambda rtz : shift(rtz,offset), cmap2)
surface2.map_geom_from_op( Trefoil )

ring = surface1+surface2
ring.set_linewidth(0)
ring.shade(.2,direction=illum).hilite(.8,direction=illum)

# 3. Construct figure, add surfaces, and plot ......................
fig = plt.figure(figsize=plt.figaspect(1))
info, infocolor = 'S3Dlib.org',  [0.1,0.15,0.15] 
text = fig.text(0.05, 0.07, info, color=infocolor, fontsize=45, fontweight='bold'  )
ax = plt.axes(projection='3d')
minmax = (-2,2)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.set_axis_off()
ax.set_proj_type('ortho')
ax.set_facecolor( [ 0, 0, 0 ] )
ax.view_init(elev,azim)

ax.add_collection3d(ring)

fig.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return ring,

def update_fig(frame):
    global ring
  
    ax.collections.remove(ring)

    offset = frame*(2*np.pi)

    surface1 = copy.copy(surf1)
    surface1.map_cmap_from_op( lambda rtz : shift(rtz,offset), cmap1)
    surface1.map_geom_from_op( Trefoil )

    surface2 = copy.copy(surf2)
    surface2.map_cmap_from_op( lambda rtz : shift(rtz,offset), cmap2)
    surface2.map_geom_from_op( Trefoil )

    ring = surface1+surface2
    ring.set_linewidth(0)
    ring.shade(.2,direction=illum).hilite(.8,direction=illum)

    ax.add_collection3d(ring)

    return ring,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0.0, 1.0, 91),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")

