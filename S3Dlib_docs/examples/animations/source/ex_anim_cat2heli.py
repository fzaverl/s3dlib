import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Helicoid Transformation Animation

# 1. Define functions to examine ....................................

def catenoid_helicoid(rtz, A) :
    r,t,z = rtz
    A = A*np.pi  #  -1 < A < 1
    cosA, sinA = np.cos(A), np.sin(A)
    U, V = t, z   
    x =  cosA * np.sinh(V) * np.sin(U) +   sinA * np.cosh(V) * np.cos(U)
    y = -cosA * np.sinh(V) * np.cos(U) +   sinA * np.cosh(V) * np.sin(U)
    Z = ( U/np.pi- 1.0 ) *cosA +  V * sinA
    return x,y,Z

def colormap_by_A(A) :
    hue = (A + 1.0) / 2.0
    return cmu.hsv_cmap_gradient( [hue,1.0,0.25], [hue,0.5,1] )

def indicator_by_A(fig, A, vOld=None) :
    symbol, blank = r'$\blacktriangleright$', r'$\blacksquare$'
    horz, vCen, vRng = 0.8, 0.5, 0.28
    vert = vCen + vRng*A
    if vOld is not None: 
        fig.text(horz,vOld,blank, ha='right', va='center', fontsize='x-large', color='w')
    fig.text(horz,vert,symbol, ha='right', va='center', fontsize='large')
    return vert

# 2. Setup and map surfaces .........................................
start_A = -1
rez = 4

surface = s3d.CylindricalSurface(rez, basetype='squ_s', cmap='hsv', antialiased=True)
surface.map_geom_from_op( lambda rtz : catenoid_helicoid(rtz,start_A), returnxyz=True )

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
cbar = plt.colorbar(surface.cBar_ScalarMappable, ax=ax, ticks=np.linspace(-1,1,5), shrink=0.6 )
cbar.set_label(r'A Parameter, $\pi$  ', rotation=270, labelpad = 15)
cbar.ax.tick_params(labelsize='small')
# surface mapping placed here to change colormap after colorbar defined.
surface.map_cmap_from_op( lambda rtz : rtz[0],colormap_by_A(start_A))
# ....
minmax = (-1.1,1.1)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.view_init(25, -28)
ax.set_axis_off()
prevIndicator = indicator_by_A(fig, start_A)

ax.add_collection3d(surface)

fig.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return surface,

def update_fig(frame):
    global surface
    global prevIndicator
    ax.collections.remove(surface)

    surface = s3d.CylindricalSurface(rez, basetype='squ_s', antialiased=True)
    surface.map_geom_from_op( lambda rtz : catenoid_helicoid(rtz,frame), returnxyz=True )
    surface.map_cmap_from_op( lambda rtz : rtz[0],colormap_by_A(frame))
    prevIndicator = indicator_by_A(fig, frame, prevIndicator)

    ax.add_collection3d(surface)

    return surface,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(-1.0, 1.0, 81),
                    init_func=init_fig, blit=False, repeat=True, interval=50)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
#ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
#print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
plt.show()
