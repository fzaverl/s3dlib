import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Cylindrical Vibration Animation

# 1. Define functions to examine ....................................

n,m = 3,3
Um, Vm, Wm = 0.2, 0.1, 0.1
cyLen = 2.5

def displacements(rtz,time) :
    r,t,z = rtz
    eit = np.cos(2*np.pi*time)
    Z = np.pi*z/2
    u = Um*np.cos(n*t)*np.cos(m*Z)*eit
    v = Vm*np.sin(n*t)*np.sin(m*Z)*eit
    w = Wm*np.cos(n*t)*np.sin(m*Z)*eit
    return [u, v, w] 

def newCoor(rtz,time) :
    r,t,z = rtz
    u, v, w = displacements(rtz,time)
    R = r + u  
    T = t + v/r  # small angle displacements:  v ~ r*dt
    Z = z + w
    return R,T,Z

# 2. Setup and map surfaces .........................................
start_time = 0
rez = 5
cboard = cmu.rgb_cmap_gradient([0.25,0.15,0], [1,.9,.75], 'cardboard' )

cylinder = s3d.CylindricalSurface(rez,basetype='tri')
cylinder.map_geom_from_op(lambda rtz : newCoor(rtz, start_time ) )
cylinder.transform(s3d.eulerRot(0,40), scale=[1,1,cyLen])
cylinder.map_cmap_from_normals( 'cardboard' ).shade(0.8)

# 3. Construct figure, add surfaces, and plot .....................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
mnmx = [-1.8,1.8]
ax.set(xlim=mnmx, ylim=mnmx, zlim=mnmx )
ax.view_init(elev=10)
ax.set_axis_off()

ax.add_collection3d(cylinder)

fig.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return cylinder,

def update_fig(frame):
    global cylinder
    ax.collections.remove(cylinder)

    cylinder = s3d.CylindricalSurface(rez,basetype='tri')
    cylinder.map_geom_from_op(lambda rtz : newCoor(rtz,frame ) )
    cylinder.transform(s3d.eulerRot(0,40), scale=[1,1,cyLen])
    cylinder.map_cmap_from_normals( 'cardboard' ).shade(0.8)

    ax.add_collection3d(cylinder)
    return cylinder,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0, 1, 40),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
