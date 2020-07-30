import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d

#.. RGB Color Space Animation

# 1. Define function to examine .....................................

def xyz2rgb(xyz) :
    # rgb color cube from -1 to 1
    x,y,z = xyz
    r = ( x + 1 ) /2.0
    g = ( y + 1 ) /2.0
    b = ( z + 1 ) /2.0  
    return r,g,b

def get_surface(Z1Z2) :
    z1, z2 = Z1Z2
    # z1 = f(1,1)    z2 = f(-1,-1)
    def topface(xyz, Z1, Z2) :
        x,y,z = xyz
        Z = ( Z2+Z1)/2.0 - (Z2-Z1)*(x+y)/4.0
        return x,y,Z
    def sideface(xyz, Z1, Z2, isX=True) :
        x,y,z = xyz
        Wmax = (Z2+Z1)/2
        Wmin = Z1
        X,Y = x,y
        if isX : 
            L,H,X = y, x, np.ones(len(x))
        else:
            L,H,Y = x,-y, np.ones(len(y))
        W = -(Wmax-Wmin)*L/2 + (Wmax+Wmin)/2
        Z = (W+1)*H/2 + (W-1)/2
        return X,Y,Z
    
    rez = 5
    ps = lambda : s3d.PlanarSurface(rez)
    top = ps().map_geom_from_op(lambda xyz : topface(xyz,z1,z2))
    front = ps().map_geom_from_op(lambda xyz : sideface(xyz,z1,z2))
    side = ps().map_geom_from_op(lambda xyz : sideface(xyz,z1,z2,False))
    surface = top+front+side
    surface.map_color_from_op( xyz2rgb )
    surface.shade(0.925, direction=[1,1,1])
    return surface

def get_frames(n) :
    negmin = -0.9999
    f2 = np.linspace(-1.0, 1.0, n)
    f1 = np.full(n,negmin)
    x = np.full(n,1.0)
    f2 = np.concatenate( ( f2, np.full(n,1.0) ) )
    f1 = np.concatenate( ( f1,np.linspace( negmin, 1.0, n) ) )
    f2 = np.concatenate( ( f2,np.linspace( 1.0, negmin, n) ) )
    f1 = np.concatenate( ( f1,np.linspace( 1.0, negmin, n) ) )
    return np.transpose([f1,f2])

# 2. Setup and map surfaces .........................................
framesize = 30
framearray = get_frames(framesize)

surface = get_surface( framearray[framesize] )

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(1),facecolor='black')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax.set_facecolor('black')
ax.set_proj_type('ortho')
ax.set_axis_off()
ax.view_init(azim=25)

ax.add_collection3d(surface)

# 4. Animation ......................................................

def init_fig():
    return surface,

def update_fig(frame):
    global surface
    global framearray

    ax.collections.remove(surface)
    surface = get_surface( framearray[frame] )
    ax.add_collection3d(surface)

    return surface,

ani = FuncAnimation(fig, update_fig, frames = range(3*framesize),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
