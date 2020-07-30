import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d

#.. HSV Color Space Animation

# 1. Define function to examine .....................................

def xyz2hsv(xyz,rot) :
    # hsv color cylinder from -1<z<1, 0<R<1
    x,y,z = xyz
    theta = np.arctan2(y,x)
    theta = np.where(theta<0,theta+2*np.pi,theta)
    h = np.mod( rot + theta/(2*np.pi), 1.0 )
    s = np.sqrt( x**2 + y**2)
    v = ( z + 1 ) /2.0  
    return h,s,v

def get_surface(Z1TSTE_a) :
    z1, ts, te, ag = Z1TSTE_a
    # ts, te - theta start, end
    def outsideface(rtz, Zm, Ts, Te, isTop=True) :
        r,t,z = rtz
        T = (Te-Ts)*t + 2*np.pi*Ts
        T = np.mod(T,2*np.pi)
        if isTop : Z = np.full(len(z),Zm)
        else:      Z = (Zm+1)*z/2 + (Zm-1)/2
        return r,T,Z
    def insideface(xyz, Zm, T, back=1) :
        x,y,z = xyz
        A = ( x + 1.0 )/2.0
        Z = back*(Zm+1)*y/2 + (Zm-1)/2
        X = A*np.cos(T*2*np.pi) 
        Y = A*np.sin(T*2*np.pi) 
        return X,Y,Z

    rez=5
    front = s3d.PlanarSurface(rez,color='C0')
    front.map_geom_from_op(lambda xyz : insideface(xyz,z1,ts))
    back = s3d.PlanarSurface(rez,color='C1')
    back.map_geom_from_op(lambda xyz : insideface(xyz,z1,te,-1))
    top = s3d.PolarSurface(rez,basetype='hex_c',color='C2')
    top.map_geom_from_op(lambda rtz : outsideface(rtz,z1,ts,te))
    side = s3d.CylindricalSurface(rez,basetype='tri_s',color='C3')
    side.map_geom_from_op(lambda rtz : outsideface(rtz,z1,ts,te,False))
    surface = front+back+top+side
    surface.map_color_from_op(lambda xyz : xyz2hsv(xyz,ag), rgb=False )
    surface.shade(0.925)
    return surface

def get_frames(n) :
    negmin = -0.9999
    ts_0 = 1/12
    te_0 = 1 + ts_0
    ts_f,te_f = 1/3, 5/6
  
    z = np.linspace(negmin, 1.0, n)
    z = np.concatenate( ( z, np.full(2*n,1.0) ) )
    z = np.concatenate( ( z, np.linspace( 1.0, negmin, n) ) )
    
    ts = np.full(n,ts_0)
    ts = np.concatenate( ( ts, np.linspace( ts_0, ts_f, n) ) )
    ts = np.concatenate( ( ts, np.full(n,ts_f) ) )
    ts = np.concatenate( ( ts, np.linspace( ts_f, ts_0, n) ) )

    te = np.full(n,te_0)
    te = np.concatenate( ( te, np.linspace( te_0, te_f, n) ) )
    te = np.concatenate( ( te, np.full(n,te_f) ) )
    te = np.concatenate( ( te, np.linspace( te_f, te_0, n) ) )

    a =  np.zeros(2*n)
    a =  np.concatenate( ( a,  np.linspace( 0,   1.0, 2*n) ) )

    return np.transpose([z,ts,te,a])

# 2. Setup and map surfaces .........................................
framesize = 30
framearray = get_frames(framesize)

surface = get_surface(framearray[0])

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(1),facecolor='black')
fig.text(0.975,0.975,str(surface), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax.set_axis_off()
ax.set_proj_type('ortho')
ax.view_init(azim=25)
ax.set_facecolor('black')

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

ani = FuncAnimation(fig, update_fig, frames = range(4*framesize),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
