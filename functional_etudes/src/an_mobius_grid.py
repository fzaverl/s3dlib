import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d

#.. influenced by M.C.Escher - Mobius Strip
#   https://mcescher.com/gallery/mathematical/

elev,azim = 0,0
light_direction = [1,.5,1 ]
illum = s3d.rtv(light_direction,elev,azim)

# 1. Define functions to examine ....................................

def sFunc(x,n,max=1,cycles=1) :
    # domain and range [0,max], linear for n=1
    def nFunc(x,n) :
        # normalized function 0-1
        a = 2*x-1
        sgna, absa = np.sign(a), np.abs(a)
        y = ( 1 + sgna*np.power(absa,n) )/2
        return y
    X = np.mod(x*cycles/max,1)
    Y = (max/cycles)*( nFunc(X,n) + np.floor(x*cycles/max))
    return Y

def twistFunction(rtz,gsize) :
    r,t,z = rtz
    twists = -1
    thickness = 0.25
    zoff = 0.3
    T = 2*t  # <-- expand the domain for two surface faces.
    w = thickness*z*(1 + 0.4*np.cos(2*T)) # <-- widen back.
    f = sFunc( T, 2, 4*np.pi, 2) # <----- non-linear twist.
    phi = 0.5*f*twists
    gap = gsize*np.ones(len(t))
    R = 1  + w*np.cos(phi) + gap*np.sin(phi)
    Z = w*np.sin(phi) - gap*np.cos(phi)
    Z = Z - zoff*np.sin(2*T) # <----------- warp flat ring.
    return R,T,Z

def sheargrid(rtz) :
    r,t,z = rtz
    shear = 0.06
    T = t - shear*z*np.sin(t)
    return r,T,z

def getAnts(offset=0) :
    #  .... ant coordinates.......
    t = np.linspace(0,2*np.pi,9,endpoint=False) + offset
    t = np.array( [t,t+.1,t+.15]).flatten()
    rtz = [np.ones(len(t)),t,np.zeros(len(t))]
    RTZ = twistFunction(rtz,0.15)
    X,Y,Z = s3d.CylindricalSurface.coor_convert(RTZ,tocart=True)
    ants=None
    for i in range(0,len(t)) :
        ant = s3d.SphericalSurface(2, color=antcolor)
        if i>=len(t)/3 : sc = 0.1
        else: sc = 0.15
        ant.transform(scale=sc,translate=[X[i],Z[i],-Y[i]]).hilite(.5).shade()
        if ants is None : ants = ant
        else : ants += ant
    return ants

# 2. Setup and map surfaces .........................................
ttlcolor,fgbgcolor  =  [0.502, 0.200, 0.278, 0.1] , [0.949,0.918,0.875]
antcolor = [0.675, 0.400, 0.298]  
rez = 8

twist = s3d.CylindricalSurface(rez, basetype='squ_s')
twist.map_color_from_image('images/grid_34_8_100.png')
twist.map_geom_from_op(sheargrid)
twist.map_geom_from_op( lambda rtz : twistFunction(rtz,0.01) )
twist.transform(s3d.eulerRot(0,90)).shade(direction=illum).hilite(.3)

surface = twist + getAnts(0.2)

# 3. Construct figures, add surface, plot ...........................
figratio = 82/50
fsize=3
fig = plt.figure(figsize=(fsize, fsize*figratio), facecolor=fgbgcolor)
text = fig.text(0.90, 0.85, 'S3Dlib.org', color=ttlcolor, ha='right',
    rotation=90, fontsize=55, fontweight='bold'  )
ax = plt.axes(projection='3d')
mm = 0.75
ynn = mm/figratio
minmax,yminmax = (-mm,mm), (-ynn,ynn)
ax.set_proj_type('ortho')
ax.set_facecolor(fgbgcolor)
ax.set(xlim=minmax, ylim=yminmax, zlim=minmax )
ax.view_init(elev,azim)
ax.set_axis_off()

ax.add_collection3d(surface)

fig.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return surface,

def update_fig(frame):
    global surface
  
    ax.collections.remove(surface)

    twistt = copy.copy(twist)
    surface = twistt + getAnts(frame)
    ax.add_collection3d(surface)

    return surface,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0.0, 2*np.pi/9, 51),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
