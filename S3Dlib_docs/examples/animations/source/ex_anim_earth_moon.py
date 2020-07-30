import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d

#.. Earth-Moon Animation

# 1. Define functions to examine ....................................
# 2. Setup and map surfaces .........................................
start_theta = 0
rez=6
lightDirection = [1,0.8,1]

moon = s3d.SphericalSurface(rez-1)
moon.map_color_from_image('data/moon.png')
moon.shade(contrast=2,direction=lightDirection)
moon.transform(scale=[.4,.4,.4],translate=[-2,1,0.5])
# .......
earth = s3d.SphericalSurface(rez)
earth.map_color_from_image('data/earth.png')
earth.map_geom_from_image('data/elevation.png',0.06)
earthcopy = copy.copy(earth)
earthcopy.transform(scale=[1.2,1.2,1.2], translate=[.5,-.25,0],rotate=s3d.eulerRot(start_theta,0) )
earthcopy.shade(contrast=1.7,direction=lightDirection)
 
# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1), facecolor=(0,0,0) )
ax = fig.add_subplot(111, projection='3d')
minmax = (-1,1)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_facecolor('black')
ax.set_axis_off()

ax.add_collection3d(moon)
ax.add_collection3d(earthcopy)

# 4. Animation ......................................................

def init_fig():
    return earthcopy,

def update_fig(frame):
    global earthcopy
    ax.collections.remove(earthcopy)

    earthcopy = copy.copy(earth)
    earthcopy.transform(scale=[1.2,1.2,1.2], translate=[.5,-.25,0],rotate=s3d.eulerRot(frame,0) )
    earthcopy.shade(contrast=1.7,direction=lightDirection)

    ax.add_collection3d(earthcopy)
    return earthcopy,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0, 360, 60),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
