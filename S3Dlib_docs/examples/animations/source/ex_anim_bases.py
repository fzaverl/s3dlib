import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Relative to the Viewer Animation

# 1. Define functions to examine ....................................

view_elev, start_time = 33, 0    #  0.0 <= time <= 1.0
illum_dir = [.5,1,.7]

# 2. Setup and map surfaces .........................................

cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')
low_color, hi_color = [0.25,0.13,0], [1,.9,.75]

surfdata = [ ['tetra', [0,-2,0] ], ['octa', [2,1,0] ],
             ['icosa', [-2,1,0] ], ['cube', [0,0,2] ], ['dodeca', [0,0,-2] ] ]

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
ax = fig.add_subplot(111, projection='3d')
minmax = (-1.8,1.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)

surface_set = []
for data in surfdata :
    surface = s3d.SphericalSurface(basetype=data[0],color=hi_color)
    surface.transform(translate=data[1])
    surface_set.append(surface)
    # text location dependent on object location
    ax.text(data[1][0], data[1][1], data[1][2] +1.3, data[0],
                horizontalalignment='center', verticalalignment='center')

# now, all object properties are dependent on azimuth and illumination direction.
azim = 360*start_time
illum = s3d.rtv(illum_dir, view_elev, azim)

surface_copies = []
for orig_surface in surface_set :
    surface = copy.copy(orig_surface)
    surface.map_cmap_from_normals( cmap='cardboard', direction=illum )
    ax.add_collection3d(surface)
    surface_copies.append(surface)

ax.set_axis_off()
ax.view_init(elev=view_elev, azim=azim)

plt.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return surface_copies,

def update_fig(frame):
    global surface_copies

    for surface in surface_copies :     
        ax.collections.remove(surface)

    azim = 360*frame
    illum = s3d.rtv(illum_dir, view_elev, azim)

    surface_copies = []
    for orig_surface in surface_set :
        surface = copy.copy(orig_surface)
        surface.map_cmap_from_normals( cmap='cardboard', direction=illum )
        ax.add_collection3d(surface)
        surface_copies.append(surface)

    ax.view_init( elev=view_elev, azim=azim )

    return surface_copies,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0.0, 1.0, 91),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
