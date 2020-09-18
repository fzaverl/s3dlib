import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import s3dlib.surface as s3d

# data source =======================================================

with open('data/powerdata.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    all_data = []
    for row in csv_reader:
        line_count += 1
        if line_count ==1 :
            heading = row
        else:  all_data.append(row)    

data = np.array(all_data).T.astype(float)

# Figure 1 -XY plot of data  ========================================
'''
x = np.linspace(0, len(data[0])-1, len(data[0]))
plt.xlabel('Frame')
plt.ylabel('Power, normalized offset')
plt.plot(x, data[0], label=heading[0], color='r')
plt.plot(x, data[1], label=heading[1], color='g')
plt.plot(x, data[2], label=heading[2], color='b')
plt.legend()
plt.show()
'''
# Figure 2 -3D plot of data  ========================================

# 1. Define functions to examine ....................................

def data3D(xyz,frameID) :
    r,t,p = s3d.SphericalSurface.coor_convert(xyz,False)
    ax,ay,az = np.abs(xyz)
    emx = 1/np.sqrt(3)
    G = (r-emx)/(1-emx)
    F = (np.ones(len(r)).astype(float) - np.cos(np.pi*G))/2
    M = ax*data[0,frameID] + ay*data[1,frameID] + az*data[2,frameID]
    R = F*M
    XYZ = s3d.SphericalSurface.coor_convert([R,t,p], True)
    return XYZ

def Octaahedron(rez) :
    t = s3d.SphericalSurface(basetype='octa')
    surface = s3d.Surface3DCollection(t.vertexCoor, t.fvIndices)
    surface.triangulate(rez)
    surface.name = 'octahedron'
    return surface

def freq_value(fig,frameID) :
    symbol = str(frameID)+' '
    horz, vert = 0.90, 0.95
    fig.text(horz,vert,symbol, ha='left', va='center', fontsize='large',
        bbox=dict(boxstyle="square", ec=(1., 1., 1.),fc=(1., 1, 1), ))
    return
    
# 2. Setup and map surfaces .........................................
frameID = 25
rez = 5

surface = Octaahedron(rez)
surface.map_color_from_op( lambda xyz : np.abs(xyz) )
surface.map_geom_from_op(lambda rtp : data3D(rtp, frameID) )
surface.shade(.5)

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
minmax,ticks = (-0.6, 0.6), [0]
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_zticks(ticks)
ax.set_xlabel(heading[0])
ax.set_ylabel(heading[1])
ax.set_zlabel(heading[2])
fig.text(.90,.95,'Frame: ', ha='right', va='center', fontsize='x-large')
freq_value(fig,frameID)

ax.add_collection3d(surface)

fig.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return surface,

def update_fig(frame):
    global surface
    ax.collections.remove(surface)

    surface = Octaahedron(rez)
    surface.map_color_from_op( lambda xyz : np.abs(xyz) )
    surface.map_geom_from_op(lambda rtp : data3D(rtp, int(frame) ) )
    surface.shade(.5)
    freq_value(fig, int(frame) )

    ax.add_collection3d(surface)
    return surface,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0, len(data[0])-1, len(data[0]) ),
        init_func=init_fig, blit=False, repeat=True, interval=50)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
#====================================================================



