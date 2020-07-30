import warnings
import copy

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from colorspacious import cspace_converter
import s3dlib.surface as s3d

#.. Lab Space Animation

# 1. Define function to examine ...............

elev, illum = 30, [1,1,1]

labRange = 50

def Cube() :  
    v = [ 
        [  0, 0, 0 ],  [  0, 1, 0 ],  [ 1 , 1, 0 ],  [  1, 0, 0 ],
        [  0, 0, 1 ],  [  0, 1, 1 ],  [ 1 , 1, 1 ],  [  1, 0, 1 ]  ]
    f = [ [0,1,2,3], [3,2,6,7], [2,1,5,6], [1,0,4,5], [0,3,7,4], [4,7,6,5] ]
    vertexCoor = np.array(v).astype(float)
    faceIndices = np.array(f)
    surface = s3d.Surface3DCollection(vertexCoor, faceIndices)
    surface.triangulate(0)
    surface.set_color([1,1,1,.1])
    surface.transform(translate=[-.5,-.5,00])
    surface.transform(scale=2.2*labRange)
    return surface

def dirBound(direction) :
    # simple binary-chop to determime limits in a Lab direction.
    def rgbErr(lab) :
        # determine out of gamut error.
        rgb = cspace_converter("CAM02-UCS", "sRGB1" )(lab)
        rgbRng = np.clip(rgb,0,1)
        errMagn = np.linalg.norm(rgb-rgbRng)
        return errMagn
    #........................................
    def toLab(normal, n) :
        # centered at [50,0,0]
        Lmax, abMax = 50, labRange
        L = Lmax*( n*normal[0] + 1.0 )
        a = n*abMax*normal[1]
        b = n*abMax*normal[2]
        return [L,a,b]
    #........................................
    normal = direction/np.linalg.norm(direction)
    maxdir = np.amax(np.abs(normal))
    nDir = normal/maxdir
    nLw, nHi = 0, 1
    maxErr, maxRGBerr, maxChops = 0.001, 0.004, 18
    inGamut = toLab(nDir,nLw)
    lastRGBgoodValue = cspace_converter("CAM02-UCS", "sRGB1" )(inGamut)
    delta_rgbErr = 0
    for i in range(maxChops) :
        nMd = (nHi+nLw)/2
        tempGamut = toLab(nDir,nMd)
        errMd = rgbErr( toLab(nDir,nMd))
        if errMd < maxErr :
            nLw = nMd
            inGamut = tempGamut
            oldRGBgoodValue = lastRGBgoodValue
            lastRGBgoodValue = cspace_converter("CAM02-UCS", "sRGB1" )(inGamut)
            delta_rgbErr = np.linalg.norm(lastRGBgoodValue-oldRGBgoodValue)
            if delta_rgbErr < maxRGBerr :  break
        else :
            nHi = nMd
        if i is maxChops-1 :
            chopErr = delta_rgbErr*255/np.sqrt(3)
            if chopErr > 2 :
                warnings.warn('BC algorithm maxed out!! {:6.2f}'.format(chopErr ))
    return inGamut

def toCoor(xyz) :
    dir = np.transpose(xyz)
    labCoor = []
    for d in dir :
        l,a,b = dirBound(d)
        l = l/50 -1
        a = a/labRange
        b = b/labRange
        labCoor.append([a,b,l])
    return np.transpose(labCoor)

def labColor(rtp) :
    x,y,z= s3d.SphericalSurface.coor_convert(rtp, tocart=True)
    L = 50*(z+1)
    a = labRange*x
    b = labRange*y
    Lab = np.transpose([L,a,b])
    rgb = cspace_converter("CAM02-UCS", "sRGB1" )(Lab)
    return np.transpose(rgb) 

# Setup surface ................................................
rez=4
start_azim = 0

surface = s3d.SphericalSurface(rez)
labCoor = toCoor(surface.vertices)
print("lab coordinate calculation completed....")
surface.map_geom_from_op( lambda rtp : labCoor , returnxyz=True)
surface.map_color_from_op(labColor)
surface.transform(scale=labRange,translate=[0,0,labRange])
surface = (surface + Cube() )

reldir = s3d.rtv(illum,elev,start_azim)
surfcopy = copy.copy(surface)
surfcopy.shade(.3,direction=reldir)
surfcopy.set_edgecolor([0,0,0,0])

# Construct figure, add surface, plot ..........................

fig = plt.figure(figsize=(5,5), facecolor='white')
ax = plt.axes(projection='3d')
adj = 1.1
xylim, zlim = (-adj*50, adj*50) , (0, adj*100)
ax.set(xlim=xylim, ylim=xylim, zlim=zlim)
ax.set_facecolor('black')
ax.set_axis_off()
ax.view_init(elev, start_azim)

ax.add_collection3d(surfcopy)

fig.tight_layout()

# 4. Animation ......................................................

def init_fig():
    return surfcopy,

def update_fig(frame):
    global surfcopy
    ax.collections.remove(surfcopy)

    azim = 360*frame
    reldir = s3d.rtv(illum,elev,azim)
    surfcopy = copy.copy(surface)
    surfcopy.shade(.3,direction=reldir)
    surfcopy.set_edgecolor([0,0,0,0])

    ax.add_collection3d(surfcopy)
    ax.view_init(elev=elev, azim=azim)

    return surfcopy,

ani = FuncAnimation(fig, update_fig, frames=np.linspace(0, 1, 72),
                    init_func=init_fig, blit=False, repeat=True, interval=100)

print(">>>>>>>>>>>>>>> Animation completed, file save proceeds")
#ani.save('ZZZ.mp4')                                   # use for movie file.
ani.save(None,writer=animation.FFMpegFileWriter())    # use for temp files.
print(">>>>>>>>>>>>>>> Save completed, screen display proceeds")
#plt.show()
print(">>>>>>>>>>>>>>> process completed")
