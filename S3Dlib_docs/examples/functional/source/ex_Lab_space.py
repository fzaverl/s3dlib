import warnings

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import s3dlib.surface as s3d
from colorspacious import cspace_converter

#.. Lab Color Space

labRange = 50

# 1. Define function to examine ...............

def dirBound(direction) :
    # simple binary-chop to determime limits in a Lab direction.
    #........................................
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
rez=5

surface = s3d.SphericalSurface(rez)
labCoor = toCoor(surface.vertices)
print("lab coordinate calculation completed....")
surface.map_geom_from_op( lambda rtp : labCoor , returnxyz=True)
surface.map_color_from_op(labColor).shade(.5)
surface.transform(translate=[0,0,1])
surface.transform(scale=labRange)

# Construct figure, add surface, plot ..........................

fig = plt.figure(figsize=(5,5))
ax = plt.axes(projection='3d')
ax.set(xlim=(-50,50), ylim=(-50,50), zlim=(0,100))
ax.set_xticks([-50,0,50])
ax.set_yticks([-50,0,50])
ax.set_zticks([0,50,100])
ax.set_xlabel('a')
ax.set_ylabel('b')
ax.set_zlabel('L')
ax.set_title( 'sRGB Lab Color Space', fontsize='x-large' )
ax.set_proj_type('ortho')

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()