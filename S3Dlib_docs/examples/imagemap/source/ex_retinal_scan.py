import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Function and Image Mapping 

# 1. Define functions to examine ....................................

def collapse( rtp) :
    r,t,p = rtp
    T = np.where(t>np.pi,t,-t)
    R = np.where(t>np.pi,r+0.03,r)
    return R,T,p

# 2. Setup and map surfaces .........................................
rez=6

cmap = cmu.binary_cmap(negColor='black' ,posColor=[0.918,0.549,0.377])
imagefilename = 'data/retinal_scan.png'
vp = [0.0,0.361,0.5,0.639]

eye = s3d.SphericalSurface(rez,basetype='octa')
eye.map_cmap_from_op( lambda rtp : rtp[1] ,cmap=cmap) 
eye.map_color_from_image(imagefilename,viewport=vp)
eye.map_geom_from_op(collapse)
eye.shade(direction=[0,1,1],contrast=0.5)
eye.transform(rotate=s3d.eulerRot(-200,30))

# 3. Construct figures, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1) )
ax = plt.axes(projection='3d')
minmax = (-0.75,0.75)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_title('Retinal  Scan')
ax.set_axis_off()

ax.add_collection3d(eye)

fig.tight_layout()
plt.show()