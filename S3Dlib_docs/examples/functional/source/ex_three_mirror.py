import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Mirrored Colormap Usage

# 1. Define functions to examine ....................................

def twistFunction(rtz,twists) :
    r,t,z = rtz
    # Note: sliced surface needed due to discontinuity @ t=0 if twists is odd
    thickness = 0.33
    w = thickness*z 
    phi = 0.5*t*twists
    R = 1 + w * np.cos(phi)
    Z = w * np.sin(phi)
    return R,t,Z

# 2 & 3. Setup surfaces and plot ....................................
rez = 4
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboardMrrd',mirrored=True)
hsvr = cmu.reversed_cmap(cmu.hue_cmap())
cmaps = [hsvr,cmu.mirrored_cmap('magma'),'cardboardMrrd']
title = ['HSV','magma_m','cardboardMrrd']

fig = plt.figure(figsize=(8,2.1))
for i in range(0,3) :
    twist = s3d.CylindricalSurface(rez, basetype='squ_s')
    twist.map_geom_from_op( lambda rtz : twistFunction(rtz,1) )
    twist.map_cmap_from_normals(cmap=cmaps[i], direction=[1,1,1])
    if i==0 : twist.shade(.5)

    ax = fig.add_subplot(1,3,i+1, projection='3d')
    ax.set(xlim=(-0.8,0.8), ylim=(-0.8,0.8), zlim=(-0.8,0.8) )
    ax.set_title(title[i])
    ax.set_axis_off()
    ax.view_init(30,-20)
    plt.colorbar(twist.cBar_ScalarMappable, ax=ax, shrink=0.6 )

    ax.add_collection3d(twist)

fig.tight_layout()
#....................................................................
redblue = cmu.hsv_cmap_gradient('b','+r')

surf = s3d.CylindricalSurface(rez, basetype='squ_s', linewidth=0)
surf.map_geom_from_op( lambda rtz : twistFunction(rtz,1) )
surf.map_cmap_from_normals(cmap=redblue, direction=[1,1,1])
surf.set_surface_alpha(0.8)

surfVect = s3d.CylindricalSurface(2, basetype='squ_s', linewidth=0)
surfVect.map_geom_from_op( lambda rtz : twistFunction(rtz,1) )
surfVect.set_surface_alpha(0.0)

fig = plt.figure(figsize=(4,3))
ax = plt.axes(projection='3d')
ax.set_title('Face Normals')
ax.set(xlim=(-0.8,0.8), ylim=(-0.8,0.8), zlim=(-0.8,0.8) )
ax.set_axis_off()
ax.view_init(30,-20)
plt.colorbar(surf.cBar_ScalarMappable, ax=ax, shrink=0.6 )

ax.add_collection3d(surfVect.facenormals(.3))
ax.add_collection3d(surf.shade(.5))

fig.tight_layout()
#....................................................................
plt.show()