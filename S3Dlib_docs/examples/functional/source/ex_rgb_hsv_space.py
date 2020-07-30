import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. RGB and HSV Color Spaces

# 1. Define functions to examine ....................................

def hsvColor(rtz) :
    r,t,z = rtz
    return t/(2*np.pi), r, z  # all values are in [0,1]

def Cylinder(rez) :
    # .....................................................
    def fold(rtz) :
        r,t,z = rtz
        zeros = np.zeros(len(z))
        ones = np.ones(len(z))
        # fold the cylinder into 3 parts..
        alpha = -2*z + 2
        alpha = 2*(ones-z)
        alpha = np.where( z <= 0.5, ones ,       alpha )
        alpha = np.where( z <= -.5, 2*(ones+z) , alpha )
        beta = ones
        beta = np.where( z <= 0.5, 2*z  , beta)
        beta = np.where( z <= -.5, -ones, beta)
        R = np.clip(alpha,0.001,1) 
        Z = beta
        return R,t,Z
    # .....................................................
    surface = s3d.CylindricalSurface(rez)    
    surface.map_geom_from_op( lambda rtz : fold(rtz) )
    surface.name = 'cylinder'
    return surface

def Cube(rez) :  
    v = [ 
        [  0, 0, 0 ],  [  0, 1, 0 ],  [ 1 , 1, 0 ],  [  1, 0, 0 ],
        [  0, 0, 1 ],  [  0, 1, 1 ],  [ 1 , 1, 1 ],  [  1, 0, 1 ]  ]
    f = [ [0,1,2,3], [3,2,6,7], [2,1,5,6], [1,0,4,5], [0,3,7,4], [4,7,6,5] ]
    vertexCoor = np.array(v).astype(float)
    faceIndices = np.array(f)
    surface = s3d.Surface3DCollection(vertexCoor, faceIndices)
    surface.transform(scale=2, translate=[-1,-1,-1])
    surface.triangulate(rez)
    surface.name = 'cube'
    return surface

# 2. Setup and map surfaces .........................................
rez = 5

rgb_cube = Cube(rez)
rgb_cube.transform(scale=.5,translate=[.5,.5,.5])
rgb_cube.map_color_from_op( lambda xyz : xyz).shade(.5)

hsv_cyl = Cylinder(rez)
hsv_cyl.transform(scale=[1,1,0.5],translate=[0,0,0.5])
hsv_cyl.map_color_from_op(hsvColor,rgb=False).shade(.5)

# 3. Construct figure, add surfaces, and plot .....................

hsv_minmax = (-1.3,1.3)
rgb_minmax, rgb_ticks = (-0.2,1.2) , [0,1]
fig = plt.figure(figsize=(6,3))

ax1 = fig.add_subplot(121, projection='3d')
ax1.set(xlim=rgb_minmax, ylim=rgb_minmax, zlim=rgb_minmax)
ax1.set_xticks(rgb_ticks)
ax1.set_yticks(rgb_ticks)
ax1.set_zticks(rgb_ticks)
ax1.set_xlabel('R', labelpad=-10)
ax1.set_ylabel('G', labelpad=-10)
ax1.set_zlabel('B', labelpad=-10)
ax1.set_title( 'RGB Color Space' )
ax1.set_proj_type('ortho')

ax2 = fig.add_subplot(122, projection='3d')
ax2.set(xlim=hsv_minmax, ylim=hsv_minmax, zlim=(0,1))
ax2.set_xticks([])
ax2.set_yticks([0])
ax2.set_zticks([0,1])
ax2.set_zlabel('Value', labelpad=-10)
ax2.set_title( 'HSV Color Space' )
ax2.set_proj_type('ortho')

ax1.add_collection3d(rgb_cube)
ax2.add_collection3d(hsv_cyl)

fig.tight_layout()
plt.show()