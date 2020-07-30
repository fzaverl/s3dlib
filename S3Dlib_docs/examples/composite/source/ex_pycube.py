import copy
import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Composite of Copies

# 2. Setup and map surface .........................................
rez = 6

top = s3d.PlanarSurface(rez, basetype='oct1')
top.map_color_from_image('data/python.png')
top.map_geom_from_image('data/python_elevation.png',0.05)

front = copy.copy(top)
side = copy.copy(top)
bottom = copy.copy(top)
backside = copy.copy(top)
top.transform(translate=[0,0,1])
backfront = copy.copy(top)
bottom.transform(rotate=s3d.eulerRot(0,180),  translate=[0,0,-1])
front.transform(rotate=s3d.eulerRot(0,90),  translate=[0,-1,0])
backfront.transform(rotate=s3d.eulerRot(180,90),  translate=[0,0,0])
side.transform(rotate=s3d.eulerRot(90,90),  translate=[1,0,0])
backside.transform(rotate=s3d.eulerRot(-90,90),  translate=[-1,0,0])

cube = (top+front+side+bottom+backfront+backside)
cube.shade().hilite(.9)

# 3. Construct figure, add surface plot ............................

fig = plt.figure(figsize=plt.figaspect(1), facecolor='black')
fig.text(0.975,0.975,str(cube), ha='right', va='top',
        fontsize='smaller', multialignment='right', color='white')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))
ax.set_facecolor('black')
ax.set_axis_off()

ax.add_collection3d(cube)

fig.tight_layout()
plt.show()