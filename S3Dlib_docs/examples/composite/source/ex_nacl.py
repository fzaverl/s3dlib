from collections import OrderedDict
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import s3dlib.surface as s3d

#.. Sub-surface Translation

atomPos = OrderedDict()

atomPos['Chloride'] = [
    [ 1, 1, 1], [-1,-1, 1], [-1, 1, 1], [ 1,-1, 1] ,
    [ 1, 1,-1], [-1,-1,-1], [-1, 1,-1], [ 1,-1,-1] ,
    [ 1, 0, 0], [ 0, 1, 0], [-1, 0, 0], [ 0,-1, 0] , [ 0, 0, 1], [ 0, 0,-1]    ]

atomPos['Sodium'] = [
    [ 1, 0, 1], [-1, 0, 1], [ 0, 1, 1], [ 0,-1, 1] ,
    [ 1, 1, 0], [-1,-1, 0], [-1, 1, 0], [ 1,-1, 0] , [ 0, 0, 0] ,  
    [ 1, 0,-1], [-1, 0,-1], [ 0, 1,-1], [ 0,-1,-1]   ]

rez = 3
sz = 0.67
Cl_color, Na_color = 'yellowgreen', 'blueviolet'
color = Cl_color
unitCell = None

for elemName, elemPos in atomPos.items() :
    for pos in elemPos :
        atom = s3d.SphericalSurface(rez, facecolor=color)
        atom.transform(scale=[sz,sz,sz] , translate=pos)
        atom.shade(0.25,direction=[1,1,1])
        if unitCell is None : unitCell = atom
        else : unitCell += atom
    sz = 1.0 - sz
    color = Na_color

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,str(unitCell), ha='right', va='top', fontsize='smaller', multialignment='right')

ax = plt.axes(projection='3d')
minmax = (-1.5,1.5)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_proj_type('ortho')
ax.set_title('NaCl Unit Cell\n')
ax.set_axis_off()
cl_patch = mpatches.Patch(color=Cl_color, label=r'$Cl^-$')
na_patch = mpatches.Patch(color=Na_color, label=r'$Na^+$')
ax.legend(handles=[cl_patch,na_patch])

ax.add_collection3d(unitCell)

plt.show()