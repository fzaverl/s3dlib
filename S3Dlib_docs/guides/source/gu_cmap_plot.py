import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from colorspacious import cspace_converter
import s3dlib.cmap_utilities as cmu

def Lstar(cmap) :
    ''' x,y arrays for ploting L* of a colormap.'''
    if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
    N = 256
    if isinstance(cmap,ListedColormap) : N = len(cmap.colors)
    rgb = cmap(range(N))[:,0:3]
    lab = cspace_converter("sRGB1", "CAM02-UCS")(rgb)        
    L = np.transpose(lab)[0]
    x = np.linspace(0, 1, N)
    return x,L
# ...................................................
vhsv, Lhsv   = Lstar('hsv')
vdflt, Ldflt = Lstar(cmu.hue_cmap())
v_80, L_80   = Lstar(cmu.hue_cmap(8))
v_05, L_05   = Lstar(cmu.hue_cmap(.5))

fig = plt.figure(figsize=(6,2.5))
ax = fig.add_subplot(1,1,1)
ax.set(xlim=(0,1), ylim=(0,100))
ax.set_ylabel('L* value')

ax.plot(v_05,L_05,    label='0.5 ',         color='g', linestyle='-.')
ax.plot(vhsv, Lhsv,   label='1.0, HSV',     color='k')
ax.plot(vdflt, Ldflt, label='1.6, default', color='r')
ax.plot(v_80,L_80,    label='8.0  ',        color='y', linestyle=':')
ax.legend()

plt.show()