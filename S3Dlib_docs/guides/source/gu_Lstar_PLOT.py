import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import s3dlib.cmap_utilities as cmu
from colorspacious import cspace_converter

# Plot L* ,lightness, for colormaps =================================

cmaps = [ cmu.reversed_cmap('binary'),'viridis', 'plasma', 'inferno', 'magma', 'cividis' ]

# ...................................................................
def Lstar(cmap) :
    ''' x,y arrays for ploting L* of a colormap.'''
    if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
    N = 256
    if isinstance(cmap,ListedColormap) : N = len(cmap.colors)
    rgb = cmap(range(N))[:,0:3]
    lab = cspace_converter("sRGB1", "CAM02-UCS")(rgb)        
    L = np.transpose(lab)[0]
    x = np.linspace(0, 1, N)
    return x,L,cmap.name

# Figure 1: x,L* line plot. -----------------------------------------
colors = ['k','r','g','b','y','m','c']
linestyles = ['-', '--', '-.', ':']
fig = plt.figure(figsize=plt.figaspect(0.6))
ax = fig.add_axes([0.1, 0.1, 0.65, 0.85])
ax.set(xlim=(0,1), ylim=(0,100))
ax.set_ylabel('L* value')

Lcol, Lstyles = len(colors), len(linestyles)
for i in range(len(cmaps)) :
    linestyle = linestyles[i%Lstyles]
    color = colors[i%Lcol]
    x,y,name = Lstar(cmaps[i])
    ax.plot(x,y, label=name, color=color, linestyle=linestyle)
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.)

# Figure 2: Visual L* gray scale. -----------------------------------
numb =len(cmaps)
w, h = 6, numb*0.25
fig, axes = plt.subplots(nrows=numb, figsize=(w,h))
fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
if numb==1 : axes = [axes]
for ax in axes :
    ax.set_xticks([])
    ax.set_yticks([])

for i in range(len(cmaps)):
    x,y,name = Lstar(cmaps[i])
    grad =  np.float32(np.vstack((y, y)))
    ax = axes[i]
    ax.imshow(grad, aspect='auto', cmap='binary_r', vmin=0., vmax=100.)
    pos = list(ax.get_position().bounds)
    x_text = pos[0] - 0.01
    y_text = pos[1] + pos[3]/2.
    fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

# Figure 3: Visual color scale. -------------------------------------
numb =len(cmaps)
w, h = 6, numb*0.25
fig, axes = plt.subplots(nrows=numb, figsize=(w,h))
fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
if numb==1 : axes = [axes]
for ax in axes :
    ax.set_xticks([])
    ax.set_yticks([])

grad = np.linspace(0, 1, 256)
grad = np.vstack((grad, grad))
for i in range(len(cmaps)):
    ax = axes[i]
    cmap=cmaps[i]
    if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
    ax.imshow(grad, aspect='auto', cmap=cmap)
    pos = list(ax.get_position().bounds)
    x_text = pos[0] - 0.01
    y_text = pos[1] + pos[3]/2.
    fig.text(x_text, y_text, cmap.name, va='center', ha='right', fontsize=10)

# ...................................................................

plt.show()