import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import s3dlib.cmap_utilities as cmu
from colorspacious import cspace_converter

# Plot L* ,lightness, for colormaps =================================

cmu.hsv_cmap_gradient( [.67,1,1,1], [.67,1,1,0], 'BuA' )
cmu.hsv_cmap_gradient( [.67,1,1,1], [.17,1,1,0], 'BuYwA' )
cmu.hsv_cmap_gradient( [0,1,1,0], [0,1,1,1], 'ARd' )
cmu.hsv_cmap_gradient( [.67,0,1,0], [.67,1,1,1], 'BuABu', mirrored=True )
cmu.stitch_cmap('BuA', 'ARd', name='BuAARd' )
cmu.rgb_cmap_gradient( [0,0,1,0], [0,0,1,0], 'clear' )
cmu.stitch_cmap('inferno', 'clear', name='inf_Clr' )
cmu.rgb_cmap_gradient([0,0,0,0],[0,0,0,1],'ABk' )
cmu.stitch_cmap('ABk', 'inferno', name='Abk_inf' )

cmaps = [ 'BuA', 'BuYwA','BuABu', 'BuAARd' ,'inf_Clr', 'Abk_inf' ]

# ===================================================================
def show_trans_cmap(ax,cmap,Lstr=None,aspect=24) :
    # ...............................................
    def over_color(lower_color, upper_color) :  
        opaque, transparent = 0.99, 0.01
        alpha_A, alpha_B = upper_color[3], lower_color[3]
        alpha_0 = alpha_A + alpha_B*( 1 - alpha_A )
        if alpha_A >= opaque :      return upper_color
        if alpha_B <= transparent : return lower_color
        if alpha_0 <= transparent : return np.array([0,0,0,0])
        C_A, C_B = np.array(upper_color[:3]), np.array(lower_color[:3])
        c0 = C_A*alpha_A + C_B*alpha_B*( 1 - alpha_A)
        C_0 = np.divide(c0,alpha_0)
        return np.append(C_0,alpha_0)
    # ...............................................
    sqrsz, sph = 4, 3    # pixels per square, squares per height
    bgrd, blck = [1.0,1.0,1.0,1.0], [0.8,0.8,0.8,1.0]
    if Lstr is not None :
        bgrd, blck = [1.0,0.7,0.7,1.0], [0.6,1.0,0.6,1.0]
    hdim = sph*sqrsz
    wdim = aspect*hdim
    i_ck = True
    X = np.empty([hdim, wdim,4])
    X[::] = np.array(bgrd)
    for i in range(hdim) :
        if i%sqrsz == 0 : i_ck = not i_ck
        j_ck = i_ck
        for j in range(wdim) :
            if j%sqrsz == 0 : j_ck = not j_ck
            bottom = X[i,j]
            if j_ck : bottom = blck
            # ---------------------
            if Lstr is not None :
                val = Lstr[int(255*j/wdim)]/100.0
                top = [val,val,val,cmap(j/wdim)[3] ]
            else :
                top = cmap(j/wdim)
            X[i,j] = over_color(bottom, top)
    ax.imshow(X, aspect='equal')
    return
# ...................................................................
def hasTransparency(cmap) :
    if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
    N = 256
    if isinstance(cmap,ListedColormap) : N = len(cmap.colors)
    alpha = cmap(range(N))[:,3:4]
    return np.any( np.where(alpha<1,True,False) )
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
    cmap = cmaps[i]
    if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
    x,y,name = Lstar(cmap)
    grad =  np.float32(np.vstack((y, y)))
    ax = axes[i]
    if hasTransparency(cmap) :
        show_trans_cmap(ax,cmap,y)
    else :
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
    if hasTransparency(cmap) :
        show_trans_cmap(ax,cmap)
    else :
        ax.imshow(grad, aspect='auto', cmap=cmap)
    pos = list(ax.get_position().bounds)
    x_text = pos[0] - 0.01
    y_text = pos[1] + pos[3]/2.
    fig.text(x_text, y_text, cmap.name, va='center', ha='right', fontsize=10)

# ===================================================================
plt.show()