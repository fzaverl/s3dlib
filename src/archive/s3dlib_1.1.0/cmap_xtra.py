# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

'''
Auxiliary functions for use with Matplotlib.

'''
import string 
import random
import numpy as np

from matplotlib import cm, colors
from matplotlib.colors import ListedColormap

from colorspacious import cspace_converter

# +---------------------------------------------------------+
# |  NOTE: depricated from 3.6:                             |
# |            matplotlib.cm.get(name)                      |
# |        to  matplotlib.colormaps[name]                   |
# |  Note: depricated from 3.6:                             |
# |            matplotlib.cm.register_cmap(name)            |
# |        to  matplotlib.colormaps.register(name)          |
# +---------------------------------------------------------+

def Lab_cmap_gradient(lowColor='k', highColor='w', name=None, mirrored=False) :
    """
    A linear-in-Lab-space Colormap, with option of registering the map.

    Parameters
    ----------
    lowColor : RGB color, optional, default: 'black'
        Color at the low end of the Colormap range.
    
    highColor : RGB color, optional, default: 'white'
        Color at the high end of the Colormap range.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will be a string of 8 random
        characters.

    mirrored : bool
        If True, colormap is divided into two linear
        segments with the lowColor at the low and high
        values, the highColor in the middle.
        
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """
    
    if name is not None :
        if not isinstance(name, str ):
            raise ValueError('Invalid name argument (str required): ' + str(name))        

    if isinstance(lowColor, str ):  lowColor = colors.to_rgb(lowColor)
    if isinstance(highColor, str ): highColor= colors.to_rgb(highColor)

    if not isinstance(lowColor,(list, tuple, np.ndarray)) :
        raise ValueError('Invalid lowColor argument: ' + str(name))
    if not isinstance(highColor,(list, tuple, np.ndarray)) :
        raise ValueError('Invalid highColor argument: ' + str(name))  

    lowColor, highColor  = list(lowColor), list(highColor)
    if len(lowColor) == 3 : lowColor.append(1.0)
    if len(highColor) == 3 : highColor.append(1.0)
    if len(lowColor) != 4 or len(highColor) != 4  :
        raise ValueError('Invalid HSVarg argument list length ')
    lowLab = cspace_converter("sRGB1", "CAM02-UCS")(lowColor[:3])
    hghLab = cspace_converter("sRGB1", "CAM02-UCS")(highColor[:3])
    
    deltaL = np.subtract(hghLab,lowLab)
    deltaA = highColor[3] - lowColor[3]
    numbSegs = 256
    x = np.linspace(0.0,1.0,num=numbSegs)
    if mirrored : x = abs( 2.0*x - 1.0 )

    # DevNote:  insignificant time spent in loop. no need to opt.
    clist = []
    for n in x :
        lab = np.add(lowLab, np.multiply(deltaL,n) )
        alf = lowColor[3] + n*deltaA
        rgb = cspace_converter("CAM02-UCS", "sRGB1")(lab)
        rgb = np.clip(rgb,0,1)
        rgba = np.append(rgb,alf)
        clist.append(rgba)  
    cmap = ListedColormap(clist)
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    cm.register_cmap(name,cmap)
    cmap.name = name
    return cmap


def Hue_Lab_gradient( color, lowL=0.0, hiL=1.0, name=None ) :
    """
    A linear-in-Lab-space Colormap with a constant hue, with option of registering the map.
    
    Parameters
    ----------
    color : RGB color(str,list,tuple) or number in the domain[0,1]
        Vaule is used to
        
    lowL : number, optional, default: 0.0
        The minimum L* lightness value in range 0 to <1

    hiL  : number, optional, default: 1.0
        The maximum L* lightness value in range >0 to 1.

    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None and color is a string, name is assigned to the string_L.
        Otherwize, if it's None, name is assigned a string of 8 random characters.
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """
    # ...............................................................
    def get_arrays(hue,N) :
        hsv_0 = np.array([hue,1.0,1.0])
        rgb_0 = np.array(colors.hsv_to_rgb(hsv_0.T)).T
        L_0,a,b = cspace_converter("sRGB1", "CAM02-UCS"  )(rgb_0.T).T
        L_0 = L_0/100.0
        L = np.linspace(0.0,1.0,N)
        endS = (1.0 -L)/(1.0 - L_0)
        sttV = L/L_0
        ones = np.ones(len(L))
        L0vals = np.full(len(L),L_0)
        H = np.full(len(L),hue)
        S = np.where( L>L0vals, endS, ones )
        V = np.where( L>L0vals, ones, sttV )
        hsv = np.array( [H,S,V] )
        rgb = np.array(colors.hsv_to_rgb(hsv.T)).T
        Lact,a,b = cspace_converter("sRGB1", "CAM02-UCS"  )(rgb.T).T
        Lact = Lact/100.0
        Lact[0]=0.0
        Lact[len(Lact)-1] = 1.000001
        return Lact,S,V
    # ...............................................................
    def linIntrp(xVal,xArr,yArr) :
        i = np.where(xArr>xVal)[0][0]
        x0,x1 = xArr[i-1], xArr[i]
        y0,y1 = yArr[i-1], yArr[i]
        m = (y1-y0)/(x1-x0)
        yVal = (xVal-x0)*m + y0
        return yVal
    # ...............................................................

    if isinstance(color, str ): 
        hue,s,v = colors.rgb_to_hsv(colors.to_rgb(color))
        if name is None : name = color+'_L'
    elif isinstance(color,(list,tuple)) :
        hue,s,v = colors.rgb_to_hsv(colors.to_rgb(color))
    else :
        hue = color
        if hue<0.0 or hue>1.0 :
            raise ValueError('Error: required that 0.0 <= hue <= 1.0. , found {}'.format(hue))

    Lval,Sval,Vval = get_arrays(hue,100)

    N=256
    L = np.linspace(lowL,hiL,N)
    H = np.full(len(L),hue)
    S = [ linIntrp(Lst,Lval,Sval) for Lst in L ]
    V = [ linIntrp(Lst,Lval,Vval) for Lst in L ]
    hsv = np.array( [H,S,V] )
    rgb = np.array(colors.hsv_to_rgb(hsv.T)).T
    rgba = np.vstack( (rgb, np.ones(len(L)) ) ).T

    cmap = ListedColormap(rgba)
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    cm.register_cmap(name,cmap)
    cmap.name = name
    return cmap


def Cmap_Lab_gradient( cmap, lowL=0.0, hiL=1.0, name=None ) :
    """
    A linear-in-Lab-space Colormap with a hue matched to cmap and maximum saturation.
    
    Parameters
    ----------
    cmap : colormap

    lowL : number, optional, default: 0.0
        The minimum L* lightness value in range 0 to <1

    hiL  : number, optional, default: 1.0
        The maximum L* lightness value in range >0 to 1.

    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None and color is a string, name is assigned to the string_L.
        Otherwize, if it's None, name is assigned a string of 8 random characters.

    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """
    # ...............................................................

    def get_arrays(cmap,lw,hi,N) :
        cols = cmap(np.linspace(0,1,N))
        hsv = colors.rgb_to_hsv(cols[:,:3])
        hue = hsv[:,0]
        hsv_0 = np.ones((N,3))
        hsv_0[:,0] = hue
        rgb_0 = np.array(colors.hsv_to_rgb(hsv_0)).T
        L_0,a,b = cspace_converter("sRGB1", "CAM02-UCS"  )(rgb_0.T).T
        L_0 = L_0/100.0
        L = np.linspace(lw,hi,N)
        endS = (1.0 -L)/(1.0 - L_0)
        sttV = L/L_0
        ones = np.ones(len(L))
        L0vals = np.full(len(L),L_0)
        H = hue
        S = np.where( L>L0vals, endS, ones )
        V = np.where( L>L0vals, ones, sttV )
        hsv = np.array( [H,S,V] )
        rgb = np.array(colors.hsv_to_rgb(hsv.T)).T
        Lact,a,b = cspace_converter("sRGB1", "CAM02-UCS"  )(rgb.T).T
        Lact = Lact/100.0
        Lact[0]=lw
        Lact[len(Lact)-1] = 1.000001*hi
        return Lact,a,b   #<<<<< linear in lab space.
        
    # ...............................................................
    def linIntrp(xVal,xArr,yArr) :
        i = np.where(xArr>xVal)[0][0]
        x0,x1 = xArr[i-1], xArr[i]
        y0,y1 = yArr[i-1], yArr[i]
        m = (y1-y0)/(x1-x0)
        yVal = (xVal-x0)*m + y0
        return yVal
    # ...............................................................
    if isinstance(cmap,str) : 
        cmap = cm.get_cmap(cmap)
    if name is None :
        name = cmap.name + '_L'

    Lval,aval,bval = get_arrays(cmap,lowL,hiL,100)

    N=256
    L = np.linspace(lowL,hiL,N)
    a = [ linIntrp(Lst,Lval,aval) for Lst in L ]
    b = [ linIntrp(Lst,Lval,bval) for Lst in L ]
    lab = np.array( [100*L,a,b]).T
    rgb = cspace_converter("CAM02-UCS", "sRGB1")(lab)
    rgba = np.vstack( (rgb.T, np.ones(len(L)) ) ).T
    rgba = np.clip(rgba, 0.0,1.0)
    cmap = ListedColormap(rgba)
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    cm.register_cmap(name,cmap)
    cmap.name = name
    return cmap


def show_cmaps( plt, cmaps, onlyColormaps=True, show=True) :
    """
    Construct and show Matplotlib figures of colormaps.
    
    Parameters
    ----------
    plt : matplotlib.pyplot

    cmaps : List of colormaps.
        List contains colormaps or strings of registered
        colormap names.

    onlyColormaps : bool {True,False}, optional, default: True
        If True, only a color scale figure is shown.  If False,
        the L* gray scale and L* line plot figures are also shown.

    show : bool {True,False}, optional, default: True
        If True, the plt method 'show()' will be executed.

    """

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
        ax.imshow(X, aspect='auto')
        return
    # ...................................................................
    def hasTransparency(cmap) :
        if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
        N = 256
        alpha = cmap(np.linspace(0, 1, N))[:,3:4]
        return np.any( np.where(alpha<1,True,False) )
    # ...................................................................
    def Lstar(cmap) :
        # x,y arrays for ploting L* of a colormap.
        if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)
        N = 256
        rgb = cmap(np.linspace(0, 1, N))[:,0:3]
        lab = cspace_converter("sRGB1", "CAM02-UCS")(rgb)        
        L = np.transpose(lab)[0]
        x = np.linspace(0, 1, N)
        return x,L,cmap.name

    # ===================================================================

    # Figure 1: Visual color scale. -------------------------------------
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


    if not onlyColormaps :
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

        # Figure 3: x,L* line plot. -----------------------------------------
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

    # ===================================================================
    if show : plt.show()
    return

