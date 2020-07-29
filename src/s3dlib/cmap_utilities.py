# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

'''
Module containing functions to create Matplotlib color maps.

'''
# Not much going on here but these functions simplify constructing
# colormaps during development.  Most custom colormaps can be made
# with just one or two lines of code which are easier to comprehend.
# Got tired of remembering which module and what function  ;-)


import string 
import random
import copy

import numpy as np

from matplotlib import cm, colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


def rgb_cmap_gradient(lowColor='k', highColor='w', name=None, mirrored=False) :
    """
    A linear-in-RGB-space Colormap, with option of registering the map.

    Parameters
    ----------
    locColor : color, optional, default: 'black'
        Color at the low end of the Colormap range.
    
    highColor : color, optional, default: 'white'
        Color at the high end of the Colormap range.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will be a string of random
        characters and the map is not registered.

    mirrored : bool
        If True, colormap is divided into two linear
        segments with the lowColor at the low and high
        values, the highColr in the middle.
        
    Returns
    -------
    LinearSegmentedColormap
        An instance of a colormap.
    """

    if not isinstance(mirrored, bool ) :
        raise ValueError('Invalid mirrored argument (bool required): ' + str(name))       
    call_name = name
    if name is None : 
        call_name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    else :
        if not isinstance(name, str ):
            raise ValueError('Invalid name argument (str required): ' + str(name))
    if not colors.is_color_like(lowColor) :
        raise ValueError('Invalid negColor argument: ' + str(lowColor))
    if not colors.is_color_like(highColor) :
        raise ValueError('Invalid posColor argument: ' + str(highColor))
    
    clist = [lowColor,highColor]
    if mirrored : clist = [lowColor,highColor,lowColor]
    cmap = LinearSegmentedColormap.from_list(call_name, clist )
    if name is not None :
        cm.register_cmap(name,cmap)
    return cmap

def _smoothHue(x,eps) :
    # ....................................
    def sgm( x, c) :
        xref= x - c
        p = 6*xref - 1
        t1 = np.power(np.abs(p),eps)
        t2 = np.sign(p)
        d = (1+t1*t2)/2
        val = c + d/3
        return val
    # ....................................
    y = sgm(x,0.0)
    y = np.where( x >0.333 , sgm(x,0.333), y )
    y = np.where( x >0.666 , sgm(x,0.666), y )
    return y

def hsv_cmap_gradient(lowHSVarg=[0,1,1], highHSVarg=[1,1,1], name=None, mirrored=False, smooth=None) :
    """
    A linear-in-HSV-space Colormap, with option of registering the map.
    
    Using hSV or hSVA values for HSV color. SVA values are in the range [0,1].
    The h values are in the range [0,2] so that Hue is mod(h,1).

    Parameters
    ----------
    lowHSVarg : 3 or 4 array or string, optional, default: [0,1,1], , 'red'
        HSV, HSVA color at the low end of the Colormap range.  If a string,
        a named color (may be proceeded by a '+')
    
    highHSVarg : 3 or 4 array or string, optional, default: [1,1,1] , '+red'
        HSV, HSVA color at the high end of the Colormap range.  If a string,
        a named color (may be proceeded by a '+')
    
    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, the map is not registered.

    mirrored : bool
        If True, colormap is divided into two linear
        segments with the lowColor at the low and high
        values, the highColr in the middle.
        
    smooth: float, optional, default: None
        Controls the ratio of the amount of CMY to RGB color.
        CMY is extended relative to the RGB hues for values
        greater than one (relative to standard HSV colormap).
        RGB is extended relative to the CMY hues for values
        less than one.  Range is [.1,10]

    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """
    # ............................................................
    def stringToHSV(sVal) :
        addOne = False
        if sVal[0] is '+':
            addOne = True
            sVal = sVal[1:]
        h,s,v = colors.rgb_to_hsv(colors.to_rgb(sVal))
        if addOne : h += 1
        return [h,s,v]    
    # ............................................................
    if not isinstance(mirrored, bool ) :
        raise ValueError('Invalid mirrored argument (boolean required): ' + str(name))       

    if name is not None :
        if not isinstance(name, str ):
            raise ValueError('Invalid name argument (str required): ' + str(name))        

    if isinstance(lowHSVarg, str ): lowHSVarg= stringToHSV(lowHSVarg)
    if isinstance(highHSVarg, str ): highHSVarg= stringToHSV(highHSVarg)

    if not isinstance(lowHSVarg,(list, tuple, np.ndarray)) :
        raise ValueError('Invalid lowHSVarg argument: ' + str(name))
    if not isinstance(highHSVarg,(list, tuple, np.ndarray)) :
        raise ValueError('Invalid highHSVarg argument: ' + str(name))  

    lowColor, highColor  = list(lowHSVarg), list(highHSVarg)
    if len(lowColor) == 3 : lowColor.append(1.0)
    if len(highColor) == 3 : highColor.append(1.0)
    if len(lowColor) != 4 or len(highColor) != 4  :
        raise ValueError('Invalid HSVarg argument list length ')

    lowColor[0] += 1.0
    highColor[0] += 1.0
    numbSegs = 256
    delta = np.subtract(highColor,lowColor)
    x = np.linspace(0.0,1.0,num=numbSegs)
    if mirrored : x = abs( 2.0*x - 1.0 )
    # DevNote:  insignificant time spent in loop. no need to opt.
    clist = []
    for n in x :
            h,s,v,a = np.add(lowColor, np.multiply(delta,n) )
            h = np.mod(h, 1)
            if smooth is not None : h = _smoothHue(h,smooth)
            r,g,b = cm.mpl.colors.hsv_to_rgb([h,s,v])
            clist.append([r,g,b,a])  
    cmap = ListedColormap(clist)
    if name is not None :
        cm.register_cmap(name,cmap)
        # name property not assigned for ListedColormap, so...
        cmap.name = name
    return cmap

def binary_cmap(negColor='b', posColor='r', name=None ) :
    """
    A two-color Colormap, with option of registering the map.

    Parameters
    ----------
    negColor : color, optional, default: 'blue'
        Color at the low end of the Colormap range.
    
    posColor : color, optional, default: 'red'
        Color at the high end of the Colormap range.
    
    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, the map is not registered.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """
    if name is not None :
        if not isinstance(name, str ):
            raise ValueError('Invalid name argument (str required): ' + str(name))        
    if not colors.is_color_like(negColor) :
        raise ValueError('Invalid negColor argument: ' + str(negColor))
    if not colors.is_color_like(posColor) :
        raise ValueError('Invalid posColor argument: ' + str(posColor))
    
    cmap = ListedColormap([negColor,posColor])
    if name is not None :
        cm.register_cmap(name,cmap)
        # name property not assigned for ListedColormap, so...
        cmap.name = name
    return cmap
  
def hue_cmap(smooth=1.6, lowHue=0, highHue=None, name=None ) :
    """
    A 'smooth-HSV' colormap, with option of registering the map.
    
    Non linear adjustment of Hue in a HSV colormap.
    For the HSV colors, S and V are 1.
    The h values are shifted based on the smoothing parameter.

    Parameters
    ----------
    smooth: float, optional, default: 1.6
        Controls the ratio of the amount of CMY to RGB color.
        CMY is extended relative to the RGB hues for values
        greater than one (relative to standard HSV colormap).
        RGB is extended relative to the CMY hues for values
        less than one.  Range is [.1,10]

    lowHue : float or string, optional, default: 0
        Hue at the low end of the Colormap range.
        Values are in the range [0,2].  
        If a string, the hue of a named color (may be
        proceeded by a '+')

    highHue : float or string, optional, default: 1
        Hue at the upper end of the Colormap range.  
        Values are in the range [0,2].  
        If a string, the hue of a named color (may be
        proceeded by a '+')

    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, the map is not registered.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """

    # ............................................................
    def stringToHue(sVal) :
        addOne = False
        if sVal[0] is '+':
            addOne = True
            sVal = sVal[1:]
        h,s,v = colors.rgb_to_hsv(colors.to_rgb(sVal))
        if addOne : h += 1
        return h    
    # ............................................................
    if (smooth<0.1 or smooth>10) :
        raise ValueError('smooth must be between 0.1 and 10. , found {}'.format(smooth))

    if isinstance(lowHue, str ): lowHue= stringToHue(lowHue)
    if highHue is None:
        highHue = lowHue+1
    else:
        if isinstance(highHue, str ):  highHue= stringToHue(highHue)

    delta = np.subtract(highHue,lowHue)
    x = np.linspace(0, 1, 256)
    x = np.add(lowHue, np.multiply(delta,x) )
    x = np.mod(x,1)       
    h = _smoothHue(x,smooth)
    one = np.ones(len(h))
    hsv = np.transpose([h,one,one])
    rgb = cm.mpl.colors.hsv_to_rgb(hsv)
    cmap = colors.ListedColormap(rgb)
    if name is not None :
        cm.register_cmap(name,cmap)
        # name property not assigned for ListedColormap, so...
        cmap.name = name
    return cmap

def stitch_cmap( *maps, bndry=None, name=None) :
    """
    Colormap composed of multiple colormaps.

    Parameters
    ----------
    *maps : colormaps
    
    bndry : scalar or list, default: None
        Boundary between colormaps.  If  a list, 
        the number of values must be one less than
        the number of maps.  The list values must
        be in assending order.  If None, input colormaps
        are evenly spaced in the returned colormap.
        The list values are in the range 0 to 1.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will be a string of random
        characters and the map is not registered.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """

    numbMaps = len(maps)
    if numbMaps == 1 : 
        raise ValueError('Only one cmap entered, need two or greater')  
    if bndry is None :
        bndry = [ i/numbMaps for i in range(1,numbMaps) ]
    else :
        numbBndry = len(bndry)
        if numbBndry != (numbMaps-1) :
            strgMsg = 'Number of boundaries, {}, is NOT equal to one less that the number of maps, {}. '
            raise ValueError(strgMsg.format(numbBndry,numbMaps))
    bndry.append(1.0)

    colormaps = []
    for cmap in maps :
        if isinstance(cmap,str) : cmap = cm.get_cmap(cmap)  
        colormaps.append(cmap)

    N = 256
    Nbndry = np.subtract(np.multiply(bndry,N),1)
    mapIndex, colorIndexStart = 0, 0
    clist = []
    for i in range(N) :
        if i > Nbndry[mapIndex] : 
            mapIndex += 1
            colorIndexStart = i
        d = 1/(Nbndry[mapIndex]-colorIndexStart)
        index = d*(i-colorIndexStart)
        clist.append( colormaps[mapIndex](index) )
    cmap = ListedColormap(clist)
    if name is not None :
        cm.register_cmap(name,cmap)
        # name property not assigned for ListedColormap, so...
        cmap.name = name
    return cmap

def alpha_cmap( cmap, alpha, name=None )  :
    """ A registered colormap with a constant alpha channel."""
    if alpha<0.0 or alpha>1.0 :
        raise ValueError('alpha must be between 0 and 1, found {}'.format(alpha))       
    
    if isinstance(cmap,str) : 
        cmap = cm.get_cmap(cmap)
    if name is None :
        name = cmap.name + '_a'

    if isinstance(cmap,LinearSegmentedColormap) :
        # LinearSegmentedColormap ...
        cdict = copy.copy(cmap._segmentdata)
        cdict['alpha'] = [ (0.0, alpha,alpha),(1.0, alpha,alpha) ]
        cm.register_cmap(name, data=cdict)
        newCmap = cm.get_cmap(name)
    else :
        # ListedColormap ............
        colorArr = copy.copy(cmap.colors) 
        colorArr = np.insert(colorArr,3,alpha,axis=1)[:,:4]
        newCmap = ListedColormap(colorArr,name)
        cm.register_cmap(name,newCmap)
    return newCmap

def mirrored_cmap( cmap, name=None, rev=False) :
    """ A registered mirrored colormap."""  
    if isinstance(cmap,str) : 
        cmap = cm.get_cmap(cmap)
    if name is None :
        name = cmap.name + '_m'
        if rev : name += 'r'

    numbSegs = 256
    x = np.linspace(0.0,1.0,num=numbSegs)
    n = abs( 2.0*x - 1.0 )
    if rev : n = 1-n
    clist =cmap(n)
    cmap = ListedColormap(clist)
   
    cm.register_cmap(name,cmap)
    # name property not assigned for ListedColormap, so...
    cmap.name = name
    return cmap

def reversed_cmap( cmap, name=None ) :
    """ A registered reversed colormap."""
    if isinstance(cmap,str) : 
        cmap = cm.get_cmap(cmap)    
    cmap = cmap.reversed()
    if name is None :
        name = cmap.name  # note: reverse returns name + '_r'
    
    cm.register_cmap(name,cmap)
    return cmap
