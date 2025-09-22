# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

'''
Module containing functions to create Matplotlib color maps.

'''
# Not much going on here but these functions simplify constructing
# colormaps during development.  Most custom colormaps can be made
# with just one or two lines of code which are easier to comprehend.
# Got tired of remembering which library modules and what methods ;-)

import string 
import random
import copy

import numpy as np

import matplotlib as mpl
from matplotlib import colors     # update for v_1.2.0
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


def rgb_cmap_gradient(lowColor='k', highColor='w', name=None, mirrored=False) :
    """
    A linear-in-RGB-space Colormap.

    Parameters
    ----------
    locColor : color, optional, default: 'black'
        Color at the low end of the Colormap range.
    
    highColor : color, optional, default: 'white'
        Color at the high end of the Colormap range.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will be a string of 8 random
        characters.

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
        mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    return cmap

def _smoothHue(lowHue,hiHue,eps,numbSegs) :
    # return an array (lenght of numbSegs) of 'smoothed' hue values from low to high
    #-----------------------------------------------------
    def _transformHue(x,eps) :
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
        third, twothirds = 1/3, 2/3
        y = np.where( x >third , sgm(x,third), y )
        y = np.where( x >twothirds , sgm(x,twothirds), y )
        return y
    #-----------------------------------------------------

    xstart = _transformHue(np.mod(lowHue,1) ,1/eps)       #... inverse function
    xend = _transformHue(np.mod(hiHue,1),1/eps)           #... inverse function
    xstart = xstart if lowHue<1.0 else (xstart+1)
    xend = xend if hiHue<1.0 else (xend+1)
    y = np.linspace(xstart, xend, numbSegs)
    y = np.mod(y,1)
    return _transformHue(y,eps)

def hsv_cmap_gradient(lowHSV=[0,1,1], hiHSV=[1,1,1], name=None, mirrored=False, smooth=None) :
    """
    A linear-in-HSV-space Colormap.
    
    Using hSV or hSVA values for HSV color. SVA values are in the range [0,1].
    The h values are in the range [0,2] so that Hue is mod(h,1).

    Parameters
    ----------
    lowHSV : 3 or 4 array or string, optional, default: [0,1,1], , 'red'
        HSV, HSVA color at the low end of the Colormap range.  If a string,
        a named color (may be proceeded by a '+')
    
    hiHSV : 3 or 4 array or string, optional, default: [1,1,1] , '+red'
        HSV, HSVA color at the high end of the Colormap range.  If a string,
        a named color (may be proceeded by a '+')
    
    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, name is assigned a string of 8 random characters.

    mirrored : bool
        If True, colormap is divided into two linear
        segments with the lowHSV at the low and high
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
        if sVal[0] == '+':
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

    if isinstance(lowHSV, str ): lowHSV= stringToHSV(lowHSV)
    if isinstance(hiHSV, str ): hiHSV= stringToHSV(hiHSV)

    if not isinstance(lowHSV,(list, tuple, np.ndarray)) :
        raise ValueError('Invalid lowHSV argument: ' + str(lowHSV))
    if not isinstance(hiHSV,(list, tuple, np.ndarray)) :
        raise ValueError('Invalid highHSVarg argument: ' + str(hiHSV))  

    lowColor, highColor  = list(lowHSV), list(hiHSV)
    if len(lowColor) == 3 : lowColor.append(1.0)
    if len(highColor) == 3 : highColor.append(1.0)
    if len(lowColor) != 4 or len(highColor) != 4  :
        raise ValueError('Invalid HSVarg argument list length ')

    numbSegs = 256
    delta = np.subtract(highColor,lowColor)
    x = np.linspace(0.0,1.0,num=numbSegs)
    if smooth is not None :
        hH = _smoothHue(lowColor[0],highColor[0],smooth,numbSegs)

    clist = []
    for i in range( len(x)):
        n = x[i]
        h,s,v,a = np.add(lowColor, np.multiply(delta,n) )
        if smooth is None :
            h = np.mod(h, 1)
        else: h = hH[i]
        r,g,b = colors.hsv_to_rgb([h,s,v])     # update for v_1.2.0
        clist.append([r,g,b,a])  
    cmap = ListedColormap(clist)
    if mirrored :
        cmap = mirrored_cmap( cmap, name=name )
        return cmap # already registered
    
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name
        
    return cmap

def binary_cmap(negColor='b', posColor='r', name=None, bndry=None ) :
    """
    A two-color Colormap.

    Parameters
    ----------
    negColor : color, optional, default: 'blue'
        Color at the low end of the Colormap range.
    
    posColor : color, optional, default: 'red'
        Color at the high end of the Colormap range.
    
    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, name is assigned a string of 8 random characters.

    bndry : float, optional, default: 0.5
        The division between the negColor and posColor,
        with range from 0.03 to 0.97

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
    
    if bndry is None :
        cmap = ListedColormap([negColor,posColor])
    else :
        if (bndry<0.03 or bndry>0.97) :
            raise ValueError('Invalid bndry argument, found {}'.format(bndry))
        colorArr = np.array([colors.to_rgba(posColor)]*256)
        division = int(bndry*256)
        colorArr[:division, :] = np.array(colors.to_rgba(negColor))
        cmap = ListedColormap(colorArr)

    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name
    return cmap

def hue_cmap(lowHue=0, hiHue=None, smooth=1.6, name=None ) :
    """
    A 'smooth-HSV' colormap.
    
    Non linear adjustment of Hue in a HSV colormap.
    For the HSV colors, S and V are 1.
    The h values are shifted based on the smoothing parameter.

    Parameters
    ----------
    lowHue : float or string, optional, default: 0
        Hue at the low end of the Colormap range.
        Values are in the range [0,2].  
        If a string, the hue of a named color (may be
        proceeded by a '+')

    hiHue : float or string, optional, default: 1
        Hue at the upper end of the Colormap range.  
        Values are in the range [0,2].  
        If a string, the hue of a named color (may be
        proceeded by a '+')

    smooth: float, optional, default: 1.6
        Controls the ratio of the amount of CMY to RGB color.
        CMY is extended relative to the RGB hues for values
        greater than one (relative to standard HSV colormap).
        RGB is extended relative to the CMY hues for values
        less than one.  Range is [.1,10]

    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, name is assigned a string of 8 random characters.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.       
    """
    # ............................................................
    def stringToHue(sVal) :
        addOne = False
        if sVal[0] == '+':
            addOne = True
            sVal = sVal[1:]
        h,_,_ = colors.rgb_to_hsv(colors.to_rgb(sVal))
        if addOne : h += 1
        return h    
    # ............................................................
    if (smooth<0.1 or smooth>10) :
        raise ValueError('smooth must be between 0.1 and 10. , found {}'.format(smooth))

    if isinstance(lowHue, str ): lowHue= stringToHue(lowHue)
    if hiHue is None:
        hiHue = lowHue+1
    else:
        if isinstance(hiHue, str ):  hiHue= stringToHue(hiHue)

    rng = 0.001
    oneLow,oneHi = 1.0-rng , 1.0+rng
    isOne = smooth>oneLow and smooth<oneHi

    if isOne :
        delta = np.subtract(hiHue,lowHue)
        x = np.linspace(0, 1, 256)
        h = np.add(lowHue, np.multiply(delta,x) )
        h = np.mod(h,1)
    else :
        h = _smoothHue(lowHue,hiHue,smooth,256)

    one = np.ones(len(h))
    hsv = np.transpose([h,one,one])
    #rgb = cm.mpl.colors.hsv_to_rgb(hsv)
    rgb = colors.hsv_to_rgb(hsv)     # update for v_1.2.0
    cmap = colors.ListedColormap(rgb)

    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name

    return cmap

def stitch_cmap( *maps, bndry=None, name=None) :
    """
    Colormap composed of multiple colormaps.

    Parameters
    ----------
    maps : str or Colormap, optional
        A Colormap instance or registered colormap name
    
    bndry : scalar or list, default: None
        Boundary between colormaps.  If  a list, 
        the number of values must be one less than
        the number of maps.  The list values must
        be in assending order.  If None, input colormaps
        are evenly spaced in the returned colormap.
        The bndry values are in the range 0+ to 1-.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will be assigned a string of 8 random
        characters.
       
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
        # allow single bndry to be input as a float
        if isinstance(bndry, float) : bndry = [bndry]
        numbBndry = len(bndry)
        if numbBndry != (numbMaps-1) :
            strgMsg = 'Number of boundaries, {}, is NOT equal to one less that the number of maps, {}. '
            raise ValueError(strgMsg.format(numbBndry,numbMaps))
    bndry.append(1.0)

    cmapList = []
    for cmap in maps :
        if isinstance(cmap,str) : cmap = mpl.colormaps[cmap]     # update for v_1.2.0
        cmapList.append(cmap)

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
        clist.append( cmapList[mapIndex](index) )
    cmap = ListedColormap(clist)
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name
    return cmap

def stitch_color(*colorArr, bndry=None, name=None ) :
    """
    Colormap composed of multiple colors.

    Parameters
    ----------
    *colorArr : colors
    
    bndry : scalar or list, default: None
        Boundary between colors.  If  a list, 
        the number of values must be one less than
        the number of maps.  The list values must
        be in assending order.  If None, input colors
        are evenly spaced in the returned colormap.
        The bndry values are in the range 0+ to 1-.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will assigned be a string of 8 random
        characters.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """
    
    numbColors = len(colorArr)
    if numbColors == 1 : 
        raise ValueError('Only one cmap entered, need two or greater')  
    if bndry is None :
        bndry = [ i/numbColors for i in range(1,numbColors) ]
    else :
        # allow single bndry to be input as a float
        if isinstance(bndry, float) : bndry = [bndry]
        numbBndry = len(bndry)
        if numbBndry != (numbColors-1) :
            strgMsg = 'Number of boundaries, {}, is NOT equal to one less that the number of maps, {}. '
            raise ValueError(strgMsg.format(numbBndry,numbColors))
    
    RGBAvals = [ colors.to_rgba(c) for c in colorArr ]
    RGBAvals = np.array(RGBAvals)

    N=256
    bIndx = [ int(N*b) for b in bndry ]
    RGBA = np.tile(RGBAvals[0],(N,1))
    for i,bI in enumerate(bIndx) : RGBA[bI:] = RGBAvals[i+1]
    cmap = ListedColormap(RGBA)
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name

    return cmap

def alpha_cmap( cmap, alpha, constant=False, name=None )  :
    """
    Set the Colormap alpha channel.

    Parameters
    ----------
    cmap : str or Colormap, optional
        A Colormap instance or registered colormap name
    
    alpha : scalar
        Value set for the color alpha channel. 
        The alpha values are in the range 0+ to 1-.

    constant : bool { True, False }, optional, False
        If False, alpha channel values are multiplied by
        alpha.  If True, all alpha channels
        are assigned to a constant alpha.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name assigned the colormap name
        with '_a' at the end characters.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
    """

    if alpha<0.0 or alpha>1.0 :
        raise ValueError('alpha must be between 0 and 1, found {}'.format(alpha))       
    
    if isinstance(cmap,str) : 
        cmap = mpl.colormaps[cmap]      # update for v_1.2.0
    if name is None :
        name = cmap.name + '_a'

    if isinstance(cmap,LinearSegmentedColormap) :
        # LinearSegmentedColormap ...
        cdict = copy.copy(cmap._segmentdata)
        cdict['alpha'] = [ (0.0, alpha,alpha),(1.0, alpha,alpha) ]
        newCmap = LinearSegmentedColormap(name,cdict)
        mpl.colormaps.register(newCmap,name=name)     # update for v_1.2.0

    else :
        # ListedColormap ............
        colorArr = copy.copy(cmap.colors)
        if not constant : 
            colorArr = np.array(colorArr)
            if colorArr.shape[1] == 4 : alpha = alpha*colorArr[:,3]
        colorArr = np.insert(colorArr,3,alpha,axis=1)[:,:4]
        newCmap = ListedColormap(colorArr,name)
        mpl.colormaps.register(newCmap,name=name)     # update for v_1.2.0
    return newCmap

def mirrored_cmap( cmap, name=None, rev=False) :
    """ A mirrored colormap.
    
    Parameters
    ----------
    cmap : str or Colormap, optional
        A Colormap instance or registered colormap name.
    
    name : str, optional
        The registered name to identify the colormap.
        If None, '_m' will be appended to the colormap name.
    
    rev : boolean {True, False}, default: False
        If True, the reversed colormap will be used.

    """  
    if isinstance(cmap,str) : 
        cmap = mpl.colormaps[cmap]     # update for v_1.2.0
    if name is None :
        name = cmap.name + '_m'
        if rev : name += 'r'

    numbSegs = 256
    x = np.linspace(0.0,1.0,num=numbSegs)
    n = 1-abs( 2.0*x - 1.0 )
    if rev : n = 1-n
    clist =cmap(n)
    cmap = ListedColormap(clist)
   
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name
    return cmap

def reversed_cmap( cmap, name=None ) :
    """
    A reversed colormap.
    
    Parameters
    ----------
    cmap : str or Colormap, optional
        A Colormap instance or registered colormap name
    
    name : str, optional
        The registered name to identify the colormap.
        If None, '_r' will be appended to the colormap name.
    
    """

    # for convenience to eliminate need for extra import.
    # Basically, same as Matplot colors.Colormap .reversed() method
    
    current_name = None

    if isinstance(cmap,str) : 
        current_name = cmap
        cmap = mpl.colormaps[cmap]     # update for v_1.2.0
    else:
        current_name = cmap.name

    # determine if revered cmap already exists.
    new_name = current_name + '_r'
    if name is not None : new_name = name
    try:
        cmap = mpl.colormaps[new_name]
        return cmap
    except : pass

    cmap = cmap.reversed()
    mpl.colormaps.register(cmap,name=new_name)     # update for v_1.2.0
    cmap.name = new_name
    return cmap

def op_cmap(operation,rgb=True,name=None):
    """
    A Colormap defined by a function argument.

    Parameters
    ----------
    operation : function object
        Function that takes one argument,
        a Numpy array of float values in the
        range from 0 to 1.
        The function returns a 3xN color value.

    rgb : bool {True, False}, optional, default: True
        By default, RGB color values are returned by the
        operation function.  If set False, the operation
        returns HSV color values.
    
    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, the operatonal function name is
        assigned to the cmap name if not a lambda
        function.  If a lamddda function, the name is
        assigned a string of 8 random characters.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
   
    """
    # ----------------------------------------
    def getFunctionName(op,xName=None) :
        # following is a copy from surface.py
        if xName is not None : return xName
        fobar = lambda x : x
        if op.__name__ != fobar.__name__ :
            name = op.__name__
        else :
            name = None
        return name
    # ----------------------------------------
    name = getFunctionName(operation,name)
    numbSegs = 256
    x = np.linspace(0.0,1.0,num=numbSegs)
    clist = operation(x)    
    clist = np.array(clist).T
    hasAlpha = clist.shape[1] == 4
    if not rgb :
        if hasAlpha :
            alp = clist[:,[3]]
            clist = clist[:,[0,1,2]]
            #clist = cm.colors.hsv_to_rgb(clist)
            clist = colors.hsv_to_rgb(clist)     # update for v_1.2.0
            clist = np.append(clist,alp,1)
        else :
            #clist = cm.colors.hsv_to_rgb(clist)      
            clist = colors.hsv_to_rgb(clist)     # update for v_1.2.0
    
    cmap = ListedColormap(clist)
    if name is None :
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    mpl.colormaps.register(cmap,name=name)     # update for v_1.2.0
    cmap.name = name
    
    return cmap

def op_alpha_cmap(cmap,operation,fit=True,name=None):
    """
    A Colormap modified by functional values of alpha.

    Parameters
    ----------
    cmap : str or Colormap, optional
        A Colormap instance or registered colormap name

    operation : function object
        Function that takes one argument,
        a Numpy array of float values in the
        range from 0 to 1.
        The function returns a Numpy array
        of float values in the
        range from 0 to 1.
    
    name : str, optional, default: None
        The registered name to identify the colormap.
        If it's None, the map is not registered if
        a lambda function is the operation argument,
        otherwise, the operatonal function name is
        assigned to the cmap name.

    fit : bool, optional, default: True
        If True, the values return by the operation
        function will be fit to the range 0 to 1.

    Returns
    -------
    ListedColormap
        An instance of a colormap.
   
    """

    # ----------------------------------------
    def getFunctionName(op,xName=None) :
        # following is a copy from surface.py
        if xName is not None : return xName
        fobar = lambda x : x
        if op.__name__ != fobar.__name__ :
            name = op.__name__
        else :
            name = None
        return name
    # ----------------------------------------
    def fitBounds(x) :
        mV,xV = np.amin(x),np.amax(x)
        return (x - mV)/(xV - mV)
    # ----------------------------------------
    name = getFunctionName(operation,name)
    if isinstance(cmap,str) : 
        cmap = mpl.colormaps[cmap]    # update for v_1.2.0
    if name is None :
        name = cmap.name + '_a'

    numbSegs = 256
    x = np.linspace(0.0,1.0,num=numbSegs)
    alphalist = operation(x)    
    if fit : alphalist = fitBounds(alphalist)
    colorList = cmap(x).T
    colorList[3:] = alphalist
    newCmap = ListedColormap(colorList.T, name=name)
    mpl.colormaps.register(newCmap,name=name)     # update for v_1.2.0
    return newCmap

def section_cmap(cmap,lowIndx,hiIndx,name=None) :
    """
    A Colormap from a section of a colormap.

    Parameters
    ----------
    cmap : str or Colormap, optional
        A Colormap instance or registered colormap name

    lowIndx : float
        Value of the start of the input colormap.
        Range is 0 to 0.97.

    hiIndx : float
        Value of the end of the input colormap.
        Range is 0.03 to 1.0.
   
    name : str, optional
        The registered name to identify the colormap.
        If None, the name will be a string of random
        characters and the map is not registered.
       
    Returns
    -------
    ListedColormap
        An instance of a colormap.
   
    """

    # ------------------------------------------------
    def opFunc(t,cmap) :      # update for v_1.2.0
        return mpl.colormaps.get_cmap(cmap)(lowIndx + (hiIndx-lowIndx)*t).T     # update for v_1.2.0
    # ------------------------------------------------

    if ( lowIndx <0.0 or lowIndx >0.97) :
        raise ValueError('Error: required that 0.0 <= lowIndx <= 0.97. , found {}'.format(lowIndx))
    if ( hiIndx <0.03 or hiIndx >1.0) :
        raise ValueError('Error: required that 0.03 <= hiIndx <= 1.0. , found {}'.format(hiIndx))
    indxRng = hiIndx-lowIndx
    if indxRng < 0.05 :    
        raise ValueError('Error: range of lowIndx to hiIndex is > 0.1 , found {}'.format(indxRng))

    if name is None : 
        name = ''.join(random.choices(string.ascii_uppercase , k = 8))
    if isinstance(cmap,str) : 
        cmap = mpl.colormaps[cmap]     # update for v_1.2.0
    op_lFunc = lambda t : opFunc(t,cmap)        # update for v_1.2.0 
    return op_cmap(op_lFunc,name=name)     # update for v_1.2.0


class DualCmap() :
    
    def __init__(self,xcmap,ycmap,kind='sum',norm=True,name=None) :
        """
        A 2D colormap.

        Parameters
        ----------
        xcmap : colormap or string of a colormap name.

        ycmap : colormap or string of a colormap name.

        kind : string {'sum','ave','ave2','srt'} or number, default: 'sum'
            Method of combining colormaps.

        norm : bool, optional, default True
            Normalize array values in __call__

        name : string,optional.
            name to identify the dual colormap.
        """
        # FutDev: possible blend color modes?
        #         Need documentation explaination for 'kind'.

        self._name = name if name is not None else ''
        if isinstance(xcmap,str) : xcmap = mpl.colormaps[xcmap]     # update for v_1.2.0
        if isinstance(ycmap,str) : ycmap = mpl.colormaps[ycmap]     # update for v_1.2.0

        self.xmap = xcmap
        self.ymap = ycmap
        self._kind = kind
        self._norm = norm
        self.pow = 0.5 
        if isinstance(kind,(int,float)): self.distfact(kind)
        return
    
    def __call__(self,x,y):
        """
        RGBA color values from x,y positions.

        Parameters
        ----------
        x,y : arrays of N floats with values in the 
            interval [0,1].  Normalization will occur
            if 'norm' is set to True in the constructor.
            Both arrays must be of the same length, N.

        Returns
        -------
        RGBA values with N shape (4,N)
        """
        if self._norm :
            x = (x-np.min(x))/(np.max(x)-np.min(x))
            y = (y-np.min(y))/(np.max(y)-np.min(y))
        color = np.add(self.xmap(x),self.ymap(y))
        if self._kind == 'ave' : color /= 2.0
        color = np.clip(color,0.0,1.0)
        if self._kind == 'ave2' : color /= 2.0
        if self._kind == 'srt' : color = np.sqrt(color)
        if self._kind == 'pwr' : color = np.power(color,self.pow)
        return color.T
    
    def __str__(self) :
        """DualCmap object string representation."""

        name = self.__class__.__name__
        kind = self._kind
        if kind == 'pwr' :
            kind = '{}, {:.1f}'.format(kind,self.pow)
        bs = ' ( {} )'.format(kind)    
        sz = ': {}, {}'.format(self.xmap.name,self.ymap.name)
        val = name+bs+sz
        return val

    @property
    def name(self) :
        """Descriptive identifier for the DualCmap object."""
        return self._name

    @name.setter
    def name(self,val) :
        self._name = str(val)
        return

    @property
    def norm(self) :
        """Boolean indicator for value normalization""" 
        return self._norm
    
    @norm.setter
    def norm(self,val) :
        self._norm = val
        return

    def distfact( self, p ) :
        """
        Set distribution factor of resulting color values.

        Parameters
        ----------
        p : float in the range [0.2, 4]
            For a value of 1, there is no redistribution.
            A value of 0.5 is equivalent to using the
            constructor argument of kind='srt'. 

        Returns
        -------
            self : DualCmap object 
        """
        if (p<0.2) or (p>4.0) :
            raise ValueError('distfact {} must be between 0.2 and 4.0'.format(str(p)))
        self.pow = p
        self._kind = 'pwr'
        return self


# =================================================================================+


def linsegfunc(v, arr) :
    """
    Linear segment function, y = f(x) with domain/range [0,1]

    Parameters
    ----------
    v : float or array
        Independent variable in the range 0 to 1.

    arr : N X 2 array
        An array of N values (x,y), with increasing values of x and, 
        y in the range 0 to 1.

    Returns
    -------
    float
        Dependent values or array of values
    
    """
    def recurs(L,R,x,V) :   # efficient enough for intended scope
        T = (L+R)//2 
        if L==(R-1) : return L
        if x[T] < V : return recurs(T,R,x,V)
        else :        return recurs(L,T,x,V)
    x,y = np.array(arr).T
    m = np.min(x)
    d = np.max(x)-np.min(x)
    x,y = (x-m)/d , np.clip(y,0.0,1.0)
    R = len(x) - 1
    slp, b = [None]*R , [None]*R
    for i in range(R) :
        j = i+1
        slp[i] = ( y[j]-y[i] ) / ( x[j]-x[i] )
        b[i]   = y[i]-slp[i]*x[i]
    if isinstance(v,(int,float)) : v=[v]
    vals = [None]*len(v)
    for i,V in enumerate(v) :
        m = recurs(0,R,x,V)
        vals[i] = slp[m]*V+ b[m]
    return vals



# =================================================================================+
# FutDev : would be nice to have a function to create a 2D colorbar which could    |
#          be added to the Matplotlib 3D Axes.                                     |
# =================================================================================+
