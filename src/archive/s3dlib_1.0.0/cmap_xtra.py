# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

'''
Auxiliary functions for use with Matplotlib.

'''
import numpy as np

from matplotlib import cm, colors
from matplotlib.colors import ListedColormap

from colorspacious import cspace_converter


def Lab_cmap_gradient(lowColor='k', highColor='w', name=None, mirrored=False) :
    """
    A linear-in-Lab-space Colormap, with option of registering the map.

    Parameters
    ----------
    locColor : RGB color, optional, default: 'black'
        Color at the low end of the Colormap range.
    
    highColor : RGB color, optional, default: 'white'
        Color at the high end of the Colormap range.
    
    name : str, optional
        The registered name to identify the colormap.
        If it's None, the name will be a string of random
        characters and the map is not registered.

    mirrored : bool
        If True, colormap is divided into two linear
        segments with the lowColor at the low and high
        values, the highColor in the middle.
        
    Returns
    -------
    LinearSegmentedColormap
        An instance of a colormap.
    """
    
    if name is not None :
        if not isinstance(name, str ):
            raise ValueError('Invalid name argument (str required): ' + str(name))        

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
    if name is not None :
        cm.register_cmap(name,cmap)
        # name property not assigned for ListedColormap, so...
        cmap.name = name
    return cmap

