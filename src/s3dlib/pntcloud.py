# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

'''
Module containing the Point3DCloud class.
'''
import warnings
import string 
import random
import numpy as np

import s3dlib.prvt_func as pvt
from s3dlib.surface import Surface3DCollection

import matplotlib as mpl
from matplotlib import cm,colors,rcParams

from scipy.interpolate import RBFInterpolator, interpn
from skimage import measure

# interpn: determine value(s) within a cloud of values     
# RBFInterpolator: determine cloud values from sample values
# measure: for marching_cubes, surface of constant value in a cloud.


_DFT_COLOR = rcParams['lines.color'] # default, surface/line color
_DFT_CMAP  = rcParams['image.cmap']  # default, colormap

_DFT_PDG = 1.0   # single sample data distribution.


class Point3DCloud() :

    def _set_geom(self,dmn) :
        """ reset the class geometric properties """
        Np1 = self._Np1
        rngX,rngY,rngZ = dmn[0],dmn[1],dmn[2]
        x = np.linspace(rngX[0],rngX[1],Np1)
        y = np.linspace(rngY[0],rngY[1],Np1)
        z = np.linspace(rngZ[0],rngZ[1],Np1)
        pnts = np.mgrid[ rngX[0]:rngX[1]:Np1*1j,
                         rngY[0]:rngY[1]:Np1*1j,
                         rngZ[0]:rngZ[1]:Np1*1j ]
        self._domain = np.array(dmn)
        self._bounds['xlim'] = dmn[0]
        self._bounds['ylim'] = dmn[1]
        self._bounds['zlim'] = dmn[2]
        self.cloud_points = pnts
        self.axis_points = (x,y,z)
        return None

    def _initColor(self,color) :
        if isinstance(color,(list, tuple, np.ndarray)) :
            if len(color) == 3 : color.append(1.0)
        elif isinstance(color,str) :
            color = colors.to_rgba(color) 
        else :    
            raise ValueError('Invalid color argument')
        self._singleColor = color      
        cArr = np.zeros( (self._numPts, 4) )
        cArr[:] = color
        self._colors =cArr
        return None

    def _domain_contains_xyz(self,coor) :
        dmn = self._domain
        low = np.prod(coor>=dmn[:,0],axis=1)
        hi  = np.prod(coor<=dmn[:,1],axis=1)
        sSize = np.sum( low*hi )
        return sSize == len(coor)

    def __init__(self,rez,name=None, **kwargs ) :
        """
        Create a rectangular point cloud in xyz coordinates.

        Parameters
        ----------
        rez : number, optional, default: 1.0
            Number of subdivisions ( 10*rez ) along the 3 coordinate axes.
            Rez values range from 0.4 to 10.0.
            
        name : string, optional, default: None.
            Descriptive identifier for the point cloud.

        Other Parameters
        ----------------

        domain : a number, list or array, default: 1
            The domain of the cloud.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.

        color : str, optional

        size : float, optional
            Marker size of cloud points.

        marker : MarkerStyle, default : 'o'.

        """
        kwargDft = { 'color':_DFT_COLOR, 'size':10, 'marker':'.', 'cmap':_DFT_CMAP,
                     'domain': 1 }
        KW = pvt.KWprocessor(kwargDft,'Point3DCloud',**kwargs)

        self.cloudName = name
        self.cloudSource = 'UnMapped'

        if rez<0.4 or rez>10 :
            raise ValueError("Error: incorrect value 'rez' passed to pntcloud:", rez)

        self._rez = rez
        Np1 = int( 10*rez )+1
        self._Np1 = Np1
        self._numPts = Np1*Np1*Np1
        # ...........
        self._domain = None
        self._bounds = {}
        self.cloud_points = None
        self.axis_points = None
        self.domain( KW.getVal('domain') )
        # ...........
        self._bounds['vlim'] = [1,1]       
        self._values = np.ones( (Np1,Np1,Np1) )
        self._sample_size = 0
        self._pdg = None  # only defined if mapped sample density

        self._cmap = KW.getVal('cmap')
        self._pntSize = KW.getVal('size')
        self._marker = KW.getVal('marker')
        #.. following only used to create a cBar_ScalarMappable if
        #   no color mapping is called, but needed if a color bar is displayed.
        #   See the cBar_ScalarMappable @property.
        self._singleColor = None
        self._colors = None
        self._initColor( KW.getVal('color') )
        self._scmpb = None # only defined if color mapped
        # FutDev: self.hasColor and self.isUniformColor ???
        #         could be more efficient when mapping color to a surface.
        self.hasVals = False
        self.isColorMapped = False
        self.geomChanged = False
        return

    def __str__(self) :
        """ String representation of the point cloud object """
        name = self.__class__.__name__
        numPoints = self._numPts
        bs = ' ( {}, {} )'.format(self._rez,self.cloudSource)
        sz = ':  points: {}'.format(numPoints)
        return name + bs + sz

    def __call__(self,xyz) :
        """ Values at coordinates in the cloud domain.

        Parameters
        ----------
        xyz : V x 3 float array-like
            Coordinates within the domain.
        
        Returns
        -------
        V float array

        """
        if not self._domain_contains_xyz(xyz) :
            raise ValueError('Found coordinates external to the cloud domain')
        apoints = np.array(self.axis_points)
        values = self._values
        result = interpn(apoints, values, xyz)
        return result

    def get_colors_at_xyz(self,xyz) :
        """ Colors at coordinates in the cloud domain.

        Parameters
        ----------
        xyz : V x 3 float array-like
            Coordinates within the domain.
        
        Returns
        -------
        V X 4 float array

        """
        # FutDev:  if uniform color, just return an array of uniform colors.
        if not self._domain_contains_xyz(xyz) :
            raise ValueError('Found coordinates external to the cloud domain')
        apoints = np.array(self.axis_points)
        if self._colors.shape[1] == 4 :  # shouldn't be needed ...but MROI
            R,G,B,A = (self._colors.T)
        else :
            R,G,B = (self._colors.T)
        N = self._Np1
        r = interpn(apoints, R.reshape(N,N,N), xyz)
        g = interpn(apoints, G.reshape(N,N,N), xyz)
        b = interpn(apoints, B.reshape(N,N,N), xyz)
        a = np.ones(r.shape)    # << note: alpha = 1 for all colors.
        rgba = np.array( [r,g,b,a]).T
        rgba = np.clip(rgba,0.0,1.0)
        return rgba

    @property
    def name(self) :
        """Descriptive identifier for the cloud."""
        temp = '' if self.cloudName is None else self.cloudName
        return temp

    @name.setter
    def name(self,val) :
        if val is not None : self.cloudName = val
        return

    @property
    def points(self) :
        """ A (3,N) array of point cloud coordinates. """
        pnts = np.reshape(self.cloud_points,(3,-1))
        return pnts

    @property
    def point_colors(self) :
        """ A (N,4) array of point cloud colors """
        return self._colors

    @property
    def cBar_ScalarMappable(self) :
        """matplotlib.cm.ScalarMappable object for color mapped point cloud values."""
        if self._scmpb is None :  # only needed if no colormapping method is called.
            clist = [ self._singleColor, self._singleColor ]
            cmName = ''.join(random.choices(string.ascii_uppercase , k = 8))
            cmap = colors.LinearSegmentedColormap.from_list(cmName, clist )
            norm = colors.Normalize(*self._bounds['vlim'] )
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            return sm  # a uniform color colormap.
        return self._scmpb

    @property
    def bounds(self) :
        """
        Dictionary of coordinate and value ranges.
        
        Each dictionary value is a 2 float array of minimum
        and maximum values.  Keys are:

        'xlim' : x-coordinate

        'ylim' : y-coordinate

        'zlim' : z-coordinate

        'vlim' : value assignments.
            
        """
        return self._bounds

    @property
    def sample_size(self) :
        """Number of data samples used for value mapping."""
        return self._sample_size

    @property
    def val_source(self) :
        """String indicating mapping type for cloud values.""" 
        return self.cloudSource

    @property
    def density_pdg(self):
        """ pdg for sample density cloud"""
        return self._pdg
        
    def domain(self, xlim, ylim=None, zlim=None) :
        """
        Reset the domain of the point cloud object.
        Initial domain is set as [ [-1,1], [-1,1], [-1,1] ]

        Parameters
        ----------
        xlim : a number, list or array, default: 1
            The domain of the cloud.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.
            If 2D array, ylim and zlim are ignored, otherwize,
            must be a number or 1D array if ylim or zlim are None.

        ylim, zlim : a number or 2D array, default: None
            Set the ylim and/or zlim.  If None, set to xlim.

        Returns
        -------
        self : Point3DCloud object

        """
        domn = pvt._interpret_domain_args(xlim,ylim,zlim)
        self._set_geom(domn)

        return self

    def get_domain(self) :
        """ 
        domain of the point cloud.
        Usefull for getting the Matplotlib 3D axis xyz limits
        when samples are used to set the domain. (dreset=True)
        """
        return self._domain

    def get_color_for_val(self,N=1) :
        """
        Arrays of cloud values, normalized values and associated colors.
        Useful for creating Matplotlib legend patch handles
        and setting focus for the add_to3d method.  Colors are selected
        from the mapped cloud colormap.

        Parameters
        ----------
        N : integer or N array of floats, default: 1 
            For a integer, number of evenly spaced values
            and colors within the cloud.  For an array, used to 
            determine the colors for a list of values within
            the cloud.  Colors are determined from the cloud colormap.

        Returns
        -------
        carr : N x 4 array of colors
        varr : N float array of values
        norm : N float array of normalized values
        """
        # FutDev: use generically, passing in values and colormap
        #         Then place in prvt_func module for use with
        #         surfaces and lines.
        minVal, maxVal = np.min(self._values), np.max(self._values)
        vRng = maxVal - minVal

        if isinstance(N,int) :
            delta = (maxVal - minVal)/(N + 1)
            varr = [ minVal + (i+1)*delta for i in range(N)]
        else :
            varr = N

        if vRng > 0 : # PICNIC check when values have not been set.
            norm = (varr-minVal)/vRng
            cmap = self._getCmap()
            carr = cmap( norm )
        else :
            norm = [1.0]*len(varr)
            carr = [self._singleColor]*len(varr)

        return carr,varr,norm

    # ==================================================
    # Properties/Methods
    # ==================================================

    def _set_matching_sample_domain(self,coor,vals,dreset) :
        """ Samples must be contained within the domain """

        # used by cloud & density mapping from samples...
        # from Scipy: 'At least 4 data points are required'
        maxNo=5
        # .. compare cloud domain with sample domain   
        if dreset :  #.. reset the cloud domain based on data ranges.
            self._sample_size = len(coor)
            dmn = np.array([np.min(coor,axis=0),np.max(coor,axis=0)]).T
            self._set_geom(dmn)
        else :       #.. remove any data point outside of the cloud domain.
            dmn = self._domain
            low = np.prod(coor>dmn[:,0],axis=1)
            hi  = np.prod(coor<dmn[:,1],axis=1)
            sSize = np.sum( low*hi )
            self._sample_size = sSize
            if sSize < maxNo :
                raise ValueError( 'Insufficient number of samples in domain, {} < {}'.format(str(sSize),str(maxNo)))
            inDomain = (low*hi).astype(bool)
            coor = coor[inDomain]
            if vals is not None : vals = vals[inDomain]
        return coor,vals

    def _vals_from_samples(self,coor,vals) :
        # normalize domain & samples in [0,1]
        domain = self._domain
        drng = domain[:,1]-domain[:,0]
        ncoor = np.divide(coor-domain[:,0],drng)
        N = self._Np1
        Nj =N*1j
        xyz = np.mgrid[ 0:1:Nj, 0:1:Nj,0:1:Nj ]    
        xyzflat = xyz.reshape(3, -1).T
        wflat = RBFInterpolator(ncoor, vals)(xyzflat)
        wgrid = wflat.reshape(N, N, N)
        return np.array(wgrid)

    def _density_from_samples(self,coor,pdg) :
        dataPts = coor
        dm = self._domain
        numbPts = len(coor)
        # weight, for future development
        weight = np.ones( dataPts.shape[0])
        A = dm.T[1] - dm.T[0]
        B = np.array( [ dm[0][0], dm[1][0], dm[2][0] ])
        C = np.empty([len(dataPts),3])
        C[:] = A
        norm_coor = (dataPts-B)/C
        Nj = self._Np1*1j
        xyz = np.mgrid[0:1:Nj,0:1:Nj, 0:1:Nj ]
        dens = np.zeros(xyz.T.shape[:3])
        hfact = 1/3  # heuristic multiplier for sensitivity distance
        fden = np.sqrt(2)/np.power(2*numbPts,1/3)
        pdg = hfact*pdg*fden
        for i,org in enumerate(norm_coor) :
            diff = (xyz.T - org)
            R = np.linalg.norm(diff,axis=3)
            dens = dens + weight[i]*np.exp(-(R/pdg)**2)
        # normalize density.
        dmin,dmax = np.min(dens) , np.max(dens)
        dens = (dens-dmin)/(dmax-dmin)
        return dens.T

    def map_geom_from_op(self,operation,name=None) :
        """
        Functional transformation of point cloud coordinates.
        
        Parameters
        ----------
        operation : function object
            Coordinate mapping function that takes one
            argument, a 3xN Numpy array of native coordinates.
            The function returns a 3xN array of coordinates.

        name : string, optional, default: None or op function name.
            Descriptive identifier for the geometry.

        Returns
        -------
        self : surface object

        """

        pnts = np.reshape(self.cloud_points,(3,-1))
        cps = self.cloud_points.shape
        newCoor = operation(pnts)
        self.cloud_points = np.reshape(newCoor,cps)
        self.geomChanged = True
        self._bounds['xlim'] = [ np.min(newCoor[0]), np.max(newCoor[0]) ]
        self._bounds['ylim'] = [ np.min(newCoor[1]), np.max(newCoor[1]) ]
        self._bounds['zlim'] = [ np.min(newCoor[2]), np.max(newCoor[2]) ]
        return self

    def map_vals_from_op(self,operation,cmap=None, name=None) :
        """
        Functional assignment of values within the point cloud.

        Parameters
        ----------
        operation : function object
            Coordinate mapping function that takes one
            argument, a 3xN Numpy array of native coordinates.
            The function returns an array of scalar values

        cmap : str or Colormap, optional, default: None
            A Colormap instance or registered colormap name.
            If not assigned (None) , colormapping will NOT be applied.

        name : string, optional, default: None or op function name.
            Descriptive identifier for the cloud values.

        Returns
        -------
        self : point cloud object

        """
        values = operation(self.cloud_points)
        self._values = values
        self._bounds['vlim'] = [ np.min(values), np.max(values) ]
        self.hasVals = True
        self.cloudSource = 'function'
        self._sample_size = 0
        if cmap is None : 
            self._scmpb = None # only defined if color mapped
        else :
            self.map_cmap_from_cloudvals(cmap)
        if name is None : 
            if self.cloudName is None : 
                self.cloudName = pvt._getFunctionName(operation)
        else : self.cloudName = name
        return self

    def map_vals_from_sampvals(self,samples,vadj=None,dreset=False,cmap=None,name=None ) :
        """
        Assign cloud values from data samples values.

        Parameters
        ----------
        samples : V X 4, float array-like 
            An array of V number of xyz vertex coordinates and
            scaler values which are assigned to each vertex.

        vadj : str {'s','c'}, optional, default: None
            Cloud value adjustment to sample value range.
            For 's', the cloud values are scaled to the minimum
            and maximum of the sample values.
            For 'c', the cloud values are clipped at the minimum
            and maximum of the sample values.

        dreset : bool { True, False }, optional, default: False
            Domain reset based on sample distribution.
            If True, the point cloud DOMAIN IS RESET
            to the bounds of the samples.
            If False, the out of the domain SAMPLES ARE EXCLUDED.

        cmap : str or Colormap, optional, default: None
            A Colormap instance or registered colormap name.
            If not assigned (None) , colormapping will NOT be applied.

        name : string, optional, default: None
            Descriptive identifier for the cloud values.

        Returns
        -------
        self : point cloud object

        """
        # ......................................................
        def adjustValues(values,samples,adj) :
            if adj is not None :
                smin, smax = np.min(samples), np.max(samples)
                key = adj[0].upper()
                if key == 'S' :    # scaled to the data
                    vmin, vmax = np.min(values), np.max(values)
                    ds, dv = smax-smin , vmax-vmin
                    svM = ds/dv
                    svB = ( vmax*smin - vmin*smax )/dv
                    values = svM*values + svB   # Y = mX + b
                if key == 'C' :    # clipped to the data
                    values = values.clip(smin,smax)
            return values
        # ......................................................
        
        coor, vals = pvt._extract_vals_from_inputCoor(samples)
        if vals is None:
            self.hasVals = False
            raise ValueError('Samples do NOT contain values for cloud mapping.' )
        coor, vals = self._set_matching_sample_domain(coor,vals,dreset)
        #.. assign cloud values from sample data
        values = self._vals_from_samples(coor,vals)
        values = adjustValues(values,vals,vadj)
        self.hasVals = True
        self._values = values
        self._bounds['vlim'] = [ np.min(values), np.max(values) ]
        self.cloudSource = 'sampleVals[{}]'.format(self._sample_size)

        if cmap is None : 
            self._scmpb = None # only defined if color mapped
        else :
            self.map_cmap_from_cloudvals(cmap)

        if name is not None : self.cloudName = name

        return self

    def map_vals_from_sampdens(self,samples,pdg=None,dreset=False,cmap=None, name=None) :
        """
        Assign cloud values from normalized density of a sample coordinate distribution.

        Parameters
        ----------
        samples : V X 3 or V X 4, float array-like 
            An array of V number of xyz vertex coordinates. If a V x 4 array
            is passed, only the first 3 are used for the coordinates.

        pdg : number, default : 1.0
            Sample coordinate distribution, normalized within the domain
            and sample size. Values range .5 to 5.0.

        dreset : bool { True, False }, optional, default: False
            Domain reset based on sample distribution.
            If True, the point cloud DOMAIN IS RESET
            to the bounds of the samples.
            If False, the out of the domain SAMPLES ARE EXCLUDED.

        cmap : str or Colormap, optional, default: None
            A Colormap instance or registered colormap name.
            If not assigned (None) , colormapping will NOT be applied.

        name : string, optional, default: None
            Descriptive identifier for the cloud values.

        Returns
        -------
        self : point cloud object

        """
        coor, vals = pvt._extract_vals_from_inputCoor(samples)
        coor, vals = self._set_matching_sample_domain(coor,vals,dreset)
        if pdg is None : pdg = _DFT_PDG
        self._pdg = pdg
        dens = self._density_from_samples(coor,self._pdg)
        self._bounds['vlim'] = [ np.min(dens), np.max(dens) ]
        self._values = dens
        self.hasVals = True
        self.cloudSource = 'sampleDens[{}, {:.2f}]'.format(self._sample_size,pdg)

        if cmap is None : 
            self._scmpb = None # only defined if color mapped
        else :
            self.map_cmap_from_cloudvals(cmap)

        return self

    # ==================================================
    # Color Properties/Methods
    # ==================================================

    def _getSubCmap(self,subBounds) :
        """ Colormap so that colors match values for cloud and subObj """
        # only called by surface methods:
        #   map_vertvals_from_cloud
        #   map_color_from_cloud
        cmap = self._getCmap()
        if isinstance(cmap,str) :  cmap = mpl.colormaps[cmap]
        Cmin,Cmax = self.bounds['vlim']
        Smin,Smax = subBounds
        low = (Smin - Cmin)/(Cmax-Cmin)
        hi = (Smax - Cmin)/(Cmax-Cmin)
        x = np.linspace(0.0,1.0,num=256)
        cmName = ''.join(random.choices(string.ascii_uppercase , k = 8))
        cmap = colors.ListedColormap( cmap(low + (hi-low)*x),cmName )
        return cmap

    def _getCmap(self) :
        cmap = self._cmap  # .. only defined if colormapped.
        if not self.isColorMapped :
            clist = [ self._singleColor, self._singleColor ]
            cmName = ''.join(random.choices(string.ascii_uppercase , k = 8))
            cmap = colors.LinearSegmentedColormap.from_list(cmName, clist )
        return cmap

    def _map_cmap_from_op_or_vals(self,values,cmap ) :
        # Point color assignment using operation or point values.
        # called from map_cmap_from_cloudvals or map_cmap_from_op

        objMap = self._cmap if cmap is None else cmap
        if isinstance(objMap,str) : 
            objMap = mpl.colormaps[objMap]
        self._cmap = objMap
        self.isColorMapped = True

        vxyz = values.flatten()
        normalizedData = (vxyz-np.min(vxyz))/(np.max(vxyz)-np.min(vxyz))
        vcolors = objMap(normalizedData)
        self._colors = vcolors 
        norm = colors.Normalize(vmin=np.min(vxyz),vmax=np.max(vxyz))
        sm = cm.ScalarMappable(cmap=objMap, norm=norm)
        sm.set_array([])
        self._scmpb = sm
        return self

    def map_cmap_from_cloudvals(self,cmap=None) :
        """
        Point color assignment using point values.

        Parameters
        ----------
        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the default Colormap is used.
            The colormap maps the point values to colors.
        
        Returns
        -------
        self : point cloud object

        """
        vals = self._values
        return self._map_cmap_from_op_or_vals(vals,cmap)

    def map_cmap_from_op(self,operation,cmap=None) :
        """
        Functional assignment of a point color from a color map.
        Also assigns point cloud values.

        Parameters
        ----------
        operation : function object
            Coordinate mapping function that takes one
            argument, a 3xN Numpy array of native coordinates.
            The function returns an array of scalar values

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the default Colormap is used.
            The colormap maps the point values to colors.
        
        Returns
        -------
        self : point cloud object

        """
        #vals = operation(self.cloud_points)
        self.map_vals_from_op(operation)
        vals = self._values
        return self._map_cmap_from_op_or_vals(vals,cmap)

    def map_color_from_op(self,operation, rgb=True, cname=None) :
        """
        Assignment of point color from a function.

        Parameters
        ----------
        operation : function object
            Function that takes one argument,
            a 3xN Numpy array of native coordinates.
            The function returns a 3xN color value.

        rgb : bool {True, False}, optional, default: True
            By default, RGB color values are returned by the
            operation function.  If set False, the operation
            returns HSV color values.

        cname : string, optional, default: None or op function name.
            Descriptive identifier for values indicated by color.

        Returns
        -------
        self : point cloud object

        """

        xyz = self.points
        colors = np.array(operation(xyz))
        hasAlpha = colors.shape[0] == 4
        if hasAlpha : 
            A = colors[3]
            colors = colors[:3]
        else :
            A = np.ones( len(colors[0]) )
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)
        RGBA = np.concatenate((RGB.T,[A]), axis=0).T
        self._colors = RGBA
        self.isSurfaceColored = True
        cname = pvt._getFunctionName(operation,cname) 
        if cname is not None : self.valuesName = cname 

        return self

    def set_color(self,color) :
        """
        Assignment of color.

        Parameters
        ----------
        color : color assigned to all points.

        Returns
        -------
        self : point cloud object

        """

        self._initColor( color )
        self._scmpb = None # only defined if color mapped

        return self

    # ==================================================
    # Object Creation Methods
    # ==================================================

    def _valsurface(self,valList,**kwargs) :
        # ......................................................
        def _get_VF(dataVol,sval,A,B) :
            M = dataVol.shape[0] - 1
            N = dataVol.shape[1] - 1
            P = dataVol.shape[2] - 1
            spacing = (1/M, 1/N, 1/P)
            level = 0.001    # FutDev: level should not be a constant.
            vol = dataVol - sval
            verts, fvIndices, _, _ = measure.marching_cubes(vol, level, spacing=spacing)
            C = np.empty([len(verts),3])
            C[:] = A
            verts = C*verts + B
            return verts,fvIndices
        # ......................................................
        dm = self._domain
        XD,YD,ZD = dm[0],dm[1],dm[2]
        dmT = dm.T
        A = dmT[1]-dmT[0]
        B=np.array( [ XD[0],YD[0],ZD[0] ])
        cpnts = self._values
        minCpts, maxCpts = np.min(cpnts),np.max(cpnts)
        valrange = maxCpts - minCpts

        # extract from kwargs....
        # ... can use 'color' or 'cmap' ...
        curMap = self._getCmap()
        kargDft = { 'color':None, 'cmap':curMap, 'name':None  }
        KW = pvt.KWprocessor(kargDft,'implsurf',**kwargs)
        name = KW.getVal('name') 
        userColor = KW.getVal('color')
        useDefColor = userColor is not None   # use color arg if color is defined.
        if useDefColor : # if single color, convert to a list.
            if not isinstance(userColor, list) : userColor=[userColor]
        cmapobj = KW.getVal('cmap')
        if isinstance(cmapobj,str) : 
            cmapobj = mpl.colormaps[cmapobj]
        kwargs = KW.filter()

        surface = None
        minVal,maxVal = valList[0], valList[0]  # initialize

        for i,fval in enumerate(valList) :
            if fval < minVal : minVal=fval
            if fval > maxVal : maxVal=fval
            verts, faces = _get_VF(cpnts,fval,A,B)

            if useDefColor :
                color = userColor[i%len(userColor) ]
            else :
                color = cmapobj( (fval-minCpts)/valrange )
            surf = Surface3DCollection( verts, faces,color=color)
            surface = surf if surface is None else surface+surf

        surface.set_cmap(self._getCmap())
        surface._bounds['vlim'] =[ minCpts, maxCpts ] 
        # now that the surface face vals are set, set the surface vertex values
        try :
            surface._map_vertvals_from_cloud(self)
        except :
            # for density clouds, the computed surface vertices may be external
            # the domain, but the surface face colors are still assigned since
            # the samples are still in the domain.
            warnings.warn('WARNING: surface vertex values are NOT assigned !')
        if self.val_source != 'function' : surface.evert()

        if name is not None: surface._geomName = name
        return surface

    def valsurf(self, *fval, **kwargs) :
        """
        Surface3DCollection defined for a constant value a the point cloud.

        Parameters
        ----------
        fval : a number
            Scalar value of surface contour within the cloud.

        Other Parameters
        ----------------
        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        color : color or list of colors, optional, default: None

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            Used for multiple surface values if color is not assigned.

        Returns
        -------
        Surface3DCollection object

        """
        surface = self._valsurface(fval, **kwargs)
        return surface

    def valsurfSet(self, numb, **kwargs) :
        """
        Surface3DCollection defined for constant values in a point cloud.

        Parameters
        ----------
        numb : an integer
            Number of constant valued contour surfaces within the domain.
            Values are evenly spaced between minimum to maximum cloud values.

        Other Parameters
        ----------------
        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        color : color or list of colors, optional, default: None

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            Used for multiple surface values if color is not assigned.

        Returns
        -------
        Surface3DCollection object

        """
        minVal, maxVal = np.min(self._values), np.max(self._values)
        delta = (maxVal - minVal)/(numb + 1)
        fval = [ minVal + (i+1)*delta for i in range(numb)]
        surface = self._valsurface(fval, **kwargs)
        return surface

    # ==================================================
    # Matplotlib 3D Visualization Properties/Methods
    # ==================================================

    def add_to3d(self,axes,*focus,**kwargs) :
        """
        Add the cloud points to the the 3D axis.

        Parameters
        ----------
        axes : (Matplotlib 3D axes object)

        focus : number 
            Focus ranges from 0 to 1 going from low to
            high cloud point normalized values, respectively. 
            The color alpha is modified from opaque at the focus
            to fully transparent.

        span : number from 0.1 to 1.0
            The span controls the width of alpha values
            from opaque to transparent over the range
            of cloud values.

        subview : a 3 X 2 float array
            Sets the view domain within the point cloud domain.

        size : float, default : 10
            Marker size of cloud points.

        marker : string, default : 'o'
            Matplotlib marker. 

        Returns
        -------
        number of cloud points displayed in view.

        """
        # ......................................................
        def fs_func(x,focus,span=None) :
            if focus<0 or focus>1 :
                return np.zeros(len(x))
            if span is None :  
                span = focus if focus>0.5 else 1.0-focus
            if span < 0.1 : span=0.1
            y = ( 1.0 + np.cos(np.pi*(x-focus)/span) )/2.0
            test = np.logical_or ( (x>focus+span) , (x<focus-span) )
            y = np.where(test, np.zeros(len(x)), y)
            return y
        # ......................................................

        kwargDft = { 'span':None, 'subview':None, 'size':self._pntSize,
                     'marker':self._marker  }
        KW = pvt.KWprocessor(kwargDft,'Point3DCloud',**kwargs)
        xyz = self.points
        c = self._colors
        s = KW.getVal('size')
        m = KW.getVal('marker')
        g = KW.getVal('span')
        clip = KW.getVal('subview')
        # clip before applying focus to reduce computation. <<<<<
        # especially when there are no points to plot. <<<<<<<<<<
        if len(focus) == 0 : focus = None
        if (focus is not None) and ( self.hasVals ) :
            if (len(focus) > 1) and (g is None) : g = 0.1
            vxyz = self._values.flatten()
            normData = (vxyz-np.min(vxyz))/(np.max(vxyz)-np.min(vxyz))
            alphas = np.zeros(len(normData))
            augColors = np.array(c)
            for f in focus :
                a = fs_func(normData,f,g)
                alphas = np.where( a>alphas, a, alphas)
            augColors[:,3] = alphas
            c= augColors


        if clip is not None :
            clip = np.array(clip)
            coor = xyz.T
            low = np.prod(coor>=clip[:,0],axis=1)
            hi  = np.prod(coor<=clip[:,1],axis=1)
            hasVals = (np.sum(low*hi)).astype(bool)
            if not hasVals : return 0
            inDomain = (low*hi).astype(bool)
            xyz = coor[inDomain].T
            c = c[inDomain]

        # ......................................................
        axes.scatter(*xyz, c=c, marker=m, s=s)
        viewPoints = len(c)
        return viewPoints

# ==================================================

