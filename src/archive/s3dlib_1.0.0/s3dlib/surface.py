# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

"""
A module for creating 3D surface objects and auxiliary functions.
 
A surface is a collection of faces with 3D vertices ordered using
the 'right hand rule' to designate the 'outer' surface normals.  All surface faces
are joined with a minimum of one adjacent face. Adjacent faces share two
common vertices.

Surface instantiation is primarily performed using the four subclasses that
have predefined surface topologies in native coordinates and triangular faces.
The surface object geometry and color are controlled through various
object methods in the 'Surface3DCollection' base class.
Base class objects are created by addition of subclass surface objects to
form a single 'composite' surface.

Primitive surface objects may be created from faces with 3, 4, or 5
edges using the  'Surface3DCollection' base class in
xyz native coordinates.

Objects are added to the mpl_toolkits.mplot3d.axes3d using
the Axes3D.add_collection3d() method.

"""

# DevNotes:
#     0  FutDev: future development notes.
#     1. Needed to used self._facecolors instead of self.get_facecolor()
#        Need to look into this.
#     2. Edgecolor not taking the alpha for facecolors.  May need a
#        direct assignment.
#     3. Tests indicate 80-95% of exec. time for larger rez is spent on the
#        call to Collection.set_color().  May need a more direct approach?

import warnings
import copy
from functools import reduce
from time import time

import numpy as np
from scipy import interpolate

import matplotlib as mpl
from matplotlib import cm, colors, image

from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

_MINRAD = 0.01    # in subclassed polar and spherical surface split basetypes.
_MAXREZ = 8       # in subclass surfaces, maximum recursive triangulation.


class Surface3DCollection(Poly3DCollection):
    """
    Base class for 3D surfaces.
    """

    _ILLUM = (1,0,1)  # defaut light direction, (cmap_normals, shade and hilite).
    _DFT_ALR = 0.25   # default axis length ratio, head size to vector magnitude.

    @staticmethod
    def triangulateBase(rez,baseVcoor,baseFaceVertexIndices, midVectFunc) :
        """
        Recursively subdivide a triangulated surface.

        Each recursion subdivides each surface face by four.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        baseVcoor : V x 3 float array
            An array of 'v' number of xyz vertex coordinates.

        baseFaceVertextIndices : F x 3 int array
            An array of 'F' number of face vertex indices.

        midVectFunc : function object
            A function that takes two xyz coordinate (list of 3)
            representing a surface face edge.  Returns one xyz coordinate
            at the bisection of the edge, mapped on the surface.

        Returns
        -------
        indexObj : a dictionary of vertex indices, for a surface
            of F faces, E edges and V vertices.
            'face' : F x 3 int array
            'edge' : E x 2 int array

        vertCoor : V x 3 float array
            An array of V number of xyz vertex coordinates.

        """

        vCoor = baseVcoor.copy()
        fvIndices = []
        evIndices = []
        #.....................................    
        def getVerticesLine(degree, leftIndex, rightIndex) :
            def recurs_edgeCoor(degree, indexOrigin, isRight, leftIndex, rightIndex, Eindex) :
                m = degree - 1
                delta = int(np.power(2,m))
                if isRight : delta = -delta
                i = indexOrigin - delta
                mid_coor = midVectFunc( vCoor[leftIndex], vCoor[rightIndex] )
                currentCoorIndex = len(vCoor)
                vCoor.append(mid_coor)
                Eindex[i] = currentCoorIndex
                if m <= 0 :
                    evIndices.append( [leftIndex, currentCoorIndex] )                    
                    evIndices.append( [currentCoorIndex, rightIndex] )                    
                    return
                recurs_edgeCoor(m,i,False,leftIndex,currentCoorIndex, Eindex)
                recurs_edgeCoor(m,i,True,currentCoorIndex,rightIndex, Eindex)
                return
            #.....................................    
            iOrigin = int(np.power(2,degree))
            Eindex = [None]*(iOrigin+1)
            Eindex[0] = leftIndex
            Eindex[iOrigin] = rightIndex
            if degree==0 :
                evIndices.append( [leftIndex, rightIndex] )
                return Eindex
            recurs_edgeCoor(degree, iOrigin, False, leftIndex, rightIndex, Eindex)
            return Eindex
        #.....................................    
        def recurs_faceIndices(A,B,C,rez,atCenter=True):
            if rez==0 :
                abc = [ A[0], B[0], C[0] ]               
                fvIndices.append( abc )
                if atCenter : evIndices.append( [ B[0], C[0] ] )
                return
            # construct the sub triangle edge vertices....
            cIndex = int( (len(A)-1)/2 )
            cShift = cIndex+1
            Ax = A[:cShift]
            xB = A[cIndex:]
            By = B[:cShift]
            yC = B[cIndex:]
            Cz = C[:cShift]
            zA = C[cIndex:]
            # construct the center sub triangle edge vertices....
            if rez==1 :
                xz = [ A[1], C[1] ]
                yx = [ B[1], A[1] ]
                zy = [ C[1], B[1] ]
            else:               
                xz = getVerticesLine(rez-1, A[cIndex], C[cIndex] )
                yx = getVerticesLine(rez-1, B[cIndex], A[cIndex] )
                zy = getVerticesLine(rez-1, C[cIndex], B[cIndex] )
            # construct the 4 sub triangles...
            recurs_faceIndices(Ax,xz,zA,rez-1)
            recurs_faceIndices(By,yx,xB,rez-1)
            recurs_faceIndices(Cz,zy,yC,rez-1)
            recurs_faceIndices(yx[::-1],zy[::-1],xz[::-1],rez-1,False)
            return
        #.....................................
        # generate the base edge indices for each face 
        baseFaceEdgeIndices = []
        edgeArrayIndexMatrix = [[None for x in range(len(vCoor))] for y in range(len(vCoor))]   
        for face in baseFaceVertexIndices :
            A = [ face[0],face[1]]
            B = [ face[1],face[2]]
            C = [ face[2],face[0]]
            edgeIndices = [A,B,C]
            baseFaceEdgeIndices.append(edgeIndices)
            for vertexIndex in edgeIndices :
                i = vertexIndex[0]
                j = vertexIndex[1]
                if edgeArrayIndexMatrix[i][j] is None:
                    edgeArrayIndexMatrix[i][j] = getVerticesLine(rez, i, j )
                    edgeArrayIndexMatrix[j][i] = edgeArrayIndexMatrix[i][j].copy()[::-1]
        # trianglulate the base triangles
        for face in baseFaceEdgeIndices :
            A = edgeArrayIndexMatrix[face[0][0]][face[0][1]]
            B = edgeArrayIndexMatrix[face[1][0]][face[1][1]]
            C = edgeArrayIndexMatrix[face[2][0]][face[2][1]]
            x = (rez > 0)
            recurs_faceIndices(A,B,C,rez, x)
        vertexCoor = np.array(vCoor) 
        indexObj = { 'face': np.array(fvIndices) , 'edge': np.array(evIndices) }
        return indexObj ,vertexCoor

    @staticmethod
    def coor_convert(xyz, tocart=True) :
        """Coordinate transformation.
        
           To be overriden by any subclass not in Cartesian coordinates.
        """
        return xyz

    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Map surface to rectangular unit coordinates (0,1)."""

        x,y,z = xyz
        a = (x+1)/2
        b = (y+1)/2
        return [a,b]

    @staticmethod
    def _sfev(rez,fevArray) :
        """Calculate number of surface faces, edges, vertices."""

        Fo, phi, Vx = fevArray
        eta = np.power(2,rez)
        F = eta*eta*Fo
        E = np.int( 3*F/2 + eta*phi )
        V = np.int( F/2 + eta*phi + Vx )
        return [F,E,V]

    @classmethod
    def _fev(cls,rez,basetype,dercls) :
        """Sets up _sfev method, based on calls from the subclass dercls."""

        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None: basetype = dercls._default_base
        basesurf = dercls._base_surfaces
        if basetype not in basesurf :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(basesurf)))
        return cls._sfev(rez, basesurf[basetype]['fevParam'] )

    @classmethod
    def _uvw_to_xyz(cls,rst,abc) :
        """Rotational transformation at a coordinate.
        
           To be overriden by any subclass not in Cartesian coordinates.
        """

        return np.transpose(rst)

    @classmethod
    def _disp_to_xyz( cls, uvw, theta=None, phi=None) :
        """Rotational transformation from rotation angles."""
        if theta is None : return uvw
        # each set of uvw has a different rotation matrix
        # that must be evaluated one at a time.
        if phi is None : phi = np.zeros( len(theta) )
        uvwT = np.transpose(uvw)
        dispVect = []
        for i in range(len(theta)) :
            rotMx = eulerRot(theta[i],phi[i], inrad=True)
            vec = np.dot( uvwT[i], rotMx )
            dispVect.append(vec)
        return dispVect

    def __init__(self, vertexCoor, faceIndices, edgeIndices=None, *args, **kwargs):
        """ 
        3D surface of connected F faces, each face with N number of vertices.
        The surface has V number of vertices and E number of edges.
        
        Parameters
        ----------
        vertexCoor : V x 3 float numpy.ndarray 
            An array of V number of xyz vertex coordinates.

        faceIndices : F x N int numpy.ndarray
            An array of F number of N face vertex indices.

        edgeIndices : E x 2 int numpy.ndarray, optional, default: None
            An array of E number of 2 edge vertex indices.
            If None is provided, the 'edges' property is not available.

        Other Parameters
        ----------------
        *args, **kwargs
            All other arguments are passed on to mpl_toolkits.mplot3d.art3d.Poly3DCollection.
            Valid argument names include: color, edgecolors, facecolors, linewidths, cmap.

        """

        self.vertexCoor = vertexCoor.copy()
        self.baseVertexCoor = vertexCoor.copy()
        self.fvIndices = faceIndices.copy()
        self.evIndices = None
        if edgeIndices is not None :
            self.evIndices = edgeIndices.copy()

        f = self.fvIndices
        v = self.vertexCoor
        super().__init__(v[f], *args, **kwargs)

        self._surfaceShapeChanged = False       # control flag for image operations      
        self._isCompositeSurface = False        # control flag for a composite object.

        self._normal_scale = 1.0                # overridden in subclass default
        
        self.restoredClip = None                # dictionary to restore from clip operation.

        # following is used for map_geom_from_image().
        self.disp_Vector = None                 # overridden in subclass.

        self.vector_field = None                # for vector fields when transformed.
        self.scale = [1,1,1]                    # for axis when transformed.
        self.translation = [0,0,0]              # for axis when transformed.
        self.rotation =  np.identity(3)         # for axis when transformed.

        self._bounds = {}
        self._set_geometric_bounds()
        self._bounds['vlim'] = [ -1.0, 1.0 ]       
    
    def __add__(self, othr) :
        
        if not isinstance(othr, Surface3DCollection) :
            raise ValueError('Add operations can only apply between Surface3DCollection objects.')

        topVerCoor = self.vertexCoor
        botVerCoor = othr.vertexCoor
        totVerCoor = np.append(topVerCoor,botVerCoor,axis=0)
        nextIndex = len(topVerCoor)

        topFVindices = self.fvIndices
        botFVindices = othr.fvIndices
        FVtail = np.add(botFVindices,nextIndex)
        totFVindices = np.append(topFVindices,FVtail,axis=0)

        topEVindices = self.evIndices
        botEVindices = othr.evIndices
        if (topEVindices is None) or (botEVindices is None) :
            totEVindices = None
        else :
            EVtail = np.add(botEVindices,nextIndex)
            totEVindices = np.append(topEVindices,EVtail,axis=0)

        obj = Surface3DCollection(totVerCoor,totFVindices,totEVindices)
        obj._isCompositeSurface = True
        
        top_onearr = np.ones( len(topFVindices) )[:, np.newaxis]
        top_orig_colors = self._facecolors
        if len(top_orig_colors) == 1 : top_orig_colors = top_orig_colors*top_onearr

        bot_onearr = np.ones( len(botFVindices) )[:, np.newaxis]
        bot_orig_colors = othr._facecolors
        if len(bot_orig_colors) == 1 : bot_orig_colors = bot_orig_colors*bot_onearr

        total_colors = np.append(top_orig_colors,bot_orig_colors,axis=0)
        obj.set_color(total_colors)        

        return obj

    def __str__(self) :
        # parent class doesn't have specific properties...
        try:     bs = ' ( {} , {} )'.format(self._rez,self._basetype)
        except:  bs = ' '
        numFaces = len(self.fvIndices)
        numVerts = len(self.vertexCoor)
        name = self.__class__.__name__
        sz = ':  faces: {},  vertices: {}'.format(numFaces,numVerts)
        val = name+bs+sz
        return val

    def _tranformVector(self, orgCoor, operation, returnxyz) :
        """Functional tranformation of coordinates."""
        xyz = np.transpose(orgCoor)
        abc = self.coor_convert(xyz)
        rst = np.array(operation(abc))
        if returnxyz : XYZ = rst
        else :         XYZ = self.coor_convert( rst , tocart=True )
        return np.transpose(XYZ)

    def _get_vectorfield_Line3DCollection(self, vector_field, color=None, width=None, alr=None) :
        """Vector3DCollection at surface vertices."""       
        if alr is None : alr = self._DFT_ALR
        if color is None : color='black'
        if width is None : width = 1.0
        self.vector_field = vector_field
        lcol = Vector3DCollection(self.vertexCoor,vector_field,alr,colors=color, linewidths=width)
        return lcol

    def _get_face_normals(self,faceCoor) :
        """Unit normals of triangular face coordinates."""
        vfT = np.transpose(faceCoor, (1,2,0) )
        vAB = np.subtract( vfT[1], vfT[0]  )
        vAC = np.subtract( vfT[2], vfT[0]  )
        vABt = np.transpose(vAB)
        vACt = np.transpose(vAC)
        cross = np.cross(vABt,vACt)
        sz = np.linalg.norm(cross,axis=1)[:,np.newaxis]
        return np.divide(cross,sz)

    def _get_face_centers(self,faceCoor) :
        """Face centers of traiangular face coordinates."""
        vfT = np.transpose(faceCoor, (1,2,0) )
        numVerts = vfT.shape[0]
        sumV = np.sum(vfT, axis=0)
        return np.transpose(sumV)/numVerts

    def _set_geometric_bounds(self) :
        """Set _bounds dictionary values from the surface vertex coordinates."""

        v = self.vertexCoor
        xlim,ylim,zlim = np.transpose( [np.amin(v ,axis=0), np.amax(v ,axis=0) ] )
        self._bounds['xlim'] = xlim
        self._bounds['ylim'] = ylim
        self._bounds['zlim'] = zlim
        rmax_xy =  np.amax( np.linalg.norm(v[:,0:2], axis=1) )
        rmin_xy =  np.amin( np.linalg.norm(v[:,0:2], axis=1) )
        self._bounds['r_xy'] = [rmin_xy, rmax_xy ] 
        rmax_xyz = np.amax( np.linalg.norm(v, axis=1) )
        rmin_xyz = np.amin( np.linalg.norm(v, axis=1) )
        self._bounds['rorg'] = [rmin_xyz, rmax_xyz ] 
        # vestigal...
        self._bounds['rlim'] = [rmax_xy, rmax_xyz ]
        return    

    def _viewportCoor(self,xyzCoor,viewport=None) :
        xyz = np.transpose(xyzCoor)
        ab = self._map3Dto2Dsquare( xyz )
        a,b = ab
        if viewport is None : return a,b,np.full(len(a),True)
        As,Bs,Ae,Be = viewport
        Aref = 1.0-As
        Bref = 1.0-Bs
        Adelta = np.mod(Ae+Aref,1.0)
        Bdelta = np.mod(Be+Bref,1.0)
        if Adelta <=0.0 : Adelta = 1.0
        if Bdelta <=0.0 : Bdelta = 1.0
        Ap = np.mod(a+Aref,1.0)
        Bp = np.mod(b+Bref,1.0)
        Av = np.divide(Ap,Adelta)
        Bv = np.divide(Bp,Bdelta)
        inViewport = Av<=1
        inViewport = np.where(Bv>1,np.full(len(inViewport),False),inViewport)
        Av = np.clip(Av,0.0,1.0)
        Bv = np.clip(Bv,0.0,1.0)
        return Av,Bv,inViewport

    @property
    def bounds(self) :
        """
        Dictionary of surface geometric and value ranges.
        
        Each dictionary value is a 2 float array of minimum
        and maximum values of the surface.  Keys are:

            'xlim' : x-coordinate
            'ylim' : y-coordinate
            'zlim' : z-coordinate
            'r_xy' : radial distance from the z axis
            'rorg' : radial distance from the origin
            'vlim' : value functional assignments.
            
        Values are assigned using datagrid and color mapping methods.

        """

        return self._bounds

    @property
    def cBar_ScalarMappable(self) :
        """matplotlib.cm.ScalarMappable object for surface values."""
        vmin, vmax = self._bounds['vlim']
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        objMap = self.get_cmap()
        sm = cm.ScalarMappable(cmap=objMap, norm=norm)
        sm.set_array([])
        return sm        
    
    @property
    def edges(self) :
        '''Line3DCollection of the triangulated surface edges.'''

        if self.evIndices is None :
            raise ValueError('No Edge Indices argument during the Surface3DCollection object instantiation.')       
        e = self.evIndices
        v = self.vertexCoor
        lcol = Line3DCollection( v[e] )
        return lcol

    @property
    def vertices(self) :
        """A 3 x N array of surface vertices."""
        v = self.vertexCoor
        m = np.transpose(v)
        return m[0],m[1],m[2]

    @property
    def facecenters(self) :
        """A 3 x N array of surface face center coordinates."""
        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        m = np.transpose(faceCenterCoor)
        return  m[0],m[1],m[2]       
    
    def triangulate(self, rez=0) :
        """
        Planar subdivision of each face into triangular faces.
        
        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the base faces
            into triangular faces, four faces per rez.
            Rez values range from 0 to 7.
            For surfaces with 4 vertices per face, each face is first
            subdivided into two faces, then recursion proceeds.
            For surfaces with 5 vertices per face, each face is first
            subdivided into three faces, then recursion proceeds.
            
        Returns
        -------
        self : surface object
        
        """

        # .............................................................
        def divideFace4(face) :
            dist_a = np.linalg.norm( self.vertexCoor[ face[0]] - self.vertexCoor[ face[2]] )
            dist_b = np.linalg.norm( self.vertexCoor[ face[1]] - self.vertexCoor[ face[3]] )
            if dist_a < dist_b :
                face1 = [ face[0], face[1], face[2] ]
                face2 = [ face[2], face[3], face[0] ]
            else :
                face1 = [ face[1], face[2], face[3] ]
                face2 = [ face[3], face[0], face[1] ]               
            return face1, face2
        # .............................................................
        def divideFace5(face) :
            diag = lambda i : self.vertexCoor[ face[i]] - self.vertexCoor[ face[(i+2)%5] ]
            fc3 = lambda i : [ face[i], face[(i+1)%5], face[(i+2)%5] ]
            fc4 = lambda i : [ face[(i+2)%5], face[(i+3)%5], face[(i+4)%5], face[(i+5)%5] ]
            dist = [None]*5
            for i in range(5) : dist[i] = np.linalg.norm( diag(i) )
            minDiag = 0
            for i in range(1,5) :
                if dist[i] < dist[minDiag] : minDiag = i
            face1 = fc3(minDiag)
            otherFace = fc4(minDiag)
            face2, face3 = divideFace4(otherFace)
            return face1, face2, face3
        # .............................................................
        def midVFun(vectA, vectB) :
            mid = np.add(vectA,vectB)
            mid = np.multiply(0.5,mid)
            return mid       
        # .............................................................
        # FutDev: major need to optimize efficiency
        subFaceIndices = []
        if self.fvIndices.shape[1] > 3 :
            for face in self.fvIndices :
                if len(face) == 3 :
                    subFaceIndices.append(face)
                elif len(face) == 4 :
                    f1, f2 = divideFace4(face)
                    subFaceIndices.append(f1)
                    subFaceIndices.append(f2)
                elif len(face) == 5 :
                    f1, f2, f3 = divideFace5(face)
                    subFaceIndices.append(f1)
                    subFaceIndices.append(f2)
                    subFaceIndices.append(f3)
            self.fvIndices = np.array(subFaceIndices)

        indexObj ,vertexCoor = self.triangulateBase(rez,list(self.vertexCoor),self.fvIndices, midVFun)
        self.fvIndices = indexObj['face']
        self.evIndices = indexObj['edge']
        self.vertexCoor = vertexCoor
        v = self.vertexCoor
        verts = v[self.fvIndices]
        self.set_verts( verts )
        self._set_geometric_bounds()
        
        c = colors.to_rgba( mpl.rcParams['patch.facecolor'] )
        colorArr = np.tile(c,(len(self.fvIndices),1))
        self.set_color( colorArr  )
        return self        

    def normalize_scale(self) :
        """Scaling normalization array and reciprocal."""
        bounds = self.bounds
        scX = ( bounds['xlim'][1] - bounds['xlim'][0]  )
        scY = ( bounds['ylim'][1] - bounds['ylim'][0]  )
        scZ = ( bounds['zlim'][1] - bounds['zlim'][0]  )
        recip = [scX,scY,scZ]
        scale = np.reciprocal(recip)
        return scale, recip

    def facenormals(self, scale=None, color=None, width=None, alr=None) :
        """
        Vector3DCollection of scaled face normals at face centers.

        Parameters
        ----------
        scale: number, optional
            If not spectified, scaled proportional to the mean edge
            length of surface base faces.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object

        """

        if alr is None : alr = self._DFT_ALR
        if color is None : color='black'
        if width is None : width = 1.0
        if scale is None : scale = self._normal_scale
        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        faceNormalCoor = self._get_face_normals(verts)
        lcol = Vector3DCollection(faceCenterCoor,faceNormalCoor*scale,alr,colors=color, linewidths=width)
        return lcol

    def dispfield_from_op(self, operation, returnxyz=False, useBase=False, scale=1, color=None, width=None, alr=None) :
        """
        Vector3DCollection of displacement of vertices to a different position.

        Parameters
        ----------
        operation : function object
            Function that takes one coordinate
            argument, a 3xN Numpy array.
            The function returns a 3xN array of vectors.
        
        returnxyz : bool { True, False }, optional, default: False
            By default, native coordinates are returned by the
            operation function.  If set True, the operation
            returns xyz Cartesian coordinates.

        useBase : bool { True, False }, optional, default: False
            When set False, vertices of the surface, prior to
            any geometric mapping or transforms, are passed to
            the operation function.  Otherwise, when True, the
            current surface vertices are passed.

        scale: number, optional, default: 1.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object

        """

        if alr is None : alr = self._DFT_ALR
        if self._isCompositeSurface :
            warnings.warn('Vector operation not available for combined shapes.')
            return
        start_coor = self.vertexCoor
        if useBase : start_coor = self.baseVertexCoor
        v = self._tranformVector(start_coor, operation, returnxyz)
        delta = v - self.vertexCoor
        delta = scale*delta
        return self._get_vectorfield_Line3DCollection(delta, color, width, alr)

    def vectorfield_from_op(self, operation, scale=1, color=None, width=None, alr=None) :
        """
        Vector3DCollection of vectors in u,v,w coordinates at surface vertices.

        Parameters
        ----------
        operation : function object
            Function that takes one coordinate
            argument, a 3xN Numpy array of native coordinates.
            The function returns a 3xN array of vectors.
        
        scale: number, optional, default: 1.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object

        """

        if alr is None : alr = self._DFT_ALR
        if self._isCompositeSurface :
            warnings.warn('Vector operation not available for combined shapes.')
            return

        xyz = np.transpose(self.vertexCoor)
        abc = self.coor_convert(xyz)
        abc = np.array(abc)
        rst = np.array(operation(abc))
        rst = scale*rst
        delta = self._uvw_to_xyz(rst,abc)
        return self._get_vectorfield_Line3DCollection(delta, color, width, alr)

    def vectorfield_to_surface(self, surface, scale=1, color=None, width=None, alr=None) :
        """
        Vector3DCollection of vectors from surface to surface vertex coordinates.

        Parameters
        ----------
        surface : surface object
            Surface that matches the calling surface vectors.
            Should be the same basetype and rez for surface subclasses.

        scale: number, optional, default: 1.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Raises
        ------
        ValueError
            Mismatched surfaces based on the number of vertices.

        Returns
        -------
        Vector3DCollection object

        """
        surLen, selfLen = len(surface.vertexCoor) , len(self.vertexCoor)
        if surLen is not selfLen :
            raise ValueError('Surfaces have unequal number of vertices, {} != {}'.format(surLen, selfLen))            
        if alr is None : alr = self._DFT_ALR
        delta =  surface.vertexCoor - self.vertexCoor
        delta = scale*delta
        return self._get_vectorfield_Line3DCollection(delta, color, width, alr)

    def map_color_from_image(self, fname, viewport=None) :
        """
        Assign image color values to surface face colors.

        This method is not available for composite surfaces.
        The entire image domain is applied to the coloring operation.

        Parameters
        ----------
        fname : str or file-like
            The image file to read: a filename, a URL or a file-like object opened
            in read-binary mode.  Currently restricted to PNG format.
        
        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the image is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is colored for the default of None.

        Returns
        -------
        self : surface object

        """
        
        if self._isCompositeSurface :
            warnings.warn('Image operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Image operation may be anomalous after shape modification.')

        if viewport is not None :
            viewport = np.array(viewport)
            if np.any( (viewport<0) | (viewport>1) ) :
                raise ValueError('viewport list values, {}, must be between 0 and 1'.format(viewport))
            if len(viewport) != 4 :
                raise ValueError('viewport list must be length 4, found {}'.format(len(viewport)))      
        verts = self.vertexCoor[ self.fvIndices ]
        
        faceCenterCoor = self._get_face_centers(verts)
        a,b,inViewport = self._viewportCoor(faceCenterCoor,viewport)
        img = image.imread(fname)
        
        height, width = np.subtract(img.shape[:2],1)
        M_index = ( (1-b)*height ).astype(int)
        N_index = ( a*width ).astype(int)
        #'''
        orig_colors = self._facecolors
        if len(orig_colors) == 1 :
            onearr = np.ones( len(faceCenterCoor) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       

        colorMap = []
        if viewport is None : colorMap = img[M_index,N_index]
        else :
            # FutDev: use numpy methods to optimize efficiency.
            for i in range( len(orig_colors) ) :
                colorMap.append(orig_colors[i])
                if inViewport[i] : colorMap[i] = img[M_index[i],N_index[i]]    

        self.set_color(colorMap)
        self.isSurfaceColored = True
        return self

    def map_color_from_op(self, operation, rgb=True) :
        """
        Assignment of face color from a function.

        Face colors are assigned from a function
        of face coordinates.

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

        Returns
        -------
        self : surface object

        """

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        xyz = np.transpose(faceCenterCoor)
        abc = self.coor_convert(xyz)
        colors = np.array(operation(abc))
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)
        self.set_color(RGB)
        self.isSurfaceColored = True
        return self

    def map_cmap_from_datagrid(self, datagrid, cmap=None, viewport=None) :
        """
        Face color assignment using a 2D datagrid surface.
        
        Datagrid values are normalized in the range 0 to 1.

        This method is not available for composite surfaces.
        The entire datagrid domain is applied to the geometric operation.

        Parameters
        ----------
        datagrid : 2D float array

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the datagrid values to colors.

        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the datagrid is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is colored for the default of None.

        Returns
        -------
        self : surface object

        """

        if self._isCompositeSurface :
            warnings.warn('Datagrid operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Datagrid operation may be anomalous after shape modification.')

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        a,b,inViewport = self._viewportCoor(faceCenterCoor,viewport)      
        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        xd = np.linspace(0, 1, data.shape[0] )
        yd = np.linspace(0, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')
        d = []
        # FutDev: use numpy methods to optimize efficiency.
        for ab in np.transpose([a,b]) :
            d.append( g(ab[0],ab[1])[0] )   
        d = np.where(inViewport,d,np.full(len(d),dmin))
        
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cmap = self.get_cmap()
        
        orig_colors = self._facecolors

        if len(orig_colors) == 1 :
            onearr = np.ones( len(faceCenterCoor) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       
        
        norm = colors.Normalize(dmin,dmax)
        temp = cmap(norm(d))
        if viewport is None :
            colorMap = temp
        else :
            colorMap = np.where(inViewport[:, np.newaxis], temp, orig_colors )

        # FutDev: set_color used instead of set_facecolor to eliminate 'gaps' between faces,
        # however, set_edgecolor required after illumination if edges are to be displayed.
        # HOWEVER !!!!  if colorMap has a alpha<1, these don't register with the edges.
        self.set_color(colorMap)
        self.isSurfaceColored = True
        self._bounds['vlim'] = [ dmin, dmax ]
        return self

    def map_cmap_from_normals(self, cmap=None, direction=None) :
        """
        Face color assignment using normals relative to illumination.

        The dot product of face normals with the illumination direction 
        is used to assign face colors from a colormap.

        Parameters
        ----------
        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the dot product values to colors.
        
        direction : list of size 3, optional, default: [1,0,1]
            A 3D vector in xyz Cartesian coordinates designating
            the direction of the illumination source.

        Returns
        -------
        self : surface object

        """

        #.....................................    
        def getColorMap(fvo,cmap, direction) :
            unitVector = lambda v : np.divide( v, np.linalg.norm(v) )
            incidentLight = unitVector( direction )
            d = np.dot(fvo,incidentLight)
            v = np.add(d,1)
            v = np.divide(v,2)
            v = np.abs(v)
            colorMap = cmap(v)
            return colorMap
        #.....................................
        if direction is None : direction = self._ILLUM
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        verts = self.vertexCoor[ self.fvIndices ]
        fvobj = self._get_face_normals(verts)

        cMap = getColorMap(fvobj, cm_colorMap, direction)
        self.set_cmap(cm_colorMap)
        self.set_color(cMap)
        return self

    def map_cmap_from_op(self, operation, cmap=None) :
        """
        Functional assignment of a color from a color map.

        Face coordinates are used to calculate a scalar 
        which is then used to assign face colors from a
        colormap.

        Parameters
        ----------
        operation : function object
            Function that takes one argument,
            a 3xN Numpy array of native coordinates.
            The function returns a Numpy array of scalar values.

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the function return values to colors.
        
        Returns
        -------
        self : surface object

        """

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cmap = self.get_cmap()

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        xyz = np.transpose(faceCenterCoor)
        # .....
        abc = self.coor_convert(xyz)
        v = np.array(operation(abc))
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())
        
        colorMap = cmap(norm(v))
        # FutDev: set_color used instead of set_facecolor to eliminate 'gaps' between faces,
        # however, set_edgecolor required after illumination if edges are to be displayed.
        # HOWEVER !!!!  if colorMap has a alpha<1, these don't register with the edges.
        self.set_color(colorMap)
        self.isSurfaceColored = True
        self._bounds['vlim'] = [ vmin, vmax ]
        return self

    def map_geom_from_datagrid(self,datagrid, scale=1.0, viewport=None ) :
        """
        Append surface vertices proportional to 2D datagrid surface.
        
        Datagrid values are normalized in the range 0 to 1.

        This method is not available for composite surfaces.
        The entire datagrid domain is applied to the geometric operation.

        Parameters
        ----------
        datagrid : 2D array

        scale : number, optional, default: 1

        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the datagrid is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is geometrically mapped for the default of None.
        
        Returns
        -------
        self : surface object

        """

        if self._isCompositeSurface :
            warnings.warn('Datagrid operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Datagrid operation may be anomalous after shape modification.')

        a,b,inViewport = self._viewportCoor(self.vertexCoor,viewport)

        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        delta = dmax-dmin
        data = (datagrid - dmin)/delta
        xd = np.linspace(0, 1, data.shape[0] )
        yd = np.linspace(0, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')

        # FutDev: optimize efficiency (?).
        d = []
        for ab in np.transpose([a,b]) :
            d.append( g(ab[0],ab[1])[0] )
        d = np.where(inViewport,d,np.zeros(len(d)))
        d = d[:, np.newaxis]
        
        # FutDev: 'fVal' only applicable for initial subclassed
        # surfaces to get surface normals at a vertex.
        # Future work: to make this generic for any surface.
        fVal = np.array([0.,0.,1.])
        if self.disp_Vector is not None :
            fVal = self.vertexCoor*self.disp_Vector

        v =  self.vertexCoor + scale*d*fVal
        verts = v[self.fvIndices]
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( verts )
        self._surfaceShapeChanged = True
        return self

    def map_geom_from_image(self, fname, scale=1.0, viewport=None, cref='v', hzero=0) :
        """
        Append surface vertices proportional to image color component values.

        For planar and polar coordinates, the z-coordinate is appended.
        For cylindrical and spherical coordinates, the r-coordinate is appended.

        This method is not available for composite surfaces.
        The entire image domain is applied to the geometric operation.

        Parameters
        ----------
        fname : str or file-like
            The image file to read: a filename, a URL or a file-like object opened
            in read-binary mode.  Currently restricted to PNG format.

        scale : number

        cref : {'r','g','b','h','s','v'}, optional, default: 'v'
            Sets which color tuple value to use for geometry mapping.
            Note: only the first character of the string is evaluated.

        hzero : number or None, optional, default: 0
            Argument is used if cref is set to 'h'.
            The hzero magnitude indicates the Hue value for a zero coordinate
            diplacement. The hzero sign (positive or negative),
            indicates the direction of increasing value in the
            appended coordinate with hue value.  Range of hzero is [-1,1].

        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the image is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is geometrically mapped for the default of None.

        Returns
        -------
        self : surface object

        """
        # ..............................................................
        def getIndex(cref) :
            cref = cref[0].lower()
            valIndex = 'rgbhsv'.find(cref)
            if valIndex < 0 :
                raise ValueError('Invalid cref argument: ' + str(cref))
            isHSV = False
            if valIndex > 2 :
                valIndex -= 3
                isHSV = True
            return valIndex, isHSV
        # ..............................................................
        if self._isCompositeSurface :
            warnings.warn('Image operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Image operation may be anomalous after shape modification.')
        
        if (hzero<-1) or (hzero>1) :
            raise ValueError('hzero must be between -1 and 1, found {}'.format(hzero))
        ho = (1.0 -hzero)
        hsgn = np.sign(hzero)
        h_ref = lambda h : np.mod( ho + hsgn*h , 1.0)
        # ...................
        a,b,inViewport = self._viewportCoor(self.vertexCoor,viewport)
        img = image.imread(fname)
        height, width = np.subtract(img.shape[:2],1)
        x_index = ( (1-b)*height ).astype(int)
        y_index = ( a*width ).astype(int)

        valIndex, isHSV = getIndex(cref)
        rgbcolor = img[x_index,y_index]
        if isHSV :
            hsvcolor = colors.rgb_to_hsv(rgbcolor)
            hsvcolor = np.transpose(hsvcolor)
            if valIndex==0 : d= h_ref(hsvcolor[0])
            else :           d = hsvcolor[valIndex]
        else :
            rgbcolor = np.transpose(rgbcolor)
            d = rgbcolor[valIndex]

        d = np.where(inViewport,d,np.zeros(len(d)))
        d = d[:, np.newaxis]
        
        # FutDev: 'fVal' only applicable for initial subclassed surfaces
        # surfaces to get surface normals at a vertex.
        # Future work: to make this generic for any surface.
        fVal = np.array([0.,0.,1.])
        if self.disp_Vector is not None :
            fVal = self.vertexCoor*self.disp_Vector

        v =  self.vertexCoor + scale*d*fVal
        verts = v[self.fvIndices]
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( verts )
        self._surfaceShapeChanged = True
        return self
    
    def map_geom_from_op(self, operation, returnxyz=False) :
        """
        Functional transformation of surface vertex coordinates.

        This method is not available for composite surfaces.

        Parameters
        ----------
        operation : function object
            Coordinate mapping function that takes one
            argument, a 3xN Numpy array of native coordinates.
            The function returns a 3xN array of coordinates.

        returnxyz : bool { True, False }, optional, default: False
            By default, native coordinates are returned by the
            operation function.  If set True, the operation
            returns xyz Cartesian coordinates.

        Returns
        -------
        self : surface object

        """

        if self._isCompositeSurface :
            warnings.warn('Map operation not available for combined shapes.')
            return self
        v = self._tranformVector(self.vertexCoor, operation, returnxyz)
        verts = v[self.fvIndices]
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( verts )
        self._surfaceShapeChanged = True
        return self
    
    def clip(self, operation, usexyz=False) :
        """
        Remove faces from the surface.

        Parameters
        ----------
        operation : function object
            Function that takes one argument, a 3xN Numpy array of coordinates.
            The function returns a bool { True, False } indicating if
            the face, at the face centered coordinate, is to be 
            retained (True), otherwise, the face is removed from the
            surface (False).

        usexyz : bool { True, False }, default: False
            If True, face centers are passed to the operation function
            in xyz coordinates.  If False, face centers are passed
            in native coordinates.

        Returns
        -------
        self : surface object

        """

        # Clip operation MUST use xyz coordinates for combined shapes.
        # since native coordinate may be undefined.
        if self._isCompositeSurface : usexyz = True

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)

        onearr = np.ones( len(self.fvIndices) )[:, np.newaxis]
        orig_colors = self._facecolors
        if len(orig_colors) == 1 : orig_colors = orig_colors*onearr

        xyz = np.transpose(faceCenterCoor)
        if usexyz is False :  xyz = self.coor_convert(xyz)
        abc = xyz
        rst = np.array(operation(abc))

        # FutDev: better way to do this loop?
        fvi = []
        fcol = []
        for i in range(len(faceCenterCoor)) :
            if rst[i] :
                fvi.append(self.fvIndices[i])
                fcol.append(orig_colors[i])

        # this permits multiple clips without changing original, so restore is available.
        if self.restoredClip is None :
            self.restoredClip = { 'faces' : self.fvIndices, 'colors' : orig_colors }

        self.fvIndices = fvi
        self.set_color(fcol)        
        v = self.vertexCoor
        self.set_verts(v[fvi])
        return self

    def restore_from_clip(self) :
        """Reset face and colors prior to any clip operation."""
        
        # FutDev: this method needs thorough testing.
        if self.restoredClip is None :
            warnings.warn('There is no clipped surface to restore.')
            return self
        fvi = self.restoredClip['faces']
        fcol = self.restoredClip['colors']
        self.restoredClip = None
        self.fvIndices = fvi
        self.set_color(fcol)        
        v = self.vertexCoor
        self.set_verts(v[fvi])
        return self

    def set_surface_alpha(self, alpha, constant=False) :
        """
        Adjust the alpha value of the surface face colors.

        Parameters
        ----------
        alpha : scalar
              Alpha is in the range 0 to 1.

        constant : bool { True, False }, optional, False
            If False, face color values are multiplied by
            alpha.  If True, all face colors alpha channels
            are assigned to alpha.

        Returns
        -------
        self : surface object

        """

        if (alpha<0.0) or (alpha>1.0) :
            raise ValueError('surface alpha must be between 0 and 1, found {}'.format(alpha))

        orig_colors = self._facecolors
        if len(orig_colors) == 1 :
            onearr = np.ones( len(self.fvIndices) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       
        if constant :
            np.put_along_axis(orig_colors, np.array([[3]]), alpha, axis=1)
            colorMap = orig_colors
        else:
            colorMap = np.multiply(orig_colors,np.array([1.0,1.0,1.0,alpha]))
        self.set_color(colorMap)       
        return self   
   
    def shade(self, depth=0, direction=None, contrast=None) :
        """
        Reduce surface HSV color Value based on face normals.
        
        The dot product of face normals with the illumination direction 
        is used to adjust face HSV color value.

        Parameters
        ----------
        depth : scalar, optional, default: 0
            Minimum color value of shaded surface faces.
            Depth value ranges from 0 to 1.

        direction : array-like, optional, default: (1,0,1)
            A xyz vector pointing to the illumination source. 

        contrast : scalar, optional, default: 1
            Shading contrast adjustment from low to high with a value
            of 1 for linear variations with the normal direction.
            Contrast value ranges from 0.1 to 3.     
        
        Returns
        -------
        self : surface object

        """

        if direction is None : direction = self._ILLUM
        if np.any( (depth<0) | (depth>1) ) :
            raise ValueError('depth values, {}, must be between 0 and 1'.format(depth))
        if contrast is not None:
            if (contrast<0.1 or contrast>3) :
                raise ValueError('contrast must be between 0.1 and 3. , found {}'.format(contrast))
    
        verts = self.vertexCoor[ self.fvIndices ]
        norms = self._get_face_normals(verts)
        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        dprod = np.dot( norms, unitDirection )

        # 0<d<1, heuristic function for 
        # normalized domain of effective dot product, d=f(dprod)
        d = 0.5*( 1 + dprod )
        if contrast is not None :
            d =  0.5*( 1 + np.sign(dprod)*np.power( np.abs(dprod) , 1/contrast) )      
        
        # extract HSV, leaving alpha unchanged.
        orig_colors = self._facecolors
        onearr = np.ones( len(norms) )[:, np.newaxis]
        if len(orig_colors) == 1 : orig_colors = orig_colors*onearr
        alphas = orig_colors[:,3:4]
        fc_less_alpha = orig_colors[:,:3]
        hsv_vals = cm.colors.rgb_to_hsv(fc_less_alpha)

        # adjust color value with multiplier.
        cvm = ((1-depth)*d + depth*np.ones(len(d)))[:, np.newaxis]
        hsv_vals[:,2] = (cvm*hsv_vals[:,2:3])[:,0]
        rgb_vals =  cm.colors.hsv_to_rgb(hsv_vals)
        shade_colors = np.concatenate((rgb_vals, alphas), axis=1)
        self.set_color(shade_colors)
        return self

    def hilite(self, height=1, direction=None, focus=None) :
        """
        Increase surface HSV color Value and reduce Saturation, based on face normals.

        The dot product of face normals with the illumination direction 
        is used to adjust face HSV color value and saturation.
        Faces having normals with a negative component in the direction of
        illumination are not affected. The reduction of color saturation is
        proportional to the increase in color value.

        Parameters
        ----------
        height : scalar, optional, default: 1
            Maximum color value of highlighted surface faces.
            Height value ranges from 0 to 1.

        direction : array-like, optional, default: (1,0,1)
            A xyz vector pointing to the illumination source. 

        focus : scalar, optional, default: 1
            Highlighting focus adjustment from low to high with a value
            of 1 for quadratic variations with the normal direction.
            Focus value ranges from 0.1 to 3.     
        
        Returns
        -------
        self : surface object

        """

        if direction is None : direction = self._ILLUM
        if np.any( (height<0) | (height>1) ) :
            raise ValueError('height values, {}, must be between 0 and 1'.format(height))
        if focus is not None:
            if (focus<0.1 or focus>3) :
                raise ValueError('focus must be between 0.1 and 3. , found {}'.format(focus))
    
        verts = self.vertexCoor[ self.fvIndices ]
        norms = self._get_face_normals(verts)
        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        dprod = np.dot( norms, unitDirection )

        # 0<d<1, heuristic function for 
        # normalized domain of effective dot product, d=f(dprod)
        d = np.where(dprod<0 , 0, dprod)       
        if focus is None : focus = 1
        d0 = np.power(focus,1/2)/2.0
        y = (d-d0)/(1-d0)
        d = np.where(y<0 , 0, y)
        d = np.power(d,1.0+2*focus)

        # extract HSV, leaving alpha unchanged.
        orig_colors = self._facecolors
        onearr = np.ones( len(norms) )[:, np.newaxis]
        if len(orig_colors) == 1 : orig_colors = orig_colors*onearr
        alphas = orig_colors[:,3:4]
        fc_less_alpha = orig_colors[:,:3]
        hsv_vals = cm.colors.rgb_to_hsv(fc_less_alpha)

        # adjust color saturation and value, linear with normalized domain.
        hd = height*d
        omhd = (np.ones( len(norms) ) - hd)
        omhd = omhd[:,np.newaxis]
        hsv_vals[:,1] = (omhd*hsv_vals[:,1:2])[:,0]
        hsv_vals[:,2] = (omhd*hsv_vals[:,2:3])[:,0]
        hsv_vals[:,2] += hd
        rgb_vals =  cm.colors.hsv_to_rgb(hsv_vals)
        hilite_colors = np.concatenate((rgb_vals, alphas), axis=1)
        self.set_color(hilite_colors)
        return self

    def transform(self, rotate=None, scale=1.0, translate=[0,0,0] )  :
        """
        Linear transformation of the surface object.

        Parameters
        ----------
        rotate :  array-like, optional, default: None
            A 3 by 3 rotation matrix.

        scale : 3D array or scalar, optional, default: 1.0
            Multipliers of the object coordinates, in
            xyz coordinates. If a single scalar, three
            multipliers are the same value.
        
        translate : array-like, optional, default: [0,0,0]
            A xyz translation vector of object origin. 

        Returns
        -------
        self : surface object

        """
        if rotate is None : rotate = np.identity(3)
        if isinstance(scale,(int,float)) :
            scale =  np.full(3,scale).astype(float)
        #.....................................    
        def isOrth(M) :
            M = np.array(M)
            tol = 0.00001
            prod = np.dot(M,M.T)
            np.fill_diagonal(prod,0.0)
            return not (np.abs(prod)>tol).any()
        #.....................................    
        trans = lambda V,S,T : np.add(  np.dot( np.multiply(V,S),rotate) ,T)
        #.....................................    
        if not isOrth(rotate) :
            warnings.warn('Transform rotate matrix is NOT Orthoginal.')

        v = trans( self.vertexCoor, scale, translate  )
        if self.vector_field is not None : 
            nScal = np.divide(scale, np.linalg.norm(scale)   )
            self.vector_field = trans( self.vector_field, nScal, [0,0,0]  )
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( v[self.fvIndices] )
        self._surfaceShapeChanged = True
        self.scale = scale
        self.translation = translate
        self.rotation = rotate
        return self

    def get_transformAxis(self, lenmult=1.0, width=1.5, colors=None) :
        """
        Line3DCollection of the 'last' transformed surface coordinate axis.
        
        Parameters
        ----------
        lenmult : scalar or 3D array, optional, default: 1
            Scalar multiplier of the three coordinate axis.

        width : scalar, optional, default: 1.5
            Line width of the coordinate axis.   

        colors : a Color or 3D array of Colors, optional, default: ['r','g','b']

        Returns
        -------
        Line3DCollection

        """
        scale = np.multiply(lenmult,self.scale)
        dDelta = np.multiply(np.identity(3).tolist(),scale).tolist()
        npdot = np.dot(dDelta,self.rotation)
        npdot2 = np.expand_dims(npdot,axis=1)
        lines = np.insert(npdot2,0,[0,0,0],axis=1)
        lines = np.add(lines,self.translation)
        if colors is None : 
            colorMap = ['r','g','b']
        else :
            colorMap = colors
        axisCol = Line3DCollection(lines, colors=colorMap, linewidths=width)
        return axisCol      

    def svd(self, data, pc=None ) :
        """
        Transform surface geometry based on a data set.
        
        Parameters
        ----------
        data : N x 3 array of x,y,z data coordinates.
        
        pc : scalar, optional, default: None
            Sets the minimum size of the surface based on the
            relative number of data.  Value range from 0 to 1, 
            with 1 for 100%.  If None, one standard deviation
            of the data sets the surface size.
            
        Returns
        -------
        disArr : array of N floats, normalized to the surface size.
        
        t : a 4D array containing surface and data properties
            [standard deviation, rotation matrix, scaling, translation]
        
        """

        # FutDev: This is a 'proof-of-concept' for data visulaizations

        if self._isCompositeSurface :
            warnings.warn('Data operation not available for combined shapes.')
            return self
        data_mean = data.mean(axis=0)
        data_centered = data - data_mean

        Ua, sa, Vta = np.linalg.svd(data_centered)
        dataAve = np.mean(sa)       
        scale = sa/dataAve
        tranData = np.divide(np.dot(data_centered,Vta.T),scale )       
        disArr=np.linalg.norm(tranData,axis=1)      
        disArrSum = np.linalg.norm(disArr)
        lendisArr = len(disArr)
        refLen = disArrSum/np.sqrt(lendisArr)        
        s=refLen
        disArr  = np.divide(disArr,refLen)
        
        if pc is not None:
            dlen = len(disArr)
            dindex = int(dlen*pc)
            if dindex >= dlen : dindex = dlen - 1
            if dindex <  0    : dindex = 0
            s = np.sort(disArr)[dindex]
            disArr = disArr/s
            s = s*refLen

        self.transform(Vta,s*scale,data_mean) 
        return disArr, [refLen,Vta,s*scale,data_mean]


class PlanarSurface(Surface3DCollection) :
    """
    Flat square 3D surface in Cartesian coordinates.

    Methods are inherited from the Surface3DCollection.
    Planar surface geometries are created in this subclass.

    """
    
    @staticmethod
    def _get_planar_surfaces_dictionary() :
        """ Constructs the dictionary of base planar surface geometries. """

        surfacesDictionary = {}
        baseVert = []
        baseVert.append([0,0,0])

        baseVert.append([1,1,0])        
        baseVert.append([-1,1,0])        
        baseVert.append([-1,-1,0])        
        baseVert.append([1,-1,0])

        baseVert.append([1,0,0])        
        baseVert.append([0,1,0])
        baseVert.append([-1,0,0])
        baseVert.append([0,-1,0])
        qVert = []
        for i in range(0,5) : qVert.append(baseVert[i])

        quadFaces=( 
            [0,1,2],[0,2,3],[0,3,4],[0,4,1] )
        octFaces1=(
            [0,5,6],[0,6,7],[0,7,8],[0,8,5],
            [5,1,6],[6,2,7],[7,3,8],[8,4,5] )
        octFaces2=(
            [0,5,1],[0,1,6],[0,6,2],[0,2,7],
            [0,7,3],[0,3,8],[0,8,4],[0,4,5] )

        Fq = len(quadFaces)
        F1 = len(octFaces1)
        F2 = len(octFaces2)
        quad =  { 'baseFaceIndices' : quadFaces , 'baseVerticesCoor' : qVert, 'fevParam' : [Fq,2,1] }
        oct1 =  { 'baseFaceIndices' : octFaces1 , 'baseVerticesCoor' : baseVert, 'fevParam' : [F1,4,1] }
        oct2 =  { 'baseFaceIndices' : octFaces2 , 'baseVerticesCoor' : baseVert, 'fevParam' : [F2,4,1] }

        surfacesDictionary['quad'] = quad
        surfacesDictionary['oct1'] = oct1
        surfacesDictionary['oct2'] = oct2

        # ----------------------------------------------------------
        mr = 0.01
        MR = mr*np.sqrt(2)
        qVertS=copy.copy( surfacesDictionary['quad']['baseVerticesCoor'] )
        qVertS[0] =   [ MR,  0, 0]
        qVertS.append([  0, MR, 0])
        qVertS.append([-MR,  0, 0])
        qVertS.append([0  ,-MR, 0])

        baseVertS = copy.copy(surfacesDictionary['oct1']['baseVerticesCoor'])
        baseVertS[0] =   [ mr, mr, 0]
        baseVertS.append([-mr, mr, 0])
        baseVertS.append([-mr,-mr, 0])
        baseVertS.append([ mr,-mr, 0])

        quadFacesS=( 
            [5,1,2],[6,2,3],[7,3,4],[0,4,1],
            [1,5,0],[2,6,5],[3,7,6],[4,0,7] )

        octFaces1S=(
            [ 0, 5, 6],[ 9, 6, 7],[10, 7, 8],[11, 8, 5],
            [ 5, 1, 6],[ 6, 2, 7],[ 7, 3, 8],[ 8, 4, 5],
            [ 5, 0,11],[ 6, 9, 0],[ 7,10, 9],[ 8,11,10] )

        octFaces2S=(
            [ 0, 5, 1],[ 0, 1, 6],[ 9, 6, 2],[ 9, 2, 7],
            [10, 7, 3],[10, 3, 8],[11, 8, 4],[11, 4, 5],
            [ 5, 0,11],[ 6, 9, 0],[ 7,10, 9],[ 8,11,10] )


        FqS = len(quadFacesS)
        F1S = len(octFaces1S)
        F2S = len(octFaces2S)
        quadS =  { 'baseFaceIndices' : quadFacesS , 'baseVerticesCoor' : qVertS, 'fevParam' : [FqS,2,1] }
        oct1S =  { 'baseFaceIndices' : octFaces1S , 'baseVerticesCoor' : baseVertS, 'fevParam' : [F1S,4,1] }
        oct2S =  { 'baseFaceIndices' : octFaces2S , 'baseVerticesCoor' : baseVertS, 'fevParam' : [F2S,4,1] }

        surfacesDictionary['quad_c'] = quadS
        surfacesDictionary['oct1_c'] = oct1S
        surfacesDictionary['oct2_c'] = oct2S

        return surfacesDictionary
    
    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point between two points in the xy plane."""

        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        return mid
    
    _base_surfaces = _get_planar_surfaces_dictionary.__func__()
    _default_base = 'quad'  # default assignment is indicated in __init Docstring.

    @classmethod
    def fev(cls,rez, basetype=None) :
        """
        Calculates the number of faces, edges, and vertices of a PlanarSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """

        return super()._fev(rez,basetype,cls)

    @staticmethod
    def meshgrid(X,Y,Z,revnormal=False) :
        """Create a PlanarSurface object based on a meshgrid."""

        vertCoor = np.reshape(np.stack((X,Y,Z), axis=2 ), (-1,3) )
        lenY, lenX = Z.shape
        mtr = np.reshape( np.arange(lenX*lenY)  , (lenY,lenX))
        A = mtr[:,0:lenX-1][0:lenY-1,:].flatten()
        B = A + 1
        C = A + lenX +1
        D = A + lenX
        if revnormal :
            t = B
            B = D
            D = t
        faceIndices = np.stack( (A,B,C,D) ).T

        obj = PlanarSurface( )
        obj.vertexCoor = vertCoor
        obj.fvIndices = faceIndices
        obj._normal_scale = _normLength( len(faceIndices), 0, 0, 4 )
        obj._set_geometric_bounds()
        obj._surfaceShapeChanged = True        

        v = obj.vertexCoor
        fvi = obj.fvIndices
        obj.set_verts(v[fvi])
        obj.triangulate()

        return obj

    def __init__(self, rez=0, basetype=None, *args, **kwargs):
        """
        Create flat square surface of 2 x 2 units.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'quad','oct1','oct2'}, optional, default: 'quad'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        *args, **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """

        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        if basetype not in self._base_surfaces :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
        baseSurfObj = self._base_surfaces[basetype]
        baseVcoor = baseSurfObj['baseVerticesCoor']
        baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
        indexObj,vertexCoor = self.triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
        faceIndices = indexObj['face']
        edgeIndices = indexObj['edge']

        super().__init__(vertexCoor, faceIndices, edgeIndices, *args, **kwargs)

        self._normal_scale = _normLength( len(faceIndices), 0, 0, 4 )
        self._rez = rez
        self._basetype = basetype
        return

    def scale_dataframe(self,X,Y,Z) :
        """
        Linear scale and translate the surface.

        Used for scaling a surface geometry based on a datagrid.

        Parameters
        ----------
        X, Y, Z : 2d arrays
            Minimum and maximum values of the arrays are used to scale
            and translate the surface from an intial domain of
            [ (-1,1), (-1,1), (0,1) ]

        Returns
        -------
        self : PlanarSurface object

        """

        Xmin, Xmax = np.amin(X), np.amax(X)
        Ymin, Ymax = np.amin(Y), np.amax(Y)
        Zmin, Zmax = np.amin(Z), np.amax(Z)
        Xscale = (Xmax-Xmin)/2
        X0 = (Xmax+Xmin)/2
        Yscale = (Ymax-Ymin)/2
        Y0 = (Ymax+Ymin)/2
        Zscale = Zmax-Zmin
        Z0 = Zmin
        self.transform(scale=[Xscale,Yscale,Zscale],translate=[X0,Y0,Z0])
        return self

class PolarSurface(Surface3DCollection) :
    """
    Flat disk 3D surface in polar coordinates.


    Methods are inherited from the Surface3DCollection.
    Polar surface geometries are created in this subclass.

    """

    @staticmethod
    def _get_polar_surfaces_dictionary() :
        """Constructs the dictionary of base polar surface geometries."""

        surfacesDictionary = {}
        soff = 0.0001
        x0 = 0.5
        y0 = np.sqrt(3)/2
    
        baseVert = []
        baseVert.append([0,0,0])
        baseVert.append([1,0,0])
        baseVert.append([x0,y0,0])
        baseVert.append([-x0,y0,0])
        baseVert.append([-1,0,0])
        baseVert.append([-x0,-y0,0])
        baseVert.append([x0,-y0,0])

        qVert = []
        qVert.append([0,0,0])
        qVert.append([1,0,0])
        qVert.append([0,1,0])
        qVert.append([-1,0,0])
        qVert.append([0,-1,0])

        baseVertS = []
        x0S, y0S = soff*y0, soff*x0
        baseVertS.append([1,-y0S,0])
        baseVertS.append([1, y0S,0])
        baseVertS.append([x0,y0,0])
        baseVertS.append([-x0,y0,0])
        baseVertS.append([-1,0,0])
        baseVertS.append([-x0,-y0,0])
        baseVertS.append([x0,-y0,0])
        baseVertS.append([ x0S,  y0S,0])
        baseVertS.append([   0, soff,0])
        baseVertS.append([-x0S,  y0S,0])
        baseVertS.append([-x0S, -y0S,0])
        baseVertS.append([   0,-soff,0])
        baseVertS.append([ x0S, -y0S,0])

        qVertS = []
        qVertS.append([1,-soff,0])
        qVertS.append([1, soff,0])
        qVertS.append([0,1,0])
        qVertS.append([-1,0,0])
        qVertS.append([0,-1,0])
        qVertS.append([ soff, soff,0])
        qVertS.append([-soff ,soff,0])
        qVertS.append([-soff,-soff,0])
        qVertS.append([ soff,-soff,0])

        baseFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,1])
        quadFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,1])

        baseFacesS=( [7,1,2], [8,2,3], [9,3,4],  [10,4,5],  [11,5,6], [12,6,0],
                     [2,8,7], [3,9,8], [4,10,9], [5,11,10], [6,12,11] )

        quadFacesS=( [5,1,2],[6,2,3],[7,3,4],[8,4,0],
                    [2,6,5],[3,7,6],[4,8,7] )

        Fs = len(quadFaces)
        Fh = len(baseFaces)
        FsS = len(quadFacesS)
        FhS = len(baseFacesS)
        square =  { 'baseFaceIndices' : quadFaces , 'baseVerticesCoor' : qVert, 'fevParam' : [Fs,2,1] }
        hexag =  { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fh,3,1] }
        squareS =  { 'baseFaceIndices' : quadFacesS , 'baseVerticesCoor' : qVertS, 'fevParam' : [FsS,4.5,1] }
        hexagS =  { 'baseFaceIndices' : baseFacesS , 'baseVerticesCoor' : baseVertS, 'fevParam' : [FhS,6.5,1] }
     
        surfacesDictionary['squ'] = square
        surfacesDictionary['hex'] = hexag
        surfacesDictionary['squ_s'] = squareS
        surfacesDictionary['hex_s'] = hexagS
        # ... 'special case' only for fev method for polar surfaces.
        surfacesDictionary['squ_c'] = { 'fevParam' : [16,5,1] }
        surfacesDictionary['hex_c'] = { 'fevParam' : [12,4,1] }

        return surfacesDictionary

    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point at the mid radial position between two points in a plane."""
        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        a = (np.linalg.norm(vectA)+np.linalg.norm(vectB))/2
        b = np.linalg.norm(mid)
        mid = np.multiply(mid,a/b)
        return mid

    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Override map surface to rectagular unit coordinates (0,1)."""
        x,y,z = xyz
        theta = np.arctan2(y,x)
        theta = np.where( theta<0, 2.0*np.pi + theta ,theta)
        a = theta/(2*np.pi)
        r = np.sqrt(x*x+y*y)
        b = 1-r
        a = np.clip(a,0.0,1.0) 
        b = np.clip(b,0.0,1.0) 
        return [a,b]

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """
        Override coordinate transformation.

        Tranformation between polar and Cartesian coordinates.
        The angular coordinate is in radians with domain 0 to 2pi.

        Parameters
        ----------
        xyz : 3xN array
            N number of vectors in either polar or cartesian coordinates,
            including the z coordinate( which remains unchanged).

        tocart : { True, False }, default: False
            If True, input is polar and output is cartesian.
            If False, input is cartesian and output is polar.

        Returns
        -------
        list: 3xN array
            Vector of N transformed coordinates. 

        """

        if tocart :
            r, theta, z = xyz
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            abc = [x,y,z]
        else:
            x,y,z = xyz
            R = np.sqrt( x*x + y*y )
            theta = np.arctan2(y,x)
            theta = np.where(theta<0,theta+2*np.pi,theta)
            abc = [R,theta,z]
        return abc

    _base_surfaces = _get_polar_surfaces_dictionary.__func__()
    _default_base = 'hex'  # default assignment is indicated in __init Docstring.

    @classmethod
    def _uvw_to_xyz( cls,uvw, rtz) :  
        """Override rotational transformation at a polar coordinate."""
        return super()._disp_to_xyz( uvw, rtz[1] )

    @classmethod
    def fev(cls,rez, basetype=None) :
        """
        Calculates number of faces, edges, and vertices of a PolarSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.
        
        """

        return super()._fev(rez,basetype,cls)

    @staticmethod
    def _flat_split_cyclinder(rez,basetype,minrad, *args, **kwargs) :
        """
        Transform cylindrical base surface to polar coordinates. 

        note :  This method provides an alternative to the 
                poler basetypes of 'hex_s' and 'squ_s'.
                The grids constructed here may provide smoother
                geometries depending on the mapping functions.
        """

        # .................................................
        def cyl2pol(rtz):
            r,t,z = rtz
            b = (1.0+minrad)/2.0
            m = -(1.0-minrad)/2.0
            R = m*z + b
            Z = 0*z
            return R,t,Z
        # .................................................
        if basetype is 'squ_c' : cylbase = 'squ_s' 
        if basetype is 'hex_c' : cylbase = 'tri_s'
        cyl_obj = CylindricalSurface(rez,cylbase, *args, **kwargs).map_geom_from_op(cyl2pol)
        baseFaceVertexIndices = cyl_obj._base_surfaces[cylbase] ['baseFaceIndices']    
        fvIndices = cyl_obj.fvIndices
        vertexCoor = cyl_obj.vertexCoor
        evIndices = cyl_obj.evIndices
        return vertexCoor, fvIndices, evIndices, baseFaceVertexIndices,

    def __init__(self, rez=0, basetype=None, minrad=None, *args, **kwargs):
        """
        Create flat disk surface of unit radius.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'squ','hex','squ_s','hex_s','squ_c','hex_c'}, optional, default: 'hex'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.

        minrad : scalar, optional, default: 0.01
            For basetypes 'squ_c' and 'hex_c, the minimum
            distance from the z-axis of any vertex coordinate.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        *args, **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """

        if minrad is None : minrad = _MINRAD
        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        # ... check for 'special case' (ugly code but needed)
        if basetype is 'squ_c' or basetype is 'hex_c' :
            vertexCoor, faceIndices, edgeIndices, bfvi = self._flat_split_cyclinder(rez,basetype,minrad, *args, **kwargs)
            baseFaceVertexIndices = bfvi
        else :
            if basetype not in self._base_surfaces :
                raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
            baseSurfObj = self._base_surfaces[basetype]
            baseVcoor = baseSurfObj['baseVerticesCoor']
            baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
            indexObj,vertexCoor = self.triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
            faceIndices = indexObj['face']
            edgeIndices = indexObj['edge']

        super().__init__(vertexCoor, faceIndices, edgeIndices, *args, **kwargs)

        self._normal_scale = _normLength( len(faceIndices), 1.2, 1.55, np.pi )
        self._rez = rez
        self._basetype = basetype
        return

class CylindricalSurface(Surface3DCollection):
    """
    Cylindrical 3D surface in cylindrical coordinates.

    Methods are inherited from the Surface3DCollection.
    Cylindrical surface geometries are created in this subclass.

    """

    @staticmethod
    def _get_cylindrical_surfaces_dictionary() :
        """Constructs the dictionary of base cylindrical surface geometries. """

        surfacesDictionary = {}
        def construct_tricylinder() :
            y0 = np.sqrt(3.0)/2.0
            x0 = 0.5   
            baseVert = []
            baseVert.append([1,0,0])
            baseVert.append([-x0,-y0,0])
            baseVert.append([-x0,y0,0])
            baseVert.append([-1,0,1])
            baseVert.append([x0,-y0,1])
            baseVert.append([x0,y0,1])
            baseVert.append([-1,0,-1])
            baseVert.append([x0,-y0,-1])
            baseVert.append([x0,y0,-1])
            baseFaces=( 
                [0,2,5],[2,1,3],[1,0,4],
                [0,5,4],[2,3,5],[1,4,3], 
                [0,8,2],[2,6,1],[1,7,0],
                [0,7,8],[2,8,6],[1,6,7])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,3,0] }

        def construct_tri2cylinder() :
            y0 = np.sqrt(3.0)/2.0
            x0 = 0.5   
            baseVert = []

            baseVert.append([-1,0,0])
            baseVert.append([x0,-y0,0])
            baseVert.append([x0,y0,0])

            baseVert.append([1,0,1])
            baseVert.append([-x0,-y0,1])
            baseVert.append([-x0,y0,1])

            baseVert.append([1,0,-1])
            baseVert.append([-x0,-y0,-1])
            baseVert.append([-x0,y0,-1])
            baseFaces=( 
                [0,4,5],[1,3,4],[2,5,3],
                [0,8,7],[1,7,6],[2,6,8],
                [0,7,4],[1,6,3],[2,8,5],
                [0,5,8],[1,4,7],[2,3,6])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,3,0] }

        def construct_tri2cylinder_sliced() :
            y0 = np.sqrt(3.0)/2.0
            x0 = 0.5   
            baseVert = []

            yof = 0.0001

            baseVert.append([-1,0,0])
            baseVert.append([x0,-y0,0])
            baseVert.append([x0,y0,0])

            baseVert.append([1,yof,1])
            baseVert.append([-x0,-y0,1])
            baseVert.append([-x0,y0,1])

            baseVert.append([1,yof,-1])
            baseVert.append([-x0,-y0,-1])
            baseVert.append([-x0,y0,-1])

            baseVert.append([1,-yof,1])
            baseVert.append([1,-yof,-1])
            
            baseFaces=( 
            #    [0,4,5],[1,3,4],[2,5,3],
            #    [0,8,7],[1,7,6],[2,6,8],
            #    [0,7,4],[1,6,3],[2,8,5],
                [0,4,5],[1,9,4],[2,5,3],
                [0,8,7],[1,7,10],[2,6,8],
                [0,7,4],[1,10,9],[2,8,5],

                [0,5,8],[1,4,7],[2,3,6])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,4,1] }

        def construct_quadcylinder() :
            xy0 = 1.0/np.sqrt(2.0)
            baseVert = []
            baseVert.append([1,0,0])
            baseVert.append([0,-1,0])
            baseVert.append([-1,0,0])
            baseVert.append([0,1,0])
            baseVert.append([xy0,xy0,1])
            baseVert.append([-xy0,xy0,1])
            baseVert.append([-xy0,-xy0,1])
            baseVert.append([xy0,-xy0,1])
            baseVert.append([xy0,xy0,-1])
            baseVert.append([-xy0,xy0,-1])
            baseVert.append([-xy0,-xy0,-1])
            baseVert.append([xy0,-xy0,-1])

            baseFaces=( 
                [1,0,7],[2,1,6],[3,2,5],[0,3,4],
                [4,3,5],[5,2,6],[6,1,7],[7,0,4],
                [0,8,3],[3,9,2],[2,10,1],[1,11,0],
                [0,11,8],[3,8,9],[2,9,10],[1,10,11] )

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,4,0] }

        def construct_quad2cylinder() :
            xy0 = 1.0/np.sqrt(2.0)
            baseVert = []

            baseVert.append([xy0,xy0,0])
            baseVert.append([-xy0,xy0,0])
            baseVert.append([-xy0,-xy0,0])
            baseVert.append([xy0,-xy0,0])

            baseVert.append([1,0,1])
            baseVert.append([0,-1,1])
            baseVert.append([-1,0,1])
            baseVert.append([0,1,1])

            baseVert.append([1,0,-1])
            baseVert.append([0,-1,-1])
            baseVert.append([-1,0,-1])
            baseVert.append([0,1,-1])


            baseFaces=( 
                [0,7,4],[1,6,7],[2,5,6],[3,4,5],
                [0,8,11],[1,11,10],[2,10,9],[3,9,8],
                [0,4,8],[1,7,11],[2,6,10],[3,5,9],
                [0,11,7],[1,10,6],[2,9,5],[3,8,4] )

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,4,0] }

        def construct_quad2cylinder_sliced() :
            xy0 = 1.0/np.sqrt(2.0)
            baseVert = []

            baseVert.append([xy0,xy0,0])
            baseVert.append([-xy0,xy0,0])
            baseVert.append([-xy0,-xy0,0])
            baseVert.append([xy0,-xy0,0])

            yof = 0.0001
            baseVert.append([1,yof,1])
            baseVert.append([0,-1,1])
            baseVert.append([-1,0,1])
            baseVert.append([0,1,1])

            baseVert.append([1,yof,-1])
            baseVert.append([0,-1,-1])
            baseVert.append([-1,0,-1])
            baseVert.append([0,1,-1])

            baseVert.append([1,-yof,1])
            baseVert.append([1,-yof,-1])

            baseFaces=( 
            #    [0,7,4],[1,6,7],[2,5,6],[3,4,5],
            #    [0,8,11],[1,11,10],[2,10,9],[3,9,8],
            #    [0,4,8],[1,7,11],[2,6,10],[3,5,9],
            #    [0,11,7],[1,10,6],[2,9,5],[3,8,4] )

                [0,7,4],[1,6,7],[2,5,6],[3,12,5],
                [0,8,11],[1,11,10],[2,10,9],[3,9,13],
                [0,4,8],[1,7,11],[2,6,10],[3,5,9],
                [0,11,7],[1,10,6],[2,9,5],[3,13,12] )


            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,5,1] }

        surfacesDictionary['tri'] = construct_tricylinder()
        surfacesDictionary['squ'] = construct_quadcylinder()
        surfacesDictionary['tri2'] = construct_tri2cylinder()
        surfacesDictionary['squ2'] = construct_quad2cylinder()
        surfacesDictionary['tri_s'] = construct_tri2cylinder_sliced()
        surfacesDictionary['squ_s'] = construct_quad2cylinder_sliced()
        return surfacesDictionary
    
    @staticmethod
    def _midVectorFun(vectA, vectB) :  
        """Mid-point between two points on a unit cylinder."""
        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        z = np.dot(mid,[0,0,1])
        v = np.subtract(mid,[0,0,z])
        xy =  v/np.linalg.norm(v)
        return np.add(xy,[0,0,z])

    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Override map surface to rectagular unit coordinates (0,1)."""
        x,y,z = xyz
        theta = np.arctan2(y,x)
        theta = np.where( theta<0, 2.0*np.pi + theta ,theta)
        a = theta/(2*np.pi)
        b =  z/2 + 0.5 
        a = np.clip(a,0.0,1.0) 
        b = np.clip(b,0.0,1.0) 
        return [a,b]

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """
        Override coordinate transformation.

        Tranformation between cylindrical and Cartesian coordinates.
        The angular coordinate is in radians with domain 0 to 2pi.

        Parameters
        ----------
        xyz : 3xN
            N number of vectors in either cylindrical or cartesian coordinates,
            including the z coordinate( which remains unchanged).

        tocart : { True, False }, default: False
            If True, input is cylindrical and output is cartesian.
            If False, input is cartesian and output is cylindrical.

        Returns
        -------
        list: 3xN array
            Vector of N transformed coordinates. 

        """

        if tocart :
            r, theta, z = xyz
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            abc = [x,y,z]
        else:
            x,y,z = xyz
            R = np.sqrt( x*x + y*y )
            theta = np.arctan2(y,x)
            theta = np.where(theta<0,theta+2*np.pi,theta)
            abc = [R,theta,z]
        return abc

    _base_surfaces = _get_cylindrical_surfaces_dictionary.__func__()
    _default_base = 'tri'  # default assignment is indicated in __init Docstring.

    @classmethod
    def _uvw_to_xyz(cls, uvw, rtz) :
        """Override rotational tranformation at a cylindrical coordinate.""" 
        return super()._disp_to_xyz( uvw, rtz[1], None )

    @classmethod
    def fev(cls,rez, basetype=None) :
        """
        Calculates number of faces, edges, and vertices of a CylindricalSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """

        return super()._fev(rez,basetype,cls)
               
    def __init__(self, rez=0, basetype=None, *args, **kwargs):
        """
        Create cylindrical surface of unit radius and length 2.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'tri','squ','tri2','squ2'}, optional, default: 'tri'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        *args, **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        """

        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        if basetype not in self._base_surfaces :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
        baseSurfObj = self._base_surfaces[basetype]
        baseVcoor = baseSurfObj['baseVerticesCoor']
        baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
        indexObj,vertexCoor = self.triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
        faceIndices = indexObj['face']
        edgeIndices = indexObj['edge']
        super().__init__(vertexCoor, faceIndices, edgeIndices, *args, **kwargs)

        self._normal_scale =  _normLength( len(faceIndices), 1.14, 1.31)
        self._rez = rez
        self._basetype = basetype
        self.disp_Vector = np.array( [1,1,0] )
        return None

class SphericalSurface(Surface3DCollection) :
    """
    Spherical 3D surface in spherical coordinates.

    Methods are inherited from the Surface3DCollection.
    Spherical surface geometries are created in this subclass.

    """

    @staticmethod
    def get_dodecahedron() :
        """Numpy arrays of dodecahedron vertices (20,3) and face vertex indices (12,5)."""
        
        v,f = SphericalSurface._get_platonic_solids_dictionary(True)
        return np.array(v), np.array(f)

    @staticmethod
    def _get_platonic_solids_dictionary(getArrays = False) :
        """Constructs the dictionary of base spherical surface geometries."""
        # ....................................................
        def get_dodecahedron() :
            
            zeta = 1/3
            ksi = np.sqrt(5)/3
            lam_squ = 8/9
            hEdge_squ = (1-ksi)/2

            Ax = 2/3
            Apx = 1/3
            Bx = np.sqrt(lam_squ-hEdge_squ)
            Cx = ksi 
            Dx = 1-Bx
            
            Apy = np.sqrt(3)/3
            By = np.sqrt(hEdge_squ)
            Cy = Apy
            Dy = By + Cy

            # -----------------------

            docTop = [ [0,0,1],
                [Ax, 0, ksi], [-Apx, Apy, ksi], [-Apx, -Apy, ksi],
                [Cx, Cy, zeta], [Dx, Dy, zeta], [-Bx, By, zeta],
                [-Bx, -By, zeta], [Dx, -Dy, zeta], [Cx, -Cy, zeta]
                ]
            docBottom = []
            for ver in docTop : docBottom.append( np.multiply(ver,[-1,1,-1]).tolist() )
            docbaseVert = np.concatenate((docTop,docBottom)).tolist()

            docbaseFaces = [
                [0,1,4,5,2], [0,2,6,7,3] , [0,3,8,9,1],
                [1,9,17,16,4], [2,5,15,14,6], [3,7,19,18,8],
                [10,12,16,17,13], [10,11,14,15,12], [10,13,18,19,11],
                [11,19,7,6,14],[13,17,9,8,18],[12,15,5,4,16]
            ]
            return docbaseVert, docbaseFaces
        def construct_tetrahedron() :
            cSqr = (np.tan(np.pi/3))*(np.tan(np.pi/3))
            zeta = (cSqr-1)/(2*cSqr)
            chi = np.sqrt(3*cSqr-1)*np.sqrt(cSqr+1)/(2*cSqr)
            x0 = chi*np.cos(np.pi/3)
            y0 = chi*np.sin(np.pi/3)
        
            baseVert = []
            baseVert.append([0,0,1])
            baseVert.append([chi,0,-zeta])
            baseVert.append([-x0,y0,-zeta])
            baseVert.append([-x0,-y0,-zeta])

            baseFaces=([0,1,2],[0,2,3],[0,3,1],[1,3,2])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_cube() :
            zeta = 1/3
            x0 = np.sqrt(2)/3
            y0 = x0*np.sqrt(3)
            chi = 2*x0
        
            cubeBaseVert = []
            cubeBaseVert.append([0,0,1])
            cubeBaseVert.append([chi,0,zeta])
            cubeBaseVert.append([-x0,y0,zeta])
            cubeBaseVert.append([-x0,-y0,zeta])
            cubeBaseVert.append([-chi,0,-zeta])
            cubeBaseVert.append([x0,-y0,-zeta])
            cubeBaseVert.append([x0,y0,-zeta])
            cubeBaseVert.append([0,0,-1])

            cubeBaseFaces=( 
                [0,1,6,2],[0,2,4,3],[0,3,5,1],
                [7,6,1,5],[7,4,2,6],[7,5,3,4] )
            
            baseFaces = []
            baseVert = cubeBaseVert.copy()

            for face in cubeBaseFaces :
                vcSum = [0,0,0]
                # get square face center
                for vI in face : vcSum = np.add(vcSum,cubeBaseVert[vI])
                vcSum = np.divide(vcSum,4)
                vcCenter = vcSum/np.linalg.norm(vcSum)
                vIndex = len(baseVert)
                baseVert.append(vcCenter)
                # subdivide cube face into 4 triangles
                for i in range(0,4) :
                    newFace = []
                    newFace.append(vIndex)
                    newFace.append(face[i])
                    newFace.append(face[ (i+1)%4])
                    baseFaces.append(newFace)

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_cubeS() :
            zeta = 1/3
            x0 = np.sqrt(2)/3
            y0 = x0*np.sqrt(3)
            chi = 2*x0

            soff = 0.01
            x0S, y0S = soff/2, soff*np.sqrt(3.0)/2.0
            w = y0S/2.
            cubeBaseVert = []
            cubeBaseVert.append([chi,-w,zeta])
            cubeBaseVert.append([chi, w,zeta])

            cubeBaseVert.append([-x0,y0,zeta])
            cubeBaseVert.append([-x0,-y0,zeta])
            cubeBaseVert.append([-chi,0,-zeta])
            cubeBaseVert.append([x0,-y0,-zeta])
            cubeBaseVert.append([x0,y0,-zeta])
            
            cubeBaseVert.append([ soff,   w, 1])                     
            cubeBaseVert.append([  x0S, y0S, 1])
            cubeBaseVert.append([ -x0S, y0S, 1])
            cubeBaseVert.append([ -soff,   0, 1])
            cubeBaseVert.append([ -x0S,-y0S, 1])
            cubeBaseVert.append([  x0S,-y0S, 1])
            cubeBaseVert.append([ soff,  -w, 1])

            cubeBaseVert.append([ soff,   w,-1])                     
            cubeBaseVert.append([  x0S, y0S,-1])
            cubeBaseVert.append([ -x0S, y0S,-1])
            cubeBaseVert.append([ -soff,   0,-1])
            cubeBaseVert.append([ -x0S,-y0S,-1])
            cubeBaseVert.append([  x0S,-y0S,-1])
            cubeBaseVert.append([ soff,  -w,-1])

            cubeBaseFaces=( 
                [ 7, 1, 6, 8], [ 8, 6, 2, 9], [ 9, 2, 4,10], [10, 4, 3,11], [11, 3, 5,12], [12, 5, 0,13], 
                [14,15, 6, 1], [15,16, 2, 6], [16,17, 4, 2], [17,18, 3, 4], [18,19, 5, 3], [19,20, 0, 5] )
            
            baseFaces = []
            baseVert = cubeBaseVert.copy()

            for face in cubeBaseFaces :
                vcSum = [0,0,0]
                # get square face center
                for vI in face : vcSum = np.add(vcSum,cubeBaseVert[vI])
                vcSum = np.divide(vcSum,4)
                vcCenter = vcSum/np.linalg.norm(vcSum)
                vIndex = len(baseVert)
                baseVert.append(vcCenter)
                # subdivide cube face into 4 triangles
                for i in range(0,4) :
                    newFace = []
                    newFace.append(vIndex)
                    newFace.append(face[i])
                    newFace.append(face[ (i+1)%4])
                    baseFaces.append(newFace)

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,8,1] }
        def construct_octahedron():
            baseVert = []
            baseVert.append([0,0,1])
            baseVert.append([1,0,0])
            baseVert.append([0,1,0])
            baseVert.append([-1,0,0])
            baseVert.append([0,-1,0])
            baseVert.append([0,0,-1])

            baseFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,1],
                [5,2,1],[5,3,2],[5,4,3],[5,1,4])

            Fo = len(baseFaces)

            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_octahedronS():
            soff = 0.01

            baseVertS = []
            x0S, y0S = soff/np.sqrt(2.0), soff/np.sqrt(2.0)
            baseVertS.append([1,-y0S,0])

            baseVertS.append([1, y0S,0])
            baseVertS.append([0,1,0])
            baseVertS.append([-1,0,0])
            baseVertS.append([0,-1,0])

            baseVertS.append([ x0S,  y0S, 1])
            baseVertS.append([-x0S,  y0S, 1])
            baseVertS.append([-x0S, -y0S, 1])
            baseVertS.append([ x0S, -y0S, 1])

            baseVertS.append([ x0S,  y0S,-1])
            baseVertS.append([-x0S,  y0S,-1])
            baseVertS.append([-x0S, -y0S,-1])
            baseVertS.append([ x0S, -y0S,-1])

            baseFacesS=( [5,1,2], [6,2,3], [7,3,4], [8,4,0],
                         [2,6,5], [3,7,6], [4,8,7],
                         [9,2,1], [10,3,2], [11,4,3], [12,0,4],
                         [2,9,10], [3,10,11], [4,11,12] )

            FoS = len(baseFacesS)

            return { 'baseFaceIndices' : baseFacesS , 'baseVerticesCoor' : baseVertS, 'fevParam' : [FoS,5,1] }
        def construct_dodecahedron() :
            docbaseVert, docbaseFaces = get_dodecahedron()
            baseFaces = []
            baseVert = docbaseVert.copy()
            
            for face in docbaseFaces :
                vcSum = [0,0,0]
                # get pentagon face center
                for vI in face : vcSum = np.add(vcSum,docbaseVert[vI])
                vcSum = np.divide(vcSum,4)
                vcCenter = vcSum/np.linalg.norm(vcSum)
                vIndex = len(baseVert)
                baseVert.append(vcCenter)
                # subdivide pentagon face into 5 triangles
                for i in range(0,5) :
                    newFace = []
                    newFace.append(vIndex)
                    newFace.append(face[i])
                    newFace.append(face[ (i+1)%5])
                    baseFaces.append(newFace)

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_icosahedron():
            cSqr = (np.tan(np.pi/5))*(np.tan(np.pi/5))
            zeta = (1-cSqr)/(2*cSqr)
            chi = np.sqrt(3*cSqr-1)*np.sqrt(cSqr+1)/(2*cSqr)
            thetaInc = 2*np.pi/5
            thetaOff = thetaInc/2

            baseVert = []
            baseVert.append([0,0,1])
            for i in range(0,5) :
                theta = i*thetaInc
                baseVert.append( [ chi*np.cos(theta), chi*np.sin(theta), zeta] )
            for i in range(0,5) :
                theta = thetaOff+i*thetaInc
                baseVert.append( [ chi*np.cos(theta), chi*np.sin(theta), -zeta] )
            baseVert.append([0,0,-1])

            baseFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,1],
                [6,2,1],[7,3,2],[8,4,3],[9,5,4],[10,1,5],
                [2,6,7],[3,7,8],[4,8,9],[5,9,10],[1,10,6],
                [11,7,6],[11,8,7],[11,9,8],[11,10,9],[11,6,10])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        # ....................................................
        if getArrays : return get_dodecahedron()
        surfacesDictionary = {}
        surfacesDictionary['tetra']  = construct_tetrahedron()
        surfacesDictionary['cube']   = construct_cube()
        surfacesDictionary['cube_s'] = construct_cubeS()
        surfacesDictionary['octa']   = construct_octahedron()
        surfacesDictionary['octa_s'] = construct_octahedronS()
        surfacesDictionary['dodeca'] = construct_dodecahedron()
        surfacesDictionary['icosa']  = construct_icosahedron()
        # ... 'special case' only for fev method for spherical surfaces.
        surfacesDictionary['octa_c'] = { 'fevParam' : [16,5,1] }
        surfacesDictionary['cube_c'] = { 'fevParam' : [12,4,1] }
        return surfacesDictionary

    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point between two points on a unit sphere."""
        mid = np.add(vectA,vectB)
        return mid/np.linalg.norm(mid)
    
    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Override map surface to rectagular unit coordinates (0,1)."""
        x,y,z = xyz
        t = np.arctan2(y,x)
        theta = np.where( t<0, 2.0*np.pi + t ,t)
        phi = np.arcsin(z)
        a = theta/(2*np.pi)
        b =  phi/np.pi + 0.5
        a = np.clip(a,0.0,1.0) 
        b = np.clip(b,0.0,1.0) 
        return [a,b]

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """
        Override coordinate transformation.

        Tranformation between spherical and Cartesian coordinates.
        The azimuthal coordinate is in radians with domain 0 to 2pi.
        The polar angle is in radians with domain 0 to pi.

        Parameters
        ----------
        xyz : 3xN array
            N number of vectors in either spherical or cartesian coordinates,

        tocart : { True, False }, default: False
            If True, input is spherical and output is cartesian.
            If False, input is cartesian and output is spherical.

        Returns
        -------
        list: 3xN array
            Vector of N transformed coordinates. 

        """

        if tocart :
            r, theta, phi = xyz
            x = r*np.sin(phi)*np.cos(theta)
            y = r*np.sin(phi)*np.sin(theta)
            z = r*np.cos(phi)
            abc = [x,y,z]
        else:
            x,y,z = xyz
            R = np.sqrt( x*x + y*y + z*z )
            phi = np.arccos(z/R)
            theta = np.arctan2(y,x)
            theta = np.where(theta<0,theta+2*np.pi,theta)
            abc = [R,theta,phi]
        return abc

    _base_surfaces = _get_platonic_solids_dictionary.__func__()
    _default_base = 'icosa'  # default assignment is indicated in __init Docstring.

    @classmethod
    def _uvw_to_xyz(cls, uvw, rtp) :
        """Override rotational tranformation at a spherical coordinate.""" 
        return super()._disp_to_xyz( uvw, rtp[1], rtp[2] )

    @classmethod
    def fev(cls, rez, basetype=None) :
        """
        Calculates number of faces, edges, and vertices of a SphericalSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """

        return super()._fev(rez,basetype,cls)

    @staticmethod
    def _spherical_split_cyclinder(rez,basetype,minrad, *args, **kwargs) :
        """
        Transform cylindrical base surface to spherical coordinates. 

        note :  This method provides an alternative to the 
                spherical basetypes of 'hex_s' and 'square_s'.
                The grids constructed here may provide smoother
                geometries depending on the mapping functions.
        """

        # .................................................
        def cyl2sph(rtz):
            r,t,z = rtz
            minphi = np.arccos(minrad)
            phi = z*minphi
            Z = np.sin(phi)
            R = np.cos(phi)
            return R,t,Z
        # .................................................
        if basetype is 'octa_c' : cylbase = 'squ_s' 
        if basetype is 'cube_c' :    cylbase = 'tri_s'
        cyl_obj = CylindricalSurface(rez,cylbase, *args, **kwargs).map_geom_from_op(cyl2sph)
        baseFaceVertexIndices = cyl_obj._base_surfaces[cylbase] ['baseFaceIndices']    
        fvIndices = cyl_obj.fvIndices
        vertexCoor = cyl_obj.vertexCoor
        evIndices = cyl_obj.evIndices
        return vertexCoor, fvIndices, evIndices, baseFaceVertexIndices,

    def __init__(self, rez=0, basetype=None, minrad=None, *args, **kwargs):
        """
        Create spherical surface of unit radius.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'tetra','octa','icosa','cube','dodeca', \
            'octa_s','cube_s','octa_c','cube_c'}, optional, default: 'icosa'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.
            All basetypes are based on the five platonic solids.

        minrad : scalar, optional, default: 0.01
            For basetypes 'octa_c' and 'cube_c, the minimum
            distance from the z-axis of any vertex coordinate.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        *args, **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """

        if minrad is None : minrad = _MINRAD
        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        # ... check for 'special case' (ugly code but needed)
        if basetype is 'octa_c' or basetype is 'cube_c' :
            vertexCoor, faceIndices, edgeIndices, bfvi = self._spherical_split_cyclinder(rez,basetype,minrad, *args, **kwargs)
            baseFaceVertexIndices = bfvi
        else :
            if basetype not in self._base_surfaces :
                raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
            baseSurfObj = self._base_surfaces[basetype]
            baseVcoor = baseSurfObj['baseVerticesCoor']
            baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
            indexObj,vertexCoor = self.triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
            faceIndices = indexObj['face']
            edgeIndices = indexObj['edge']
        
        super().__init__(vertexCoor, faceIndices, edgeIndices, *args, **kwargs)

        self._normal_scale = _normLength( len(faceIndices), 0.61, 1.52)
        self._rez = rez
        self._basetype = basetype
        self.disp_Vector = np.array( [1,1,1] )
        return None



class Vector3DCollection(Line3DCollection) :
    """
    Collection of 3D vectors represented as arrows.

    """

    def __init__(self, location, vect , alr=None, *args, **kwargs) :
        """
        Create a collection of 3D vectors.

        Parameters
        ----------
        location : N x 3 float array
            Cartesian coordinate location (tails) for N number of vectors.

        vect : N x 3 float array
            N number of vectors in Cartesian coordinates.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Other Parameters
        ----------------
        *args, **kwargs
            All other parameters are passed on to mpl_toolkits.mplot3d.art3d.Line3DCollection.
            Valid keywords include: colors, linewidths.
        
        """
        if alr is None : alr = 0.25
        lines = self._quiverLines(location,vect,alr)
        super().__init__(lines, *args, **kwargs)

    def __str__(self) :
        # note: each vector representation is composed of 3 lines
        numbVect = int(len(self._segments3d)/3)
        name = self.__class__.__name__
        val = name + ':  vectors: {}'.format(numbVect)
        return val

    def _quiverLines( self, vC, vN, alr) :
        # ..................................................
        def quiver( *args, length=1, arrow_length_ratio=alr, pivot='tail', normalize=False):
            # taken directly from mpl_toolkits.mplot3d.axes3d source 
            # and removed axis reference.  Acknowledgement:
            # Copyright (c) 2012-2013 Matplotlib Development Team; All Rights Reserved
            # ...............................................
            def calc_arrow(uvw, angle=15):
                """
                To calculate the arrow head. uvw should be a unit vector.
                We normalize it here:
                """
                # get unit direction vector perpendicular to (u,v,w)
                norm = np.linalg.norm(uvw[:2])
                if norm > 0:
                    x = uvw[1] / norm
                    y = -uvw[0] / norm
                else:
                    x, y = 0, 1

                # compute the two arrowhead direction unit vectors
                ra = np.radians(angle)
                c = np.cos(ra)
                s = np.sin(ra)

                # construct the rotation matrices
                Rpos = np.array([[c+(x**2)*(1-c), x*y*(1-c), y*s],
                                [y*x*(1-c), c+(y**2)*(1-c), -x*s],
                                [-y*s, x*s, c]])
                # opposite rotation negates all the sin terms
                Rneg = Rpos.copy()
                Rneg[[0, 1, 2, 2], [2, 2, 0, 1]] = \
                    -Rneg[[0, 1, 2, 2], [2, 2, 0, 1]]

                # multiply them to get the rotated vector
                return Rpos.dot(uvw), Rneg.dot(uvw)
            # ...............................................
            argi = 6
            if len(args) < argi:
                raise ValueError('Wrong number of arguments. Expected %d got %d' %
                                (argi, len(args)))
            input_args = args[:argi]
            masks = [k.mask for k in input_args
                    if isinstance(k, np.ma.MaskedArray)]
            bcast = np.broadcast_arrays(*input_args, *masks)
            input_args = bcast[:argi]
            masks = bcast[argi:]
            if masks:
                mask = reduce(np.logical_or, masks)
                input_args = [np.ma.array(k, mask=mask).compressed()
                            for k in input_args]
            else:
                input_args = [np.ravel(k) for k in input_args]

            if any(len(v) == 0 for v in input_args):
                return []
            shaft_dt = np.array([0, length])
            arrow_dt = shaft_dt * arrow_length_ratio
            if pivot == 'tail':
                shaft_dt -= length
            elif pivot == 'middle':
                shaft_dt -= length/2.
            elif pivot != 'tip':
                raise ValueError('Invalid pivot argument: ' + str(pivot))
            XYZ = np.column_stack(input_args[:3])
            UVW = np.column_stack(input_args[3:argi]).astype(float)
            norm = np.linalg.norm(UVW, axis=1)
            mask = norm > 0
            XYZ = XYZ[mask]
            if normalize:
                UVW = UVW[mask] / norm[mask].reshape((-1, 1))
            else:
                UVW = UVW[mask]
            if len(XYZ) > 0:
                shafts = (XYZ - np.multiply.outer(shaft_dt, UVW)).swapaxes(0, 1)
                head_dirs = np.array([calc_arrow(d) for d in UVW])
                heads = shafts[:, :1] - np.multiply.outer(arrow_dt, head_dirs)
                heads.shape = (len(arrow_dt), -1, 3)
                heads = heads.swapaxes(0, 1)
                lines = [*shafts, *heads]
            else:
                lines = []
            return lines
        # ..................................................
        x,y,z = np.transpose( vC )
        u,v,w = np.transpose( vN )
        lines = quiver(x,y,z,u,v,w)
        return lines




def _normLength( N, exp, fac, limitarea=None) :
    """
    Default face normal vector length.
    
    Length proportional to the mean triangulated face edge length.
    """
    if limitarea is None : limitarea = 4*np.pi
    # aratio - area ratio as a function of the number of faces, N.
    # heuristic values, exp and fac, are adjustment for small rez surfaces,
    # with aratio limit to 1 for large N. (continuous smooth surface)
    aratio = 1.0 - fac*np.power(N,-exp)
    sc = 0.5*(1.520)  # scale of normalized edge length
    nLen = sc*np.sqrt(aratio*limitarea/N)
    return nLen

def eulerRot(theta, phi, psi=0, useXconv=True, inrad=False) :
    """Rotation matrix from Euler angles, X or Y convention."""
    # X-convention (zxz') and Y-convention (zyz') rotations.
    tht_Rad = theta
    phi_Rad =  phi
    psi_Rad =  psi
    if not inrad :
        tht_Rad = theta*np.pi /180.0
        phi_Rad =  phi*np.pi /180.0
        psi_Rad =  psi*np.pi /180.0

    cT,sT = np.cos(tht_Rad),  np.sin(tht_Rad)
    cP,sP = np.cos(phi_Rad),  np.sin(phi_Rad)
    cS,sS = np.cos(psi_Rad),  np.sin(psi_Rad)
    T_t = [ [ cT, sT, 0 ], [-sT, cT, 0 ], [ 0 , 0 , 1 ]  ]
    T_X = [ [ 1 , 0 , 0 ], [ 0 , cP, sP], [ 0 ,-sP, cP]  ]
    T_Y = [ [ cP, 0 ,-sP], [ 0 , 1 , 0 ], [ sP, 0 , cP]  ]
    T_s = [ [ cS, sS, 0 ], [-sS, cS, 0 ], [ 0 , 0 , 1 ]  ]
    if useXconv : T_p = T_X
    else:         T_p = T_Y
    M = np.matmul( T_p, T_t )
    if psi != 0 : M = np.matmul( T_s, M )
    return M

def setupAxis(axes, length=1.5, width=2, color='black', offset=0.0, negaxis=True) :
    """Set XYZ origin coordinate axis arrows."""
    seglength = length-offset
    endCoor = np.multiply(np.identity(3),seglength)
    strCoor = np.multiply(np.identity(3),offset)
    alr = 0.2/(seglength)
    x,y,z = np.transpose( strCoor )
    u,v,w = np.transpose( endCoor )
    axes.quiver( x,y,z,u,v,w, arrow_length_ratio=alr, color=color, linewidth=width)
    if negaxis :
        endCoor = np.multiply(np.identity(3),-seglength)
        strCoor = np.multiply(np.identity(3),-offset)
        x,y,z = np.transpose( strCoor )
        u,v,w = np.transpose( endCoor )
        axes.quiver( x,y,z,u,v,w, arrow_length_ratio=0, color=color)
    pos = length + 0.15

    axes.text(pos,0,0,'X', color = color, fontweight='bold', fontsize='large',
            horizontalalignment='center', verticalalignment='center')
    axes.text(0,pos,0,'Y', color = color, fontweight='bold', fontsize='large',
            horizontalalignment='center', verticalalignment='center')
    axes.text(0,0,pos,'Z', color = color, fontweight='bold', fontsize='large',
            horizontalalignment='center', verticalalignment='center')
    return      

def standardAxis( axes, **args ) :
    """Set XYZ origin coordinate axis arrows, viewing angle, projection type and clear planes."""
    setupAxis(axes, **args)
    minmax = (-1,1)
    axes.set(xlim=minmax, ylim=minmax, zlim=minmax)
    axes.view_init(30, 30)
    axes.set_proj_type('ortho')
    axes.set_axis_off()
    return

def rtv(direction,elev=None,azim=None) :
    """Transform direction vector 'relative to view' axis viewing angles."""
    if elev is None : elev = 30
    if azim is None : azim = -60
    azim_Rad = azim*np.pi /180.0
    elev_Rad = elev*np.pi /180.0
    cA,sA = np.cos(azim_Rad),  np.sin(azim_Rad)
    cE,sE = np.cos(elev_Rad),  np.sin(elev_Rad)
    TM = [ [ cE*cA, cE*sA, sE ],  [ -sA, cA, 0.0 ],  [  -sE*cA, -sE*sA, cE ] ]
    reldir = np.dot(direction,TM)
    return reldir

def elev_azim_2vector(elev, azim, inrad=False) :
    """Unit vector in the terms of angular position."""
    tht_Rad = azim
    phi_Rad =  elev
    if not inrad :
        tht_Rad = azim*np.pi /180.0
        phi_Rad = elev*np.pi /180.0
    k = np.cos(phi_Rad)
    r = np.sin(phi_Rad)
    i = r*np.cos(tht_Rad)
    j = r*np.sin(tht_Rad)
    return [i,j,k]

