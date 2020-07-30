.. _wireframe:


***************************************************
Functional, Datagrid and Meshgrid wireframe plots
***************************************************

This is a comparison among wireframe plots constructed from either functional, Datagrid or Meshgrid
methods for 3D visualization of data values at xy positions.


Surface geometry is that taken from the matplotlib 
`3D wireframe plot <https://matplotlib.org/3.1.1/gallery/mplot3d/wire3d.html#sphx-glr-gallery-mplot3d-wire3d-py>`_
example.


.. image:: images/wireframe.png
   :class: sphx-glr-single-img



============================================================================================


**Functional**

* Pros

  * Most useful when data can be represented functionally.

  * No need for an additional data arrays.

  * Function in native coordinates for planar, polar, cylindrical or spherical coordinates.

  * Images may be premapped to the surface.

* Cons

  * The starting normalized coordinates require scaling.

  * Only predefined grids may be used, controlled only by rez.

============================================================================================



**Datagrid**

* Pros

  * Data is automatically normalized and only 2D array is needed.

  * Data surface is smoothed at higher rez.

  * May be used for planar, polar, cylindrical or spherical objects.

  * May be placed within in a surface viewport.

  * Images may be premapped to the surface.

* Cons

  * The starting normalized coordinates may require scaling.

  * Only predefined grids may be used, controlled only by rez.

============================================================================================



**Meshgrid**

* Pros

  * Grid resolution is set by the resolution of the data.

  * Grid scaling is automatic.

* Cons

  * Only PlanarSurface objects can be used in xyz coordinates.

  * Image mapping to the surface is not applicable.

  * Surface coloring, shading and highlighting should be executed before and after normalization scaling.

============================================================================================


.. note::
   These plots add the surface property *edges* to the axis3d, not the surface.  For wireframe
   plots with variations in edge color, a surface with transparent faces should be used as
   demonstrated in the :ref:`surface_edges` example.


An additional comparison between datagrids and meshgrids is in the :ref:`data_mesh` example.


.. literalinclude:: source/ex_wireframe.py
   :language: python




