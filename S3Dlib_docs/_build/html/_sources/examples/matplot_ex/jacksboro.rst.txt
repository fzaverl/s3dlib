.. _jacksboro:


*******************
Datagrid Geometry
*******************

The data used to construct the surface geometry is that used for the Matplotlib 
`Custom hillshading in a 3D surface plot <https://matplotlib.org/3.1.1/gallery/mplot3d/custom_shaded_3d_surface.html#sphx-glr-gallery-mplot3d-custom-shaded-3d-surface-py>`_
example.


.. image:: images/jacksboro.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_jacksboro.py
   :language: python
   :emphasize-lines: 13


This example uses the datagrid to construct the normalized planar surface. A color map is
applied to the face normals to visualize the 3D surface.  Note that the array order of the originating
datagrid is reversed, shown in the highlighted line, to be consistent with the datagrid coordinate
orientation definition.  This is further discussed in the :ref:`image_mapping` guide.

.. _vert_color:

Vertical Colorization
-------------------------

The referenced matplotlib example uses a cmap applied to the vertical location of the surface.  This can also
be easily constructed since the surface geometry is already applied.  A simple lambda function
can be used, instead of defining a separate functional operation.
The resulting plot is shown below.


.. image:: images/jacksboro_z.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_jacksboro_z.py
   :language: python
   :emphasize-lines: 30,31

Shading and highlighting is applied to the surface using the light direction shown in the highlighted line.
This provided a more direct comparison to the Matplotlib example which uses hillshading.  Surface
scaling using the dataframe arrays also resized and positioned the surface for displaying data units.
 
