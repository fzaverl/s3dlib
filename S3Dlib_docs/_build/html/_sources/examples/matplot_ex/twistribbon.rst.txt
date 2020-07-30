.. _twistribbon:

.. role::  raw-html(raw)
    :format: html





************************
Cylindrical Coordinates
************************

Surface Normal Color Mapping
====================================================================================================

Surface geometry is that taken from the Matplotlib 
`More triangular 3D surfaces <https://matplotlib.org/3.1.1/gallery/mplot3d/trisurf3d_2.html#sphx-glr-gallery-mplot3d-trisurf3d-2-py>`_
example.



.. image:: images/twistribbon.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_twistribbon.py
   :language: python
   :emphasize-lines: 37,42


.. note::
   A sliced cylindrical surface, basetype='squ_s', was used for the twist surface because the surface
   twists and joins together at :raw-html:`&theta; = 0 amd 2&pi;`.  Also, for this reason, the mirrored
   color map is applied since the normals at these two locations are in opposite directions.  As noted
   in the code comments, this condition applies for cases when the parameter 'twists' is odd.
   A further description of mirrored colormap usage is given in the  :ref:`cyclic_mirror` example.


Color Mapping in the Z direction
====================================================================================================

The referenced Matplotlib example uses a cmap applied to the vertical location of the surface.  This can also
be easily constructed by uncommenting the highlighted lines to apply a surface cmap.  A simple lambda function
can be used instead of defining a separate functional operation  since the surface geometry is already applied.
The resulting plot is shown below.  When this is done, the shade operation is no longer needed (but can be
applied after the mapping operation to highlight the surface geometry)

.. image:: images/twistribbon_z.png
   :class: sphx-glr-single-img


Parametric Surfaces
====================================================================================================

The *twist* surface was constructed using an additional parameter in the geometry mapping function.
Multiple surfaces, varying with the 'twists' parameter are shown below:

.. image:: images/twist_multi.png
   :class: sphx-glr-single-img

Surface normal color mapping is particularly useful for complex shaped surfaces.
The comprehension of the higher twist surfaces is more difficult using a z color mapping.
Using same twist function, steps 2 and 3 were replaced by the following code to produce the
above figure.


.. literalinclude:: source/ex_twist_multi.py
   :language: python
   :lines: 20-37

