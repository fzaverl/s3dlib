.. _face_edge:


************************
Face and Edge Colors
************************

Solid Colors
========================================================================

The colors of surface faces and edges are setup during the object
instantiation using the named argument assignments of *color*, *facecolor*, and *edgecolor*.
With no assignments,
the default color for faces and edges is the Matplotlib color
'C0'.  The assignment of three arguments can be made in several
combinations and, at first, slightly confusing.
From the Matplotlib *Poly3DCollection* code comments::

    '''
    Note that this class does a bit of magic with the _facecolors
    and _edgecolors properties.
    '''

The logic for how colors are assigned is:

1. If color is not assigned, the color will be the 'C0' default color (blue).
2. If facecolor is not assigned, facecolor will be the color.
3. If edgecolor is not assigned, edgecolor will be the facecolor only if the color is **not** assigned, otherwise it will be the color.

Using three different colors for the three arguments, the eight combinations for assignment
to a SphericalSurface are shown below.  The most confusing is when edgecolor is not assigned, in which
case it will take the value of the facecolor or color argument value ( *F* or *CF* ).  In general, either
use just facecolor or the combination of both facecolor and edgecolor to avoid unexpected results.


.. image:: images/CFEcolor.png
   :class: sphx-glr-single-img


When the facecolors are changed after instantiation, either by colormap operations, colormapping normals
or shading, the edgecolors will be reassigned.  This is shown for the shading case where shading has
been applied to the surface. If edges are to be shown, use the object method after shade::

    surface.set_edgecolor(edgecolor)

The script to generate the above plots is given below.


.. literalinclude:: source/gu_CFEcolor.py
   :language: python


Transparency
========================================================================

When the facecolor has an alpha channel value less than one, the edgecolor is not assigned
the alpha and retains the full opacity.  This effect is seen in the following plot for the
first surface on the left.  Even when the edgecolor is assigned transparent and shading applied,
only the facecolor RGB channels are assigned to the edges.  The edge alpha is set to one.
This case is shown for the middle surface.

**To hide the edges** when using facecolor with a transparency, assign the *linewidth* to zero,
as shown for the right surface.

.. image:: images/CFEcolor2.png
   :class: sphx-glr-single-img

The script fragment to generate the above plot is given below.

.. literalinclude:: source/gu_CFEcolor2.py
   :language: python
   :lines: 11-16


**To hide the faces** but still apply shading or a color map to the edges, use the
object method::

    surface.set_facecolor(transparent_color)

where *transparent_color* is any color with the alpha channel set to zero.  


.. image:: images/CFEcolor3.png
   :class: sphx-glr-single-img

.. note::

   As described in the :ref:`helloworld` tutorial, edges can be displayed using
   the object property *edges*.  However, **only** a solid color can be assigned to
   this property with the object constructor.  Color variations are not available.


The script fragment to generate the above plot is given below.

.. literalinclude:: source/gu_CFEcolor2.py
   :language: python
   :lines: 19-23





