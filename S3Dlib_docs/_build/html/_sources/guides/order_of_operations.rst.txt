.. _order_of_operations:


*************************
OOPs, Order of Operations
*************************

In general, three types of operations on a surface object affect the final visualization.

* geometric
* color
* transform

These may be applied in any sequence but the result is dependent on the sequence.  There is no
sequence which is correct but dependent on the intended result.
In the :ref:`example` section, various sequences were used to establish the final surface.

Color, Geometry, Rotation
-------------------------------------

The following figure illustrates the different permutations of sequence for identical operations.

.. image:: images/order_ops.png
   :class: sphx-glr-single-img

For these operations, if a 'donut' geometric shape is the final desired result, geometric
transforms should preceed any transforms.  This includes not only rotations, but scaling
and translations.
This is probably true in most cases but not a given rule.

In each of these figures, shading was applied to the final object.  Shading, in a sense, is
an operation that is both controlling color but dependent on orientation relative to the
relative illumination direction. In the :ref:`anim_bases` animation example, this relation between
orientation and view was considered.


Scale, Rotation
-------------------------------------

The order of scaling and rotating a surface affects the final configuration of the surface.
In the following example, a cube is transformed by scaling and rotating in different sequences
by separate transform method operations.
The lower set demonstrates the effect of transform order.  Included in the surface plots are
shown the transformed x,y,z axes colored red, green, and blue respectively.  The transformed
axes reflect the coordinates of the last transformation. As a result, the lower right axes are aligned
with the original axes and the lower left axes all have the same length.

.. image:: images/order_trans.png
   :class: sphx-glr-single-img

When rotations and scaling are performed using one transform method call, the order of operation
is to scale, then rotate.
A single 'rotational' transformation matrix can be computed by multiplying the rotation matrix 
by the scaling array.  The resulting rotational matrix will not be orthogonal but may be
used as a *rotate* parameter in the transform method without scaling.  As seen in the following
plots, the transform effect is the same as using the order of rotations then scaling.  Also
note in the right plot, the axes have been skewed since the rotation is not orthogonal.


.. image:: images/order_trans2.png
   :class: sphx-glr-single-img





Python Script
-------------------------------------

The script used to create plots for comparing the order of operation for color, geometry
and rotation is given below

.. literalinclude:: source/gu_order_ops.py
   :language: python

The script used to create plots for comparing the order of operation for rotation and scaling
is given below

.. literalinclude:: source/gu_order_trans.py
   :language: python




