.. _knot:


****************************
Multiple Geometric Maps
****************************

.. image:: images/knot.png
   :class: sphx-glr-single-img

This example demonstrates geometric mapping from one surface to another in sequence.  First,
a cylindrical surface is mapped to a torus using the function, *torusFunc*. This is identical
to the :ref:`cyclic_mirror` example but using a value of 0.2 for the *ratio* parameter.
This surface geometry is then further mapped using the *knot* function to produce the
final geometry.

.. literalinclude:: source/ex_knot.py
   :language: python
   :emphasize-lines: 29,30

