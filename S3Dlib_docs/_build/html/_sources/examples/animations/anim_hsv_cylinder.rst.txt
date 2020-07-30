.. _anim_hsv_cylinder:


****************************
HSV Mapping
****************************

.. image:: images/anim_hsv_cylinder.png
   :class: sphx-glr-single-img


Animation control:

===========================      ============================================
Visualization                    Frame Value
===========================      ============================================
Surface geometry                 sectioning **parameter** per frame
Surface position                 fixed to the coordinate axis
Surface color                    constant
Shading and highlighting         fixed to the coordinate axis
Axis coordinate                  constant 
===========================      ============================================


This is a composite surface of Planar, Polar and Cylindrical surfaces.
As a result, a single color mapping function is used instead of defining
three separate functions which would be required if treated separately,
one for each surface native coordinate system.
  
.. note::

   Any composite surface uses a xyz coordinate system. 



.. literalinclude:: source/ex_anim_hsv_cylinder.py
   :language: python
   :emphasize-lines: 107

