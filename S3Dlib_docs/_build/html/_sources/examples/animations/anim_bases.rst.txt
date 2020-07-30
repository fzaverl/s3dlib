.. _anim_bases:


********************************************************
View_init() Azim Reset and Shading
********************************************************

.. image:: images/anim_bases.png
   :class: sphx-glr-single-img

Animation control:

===========================      ============================================
Visualization                    Frame Value
===========================      ============================================
Surface geometry                 constant 
Surface position                 fixed to the coordinate axis
Surface color                    constant
Shading and highlighting         illumination **direction** per frame
Axis coordinate                  **azim** per frame using *view_init* 
===========================      ============================================

Changing the axis coordinate *and* shading/highlighting results in 'stationary viewer' perception.
The view is rotated about multiple surface objects.



.. literalinclude:: source/ex_anim_bases.py
   :language: python
   :emphasize-lines: 73,77

