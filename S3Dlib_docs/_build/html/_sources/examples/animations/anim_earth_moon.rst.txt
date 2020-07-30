.. _anim_earth_moon:


****************************
Transform and Shade
****************************

.. image:: images/anim_earth_moon.png
   :class: sphx-glr-single-img


Based on the static plot from the :ref:`earth_shaded` example. Animation control:

===========================      ============================================
Visualization                    Frame Value
===========================      ============================================
Surface geometry                 constant 
Surface position                 **rotation** transform per frame
Surface color                    constant
Shading and highlighting         fixed to the coordinate axis
Axis coordinate                  constant
===========================      ============================================


Copies of the geometric and colored surface were used so that those two operations need not
be repeated for each frame, only the rotation and shading.


.. literalinclude:: source/ex_anim_earth_moon.py
   :language: python
   :emphasize-lines: 50

