.. _anim_cat2heli:


****************************
Parametric Animation
****************************



.. image:: images/anim_cat2heli.png
   :class: sphx-glr-single-img

Based on the static plot from the :ref:`cat2heli_disp` example. Animation control:

===========================      ============================================
Visualization                    Frame Value
===========================      ============================================
Surface geometry                 functional **parameter** per frame
Surface position                 fixed to the coordinate axis
Surface color                    **color** per frame
Shading and highlighting         fixed to the coordinate axis
Axis coordinate                  constant 
===========================      ============================================






.. literalinclude:: source/ex_anim_cat2heli.py
   :language: python
   :emphasize-lines: 73-75

