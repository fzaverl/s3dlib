.. _anim_retinal_scan:


****************************
Figure View_init() Reset
****************************

.. image:: images/anim_retinal_scan.png
   :class: sphx-glr-single-img

Based on the static plot from the :ref:`retinal_scan` example. Animation control:

===========================      ============================================
Visualization                    Frame Value
===========================      ============================================
Surface geometry                 constant 
Surface position                 fixed to the coordinate axis
Surface color                    constant
Shading and highlighting         fixed to the coordinate axis
Axis coordinate                  **azim** per frame using *view_init* 
===========================      ============================================

Changing only the axis coordinate results in 'stationary object' perception.

.. literalinclude:: source/ex_anim_retinal_scan.py
   :language: python
   :emphasize-lines: 59

