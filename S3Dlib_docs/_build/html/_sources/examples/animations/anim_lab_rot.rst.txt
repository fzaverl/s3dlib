.. _anim_lab_rot:


********************************************************
Lab Mapping
********************************************************

.. image:: images/anim_lab_rot.png
   :class: sphx-glr-single-img

Based on the static plot from the :ref:`Lab_space` example. Animation control:


===========================      ============================================
Visualization                    Frame Value
===========================      ============================================
Surface geometry                 constant 
Surface position                 fixed to the coordinate axis
Surface color                    constant
Shading and highlighting         illumination **direction** per frame
Axis coordinate                  **azim** per frame using *view_init* 
===========================      ============================================

To minimize execution time, the Lab surface object was constructed once,
then copied for each frame construction.
Changing the axis coordinate *and* shading/highlighting results in 'stationary viewer' perception.

.. literalinclude:: source/ex_anim_lab_rot.py
   :language: python
   :emphasize-lines: 142,146

