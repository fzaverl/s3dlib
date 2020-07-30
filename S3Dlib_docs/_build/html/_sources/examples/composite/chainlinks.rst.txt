.. _chainlinks:


*******************
Surface Addition
*******************

.. image:: images/chainlinks.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_chainlinks.py
   :language: python
   :emphasize-lines: 32,42


Composite surfaces are constructed by simply adding surface objects together to create a
single object, which is then added to the axis.  


.. warning:: 

   Individual surface objects may be added to the axis, 
   however the z-order is **not** preserved among the face polygons of the different
   objects.  The average z-order of face polygons of each object is used as a reference
   to position the individual polygon faces, not referenced to other object polygon faces.
   The effect of adding the individual objects to the axis instead of the composite is
   exemplified in the figure below.


.. image:: images/chainlinks_separate.png
   :class: sphx-glr-single-img
