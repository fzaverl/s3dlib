.. _screw:

.. role::  raw-html(raw)
    :format: html


**************************************
Sliced Polar Surface
**************************************

.. image:: images/screw.png
   :class: sphx-glr-single-img


The function in this example is not cyclic with :raw-html:`&theta;`, e.g. 
:raw-html:`f( &theta;=0 ) &ne; f( &theta;=2&pi; )`.
Therefore, a PolarSurface object was used with a basetype *hex_s* which
is not continuous at 0 and :raw-html:`2&pi;`.
Evaluation was made in the domain of 
:raw-html:`-3&pi; &le; &theta; &le; 3&pi;`,
as seen in the highlighted lines with k=3.


.. literalinclude:: source/ex_screw.py
   :language: python
   :emphasize-lines: 11,22

