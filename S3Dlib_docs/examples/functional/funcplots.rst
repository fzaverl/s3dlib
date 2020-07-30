.. _funcplots:


*****************************
Function Plots, z = f(x,y) 
*****************************

.. image:: images/funcplots.png
   :class: sphx-glr-single-img

These examples show the basic functional plotting using:

*one code statement to create a surface object, one code statement to create the
surface geometry and one code statement to color the surface*. Finally, *one code statement
to add that object to a Matplotlib 3D axis.*

Also, several functions were used to demonstrate:

*a functional relationship should look like a function*

The functions are given on the Wikipedia page
`Test functions <https://en.wikipedia.org/wiki/Test_functions_for_optimization#Test_functions_for_single-objective_optimization>`_
for optimization.


.. literalinclude:: source/ex_funcplots.py
   :language: python
   :emphasize-lines: 58-60,71

