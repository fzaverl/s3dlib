.. s3dlibdoc documentation master file

Introduction
========================================


S3Dlib is a Matplotlib-based API consisting of a set of 3D plotting classes
and functions to be used in conjunction with the
`Matplotlib <https://matplotlib.org>`_
plotting library.

.. raw:: html

    <div style='text-align:center;'><div style='display:inline-block;'>

.. image:: images/figure_1.png

.. image:: images/figure_2.png

.. image:: images/figure_3.png

.. image:: images/figure_4.png

.. image:: images/figure_5.png

.. raw:: html

    </div></div>


Now, S3Dlib makes hard things easy.

Although Matplotlib provides 3D capability, ease of use is limited to simple projections
to the z-coordinate.  Complex 3D surfaces require considerable effort to construct grids for 
such surfaces along with functional mapping.

S3Dlib takes an object oriented approach to 3D surface construction.
With predefined three-dimensional surface grids in native coordinate systems, S3Dlib makes
plotting easy for more advanced functional relations in 3D.  S3Dlib classes are directly derived from
Matplotlib classes. As a result, all the object properties provided by that library are available when using S3Dlib.
Once created, S3Dlib surface objects are directly added to the Matplotlib Axis3D.


Documentation
========================================

This is the documentation for S3Dlib version 1.0.

The basic procedure for quickly seeing how S3Dlib can generate plots is 
demonstrated in the  :ref:`helloworld` tutorial.  
Get started by going over the :ref:`tutorial` with more detailed explanation in
the :ref:`guides` pages.

Check out the :ref:`example` gallery to see various S3Dlib applications.
Numerous examples also provide explanations on usage and various features. 

Developer's Notes
========================================

Unlike numerous 3D packages, S3Dlib is not a substitute for, but a complement to Matplotlib
for rendering 3D surfaces.

The main objective of developing S3Dlib was that, given a function in native coordinates,
it should only take one code statement to create a surface object, one
code statement to create the surface geometry, and one code statement to color the surface.
Finally, one code statement to add that object to a Matplotlib 3D axis.
Any further annotation would be achieved using Matplotlib.  

A secondary objective was that the function used to define the surface geometry
is easily comprehensible by examining the code, and not obscured in the grid
creation algorithm.  In other words, a functional relationship should look like a function
in native coordinates.
 
As shown in the various examples, these objectives have been reasonably achieved. 


Any comments regarding computational errors or code improvements, and additional capabilities
would be encouraged and appreciated.  In particular, there are numerous ways of achieving
the same thing using Numpy arrays, so optimization suggestions would be very helpful.
Comments based on best practices would also be helpful, including documentation improvements
(please provide references if available).
Any examples using S3Dlib would also be welcome.  Click
`here <mailto:comments@s3dlib.org?subject=Comments-S3Dlib>`_ 
to send an email for your comments.

This documentation was built with a similar look and feel of the Matplotlib documentation
for several reasons.  First, familiarity with Matplotlib documentation makes this documentation
easy to use for the intended audience.  Second, it's simply just a good design.  Finally, using the expertise of other
developers definitely helps getting the documentation up and running in the shortest amount of time.

Hope this package is helpful.

.. toctree::
   :hidden:
   
   inst_index.rst
   doc_index.rst


