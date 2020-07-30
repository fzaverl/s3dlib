.. _guides:


***************
Guides
***************

Matplotlib opened the door to an object oriented approach for 3D surface construction by providing the two
Axes3D methods::

    Axes3D.add_collection3d(object)
    Axes3D.collections.remove(object)

where an object can be a Poly3DCollection or Line3DCollection.  S3Dlib classes are inherited from these two,
allowing an object oriented approach having:

* Encapsulation - The internal coordinate structure of surfaces are hidden and generated during
  instantiation.  Surface objects use methods which then define geometry, orientation and color. 
  Object properties, if needed, are accessible through methods.

* Polymorphism - A base class, inherited from the Matplotlib classes, provides the
  methods used by the subclasses.  These subclass objects are then polymorphically added to the axe3D.
  
* Composition - Being derived from the base class, the subclass objects of different classes may be added
  together to form a complex compound surface object with methods available from the base class.

Most of the S3Dlib surface object methods return the object.  This provides concise expressions through the use of
method chaining. Being an object, duplication of objects is simply applied using the copy method.

One of the most significant features of using an object is the creation of animations.  Objects
can be removed, changed and added back from the axis between frames.  The axis view may also be
changed between frames without affecting the object it contains.  In addition, object methods may
be called between frames and then the object can be added back to the axis.


The S3Dlib methodology is
to utilize one of the predefined surface topologies and to apply functional operations
in native coordinates.  Surface coloring and lighting are independently applied
for 3D visualization and/or to indicate additional functional values.

This approach does have its downsize; sacrificing grid optimization for rapid
development.  This is particularly dependent when the resulting surface has regions of high curvature.
However, with increasing surface resolution through simple script parameter assignments, this drawback
can usually be eliminated with moderate increases in execution time.


.. note::

   When using S3Dlib, keep in mind the two concepts:

   * Surface coordinates are normalized.

   * Any function definition must operate on coordinate numpy arrays.

   The numerous examples are provided to demonstrates these.



.. toctree::
   :maxdepth: 2
   
   base_surfaces
   color_maps
   shading
   face_edge
   image_mapping
   vectorfields
   order_of_operations
   orientation

