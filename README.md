
# S3Dlib

### Python classes to create 3D surface objects rendered in *Matplotlib*

![S3Dlib logo](s3dliblogo.png)

Detailed documentation and examples are provided at [s3dlib.org](https://s3dlib.org)

---

A 3D surface object is a collection of faces with vertices ordered using
the 'right hand rule' to designate the 'outer' surface normals.  All surface faces
are joined with a minimum of one adjacent face. Adjacent faces share two
common vertices.

The surface object geomentry and color are controlled through various
object methods in the 'Surface3DCollection' base class.
Surface instantiation is performed using the four subclasses that
have predefined surface topologies in native coordinates.
Base class objects may be created by addition of subclass objects to
form a single 'composite' surface.

Objects are added to the mpl_toolkits.mplot3d.Art3d using
the Axis3d.add_collection3d() method.

Included is a module containing functions to create custom Matplotlib color maps.

