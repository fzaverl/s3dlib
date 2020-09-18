## Tranformation function studies

This directory contains Python scripts and companion
figures using the [S3Dlib](https://s3dlib.org) package.
The studies are used to note limitations and functional improvements to the package.

---

![Knot](images/an_escher1.png)

> ### an_escher1.py  
>
> Inspired by M.C.Escher - [Knots](https://mcescher.com/gallery/mathematical/)  
> The knot function is based on a trefoil knot, described in [Wikipedia](https://en.wikipedia.org/wiki/Trefoil_knot)  
>
> Since the ring cross-section is square, twists can be multiples of 0.25.
> The right-handed and left-handed varients are controlled by the sign for Z in the Trefoil function.
> The clock and counter-clockwise rotation of the twist in the knot is controlled by the sign
> for t0 in the twistFunction function.  
>
> The following default values were set to match the Escher figure.
> - number of twists, 0.75
> - twist offset, 0.25
> - numerical signs used in the trefoil and twist functions.
> - colors and illumination direction.

---

![Knot2](images/escher2.png)

> ### escher2.py  
>
> Inspired by M.C.Escher - [Knots](https://mcescher.com/gallery/mathematical/)  
> The knot function is based on a trefoil knot, described in [Wikipedia](https://en.wikipedia.org/wiki/Trefoil_knot)  
>
> This example differs from the previous case in that, for this surface, both the outer and inner surfaces are exposed.
> As a result, illumination rendering is not as consistent. Also, twists needed to be multiples of 0.5 for continuous
> surface coloring.  
>
> For this case, a compound surface was constructed. The two initial surface objects were geometrically
> mapped prior to being added together. In the current version, geometric mapping is not permitted for
> compound surfaces. However, it would have been effective to combine the objects prior to mapping
> since both sub-surfaces were of the same type of sub-class object.  
>
> Dev notes for future versions:
> - mirrored shading and highlighting (effectively having two opposing light sources)
> - allow geometric transforms for compound surfaces, iff sub-surfaces are of the same derived class.
>   In such cases, the native coordinates are known and transforms can be made using the sub-class coor transform. 

---

![Knot3](images/escher3.png)

> ### escher3.py  
>
> Inspired by M.C.Escher - [Knots](https://mcescher.com/gallery/mathematical/)  
> Third in the set using a torus geometry.  
>
> As in the previous case, a compound surface was constructed. Again, the same restriction on
> compound surfaces requiring the work around. Also for surface with alpha set, the edges are shown
> if the linewidth is not set to zero, hence the need for this assignment.  
>
> Dev notes for future versions:
> - evaluate the issue of requiring linewidth to be set to zero for faces with non-unity alphas. 

---


![Knot2](images/an_escher4.png)

> ### an_escher4.py  
>
> Exploration of creating visualization of movement without changes in orientation or geometry.  
>
> The s3dlib.cmap_utilites were used to create custom stiched color maps with gradients in the alpha
> channels. Using the geometries in the previous case, the maps were simply applied
> to varying angular positions of the geometry.  
>
> Demonstration of using a:
> - color map with changes in alpha.
> - compound color map with regions having alpha set to zero.
> - color map applied to various coordinates of the surface geometry.  
>
> Dev notes for future versions:
> - display a Matplotlib colorbar which is semi-transparent or transparent. 

---


![power](images/an_power.png)

> ### an_power.py  
>
> Exploration of an alternative visualization for three line plots over a periodic interval.
> The following 2D plots of the data used in the above plot are shown below. 
>
> ![powerplot](images/power_plot.png)
>
> This data was obtained by modifying the
> [powerdata](https://picower.mit.edu/news/information-flows-through-brains-heirarchy-higher-regions-use-higher-frequency-waves) 
> to create a period over the frames range with a value range of 0 to 1.
> >
> Notes:
> - Only problem during code development was noting that the frame argument
>   passed into the update_fig() function was of type float.  This argument value
>   is an index array and needed to be passed as an int to the data3D() method.  
> - This 3D method of visualization does not provide any new conceptional information 
>   for this particular data set.  However, this visulization method may be useful for other
>   types of data sets, particularly when the frames represent a cyclic time sequence.
> - The alternative 3D object can be obtained by xyz scaling ( ellipsoid ) or modidying the
>   expression for F in the data3D method ( power of F results in sharpening or smoothing ).
> - The simplest coloring of the surface was obtained by directly using RGB space.
> - More complex coloring could define a function to pass as an argument to map_geom_from_op().  

---









