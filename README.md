yangtiao
================
The package can be used to denoise image by block-B-spline.


Denoise Image:
------------------------------
when  get a noisy image:

![](README_files/2.jpeg)

Firstly, search jump points,and interpolate between the points:

![](README_files/3.jpeg)

Then, depend on the jump points, segment the image to make each block a continuous image:

![](README_files/4.jpeg)

Finally, The B-spline is used to fit the surface in each block:

![](README_files/5.jpeg)

The following image is the original image:

![](README_files/1.jpeg)


The R package also can be use to calculate the value of B-spline basis，and provide 2D(`sin(x)` and `circular`) and 3D(`two dimensional density surface` and `sphere`) examples in the `demo` .

Some examples' result diagram:
------------------------------

![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-1-1.png)
![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-1-2.png)
![](README_files/density.png)
![](README_files/sphere.png)
