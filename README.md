## Naive Surface Nets in C#

[screen1]: https://github.com/TomaszFoster/NaiveSurfaceNets/blob/master/screenshot/screen1.png "Initial Setup"
[screen2]: https://github.com/TomaszFoster/NaiveSurfaceNets/blob/master/screenshot/screen2.png "Running"
[screen3]: https://github.com/TomaszFoster/NaiveSurfaceNets/blob/master/screenshot/screen3.png "Closeup"

This is an implementation of Naive Surface nets adapted from S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech Report.

I relied on Gibsons paper, the incredibly well documented [Mikola Lysenkos implementation](https://github.com/mikolalysenko/mikolalysenko.github.com/blob/master/Isosurface/js/surfacenets.js) and his [explanation](https://0fps.net/2012/07/12/smooth-voxel-terrain-part-2/), as well as considerable help from my friend and colleague [Andrew Gotow](https://github.com/andrewgotow) who also wrote the shader for this visualization.

Please note: This code was written with modularity in mind, separating each task into it's own function. There are a LOT of performance improvements that can be made by not having multiple for loops, however, I wrote this for clarity so that you can better understand the concepts involved.

In this code, I add several child game objects and create a mesh that maintains a certain distance from all of these objects. These objects can have either a box, sphere, or capsule collider. You can also use lines (capsule collider with zero radius) and points (sphere with zero radius). This is just an example but I'm sure you can add your own implementation

We'll start by creating a sample space, then take samples at each voxel in that space based on the sample space size and our resolution,
  
You can take this a step further and create a different distance for each object. This can be very useful in visualizing electric fields, radiation visualization, etc.

In these screens you can see the mesh after it's rendered. Also, I have Shaded Wireframe turned on so that you can better see the triangulation
![alt text][screen1]
![alt text][screen2]

Here's a closeup of the triangle mesh. Nice and smooth tessellation.
![alt text][screen3]