## Naive Surface Nets in C#

This is an implementation of Naive Surface nets adapted from S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech Report.

I relied on Gibsons paper, [Mikola Lysenkos implementation](https://github.com/mikolalysenko/mikolalysenko.github.com/blob/master/Isosurface/js/surfacenets.js), as well as considerable help from my friend and colleague [Andrew Gotow](https://github.com/andrewgotow).

Please note: This code was written with modularity in mind, separating each task into it's own function. There are a LOT of performance improvements that can be made by not having multiple for loops, however, I wrote this for clarity so that you can better understand the concepts involved.

In this code, I add several child game objects and create a mesh that maintains a certain distance from all of these objects. These objects can have either a box, sphere, or capsule collider. You can also use lines (capsule collider with zero radius) and points (sphere with zero radius). This is just an example but I'm sure you can add your own implementation

We'll start by creating a sample space, then take samples at each voxel in that space based on the sample space size and our resolution,
  
You can take this a step further and create a different distance for each object. This can be very useful in visualizing electric fields, radiation visualization, etc.

[screen1]: https://github.com/tomaszfoster/screenshots/screen1.png "Screenshot 1"

![alt text][screen1]