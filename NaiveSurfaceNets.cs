// The MIT License (MIT)
//
// Copyright (c) 2017-2018 Tomasz Foster
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

/**
 * SurfaceNets in C# for Unity3D
 *
 * Written by Tomasz Foster (C) 2017
 * https://www.github.com/tomaszfoster
 *
 * MIT License
 *
 * Based on: S.F. Gibson, "Constrained Elastic Surface Nets". (1998) MERL Tech Report.
 * Based on: Mikola Lysenko, "SurfaceNets in JavaScript" 
 * https://github.com/mikolalysenko/mikolalysenko.github.com/blob/master/Isosurface/js/surfacenets.js
 *
 * Please note: This code was written with modularity in mind, separating each task
 * into it's own function. There are a LOT of performance improvements that can be made
 * by not having multiple for loops, however, I wrote this for clarity so that you can
 * better understand the concepts involved.
 * 
 * In this code, I add several child game objects and create a mesh that maintains a certain distance
 * from all of these objects. These objects can have either a box, sphere, or capsule collider. You can also
 * use lines (capsule collider with zero radius) and points (sphere with zero radius). This is just an example
 * but I'm sure you can add your own implementation
 * 
 * We'll start by creating a sample space, then take samples at each voxel in that space
 * based on the sample space size and our resolution,
 *   
 * You can take this a step further and create a different distance for each object.
 * This can be very useful in visualizing electric fields, radiation visualization, etc.
 *
 **/

// some prerequisites to set up our Unity environment
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;
using System;

public class GenerateSurfaceNets : MonoBehaviour {

	// how farm from the surface we want to create our mesh
	public float threshold = 0.25f;
	// the resolution
	public float samplesPerMeter = 1.0f;
	// we set up a sample space. we don't want to sample a huge area to create our mesh
	public Bounds sampleSpace;
	public GameObject[] fieldObjects;
	public Material material;

	// Internal buffer
	private int[] buffer;

	public Vector3 sampleResolution {
		get {
			var samples = sampleSpace.size * samplesPerMeter;
			return new Vector3( Mathf.RoundToInt(samples.x), Mathf.RoundToInt(samples.y), Mathf.RoundToInt(samples.z) );
		}
	}

	/*
	 * Give each voxel a position, an index, a corner mask, and an edge mask.
	 * isOnSurface basically asks if our mesh goes through this voxel
	 */
	public struct Voxel {
		public int vertexIndex;
		public Vector3 vertexPosition;
		public bool isOnSurface;
		public int voxelEdgeMask;
		public int cornerMask;
	}

	public struct sample {
		public Vector3 point;
		public float value;
	}

	// create an array of voxels that'll fill our sample space
	private Voxel[,,] voxels;
	private List<Vector3> vertices;
	private List<Vector3> newVertices;
	private List<int> triangles;

	// Define the edge table and cube edges with a set size
	private int[] edgeTable = new int[256];
	private int[] cubeEdges = new int[24];

	/* 
	 * Set up the 8 corners of each voxel like this:
	 * Very important that this is set up the same way as the cube edges
	 * 
	 *  y         z
	 *  ^        /     
	 *  |
	 *    6----7
	 *   /|   /|
	 *  4----5 |
	 *  | 2--|-3
	 *  |/   |/
	 *  0----1   --> x
	 * 
	 */

	private Vector3[] voxelCornerOffsets = new Vector3[8] {
		new Vector3(0,0,0),		// 0
		new Vector3(1,0,0), 	// 1
		new Vector3(0,1,0), 	// 2
		new Vector3(1,1,0), 	// 3
		new Vector3(0,0,1), 	// 4
		new Vector3(1,0,1), 	// 5
		new Vector3(0,1,1), 	// 6
		new Vector3(1,1,1)  	// 7
	};

	/*
	 * convert coordinates from the voxel index to world space position
	 * i.e. convert voxel[u,v,w] to a Vector3 with an x,y,z in our world space coordinate plane 
	 */
	private Vector3 GetWorldSpaceSamplePosition (int x, int y, int z) {
		return sampleSpace.min + new Vector3 (
			x * sampleSpace.size.x / sampleResolution.x,
			y * sampleSpace.size.y / sampleResolution.y,
			z * sampleSpace.size.z / sampleResolution.z);
	}

	/*
	 * Build an intersection table. This is a 2^(cube config) -> 2^(edge config) map
	 * There is only one entry for each possible cube configuration
	 * and the output is a 12-bit vector enumerating all edges
	 * crossing the 0-level
	 */

	void GenerateIntersectionTable() {

		for(int i=0; i<256; ++i) {
			int em = 0;
			for(int j=0; j<24; j+=2) {
				var a = Convert.ToBoolean(i & (1 << cubeEdges [j]));
				var b = Convert.ToBoolean(i & (1 << cubeEdges[j+1]));
				em |= a != b ? (1 << (j >> 1)) : 0;
			}
			edgeTable[i] = em;
		}
	}

	/* 
	 * Utility function to build a table of possible edges for a cube with each
	 * pair of points representing one edge i.e. [0,1,0,2,0,4,...] would be the 
	 * edges from points 0 to 1, 0 to 2, and 0 to 4 respectively:
	 * 
	 *  y         z
	 *  ^        /     
	 *  |
	 *    6----7
	 *   /|   /|
	 *  4----5 |
	 *  | 2--|-3
	 *  |/   |/
	 *  0----1   --> x
	 * 
	 */

	void GenerateCubeEdgesTable () {
		int k = 0;
		for( int i=0; i<8; ++i) {
			for( int j=1; j<=4; j<<=1) {
				int p = i^j;
				if(i <= p) {
					cubeEdges[k++] = i;
					cubeEdges[k++] = p;
				}
			}
		}
	}

	/*
	 * Not necessary but let's draw a wire cube mesh around our sample space just
	 * to help us visualize it
	 */

	void OnDrawGizmos() {
		Gizmos.color = Color.blue;
		Gizmos.DrawWireCube (sampleSpace.center, sampleSpace.size);
	}

	// Set up the size of each list of voxels, verticies, triangles
	void Awake() {
		voxels = new Voxel[(int)sampleResolution.x, (int)sampleResolution.y, (int)sampleResolution.z];
		vertices = new List<Vector3> ();
		newVertices = new List<Vector3> ();
		triangles = new List<int> ();
	}

	/** 
	 * this doens't have to implement at run time, you can make 
	 * a custom editor for Unity but we won't get into that here
	 **/
	void Start () {
		GenerateCubeEdgesTable ();
		GenerateIntersectionTable ();
		CalculateVertexPositions ();
		GenerateFieldMesh ();
	}

	/*
	 * SampleField: Takes in a vector3 from a voxel and looks at every child object
	 * and returns the smallest distance from each possible object.
	 */
	 
	float SampleField (Vector3 position) {

		// set the default final result to infinity
		float finalResult = Mathf.Infinity;

		for (int i = 0; i < fieldObjects.Length; i++) {

			// get the type of collider on the object
			Collider collider = fieldObjects[i].GetComponent<Collider> ();

			if (collider is SphereCollider) {
				SphereCollider sphereCollider = (SphereCollider) collider;
				finalResult = Mathf.Min (
					finalResult, 
					SphericalField (
						sphereCollider.center, 
						sphereCollider.radius, 
						position, 
						fieldObjects [i].transform 
					)
				);
			} else if (collider is CapsuleCollider) {
				CapsuleCollider capsuleCollider = (CapsuleCollider) collider;
				finalResult = Mathf.Min (
					finalResult, 
					CapsuleField (
						capsuleCollider.center, 
						capsuleCollider.radius, 
						position, 
						fieldObjects [i].transform, 
						capsuleCollider.height
					)
				);
			} else if (collider is BoxCollider) {
				BoxCollider boxCollider = (BoxCollider)collider;
				finalResult = Mathf.Min (
					finalResult,
					BoxField(
						position,
						boxCollider.center,
						boxCollider.size,
						fieldObjects[i].transform
					)
				);
			}
		}
		return finalResult;

	}

	// return value for a box field
	float BoxField( Vector3 samplePosition, Vector3 centerOfBox, Vector3 boxSize, Transform boxTransform ) {
		
		Vector3 s = boxTransform.InverseTransformPoint (samplePosition);
		s -= centerOfBox;

		Vector3 d = new Vector3( Mathf.Abs(s.x), Mathf.Abs(s.y), Mathf.Abs(s.z) ) - (boxSize * 0.5f);
		Vector3 maxD = new Vector3 (Mathf.Max (d.x, 0.0f), Mathf.Max (d.y, 0.0f), Mathf.Max (d.z, 0.0f));
		return Mathf.Min ( Mathf.Max( d.x, Mathf.Max(d.y, d.z) ), 0.0f ) + Vector3.SqrMagnitude( maxD );
	}

	// Return value for a spherical field
	float SphericalField( Vector3 centerOfSphere, float sphereRadius, Vector3 samplePosition, Transform sphereTransform ) {
		Vector3 s = sphereTransform.InverseTransformPoint (samplePosition);
		return (s - centerOfSphere).magnitude - sphereRadius;
	}

	// Return value for a capsule field
	float CapsuleField( Vector3 centerOfCapsule, float capsuleRadius, Vector3 samplePosition, Transform capsuleTransform, float capsuleHeight ) {
		Vector3 s = capsuleTransform.InverseTransformPoint (samplePosition);
		Vector3 p1 = centerOfCapsule + Vector3.up * (capsuleHeight * 0.5f - capsuleRadius);
		Vector3 p2 = centerOfCapsule - Vector3.up * (capsuleHeight * 0.5f - capsuleRadius);
		float t = Vector3.Dot ((s - p1), (p2 - p1).normalized);
		t = Mathf.Clamp01 (t);
		Vector3 closestPoint = p1 + (p2 - p1).normalized * t;
		return (s - closestPoint).magnitude - capsuleRadius;
	}


	void GenerateFieldMesh() {

		Mesh mesh = new Mesh ();
		mesh.name = "field mesh";

		gameObject.AddComponent<MeshFilter>();
		gameObject.AddComponent<MeshRenderer>();
		ComputeMesh ();
		mesh.Clear ();
		mesh.SetVertices (newVertices);
		mesh.SetTriangles (triangles, 0);
		mesh.RecalculateNormals ();

		GetComponent<MeshFilter> ().mesh = mesh;
		gameObject.GetComponent<Renderer>().material = material;

	}

	/** 
	 * calculate the postion of each vertex. this sets up our 3 dimensional grid of voxels
	 * while also sampling each voxel
	 **/

	void CalculateVertexPositions () {
		for (int x = 0; x < voxels.GetLength(0)-1; x++ ) {
			for (int y = 0; y < voxels.GetLength(1)-1; y++ ) {
				for (int z = 0; z < voxels.GetLength(2)-1; z++ ) {

					// default values.
					voxels [x, y, z].isOnSurface = false;
					voxels [x, y, z].voxelEdgeMask = 0;
					voxels [x, y, z].vertexIndex = -1;
					voxels [x, y, z].cornerMask = 0;

					int cornerMask = 0;

					// sample the 8 corners for the voxel and create a corner mask
					sample[] samples = new sample[8];
					for (int i = 0; i < 8; i++) {
						var offset = voxelCornerOffsets [i];
						var pos = GetWorldSpaceSamplePosition (x + (int)offset.x, y + (int)offset.y, z + (int)offset.z);
						var sample = SampleField (pos);
						samples [i].value = sample;
						samples [i].point = pos;
						cornerMask |= ((sample > threshold) ? (1 << i) : 0);
					}

					//Check for early termination if cell does not intersect boundary
					if(cornerMask == 0 || cornerMask == 0xff) {
						continue;
					}

					// get edgemask from table using our corner mask
					int edgeMask = edgeTable [cornerMask];
					int edgeCrossings = 0;
					var vertPos = Vector3.zero;

					for(int i=0; i < 12; ++i) {

						//Use edge mask to check if it is crossed
						if( !( (edgeMask & (1<<i)) > 0 ) ) {
							continue;
						}

						//If it did, increment number of edge crossings
						++edgeCrossings;

						//Now find the point of intersection
						int e0 = cubeEdges[ i<<1 ];
						int e1 = cubeEdges[(i<<1)+1];
						float g0 = samples[e0].value;
						float g1 = samples[e1].value;
						float t  = (threshold - g0 ) / (g1 - g0);

						vertPos += Vector3.Lerp (samples [e0].point, samples [e1].point, t);
					}
					vertPos /= edgeCrossings;

					voxels [x, y, z].vertexPosition = vertPos;
					voxels [x, y, z].isOnSurface = true;
					voxels [x, y, z].voxelEdgeMask = edgeMask;
					voxels [x, y, z].vertexIndex = vertices.Count;
					voxels [x, y, z].cornerMask = cornerMask;
					vertices.Add (vertPos);
				}
			}
		}
	}


	public void ComputeMesh()
	{
		// set the size of our buffer
		buffer = new int[4096];

		// get the width, height, and depth of the sample space for our nested for loops
		int width = voxels.GetUpperBound(0) + 1;
		int height = voxels.GetUpperBound(1) + 1;
		int depth = voxels.GetUpperBound(2) + 1;

		int n = 0;
		int[] pos = new int[3];
		int[] R = new int[]{1, width + 1, (width + 1) * (height + 1)};
		float[] grid = new float[8];
		int bufferNumber = 1;

		// resize the buffer if it's not big enough
		if (R[2] * 2 > buffer.Length)
			buffer = new int[R[2] * 2];

		for (pos[2] = 0; pos[2] < depth - 1; pos[2]++, n += width, bufferNumber ^= 1, R[2] = -R[2] )
		{
			var bufferIndex = 1 + (width + 1) * (1 + bufferNumber * (height + 1));

			for(pos[1] = 0; pos[1] < height - 1; pos[1]++, n++, bufferIndex += 2)
			{
				for (pos[0] = 0; pos[0] < width - 1; pos[0]++, n++, bufferIndex++)
				{
					// get the corner mask we calculated earlier
					var mask = voxels[pos[0], pos[1], pos[2]].cornerMask;

					// Early Termination Check
					if (mask == 0 || mask == 0xff)
					{
						continue;
					}

					// get edge mask
					var edgeMask = edgeTable[mask];

					var vertex = new Vector3();
					var edgeIndex = 0;

					//For Every Cube Edge
					for (var i = 0; i < 12; i++)
					{
						//Use Edge Mask to Check if Crossed
						if (!Convert.ToBoolean(edgeMask & (1 << i)))
						{
							continue;
						}

						//If So, Increment Edge Crossing #
						edgeIndex++;

						//Find Intersection Point
						var e0 = cubeEdges[i << 1];
						var e1 = cubeEdges[(i << 1) + 1];
						var g0 = grid[e0];
						var g1 = grid[e1];
						var t = g0 - g1;
						if (Math.Abs(t) > 1e-16)
							t = g0 / t;
						else
							continue;

						//Interpolate Vertices, Add Intersections
						for(int j = 0, k = 1; j < 3; j++, k <<=1)
						{
							var a = e0 & k;
							var b = e1 & k;
							if (a != b)
								vertex[j] += Convert.ToBoolean(a) ? 1f - t : t;
							else
								vertex[j] += Convert.ToBoolean(a) ? 1f : 0;
						}
					}

					//Average Edge Intersections, Add to Coordinate
					var s = 1f / edgeIndex;
					for(var i = 0; i < 3; i++)
					{
						vertex[i] = pos[i] + s * vertex[i];
					}
					vertex = voxels[pos[0], pos[1], pos[2]].vertexPosition;

					//Add Vertex to Buffer, Store Pointer to Vertex Index
					buffer[bufferIndex] = newVertices.Count;
					newVertices.Add(vertex);

					//Add Faces (Loop Over 3 Base Components)
					for (var i = 0; i < 3; i ++)
					{
						//First 3 Entries Indicate Crossings on Edge
						if(!Convert.ToBoolean(edgeMask & (1 << i)))
						{
							continue;
						}

						//i - Axes, iu, iv - Ortho Axes
						var iu = (i + 1) % 3;
						var iv = (i + 2) % 3;

						//Skip if on Boundary
						if (pos[iu] == 0 || pos[iv] == 0)
							continue;

						//Otherwise, Look Up Adjacent Edges in Buffer
						var du = R[iu];
						var dv = R[iv];

						//Flip Orientation Depending on Corner Sign
						if (Convert.ToBoolean(mask & 1))
						{
							triangles.Add(buffer[bufferIndex]);
							triangles.Add(buffer[bufferIndex - du - dv]);
							triangles.Add(buffer[bufferIndex - du]);
							triangles.Add(buffer[bufferIndex]);
							triangles.Add(buffer[bufferIndex - dv]);
							triangles.Add(buffer[bufferIndex - du - dv]);
						}
						else
						{
							triangles.Add(buffer[bufferIndex]);
							triangles.Add(buffer[bufferIndex - du - dv]);
							triangles.Add(buffer[bufferIndex - dv]);
							triangles.Add(buffer[bufferIndex]);
							triangles.Add(buffer[bufferIndex - du]);
							triangles.Add(buffer[bufferIndex - du - dv]);
						}
					}

				}
			}
		}

	}

}