#ifndef __VOXEL_H__
#define __VOXEL_H__

// requires a point structure to store the triangle's vertices
#include "Vertex.h"

// requires a triangle structure to store the triangles intersecting the voxel
#include "Triangle.h"

#include <vector>

// class representing a voxel
template<typename T>
class Voxel
{
public:
	//////////////////////////////////////////////////////////////////////////////
	// attributes
	//////////////////////////////////////////////////////////////////////////////
	
	std::vector<Vertex<T>> _vertices;	// list of vertices inside the voxel
	std::vector<Triangle<T>*> _triangles;	// list of trianges intersecting the voxel

	//////////////////////////////////////////////////////////////////////////////
	// methods
	//////////////////////////////////////////////////////////////////////////////
	
	// constructor (empty by default)
	Voxel();

	// adding a pointer to a vertex into the voxel
	void addVertex(Vertex<T> v);

	// adding a triangle to the list of triangles intersecting the voxel
	void addTriangle(Triangle<T> *t);
};

#include "Voxel.hxx"

#endif

