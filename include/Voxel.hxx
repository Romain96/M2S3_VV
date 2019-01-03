#ifndef __VOXEL_HXX__
#define __VOXEL_HXX__

#include <iostream>
#include "Voxel.h"

// constructor
template<typename T>
Voxel<T>::Voxel()
{
	// nothing
}

// adds a vertex to the list of vertices belonging to the vertex
template<typename T>
void Voxel<T>::addVertex(Vertex<T> v)
{
	_vertices.push_back(v);
}

// adds a triangle to the list of triangles intersecting the voxel 
// triangles pointer to triangles in fact
template<typename T>
void Voxel<T>::addTriangle(Triangle<T> *t)
{
	_triangles.push_back(t);
}

#endif

