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
// vertices are in fact pointers to a list
template<typename T>
void Voxel<T>::addVertex(vec3<T> *v)
{
	_vertices.push_back(v);
}

#endif

