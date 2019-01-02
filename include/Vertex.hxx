#ifndef __VERTEX_HXX__
#define __VERTEX_HXX__

#include "Vertex.h"

// unique ID initialization
template<typename T> 
unsigned long long Vertex<T>::_idgen = 0;

// default constructor
template<typename T>
Vertex<T>::Vertex()
{
	// unique ID
	_id = _idgen++;
}

// parametrized constructor
template<typename T>
Vertex<T>::Vertex(T x, T y, T z) :
	_coord(vec3<T>(x,y,z))
{
	// unique ID
	_id = _idgen++;
}

#endif

