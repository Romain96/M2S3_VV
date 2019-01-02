#ifndef __VERTEX_HXX__
#define __VERTEX_HXX__

#include <iostream>
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

template<typename T>
bool Vertex<T>::operator== (const Vertex<T> &v)
{
	std::cout << "custom operator : HELLO WORLD" << std::endl;
	return this->_id == v._id;
}

#endif

