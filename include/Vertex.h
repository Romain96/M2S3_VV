#ifndef __VERTEX_H__
#define __VERTEX_H__

// needing a 3D vector class
#include "gmath/vec3.h"

// class representing a vertex : a unique ID and  X Y Z coordinates
template<typename T>
class Vertex
{
public:
	//////////////////////////////////////////////////////////////////////////////
	// attributes
	//////////////////////////////////////////////////////////////////////////////
	
	// unique ID cor comparisons
	static unsigned long long _idgen;	// should be enough...(>= 2^64)
	unsigned long long _id;	// unique ID

	// X Y Z coordinates
	vec3<T> _coord;

	//////////////////////////////////////////////////////////////////////////////
	// methods
	//////////////////////////////////////////////////////////////////////////////	

	// constructor
	Vertex();
	Vertex(T x, T y, T z);

	bool operator== (const Vertex<T> &v);
};

#include "Vertex.hxx"

#endif

