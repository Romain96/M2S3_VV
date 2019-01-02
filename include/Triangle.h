#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

// needing point structure for triangle vertices
#include "Triangle.h"

// needing a vertex structure
#include "Vertex.h"

// triangle class
template<typename T>
class Triangle
{
public:
	//////////////////////////////////////////////////////////////////////////////
	// attributes
	//////////////////////////////////////////////////////////////////////////////
	
	// vertices
	Vertex<T> _a;
	Vertex<T> _b;
	Vertex<T> _c;

	// lenghts - TODO
	
	T _A;
	T _B;
	T _C;

	// bilinear interpolation coefficients - TODO

	T _Ka;
	T _Kb;
	T _Kc;

	//////////////////////////////////////////////////////////////////////////////
	// methods
	//////////////////////////////////////////////////////////////////////////////

	// constructors
	Triangle();
	Triangle(Vertex<T> a, Vertex<T> b, Vertex<T> c);

	// methods for intersection... - TODO
};

#include "Triangle.hxx"

#endif

