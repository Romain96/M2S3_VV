#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

// nedding point structure for triangle vertices
#include "gmath/vec3.h"

// triangle class
template<typename T>
class Triangle
{
public:
	//////////////////////////////////////////////////////////////////////////////
	// attributes
	//////////////////////////////////////////////////////////////////////////////
	
	// vertices
	vec3<T> _a;
	vec3<T> _b;
	vec3<T> _c;

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
	Triangle(vec3<T> a, vec3<T> b, vec3<T> c);

	// methods for intersection... - TODO
};

#include "Triangle.hxx"

#endif

