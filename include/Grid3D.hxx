#ifndef __GRID3D_HXX__
#define __GRID3D_HXX__

#include <iostream>
#include "Grid3D.h"

// parametrized constructor	
template<typename T>
Grid3D<T>::Grid3D(T x, T y, T z, unsigned int n) :
	_x(x),
	_y(y),
	_z(z),
	_n(n)

{
	// debug
	std::cout << "[Grid3D::Grid3D] : 3D grid initialised with parameters x=" << x << " y=" << y << " z=" << z << " n=" << n << std::endl;
}

// iterates on all triangles containes in "triangles"
// for each of the three vertices if a vertex
// is the "target" then pointing to "goal"
template<typename T>
void Grid3D<T>::__pointerReconstruction(std::vector<Triangle<T>> &triangles, Vertex<T> target, Vertex<T> goal)
{
	
	for (typename std::vector<Triangle<T>>::iterator it = triangles.begin(); it != triangles.end(); it++)
	{
		// comparison uses the custom comparator
		// first vertex
		if ((*it)._a == target)
		{
			(*it)->_a == goal;
		}

		// second vertex
		if ((*it)._b == target)
		{
			(*it)->_b == goal;
		}

		// third vertex
		if ((*it)._c == target)
		{
			(*it)->_c == goal;
		}
	}
}

#endif
