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
	// reserving capacity
	_voxels.reserve(n*n*n);
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
			(*it)._a == goal;
		}

		// second vertex
		if ((*it)._b == target)
		{
			(*it)._b == goal;
		}

		// third vertex
		if ((*it)._c == target)
		{
			(*it)._c == goal;
		}
	}
}

// voxelization
template<typename T>
void Grid3D<T>::voxelize(std::vector<Vertex<T>> &vertices, std::vector<Triangle<T>> &triangles)
{
	// iterating on all vertices
	for (typename std::vector<Vertex<T>>::iterator it = vertices.begin(); it != vertices.end(); it++)
	{
		// finding exact voxel coordinate that hold the vertex
		int x = static_cast<int>((*it)._coord.X() + this->_x/2.0);
		int y = static_cast<int>((*it)._coord.Y() + this->_y/2.0); 
		int z = static_cast<int>((*it)._coord.Z() + this->_z/2.0);
		
		// serializing
		int key = x*_n*_n + y*_n + z;

		// inserting vertex at key position
		_voxels[key].addVertex(*it);

		// making triangles containing this vertex pointing on the newly created one in the voxel
		__pointerReconstruction(triangles, (*it), _voxels[key]._vertices.back());
	}
}

#endif
