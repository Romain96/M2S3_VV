#ifndef __GRID3D_H__
#define __GRID3D_H__

// needing point structure for vertices
#include "gmath/vec3.h"

// needing triangle structure
#include "Triangle.h"

#include <vector>

// interface for 3D grid based space partitionning class
template<typename T>
class Grid3D
{
public:
	//////////////////////////////////////////////////////////////////////////////
	// attributes
	//////////////////////////////////////////////////////////////////////////////
	T _x;	// length on the X axis
	T _y;	// length on the Y axis
	T _z;	// length on the Z axis

	unsigned int _n;	// subdivision factor (size of box is (x/n, y/n, z/n))

	//////////////////////////////////////////////////////////////////////////////
	//  methods
	//////////////////////////////////////////////////////////////////////////////
	
	// constructors
	Grid3D() = delete;
	Grid3D(T x, T y, T z, unsigned int n);

	// used in voxelization
	void __pointerReconstruction(std::vector<Triangle<T>> &triangles, Vertex<T> target, Vertex<T> goal);

	// voxelization
	void voxelize(std::vector<vec3<T>> &vertices, std::vector<Triangle<T>> &triangles);
};

#include "Grid3D.hxx"

#endif

