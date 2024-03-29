#ifndef __GRID3D_H__
#define __GRID3D_H__

// needing point structure for vertices
#include "Vertex.h"

// needing triangle structure
#include "Triangle.h"

// needing a voxel structure
#include "Voxel.h"

// needing the intersection
#include "gmath/intersection3.h"

#include <map>

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

	std::map<int, Voxel<T>> _voxels;	// the actual grid itself

	//////////////////////////////////////////////////////////////////////////////
	//  methods
	//////////////////////////////////////////////////////////////////////////////
	
	// constructors
	Grid3D() = delete;
	Grid3D(T x, T y, T z, unsigned int n);

	// used in voxelization : repair triangle to vertices pointers
	void __pointerReconstruction(std::vector<Triangle<T>> &triangles, Vertex<T> target, Vertex<T> goal);

	// used in voxelization : find all triangles intersecting a voxel
	void __findTriangleVoxelIntersection(std::vector<Triangle<T>> &triangles);

	// voxelization
	void voxelize(std::vector<Vertex<T>> &vertices, std::vector<Triangle<T>> &triangles);
};

#include "Grid3D.hxx"

#endif

