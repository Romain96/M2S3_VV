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
	// nothing
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

// triangle voxel intersection
// iterates on all voxels 
// for each voxel find all triangles intersecting it
// adding theses triangles to the voxel
template<typename T>
void Grid3D<T>::__findTriangleVoxelIntersection(std::vector<Triangle<T>> &triangles)
{
	// translation steps
	double dx = this->_x/this->_n;
	double dy = this->_y/this->_n;
	double dz = this->_z/this->_n;

	// iterating on all voxels
	for (double baseX = -this->_x/2.0; baseX < this->_x/2.0; baseX += dx)
	{
		for (double baseY = -this->_y/2.0; baseY < this->_y/2.0; baseY += dy)
		{
			for (double baseZ = -this->_z/2.0; baseZ < this->_z/2.0; baseZ += dz)
			{
				// computing the 8 vertices delimitating the voxel
				vec3<double> pos0(baseX, baseY, baseZ);
				vec3<double> posX(baseX + dx, baseY, baseZ);
				vec3<double> posY(baseX, baseY + dy, baseZ);
				vec3<double> posXY(baseX + dx, baseY + dy, baseZ);
				vec3<double> posZ(baseX, baseY, baseZ + dz);
				vec3<double> posZX(baseX + dx, baseY, baseZ + dz);
				vec3<double> posZY(baseX, baseY + dy, baseZ + dz);
				vec3<double> posZXY(baseX + dx, baseY + dy, baseZ + dz);

				// iterating on all triangles
				for (typename std::vector<Triangle<T>>::iterator it = triangles.begin(); it != triangles.end(); it++)
				{
					// creating the plane associated with the triangle
					double x1 = (*it)._b._coord.X() - (*it)._a._coord.X();
					double y1 = (*it)._b._coord.Y() - (*it)._a._coord.Y();
					double z1 = (*it)._b._coord.Z() - (*it)._a._coord.Z();

					double x2 = (*it)._c._coord.X() - (*it)._a._coord.X();
					double y2 = (*it)._c._coord.Y() - (*it)._a._coord.Y();
					double z2 = (*it)._c._coord.Z() - (*it)._a._coord.Z();

					vec3<double> axis1(x1, y1, z1);
					vec3<double> axis2(x2, y2 ,z2);

					plan3<double> tri((*it)._a._coord, axis1, axis2);

					// testing intersection between this plane and the 12 edges of the voxel
					int n = 0;
					intersection3<double> inter;

					// pos0 - posX
					if (inter.intersect(pos0, posX, tri, 0.0))
						n++;
					// pos0 - posY
					if (inter.intersect(pos0, posY, tri, 0.0))
						n++;
					// posX - posXY
					if (inter.intersect(posX, posX, tri, 0.0))
						n++;
					// posY - posXY
					if (inter.intersect(posY, posXY, tri, 0.0))
						n++;
					
					// posZ - posZX
					if (inter.intersect(posZ, posZX, tri, 0.0))
						n++;
					// posZ - posZY
					if (inter.intersect(posZ, posZY, tri, 0.0))
						n++;
					// posZX - posZXY
					if (inter.intersect(posZX, posZXY, tri, 0.0))
						n++;
					// posZY - posZXY
					if (inter.intersect(posZY, posZXY, tri, 0.0))
						n++;

					// pos0 - posZ
					if (inter.intersect(pos0, posZ, tri, 0.0))
						n++;
					// posX - posZX
					if (inter.intersect(posX, posZX, tri, 0.0))
						n++;
					// posY - posZY
					if (inter.intersect(posY, posZY, tri, 0.0))
						n++;
					// posXY - posZXY
					if (inter.intersect(posXY, posZXY, tri, 0.0))
						n++;

					std::cout << "n = " << n << std::endl;

					// if number of intersection is at least 3 then the plane intersect the voxel
					if (n >= 3)
					{
						// adding the triangle to the list of triangles intersected by the voxel
						int x = static_cast<int>(baseX + this->_x/2.0);
						int y = static_cast<int>(baseY + this->_y/2.0);
						int z = static_cast<int>(baseZ + this->_z/2.0);
						
						int key = x*_n*_n + y*_n + z;

						_voxels[key].addTriangle(&(*it));
					}

				}
			}
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
		int x, y, z;
		
		if ((*it)._coord.X() >= 0.0)
			x = static_cast<int>((*it)._coord.X() + this->_n/2.0);
		else
			x = static_cast<int>(-this->_x - this->_n/2 - (*it)._coord.X());

		if ((*it)._coord.Y() >= 0.0)
			y = static_cast<int>((*it)._coord.Y() + this->_n/2.0);
		else
			y = static_cast<int>(-this->_y - this->_n/2 - (*it)._coord.Y());

		if ((*it)._coord.Z() >= 0.0)
			z = static_cast<int>((*it)._coord.Z() + this->_n/2.0);
		else
			z = static_cast<int>(-this->_z - this->_n/2 - (*it)._coord.Z());

		// serializing
		int key = x*_n*_n + y*_n + z;

		std::cout << "XYZ = " << x << "," << y << "," << z << " key = " << key << std::endl;

		// inserting vertex at key position
		_voxels[key].addVertex(*it);

		// making triangles containing this vertex pointing on the newly created one in the voxel
		__pointerReconstruction(triangles, (*it), _voxels[key]._vertices.back());
	};

	// adding all triangles intersecting voxels as a list in the intersected voxel(s)
	__findTriangleVoxelIntersection(triangles);
}

#endif
