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


#endif
