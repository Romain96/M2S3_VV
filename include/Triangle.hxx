#ifndef __TRIANGLE_HXX__
#define __TRIANGLE_HXX__

#include "Triangle.h"

// default constructor
template<typename T>
Triangle<T>::Triangle() :
	_a(nullptr),
	_b(nullptr),
	_c(nullptr),
	_A(T(0)),
	_B(T(0)),
	_C(T(0)),
	_Ka(T(0)),
	_Kb(T(0)),
	_Kc(T(0))
{
	// nothing
}

// parametrized constructor
template<typename T>
Triangle<T>::Triangle(Vertex<T> a, Vertex<T> b, Vertex<T> c) : 
	_a(a),
	_b(b),
	_c(c),
	_Ka(T(0)),
	_Kb(T(0)),
	_Kc(T(0))
{
	// computing lengths
	_A = static_cast<T>((_c._coord.X()-_b._coord.X())*(_c._coord.X()-_b._coord.X()) + (_c._coord.Y()-_b._coord.Y())*(_c._coord.Y()-_b._coord.Y()) + (_c._coord.Z()-_b._coord.Z())*(_c._coord.Z()-_b._coord.Z()));
	_B = static_cast<T>((_c._coord.X()-_a._coord.X())*(_c._coord.X()-_a._coord.X()) + (_c._coord.Y()-_a._coord.Y())*(_c._coord.Y()-_a._coord.Y()) + (_c._coord.Z()-_a._coord.Z())*(_c._coord.Z()-_a._coord.Z()));
	_C = static_cast<T>((_b._coord.X()-_a._coord.X())*(_b._coord.X()-_a._coord.X()) + (_b._coord.Y()-_a._coord.Y())*(_b._coord.Y()-_a._coord.Y()) + (_b._coord.Z()-_a._coord.Z())*(_b._coord.Z()-_a._coord.Z()));
}

#endif

