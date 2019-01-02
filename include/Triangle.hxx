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
Triangle<T>::Triangle(vec3<T> a, vec3<T> b, vec3<T> c) : 
	_a(a),
	_b(b),
	_c(c),
	_Ka(T(0)),
	_Kb(T(0)),
	_Kc(T(0))
{
	// computing lengths
	_A = static_cast<T>((_c.X()-_b.X())*(_c.X()-_b.X()) + (_c.Y()-_b.Y())*(_c.Y()-_b.Y()) + (_c.Z()-_b.Z())*(_c.Z()-_b.Z()));
	_B = static_cast<T>((_c.X()-_a.X())*(_c.X()-_a.X()) + (_c.Y()-_a.Y())*(_c.Y()-_a.Y()) + (_c.Z()-_a.Z())*(_c.Z()-_a.Z()));
	_C = static_cast<T>((_b.X()-_a.X())*(_b.X()-_a.X()) + (_b.Y()-_a.Y())*(_b.Y()-_a.Y()) + (_b.Z()-_a.Z())*(_b.Z()-_a.Z()));
}

#endif

