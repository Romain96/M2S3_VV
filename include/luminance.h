// luminance.h: interface for the intens class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_LUMINANCE_H__16DDDDB6_B9E0_4BE8_A755_6F2550C60FA5__INCLUDED_)
#define AFX_LUMINANCE_H__16DDDDB6_B9E0_4BE8_A755_6F2550C60FA5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gmath/vec3.h"

// default luminance is simply triplet RGB
class luminance
{
	vec3<float>		v;
public:
	luminance() : v() { }
	luminance(float r, float g, float b) : v(r,g,b) { }
	luminance(float l) : v(l) { }
	float RED() const { return v.X(); }
	float GREEN() const { return v.Y(); }
	float BLUE() const { return v.Z(); }
	float GREY() const { return ((v.X()+v.Y()+v.Z())/3.0f); }

	void clamp(float min, float max) 
	{
		float r = v.X(); if (r<min) r=min; else if (r>max) r=max;
		float g = v.Y(); if (g<min) g=min; else if (g>max) g=max;
		float b = v.Z(); if (b<min) b=min; else if (b>max) b=max;
		*this = luminance(r, g, b);
	}
	void gammaNormalized(const luminance &max, float scal, float power) 
	{
		float r = v.X()/max.v.X(); if (r<0.0) r=0.0; else if (r>1.0) r=1.0;
		float g = v.Y()/max.v.Y(); if (g<0.0) g=0.0; else if (g>1.0) g=1.0;
		float b = v.Z()/max.v.Z(); if (b<0.0) b=0.0; else if (b>1.0) b=1.0;
		r = scal*pow(r,power);
		g = scal*pow(g,power);
		b = scal*pow(b,power);
		*this = luminance(r, g, b);
	}
	void gamma(float scal, float power) 
	{
		float r = v.X()/scal; if (r<0.0) r=0.0; else if (r>1.0) r=1.0;
		float g = v.Y()/scal; if (g<0.0) g=0.0; else if (g>1.0) g=1.0;
		float b = v.Z()/scal; if (b<0.0) b=0.0; else if (b>1.0) b=1.0;
		r = scal*pow(r,power);
		g = scal*pow(g,power);
		b = scal*pow(b,power);
		*this = luminance(r, g, b);
	}
	void blend(const luminance &c1, const luminance &c2, float alpha) 
	{
		float r = c1.RED()*alpha+c2.RED()*(1.0-alpha);
		float g = c1.GREEN()*alpha+c2.GREEN()*(1.0-alpha);
		float b = c1.BLUE()*alpha+c2.BLUE()*(1.0-alpha);
		*this = luminance(r, g, b);
	}
	void add(const luminance &a, const luminance &b) { v.add(a.v,b.v); }
	luminance operator+(const luminance &b) const { luminance x=*this; x.v.operator+=(b.v); return x; }
	void operator+=(const luminance &b) { v.operator+=(b.v); }
	void sub(const luminance &a, const luminance &b) { v.sub(a.v,b.v); }
	luminance operator-(const luminance &b) const { luminance x=*this; x.v.operator-=(b.v); return x; }
	void operator-=(const luminance &b) { v.operator-=(b.v); }
	void mult(const luminance &a, const luminance &b) { v.mult(a.v,b.v); }
	luminance operator*(const luminance &b) const { luminance x=*this; x.v.operator*=(b.v); return x; }

	void scale(const luminance &a, float k) { v.scale(a.v,k); }
	template <class X> void scale(X k) { v.scale((float)k); }

	void keepMin(const luminance &v1, const luminance &v2) { v.keepMin(v1.v, v2.v); }
	void keepMax(const luminance &v1, const luminance &v2) { v.keepMax(v1.v, v2.v); }

};

#endif // !defined(AFX_INTENS_H__16DDDDB6_B9E0_4BE8_A755_6F2550C60FA5__INCLUDED_)
