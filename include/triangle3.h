// triangle3.h: interface for the triangle3 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TRIANGLE3_H__6CF19E40_72AA_41BC_857A_0E338EA12B93__INCLUDED_)
#define AFX_TRIANGLE3_H__6CF19E40_72AA_41BC_857A_0E338EA12B93__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gmath/intersection3.h"

const float TANG_ERROR=0.001;

/////////////////////////////////////
/////////////////////////////////////
template <class T> class triangle3
{
protected:
	vec3<T>				pt, a, b, c;
	T					ka, kb, kc;		// coefficients for bilinear interpolation
public:
	triangle3<T>(const vec3<T> &x, const vec3<T> &y, const vec3<T> &z):a(x),b(y),c(z) { ka=T(0); kb=T(0); kc=T(0); }
	triangle3<T>() { ka=T(0); kb=T(0); kc=T(0); }

	void setVertices(const vec3<T> &x, const vec3<T> &y, const vec3<T> &z) { a=x; b=y; c=z; }

	int nVertices() const { return 3; }
	T Ka() const { return ka; }
	T Kb() const { return kb; }
	T Kc() const { return kc; }
	vec3<T> operator[](int i) const
	{
		if (i<0 || i>2) { printf("out of range in triangle3!"); return vec3<T>(T(0)); }
		if (i==0) return a;
		else if (i==1) return b;
		return c;
	}
	vec3<T> getVertex(int i) const
	{
		if (i<0 || i>2) { printf("out of range in hvTriangle!"); return vec3<T>(T(0)); }
		if (i==0) return a;
		else if (i==1) return b;
		return c;
	}

	vec3<T> bary() const
	{
		T bary=a;
		bary += b;
		bary += c;
		bary.divScale(T(3.0));
		return bary;
	}
	operator vec3<T>() const 
	{
		return this->bary();
	}
	vec3<T> inPoint() const { return pt; }

	vec3<T> normal() const
	{
		vec3<T> u; u.PVec(a, b);
		vec3<T> v; v.PVec(a, c);
		vec3<T> w; w.cross(u,v);
		T n=w.norm();
		if (n!=T(0)) w.normalize(n);
		return w;
	}

	T area() const
	{
		vec3<T> u; u.PVec(a, b);
		vec3<T> v; v.PVec(a, c);
		vec3<T> w; w.cross(u,v);
		T n=w.norm();
		return n/2.0;
	}
	operator T() const
	{
		return this->area();
	}

	operator plan3<T>() const
	{
		vec3<T> u; u.PVec(a, b);
		vec3<T> v; v.PVec(a, c);
		return plan3<T>(a, u, v);
	}

	int intersect(const line3<T> &li, T err, intersection3<T> &i)
	{
		vec3<T> u,v, aa;
		u.PVec(a,b); 
		v.PVec(a,c); 
		aa.cross(u,v);
		// if triangle is degenerated, no intersection
		if (aa.isNull()) return 0;
		T daa= aa.norm();
		aa.normalize(daa);

		// intersection with plane
		plan3<T> pl(a, aa);
		if (!i.intersect(li, pl, TANG_ERROR)) return 0;
		if (!isInside(i[0], err)) return 0;
		return 1;
	}
	virtual int intersectRay(const line3<T> &li, T err, intersection3<T> &i)
	{
		vec3<T> u,v, aa;
		u.PVec(a,b); 
		v.PVec(a,c); 
		aa.cross(u,v);
		// if triangle is degenerated, no intersection
		if (aa.isNull()) return 0;
		T daa= aa.norm();
		aa.normalize(daa);

		// intersection with plane
		plan3<T> pl(a, aa);
		if (!i.intersectRay(li, pl, TANG_ERROR)) return 0;
		if (!isInside(i[0], err)) return 0;
		return 1;
	}

	bool isInside(const vec3<T> &x, T err)
	{
		vec3<T> u,v, aa;
		u.PVec(a,b); v.PVec(a,c); aa.cross(u,v);
		if (aa.isNull()) return false;
		T daa= aa.norm();
		aa.normalize(daa);
		return isInside(x, err, aa);
	}
	bool isInside(const vec3<T> &x, T err, const vec3<T> &aa)
	{
		vec3<T> u,v, a1, a2, a3;
		T da1, da2, da3, daa;

		u.PVec(x,b); v.PVec(x,c); a1.cross(u,v);
		u.PVec(x,c); v.PVec(x,a); a2.cross(u,v);
		u.PVec(x,a); v.PVec(x,b); a3.cross(u,v);
		// if point is close to border, consider as inside
		da1 = a1.dot(aa);
		da2 = a2.dot(aa);
		da3 = a3.dot(aa);
		//if (da1< err && da1>T(0)) da1 = -da1;
		//if (da2< err && da2>T(0)) da2 = -da2;
		//if (da3< err && da3>T(0)) da3 = -da3;
		if (da1> -err && da1<T(0) ) da1= -da1;
		if (da2> -err && da2<T(0) ) da2= -da2;
		if (da3> -err && da3<T(0) ) da3= -da3;
		if (da1<T(0) || da2<T(0) || da3<T(0)) return false;
		da1 = a1.norm(); da2 = a2.norm(); da3 = a3.norm();

		daa = da1+da2+da3;
		ka=da1/daa; 
		kb=da2/daa; 
		kc=da3/daa;
		pt = x;
		return true;
	}
	
	// can be used only after computing an intersection and only if the intersection was successful
	T interpolate(T x[3]) const
	{
		return ka*x[0]+kb*x[1]+kc*x[2];
	}
	vec3<T> interpolate(const vec3<T> x[3]) const 
	{
		vec3<T> res(T(0));
		vec3<T> sv = x[0]; sv.scale(ka);
		res += sv;
		sv = x[1]; sv.scale(kb);
		res += sv;
		sv = x[2]; sv.scale(kc);
		res += sv;
		return res;
	}
};



#endif // !defined(AFX_TRIANGLE3_H__6CF19E40_72AA_41BC_857A_0E338EA12B93__INCLUDED_)
