// box3.h: interface for the box3 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BOX3_H__6BB85A0D_9425_46A7_B2C5_039BA5A6227C__INCLUDED_)
#define AFX_BOX3_H__6BB85A0D_9425_46A7_B2C5_039BA5A6227C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "gmath/vec3.h"

template <class T> class box3  
{
protected:
	vec3<T>		min, max;
public:
	// constructor
	box3<T>() { }
	box3<T>(const vec3<T> &v) { min=v; max=v;  }
	box3<T>(const vec3<T> &m, const vec3<T> &n) { min=m; max=n; }

	// size of box: width, height, length, given by resulting 3-D vector
	vec3<T> size() const { vec3<T> v; v.PVec(min, max); return v; }
	operator vec3<T>() { vec3<T> v; v.PVec(min, max); return v; }
	T sizeX() const { return max.X()-min.X(); }
	T sizeY() const { return max.Y()-min.Y(); }
	T sizeZ() const { return max.Z()-min.Z(); }

	// get min and max
	vec3<T> minPoint() const { return min; }
	vec3<T> maxPoint() const { return max; }

	// Adds a point to the box -> resize the box to include this point
	void addPoint(const vec3<T> &p)
	{
		min.keepMin(min, p);
		max.keepMax(max, p);
	}

	void enlarge(double dist)
	{
		vec3<T> s = size();
		s *= dist;
		max+=s; min-=s;
	}

	void getCorners(vec3<T> p[8]) const
	{
		vec3<T> dir; dir.PVec(min, max);
		min.boxPoints(dir, T(1), p);
	}


	vec3<T> middle() const
	{
		vec3<T> p; p.add(min, max);
		p.scale(0.5);
		return p;
	}

	vec3<T> randPointInside() const
	{
		vec3<T> rnd; rnd.random();
		vec3<T> dir; dir.PVec(min, max);
		dir.mult(dir,rnd);
		dir.add(dir, min);
		return dir;
	}
	bool operator==(const box3<T> &bb) const
		{
		return (bb.min==min) && (bb.max==max) ;
		}

	// Test if point p is inside the box defined by min-(err,err,err) and max+(err,err,err)
	bool isInside(const vec3<T> &p, T err = T(0)) const
	{
		if (p.X()<min.X()-err || p.X()>max.X()+err) return false;
		if (p.Y()<min.Y()-err || p.Y()>max.Y()+err) return false;
		if (p.Z()<min.Z()-err || p.Z()>max.Z()+err) return false;
		return true;
	}

	bool isOverlapping(const box3<T> &b) const
	{
		if (max.X() <= b.min.X() || b.max.X() <= min.X()) return false;
		if (max.Y() <= b.min.Y() || b.max.Y() <= min.Y()) return false;
		if (max.Z() <= b.min.Z() || b.max.Z() <= min.Z()) return false;
		return true;
	}

	bool hasFlatX(T err) const
	{
		vec3<T> v; v.PVec(min, max);
		if (v.X()<=err) return true;
		return false;
	}

	bool hasFlatY(T err) const
	{
		vec3<T> v; v.PVec(min, max);
		if (v.Y()<=err) return true;
		return false;
	}

	bool hasFlatZ(T err) const
	{
		vec3<T> v; v.PVec(min, max);
		if (v.Z()<=err) return true;
		return false;
	}

	void splitt2(box3<T> &ba, box3<T> &bb) const
	{
		vec3<T> v=size();
		if (v.X()>v.Y() && v.X()>v.Z() )
		{
			ba = box3<T> ( min, vec3<T>( min.X()+v.X()/2, max.Y(), max.Z()) );
			bb = box3<T> ( vec3<T>( min.X()+v.X()/2, min.Y(), min.Z()), max );
		}
		else if ( v.Y() > v.Z() )
		{
			ba = box3<T> ( min, vec3<T>( max.X(), min.Y()+v.Y()/2, max.Z() ) );
			bb = box3<T> ( vec3<T>( min.X(), min.Y()+v.Y()/2, min.Z() ), max );
		}
		else
		{
			ba = box3<T> ( min, vec3<T>( max.X(), max.Y(), min.Z()+v.Z()/2 ) );
			bb = box3<T> ( vec3<T>( min.X(), min.Y(), min.Z()+v.Z()/2 ), max );
		}
	}

	void splitt8(box3<T> bf[8]) const
	{
		vec3<T> v=size();

		bf[0].min = min;
		bf[0].max = vec3<T>(min.X()+v.X()/2, min.Y()+v.Y()/2, min.Z()+v.Z()/2);

		bf[1].min = vec3<T>(min.X(),min.Y()+v.Y()/2,min.Z());
		bf[1].max = vec3<T>(min.X()+v.X()/2,max.Y(),min.Z()+v.Z()/2);

		bf[2].min = vec3<T>(min.X()+v.X()/2,min.Y(),min.Z());
		bf[2].max = vec3<T>(max.X(),min.Y()+v.Y()/2,min.Z()+v.Z()/2);

		bf[3].min = vec3<T>(min.X()+v.X()/2,min.Y()+v.Y()/2,min.Z());
		bf[3].max = vec3<T>(max.X(),max.Y(),min.Z()+v.Z()/2);

		bf[4].min = vec3<T>(min.X(),min.Y(),min.Z()+v.Z()/2);
		bf[4].max = vec3<T>(min.X()+v.X()/2,min.Y()+v.Y()/2,max.Z());

		bf[5].min = vec3<T>(min.X(),min.Y()+v.Y()/2,min.Z()+v.Z()/2);
		bf[5].max = vec3<T>(min.X()+v.X()/2,max.Y(),max.Z());

		bf[6].min = vec3<T>(min.X()+v.X()/2,min.Y(),min.Z()+v.Z()/2);
		bf[6].max = vec3<T>(max.X(),min.Y()+v.Y()/2,max.Z());

		bf[7].min = vec3<T>(min.X()+v.X()/2,min.Y()+v.Y()/2,min.Z()+v.Z()/2);
		bf[7].max = max;
	}

	void splitt(box3<T> bf[8]) const
	{
		splitt8(bf);
	}

	void merge(const box3<T> &b1, const box3<T> &b2)
	{
		min.keepMin(b1.min, b2.min);
		max.keepMax(b1.max, b2.max);
	}

};

#endif // !defined(AFX_BOX3_H__6BB85A0D_9425_46A7_B2C5_039BA5A6227C__INCLUDED_)
