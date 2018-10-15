// phongbrdf.h: interface for the Phong-like brdf class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PHONGBRDF_H__12C54D2A_555B_41AD_A623_A0E629CFC1A2__INCLUDED_)
#define AFX_PHONGBRDF_H__12C54D2A_555B_41AD_A623_A0E629CFC1A2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "brdf.h"

// brdf corresponding to the Phong model
class phongbrdf: public brdf  
{
protected:
	luminance     	Kd;		/* diffus. 0<=Kd<=1 */
	luminance     	Ks;		/* Specul. 0<=Ks<=1 */
	float			n;		/* power of cosinus for specular lobe*/

	float			minangle;
public:
	phongbrdf():brdf(), Kd(0.5), Ks(0.0),n(0.0) { }
	phongbrdf(char *s, const luminance &diff, const luminance &spec, float phongn):brdf(s) 
	{ 
		Kd=diff; Ks=spec; n=phongn; 
	} 

	virtual bool isSeparable() const { return true; }
	virtual bool isSpecular() const { return Ks.GREY()>0.0f; }
	virtual float minSpecular(float minerr) const
	{
		return acos((float)pow((double)minerr, 1.0/(double)n));
	}


	virtual float phongDiffuseRED() const { return Kd.RED(); }
	virtual float phongDiffuseGREEN() const { return Kd.GREEN(); }
	virtual float phongDiffuseBLUE() const { return Kd.BLUE(); }
	virtual float phongSpecRED() const { return Ks.RED(); }
	virtual float phongSpecGREEN() const { return Ks.GREEN(); }
	virtual float phongSpecBLUE() const { return Ks.BLUE(); }
	virtual float phongSpecN() const { return n; }

	void apply(const luminance &itin, const vec3<float> &normal, const vec3<float> &lightdir, 
		const vec3<float> &eyedir, const vec3<float> &u, const vec3<float> &v, void *param, luminance &itout)
	{
		applyDiffuse(itin,normal,lightdir,eyedir,u,v,param,itout);
		luminance spec;
		applySpecular(itin,normal,lightdir,eyedir,u,v,param,spec);
		itout+=spec;
	}
	void applyDiffuse(const luminance &itin, const vec3<float> &normal, const vec3<float> &lightdir, 
		const vec3<float> &eyedir, const vec3<float> &u, const vec3<float> &v, void *param, luminance &itout)
	{
		itout=itin*Kd;
	}
	void applySpecular(const luminance &itin, const vec3<float> &normal, const vec3<float> &lightdir, 
		const vec3<float> &eyedir, const vec3<float> &u, const vec3<float> &v, void *param, luminance &itout)
	{
		float beta;
		itout=itin*Ks;
		vec3<float> r; r.reflection(eyedir,normal);
		beta = (float)r.dot(lightdir);
		if (beta > 0.0)
			{
			beta = (float)pow(beta,n);
			//printf("beta=%g, %g,%g,%g - %g,%g,%g\n", beta, r.X(), r.Y(), r.Z(), lightdir.X(), lightdir.Y(), lightdir.Z());
			itout.scale(beta);
			}
		else itout=luminance();
	}

};

#endif // !defined(AFX_BRDF_H__12C54D2A_555B_41AD_A623_A0E629CFC1A2__INCLUDED_)
