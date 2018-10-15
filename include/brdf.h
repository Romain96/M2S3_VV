// brdf.h: interface for the brdf class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BRDF_H__12C54D2A_555B_41AD_A623_A0E629CFC1A2__INCLUDED_)
#define AFX_BRDF_H__12C54D2A_555B_41AD_A623_A0E629CFC1A2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "luminance.h"

// virtual class 
class brdf  
{
protected:
	char		 	name[64];

	luminance	  	source;		/* self-emission */

	luminance		Kr;		/* perfect Mirror reflection */

	luminance		Kt;		/* perfect Transparency */
	float			refract; 	/* Refraction indice (transparency) */

public:
	brdf():source(), Kr(), Kt(), refract(1.0f) { name[0]='\0'; }
	brdf(char *n):source(), Kr(), Kt(), refract(1.0f) { strcpy(name, n); }

	void setName(char *n) { strcpy(name,n); }

	bool isSource() const { return source.GREY()>0.0f; }
	void setSelf(const luminance &li) { source=li; }
	void getSelf(luminance &li) const { li=source; }

	virtual bool isSeparable() const =0;
	virtual bool isMirror() const { return Kr.GREY()>0.0f; }
	virtual bool isTransparent() const { return Kt.GREY()>0.0f; }
	virtual bool isSpecular() const =0;
	virtual float minSpecular(float minerr) const =0;

	virtual void setMirrorReflection(const luminance &r)  { Kr=r; } 
	virtual void setTransparency(const luminance &t)  { Kt=t; } 
	virtual void setRefractionIndex(float ii)  { refract=ii; } // must be greater than 1
	virtual float getRefractionIndex() const { return refract; }

	// for interactive display using Phong model
	virtual float phongDiffuseRED() const =0;
	virtual float phongDiffuseGREEN() const =0;
	virtual float phongDiffuseBLUE() const =0;
	virtual float phongSpecRED() const =0;
	virtual float phongSpecGREEN() const =0;
	virtual float phongSpecBLUE() const =0;
	virtual float phongSpecN() const =0;
	
	virtual void applyMirror(const luminance &itin, void *param, luminance &itout)
	{
		itout=itin*Kr;
	}
	virtual void applyTransparency(const luminance &itin, void *param, luminance &itout)
	{
		itout=itin*Kt;
	}

	virtual void apply(const luminance &itin, const vec3<float> &normal, const vec3<float> &lightdir, 
		const vec3<float> &eyedir, const vec3<float> &u, const vec3<float> &v, void *param, luminance &itout) =0;
	virtual void applyDiffuse(const luminance &itin, const vec3<float> &normal, const vec3<float> &lightdir, 
		const vec3<float> &eyedir, const vec3<float> &u, const vec3<float> &v, void *param, luminance &itout) =0;
	virtual void applySpecular(const luminance &itin, const vec3<float> &normal, const vec3<float> &lightdir, 
		const vec3<float> &eyedir, const vec3<float> &u, const vec3<float> &v, void *param, luminance &itout) =0;

};

#endif // !defined(AFX_BRDF_H__12C54D2A_555B_41AD_A623_A0E629CFC1A2__INCLUDED_)
