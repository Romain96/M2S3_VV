// scene_stochastic.h
//
// scene defines a set of 3D meshes
// Main operation is: stochastic raytracing
// 
// By JMD 10/8/06
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SCENE_STOCH_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_)
#define AFX_SCENE_STOCH_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <GL/glew.h>
#include <GL/glut.h>

#include "scene.h"

const int MAX_SUBPIXELS = 500;
const float LUMI_VAR = 0.001;

///////////////////////////////////////////////
template <class T> class scene_stochastic: public scene<T>
{
protected:
	// surface lighting
	int						nslights;
	vec3<T>					slightsp[MAX_LIGHTS], slightsU[MAX_LIGHTS], slightsV[MAX_LIGHTS];
	luminance				slightit[MAX_LIGHTS];

	// stochastic ray tracing data
	int			nrays;			// nrays * nrays subpixels
	int			nspec;			// nspec*nspec specular sub-rays
	int			nli;			// nli*nli shadow rays
	T			mirr_angle;

public:
	// constructors
	scene_stochastic<T>(int resX, int resY): scene<T>(resX, resY) 
	{ 
		nslights=0;
		nrays = 3;
		mirr_angle=PI/8.0;
		nspec=2;
		nli=2;
	}
	void setRays(int n, int ns, int nl) { nrays=n; nspec=ns; nli = nl; } 
	void SetMirrBlurr(T mm) { mirr_angle=mm; }

	void pushAreaLight(const vec3<T> pos, const vec3<T> uu, const vec3<T> vv, const luminance &it)
	{
		if (nslights>=MAX_LIGHTS) return;
		slightsp[nslights]=pos; slightsU[nslights]=uu; slightsV[nslights]=vv;
		slightit[nslights]=it; 
		nslights++;
	}


// Stochastic RAYTRACING
////////////////////////////////////////////////////////////
	virtual void raytracing()
	{
			int subip2[4]={1, 2, 1, 2 }, subjp2[4]={1, 1, 2, 2 };
			int subip3[9]={1,3,1,2,3,3,2,2,1};
			int subjp3[9]={1,2,3,2,1,3,1,3,2};
			int subip4[16]={1,4,1,3,4,2,3,4,2,3,1,4,3,2,1,2};
			int subjp4[16]={1,4,4,2,1,2,4,3,3,1,2,2,3,1,3,4};
			int subip5[25]={2,5,2,5,3,1,1,5,3,4,4,2,4,1,3,2,4,5,4,2,1,3,5,3,1};
			int subjp5[25]={2,1,5,4,3,1,4,3,1,5,2,4,3,2,4,3,1,2,4,1,5,2,5,5,3};	

			int i,j, k;
			vec3<float> dir,up; 
			dir.PVec(scene<T>::eyepos, scene<T>::lookpos);
			dir.normalize(dir.norm());
			up.orthogonalInPlan(dir,vec3<float>(0.0,0.0,1.0));
			perspective<float> persp(scene<T>::eyepos,scene<T>::lookpos, up, 0.84f, (float)(0.84f*scene<T>::H/scene<T>::W), 1.0f, 50.0f);
			
			unsigned char *picture=(unsigned char *)malloc(sizeof(unsigned char)*scene<T>::W*scene<T>::H*3);
			for (i=0; i<scene<T>::W*scene<T>::H*3; i++) picture[i]=0;
			
			// for each pixel of the screen
			for (j=0; j<scene<T>::H; j++) 
			{
					for (i=0; i<scene<T>::W; i++) 
					{
						luminance sum(0.0);
						float colr[MAX_SUBPIXELS], colg[MAX_SUBPIXELS], colb[MAX_SUBPIXELS];
						float sumr=0.0, sumg=0.0, sumb=0.0;
						float sub = 1.0/nrays;
						bool cont=true;
						for (k=0; k<nrays*nrays && cont; k++)
						{
						float px, py;
						switch(nrays)
							{
							case 1: { px = (float)rand()/(float)RAND_MAX*sub; py = (float)rand()/(float)RAND_MAX*sub; break; }
							case 2: { px = (subip2[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp2[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }
							case 3: { px = (subip3[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp3[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }
							case 4: { px = (subip4[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp4[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }
							case 5: { px = (subip5[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp5[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }	
							default: { px = (subip5[k%25]-1+(float)rand()/(float)RAND_MAX)/5.0; py=(subjp5[k%25]-1+(float)rand()/(float)RAND_MAX)/5.0; break; }
							}
						line3<float> dd = persp.reverse(vec3<float>(1.0-2.0*((float)i+px)/(float)scene<T>::W, -1.0+2.0*((float)j+py)/(float)scene<T>::H, 0.0));
						luminance col = scene_stochastic::primaryRay(dd, k);
						col.clamp(0.0f, 1.0f);
						sum += col;
						colr[k]=col.RED(); colg[k]=col.GREEN(); colb[k]=col.BLUE();
						sumr += colr[k]; sumg += colg[k]; sumb += colb[k];
						if (k>=5) 
							{
								float varr=0.0, varg=0.0, varb=0.0;
								for (int ii=0; ii<k; ii++)
									{
										varr+= (colr[ii]-sumr/(float)(k+1))*(colr[ii]-sumr/(float)(k+1));
										varg+= (colg[ii]-sumg/(float)(k+1))*(colg[ii]-sumg/(float)(k+1));
										varb+= (colb[ii]-sumb/(float)(k+1))*(colb[ii]-sumb/(float)(k+1));
									}
								varr /= (float)(k+1); varg /= (float)(k+1); varb /= (float)(k+1);
								if (sqrt(varr)<LUMI_VAR && sqrt(varg)<LUMI_VAR && sqrt(varb)<LUMI_VAR) cont=false;
							}
						}
						picture[3*(i+j*scene<T>::W)]=(unsigned char)(sumr/(float)(k+1)*255.0);
						picture[3*(i+j*scene<T>::W)+1]=(unsigned char)(sumg/(float)(k+1)*255.0);
						picture[3*(i+j*scene<T>::W)+2]=(unsigned char)(sumb/(float)(k+1)*255.0);
					}
				// update drawing every 10th line
				if (j%5==0 || j==scene<T>::H-1)
					{
					glUseProgram(0);
					glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT );
					glViewport(0, 0, scene<T>::W, scene<T>::H);			
					glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
					glWindowPos2i(0, 0);
					scene<T>::CHECK_GL_ERROR();
					glDrawPixels(scene<T>::W, scene<T>::H, GL_RGB, GL_UNSIGNED_BYTE, picture);
					glutSwapBuffers(); 
					printf("line %d\n", j);
					}
			}
	}
	virtual void trace(int i, int j)
	{
			int subip2[4]={1, 2, 1, 2 }, subjp2[4]={1, 1, 2, 2 };
			int subip3[9]={1,3,1,2,3,3,2,2,1};
			int subjp3[9]={1,2,3,2,1,3,1,3,2};
			int subip4[16]={1,4,1,3,4,2,3,4,2,3,1,4,3,2,1,2};
			int subjp4[16]={1,4,4,2,1,2,4,3,3,1,2,2,3,1,3,4};
			int subip5[25]={2,5,2,5,3,1,1,5,3,4,4,2,4,1,3,2,4,5,4,2,1,3,5,3,1};
			int subjp5[25]={2,1,5,4,3,1,4,3,1,5,2,4,3,2,4,3,1,2,4,1,5,2,5,5,3};	
			vec3<float> dir,up; 
			dir.PVec(scene<T>::eyepos, scene<T>::lookpos);
			dir.normalize(dir.norm());
			up.orthogonalInPlan(dir,vec3<float>(0.0,0.0,1.0));
			perspective<float> persp(scene<T>::eyepos,scene<T>::lookpos, up, 0.84f, (float)(0.84f*scene<T>::H/scene<T>::W), 1.0f, 50.0f);
			unsigned char *picture=(unsigned char *)malloc(sizeof(unsigned char)*scene<T>::W*scene<T>::H*3);
			float colr[MAX_SUBPIXELS], colg[MAX_SUBPIXELS], colb[MAX_SUBPIXELS];
			float sumr=0.0, sumg=0.0, sumb=0.0;
			float sub = 1.0/nrays;
			bool cont=true;
			int k;
			luminance sum(0.0);
			for (k=0; k<nrays*nrays && cont; k++)
			{
			float px, py;
			switch(nrays)
				{
				case 1: { px = (float)rand()/(float)RAND_MAX*sub; py = (float)rand()/(float)RAND_MAX*sub; break; }
				case 2: { px = (subip2[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp2[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }
				case 3: { px = (subip3[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp3[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }
				case 4: { px = (subip4[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp4[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }
				case 5: { px = (subip5[k]-1+(float)rand()/(float)RAND_MAX)*sub; py=(subjp5[k]-1+(float)rand()/(float)RAND_MAX)*sub; break; }	
				default: { px = (subip5[k%25]-1+(float)rand()/(float)RAND_MAX)/5.0; py=(subjp5[k%25]-1+(float)rand()/(float)RAND_MAX)/5.0; break; }
				}
				line3<float> dd = persp.reverse(vec3<float>(1.0-2.0*((float)i+px)/(float)scene<T>::W, -1.0+2.0*((float)j+py)/(float)scene<T>::H, 0.0));
				printf("PIXEL %d,%d, sub %d\n", i,j, k);
				luminance col = scene_stochastic::primaryRay(dd, k);
				printf("color %d = %g,%g,%g\n", k, col.RED(), col.GREEN(), col.BLUE());
				sum += col;
				colr[k]=col.RED(); colg[k]=col.GREEN(); colb[k]=col.BLUE();
				sumr += colr[k]; sumg += colg[k]; sumb += colb[k];
				if (k>=5) 
					{
						float varr=0.0, varg=0.0, varb=0.0;
						for (int ii=0; ii<k; ii++)
							{
								varr+= (colr[ii]-sumr/(float)(k+1))*(colr[ii]-sumr/(float)(k+1));
								varg+= (colg[ii]-sumg/(float)(k+1))*(colg[ii]-sumg/(float)(k+1));
								varb+= (colb[ii]-sumb/(float)(k+1))*(colb[ii]-sumb/(float)(k+1));
							}
						varr /= (float)(k+1); varg /= (float)(k+1); varb /= (float)(k+1);
						if (sqrt(varr)<LUMI_VAR && sqrt(varg)<LUMI_VAR && sqrt(varb)<LUMI_VAR) cont=false;
					}
			}
			printf("FINAL COLOR: %g,%g,%g\n", sumr/(float)(k+1)*255.0, sumg/(float)(k+1)*255.0,sumb/(float)(k+1)*255.0);
	}	

	luminance primaryRay(line3<float> dd, int isub)
	{
		int count;
		scenepoint<T>  ip;
		if (scene<T>::verbose) scene<T>::ncbox=0;
		int first=visibility(scene<T>::bbhscene, dd, scene<T>::CLOSEST, count, &ip, 1, 0, 0, 0);
		if (scene<T>::verbose)  printf("%d box intersections\n", scene<T>::ncbox);
		// no intersection 
		if (first==0) return scene<T>::back;
		if (verbose)
		{
			rayorig[nbrays*6+0]=dd.origin().X(); rayorig[nbrays*6+1]=dd.origin().Y(); rayorig[nbrays*6+2]=dd.origin().Z();
			rayorig[nbrays*6+3]=ip.getPoint().X(); rayorig[nbrays*6+4]=ip.getPoint().Y(); rayorig[nbrays*6+5]=ip.getPoint().Z();
			nbrays++;
		}
		// compute shading for closest point
		return scene_stochastic::localIllumination(dd, &ip, 0, isub);
	}
	luminance secondaryRay(scenepoint<T> *start, line3<float> dd, int depth, int isub)
	{
		int count;
		scenepoint<T>  ip;
		bool verb=scene<T>::verbose;
		if (verb&&depth>0) scene<T>::verbose=false;
		int first=visibility(scene<T>::bbhscene, dd, scene<T>::CLOSEST, count, &ip, 1, start, 0, 0);
		if (scene<T>::verbose)  printf("%d box intersections\n", scene<T>::ncbox);
		if (verb&&depth>0) scene<T>::verbose=true;
		// no intersection 
		if (first==0) return scene<T>::back;
		if (verbose)
		{
			rayorig[nbrays*6+0]=dd.origin().X(); rayorig[nbrays*6+1]=dd.origin().Y(); rayorig[nbrays*6+2]=dd.origin().Z();
			rayorig[nbrays*6+3]=ip.getPoint().X(); rayorig[nbrays*6+4]=ip.getPoint().Y(); rayorig[nbrays*6+5]=ip.getPoint().Z();
			nbrays++;
		}
		// compute shading for closest point
		return scene_stochastic::localIllumination(dd, &ip, depth, isub);
	}


	luminance localIllumination(line3<float> dd, scenepoint<T> *ip, int depth, int isub)
	{
		bool verb=scene<T>::verbose;
		if (scene<T>::verbose) 
		{ 
			mesh<T> *m = scene<T>::mm[scene<T>::meshid[ip->getObj()]];
			printf("Local illum: %d(%d), point %g,%g,%g\n", ip->getObj(), ip->getFace(), ip->getPoint().X(), ip->getPoint().Y(),ip->getPoint().Z());  
			triangle3<T> tri = ip->getTriangle();
			printf("on triangle: (%g,%g,%g)-(%g,%g,%g)-(%g,%g,%g)\n", tri.getVertex(0).X(), tri.getVertex(0).Y(), tri.getVertex(0).Z(),
				tri.getVertex(1).X(), tri.getVertex(1).Y(), tri.getVertex(1).Z(),
				tri.getVertex(2).X(), tri.getVertex(2).Y(), tri.getVertex(2).Z() );
			unsigned int nf=m->nAdjacentFaces(ip->getFace());
			printf("%d neighbors are: ",nf); 
			unsigned int *nn=m->getAdjacentFaces(ip->getFace());
			for (int i=0; i<nf; i++) printf("%d ", nn[i]);
			printf("\n");
		}
		brdf *r = scene<T>::refl[scene<T>::brdfid[ip->getObj()]];
		vec3<T> obsdir = dd.direction(); obsdir.reverse();
		vec3<T> no = getOrientedNormal(obsdir, ip);
		vec3<T> uu; uu.orthogonal(no);
		vec3<T> vv; vv.cross(no, uu); vv.normalize(vv.norm());

		int i,j;
		luminance col(0.0f);
		for (i=0; i<scene<T>::nlights; i++)
		{
			vec3<T> ll=vec3<T>(scene<T>::lights[i].X(), scene<T>::lights[i].Y(), scene<T>::lights[i].Z());
			luminance itin=scene<T>::lightit[i];
			luminance res(0.0f);
			vec3<T> dir;
			if (scene<T>::lights[i].W()!=0.0) dir.PVec(ip->getPoint(), ll); else dir=ll;
			dir.normalize(dir.norm());
			line3<T> dd(ip->getPoint(), dir);
			scenepoint<T> shadow;
			scenepoint<T> end(ll);
			int count;

			float cc = dir.dot(no);
			if (scene<T>::verbose) { printf("Light %d (%g,%g,%g) Normal %g,%g,%g  dot: %g\n", i, itin.RED(), itin.GREEN(), itin.BLUE(), no.X(), no.Y(),no.Z(), cc);  }
			if (cc>0.0f)
			{
				itin.scale(cc);
				if (scene<T>::verbose)
				{
					srayorig[snbrays*6+0]=dd.origin().X(); srayorig[snbrays*6+1]=dd.origin().Y(); srayorig[snbrays*6+2]=dd.origin().Z();
					srayorig[snbrays*6+3]=ll.X(); srayorig[snbrays*6+4]=ll.Y(); srayorig[snbrays*6+5]=ll.Z();
					snbrays++;
				}
				//if (true)
				if (verb) scene<T>::verbose=false;
				if (visibility(scene<T>::bbhscene, dd, scene<T>::OCCLUDER, count, &shadow, 1, ip, scene<T>::lights[i].W()==0.0?0:&end, 0)==0)
				{
					r->apply(itin, no, dir, obsdir, uu, vv, 0, res);
					col+=res;
				}
				else if (scene<T>::verbose) 
				{
					mesh<T> *m = scene<T>::mm[scene<T>::meshid[shadow.getObj()]];
					printf("shadow from %d(%d) at %g,%g,%g, dist %g\n", shadow.getObj(), shadow.getFace(), shadow.getPoint().X(), shadow.getPoint().Y(), shadow.getPoint().Z(), shadow.distance());
					triangle3<T> tri = ip->getTriangle();
					printf("on triangle: (%g,%g,%g)-(%g,%g,%g)-(%g,%g,%g)\n", tri.getVertex(0).X(), tri.getVertex(0).Y(), tri.getVertex(0).Z(),
						tri.getVertex(1).X(), tri.getVertex(1).Y(), tri.getVertex(1).Z(),
						tri.getVertex(2).X(), tri.getVertex(2).Y(), tri.getVertex(2).Z() );
					unsigned int nf=m->nAdjacentFaces(shadow.getFace());
					printf("%d neighbors are: ",nf); 
					unsigned int *nn=m->getAdjacentFaces(shadow.getFace());
					for (int i=0; i<nf; i++) printf("%d ", nn[i]);
					printf("\n");
				}
				if (verb) scene<T>::verbose=true;
			}
		}
		for (i=0; i<nslights; i++)
		{
			for (j=0; j<nli*nli; j++)
			{
				float dx = float(isub%nrays)/(float)nrays+float(j%nli)/(float)nrays/(float)nli+(float)rand()/(float)RAND_MAX/(float)nrays/(float)nli;
				float dy = float(isub/nrays)/(float)nrays+float(j/nli)/(float)nrays/(float)nli+(float)rand()/(float)RAND_MAX/(float)nrays/(float)nli;
				vec3<T> ll=vec3<T>(slightsp[i].X()+T(dx)*slightsU[i].X()+T(dy)*slightsV[i].X(), slightsp[i].Y()+T(dx)*slightsU[i].Y()+T(dy)*slightsV[i].Y(), slightsp[i].Z()+T(dx)*slightsU[i].Z()+T(dy)*slightsV[i].Z());
				luminance itin=slightit[i];
				itin.scale(1.0/(T)(nli*nli));
				luminance res(0.0f);
				vec3<T> dir;
				dir.PVec(ip->getPoint(), ll);
				dir.normalize(dir.norm());
				line3<T> dd(ip->getPoint(), dir);
				scenepoint<T> shadow;
				scenepoint<T> end(ll);
				int count;

				float cc = dir.dot(no);
				if (scene<T>::verbose) { printf("ALight %d/%d (%g,%g,%g) Normal %g,%g,%g  dot: %g\n", i, j, itin.RED(), itin.GREEN(), itin.BLUE(), no.X(), no.Y(),no.Z(), cc);  }
				if (cc>0.0f)
				{
					itin.scale(cc);
					if (scene<T>::verbose)
					{
						srayorig[snbrays*6+0]=dd.origin().X(); srayorig[snbrays*6+1]=dd.origin().Y(); srayorig[snbrays*6+2]=dd.origin().Z();
						srayorig[snbrays*6+3]=ll.X(); srayorig[snbrays*6+4]=ll.Y(); srayorig[snbrays*6+5]=ll.Z();
						snbrays++;
					}
					//if (true)
					if (verb) scene<T>::verbose=false;
					if (visibility(scene<T>::bbhscene, dd, scene<T>::OCCLUDER, count, &shadow, 1, ip, scene<T>::lights[i].W()==0.0?0:&end, 0)==0)
					{
						r->apply(itin, no, dir, obsdir, uu, vv, 0, res);
						col+=res;
					}
					else if (scene<T>::verbose) 
					{
						mesh<T> *m = scene<T>::mm[scene<T>::meshid[shadow.getObj()]];
						printf("shadow from %d(%d) at %g,%g,%g, dist %g\n", shadow.getObj(), shadow.getFace(), shadow.getPoint().X(), shadow.getPoint().Y(), shadow.getPoint().Z(), shadow.distance());
						triangle3<T> tri = ip->getTriangle();
						printf("on triangle: (%g,%g,%g)-(%g,%g,%g)-(%g,%g,%g)\n", tri.getVertex(0).X(), tri.getVertex(0).Y(), tri.getVertex(0).Z(),
							tri.getVertex(1).X(), tri.getVertex(1).Y(), tri.getVertex(1).Z(),
							tri.getVertex(2).X(), tri.getVertex(2).Y(), tri.getVertex(2).Z() );
						unsigned int nf=m->nAdjacentFaces(shadow.getFace());
						printf("%d neighbors are: ",nf); 
						unsigned int *nn=m->getAdjacentFaces(shadow.getFace());
						for (int i=0; i<nf; i++) printf("%d ", nn[i]);
						printf("\n");
					}
					if (verb) scene<T>::verbose=true;
				}
			}
		}
		if (r->isMirror() && depth<MAX_BOUNCE)
		{
			vec3<T> rdd; rdd.reflection(obsdir, no);
			vec3<T> ruu; ruu.orthogonal(rdd); ruu.normalize(ruu.norm());
			vec3<T> rvv; rvv.cross(rdd,ruu); rvv.normalize(rvv.norm());
			if (depth==0) for (j=0; j<nspec*nspec; j++)
			{
				vec3<T> rdir; rdir.randomPolDir(j,nspec*nspec,0,1,rdd,ruu,rvv,0.0, mirr_angle);
				line3<T> ldd(ip->getPoint(), rdir);
				luminance res, refl=scene_stochastic::secondaryRay(ip, ldd, depth+1, isub);
				r->applyMirror(refl, 0, res);
				res.scale(1.0/T(nspec*nspec));
				if (scene<T>::verbose)
				{
					printf("reflexion %d: %g,%g,%g, col=%g,%g,%g\n", j, rdir.X(), rdir.Y(), rdir.Z(), res.RED(), res.GREEN(), res.BLUE());
				}
				col += res;
			}
			else
			{
				vec3<T> rdir; rdir.reflection(obsdir, no);
				line3<T> rdd(ip->getPoint(), rdir);
				luminance res, refl=scene_stochastic::secondaryRay(ip, rdd, depth+1, isub);
				r->applyMirror(refl, 0, res);
				col += res;
			}
		}
		if (r->isTransparent() && depth<MAX_BOUNCE)
		{
			vec3<T> rdir; 
			if (!rdir.refraction(obsdir, no, (double)r->getRefractionIndex()))
			{
				printf("error index\n");
				//rdir.orthogonalInPlan(no, obsdir);
				//line3<T> rdd(ip->getPoint(), rdir);
				//luminance res, refl=secondaryRay(ip, rdd, depth+1);
				//r->applyTransparency(refl, 0, res);
				//col += res;
			}
			else
			{
				line3<T> rdd(ip->getPoint(), rdir);
				scenepoint<T>  out;
				int objid = ip->getObj();
				mesh<T> *m=scene<T>::mm[scene<T>::meshid[objid]];
				transform<T> tr = scene<T>::pos[objid];
				int count;
				if (intersectRayObj(scene<T>::bbhobj[objid], objid, m, tr, rdd, scene<T>::CLOSEST, count, &out, 1, ip, 0, false, 0)>0)
				{
					rdir.reverse();
					vec3<T> noout = getOrientedNormal(rdir, &out);
					vec3<T> uuout; uuout.orthogonal(noout);
					vec3<T> vvout; vvout.cross(noout, uuout); vvout.normalize(vvout.norm());
					vec3<T> odir; 
					if (!odir.refraction(rdir, noout, 1.0/(double)r->getRefractionIndex()))
					{
						rdir=obsdir; rdir.reverse();
						line3<T> rdd(ip->getPoint(), rdir);
						luminance res, refl=scene_stochastic::secondaryRay(ip, rdd, depth+1, isub);
						r->applyTransparency(refl, 0, res);
						col += res;
						//col = luminance(0.0f);
					}
					else
					{
						line3<T> rddout(out.getPoint(), odir);
						luminance res, refl=scene_stochastic::secondaryRay(&out, rddout, depth+1, isub);
						r->applyTransparency(refl, 0, res);
						col += res;
					}
				}
				else
				{
					obsdir.reverse();
					line3<T> rdd(ip->getPoint(), obsdir);
					luminance res, refl=scene_stochastic::secondaryRay(ip, rdd, depth+1, isub);
					r->applyTransparency(refl, 0, res);
					col += res;
				}
			}
		}
		return col;
	}


};

#endif // !defined(AFX_LINE3_H__EC87CD56_6E08_4390_B876_CEC6D44EEA86__INCLUDED_)
