// mesh.h: interface for the mesh class.
//
// defines a polyhedron which is a list of connected triangles
// By JMD 8/8/04
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MESH_H__A6721A04_8090_4DA7_A7BB_8AAB841F9C11__INCLUDED_)
#define AFX_MESH_H__A6721A04_8090_4DA7_A7BB_8AAB841F9C11__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string.h>
#include <math.h>
#include "gmath/vec3.h"
#include "gmath/intersection3.h"
#include "triangle3.h"

const unsigned int MAX_INCIDENT_CELLS=128;
const unsigned int MAX_NEIGHBOR_CELLS=64;
const unsigned int MAX_VERTEX_ATTRIBUTES=128;

// NOTE: mesh is composed of triangles only!
template <class T> class mesh 
{
protected:
	unsigned int ns, nf;		// ns=num vertices, nf = num faces
	unsigned int *lface;		// face is a triplet of vertex indices 
	T			 *vertices;		// list of vertices (triplets x,y,z)
	T			 *normals;		// list of normals (triplets x,y,z)

	int			nattr;			// how many attributes
	T			*attribute[MAX_VERTEX_ATTRIBUTES];  // attributes such as texture coordinates
	int			sizeattr[MAX_VERTEX_ATTRIBUTES];  // dimension of an attribute = 1, 2, ...

	vec3<T>		min,max;		// bounding box

	unsigned int	**incidentfaces;	// for each vertex list of incident faces
	unsigned int	*nincidentfaces;	// for each vertex amount of incident faces
	unsigned int	**neighborfaces;	// for each face, list of neighboring faces (sharing at least one common vertex)
	unsigned int	*nneighborfaces;

public:
	//////// constructors
	mesh<T>() { ns=0; nf=0; lface=0; vertices=0; normals=0; nattr=0; incidentfaces=0; nincidentfaces=0; neighborfaces=0; nneighborfaces=0; }
	mesh<T>(int n, int f, unsigned int *lf, T *vv, T *nn=0, T *att=0, int sizeatt=0)
	{
		ns=n; nf=f; 
		lface=lf;
		vertices=vv;
		normals =nn;
		nattr=0;
		if (att!=0 && sizeatt>0)
		{
			nattr=1; attribute[0]=att; sizeattr[0]=sizeatt;
		}
		this->computeTopology();
		this->orientAllFaces(0);
		if (nn==0) this->computePhongNormals();
		updateBox();
	}
	void pushAttribute(T *att, int asize) { attribute[nattr]=att; sizeattr[nattr]=asize; nattr++; } 
	void updateBox()
	{
		int i;
		for (i=0; i<ns; i++)
		{
			vec3<T> np(vertices[i*3],vertices[i*3+1],vertices[i*3+2]);
			if (i==0) { min=np; max=min; }
			else
			{
				min.keepMin(min, np);
				max.keepMax(max, np);
			}
		}
	}
	// some parametric objects
	void makeTorus(T gr, T pr, int na, int nb)
	{
		int i,j;
		ns=0;
		vertices = (T *)malloc(sizeof(T)*3*(na+1)*(nb+1));
		normals = (T *)malloc(sizeof(T)*3*(na+1)*(nb+1));
		attribute[0]= (T *)malloc(sizeof(T)*2*(na+1)*(nb+1));
		nattr=1; sizeattr[0]=2;

		for (i=0; i<=na; i++)
		{
			T alpha = (T)i*2.0*3.141592653589793/(T)na;
			for (j=0; j<=nb; j++)
			{
				T beta = (T)j*2.0*3.141592653589793/(T)nb;
				vertices[ns*3]=(pr*cos(beta)+gr)*cos(alpha);
				vertices[ns*3+1]=(pr*cos(beta)+gr)*sin(alpha);
				vertices[ns*3+2]=pr*sin(beta);
				normals[ns*3]=cos(beta)*cos(alpha);
				normals[ns*3+1]=cos(beta)*sin(alpha);
				normals[ns*3+2]=sin(beta);
				attribute[0][ns*2]=
				attribute[0][ns*2+1]=
				ns++;
			}
		}
		lface = (unsigned int *)malloc(sizeof(unsigned int)*3*2*na*nb);
		nf = 0;
		for (i=0; i<na; i++)
		{
			for (j=0; j<nb; j++)
			{
				lface[nf*3*2]=i*(nb)+j; lface[nf*3*2+1]=(i+1<na?i+1:0)*(nb)+j; lface[nf*3*2+2]=i*(nb)+(j+1<nb?j+1:0);
				lface[nf*3*2+3]=(i+1<na?i+1:0)*(nb)+j; lface[nf*3*2+4]=i*(nb)+(j+1<nb?j+1:0); lface[nf*3*2+5]=(i+1<na?i+1:0)*(nb)+(j+1<nb?j+1:0);
				nf++;
			}
		}
		nf *= 2;
		this->computeTopology();
		this->orientAllFaces(0);
		//this->computePhongNormals();
		updateBox();
	}

	unsigned int nFaces() const { return nf; }
	unsigned int nVertices() const { return ns; }
	T *getVerticesData() { return vertices; }
	T *getNormalsData() { return normals; }
	unsigned int *getFacesData() { return lface; }
	vec3<T> getMin() const { return min; }
	vec3<T> getMax() const { return max; }

	void getFaceVertexId(unsigned int iface, unsigned int &a, unsigned int &b, unsigned int &c) const { a=lface[3*iface]; b=lface[3*iface+1]; c=lface[3*iface+2]; } 
	unsigned int getFaceVertexId(unsigned int iface, unsigned int i) const { return lface[3*iface+i]; } 
	int faceContainsVertex(unsigned int iface, int vertexid) const // -1 means not in cell
	{
		if (lface[3*iface]==vertexid) return 0;
		if (lface[3*iface+1]==vertexid) return 1;
		if (lface[3*iface+2]==vertexid) return 2;
		return -1;
	}
	vec3<T> getVertex(unsigned int id) const
	{
		return vec3<T>(vertices[id*3],vertices[id*3+1],vertices[id*3+2]);
	}
	vec3<T> getFaceVertex(unsigned int iface, unsigned int i) const // i in [0,2], because only triangles
	{
		unsigned int id = getFaceVertexId(iface,i);
		return vec3<T>(vertices[id*3],vertices[id*3+1],vertices[id*3+2]);
	}
	vec3<T> getFaceBary(unsigned int iface) const
	{
		vec3<T> vv = getFaceVertex(iface, 0);
		vv += getFaceVertex(iface, 1);
		vv += getFaceVertex(iface, 2);
		vv.scale(T(1.0/3.0));
		return vv;
	}
	// computes the normal of a cell by using cross product
	vec3<T> normalFace(int faceid) const
	{
		vec3<T> uu, vv;
		uu.PVec(getFaceVertex(faceid, 0), getFaceVertex(faceid, 1));
		vv.PVec(getFaceVertex(faceid, 0), getFaceVertex(faceid, 2));
		vec3<T> res; res.cross(uu,vv);
		T no = res.norm();
		if (no==T(0)) return vec3<T>();
		res.normalize(no);
		return res;
	}
	vec3<T> getFaceVertexNormal(unsigned int iface, unsigned int i) const // i in [0,2]
	{
		unsigned int id = getFaceVertexId(iface,i);
		return vec3<T>(normals[id*3],normals[id*3+1],normals[id*3+2]);
	}

	// Topological operations
	int nIncidentFaces(unsigned int vertexid)  
	{ 
		return nincidentfaces[vertexid]; 
	}
	unsigned int *getIncidentFaces(unsigned int vertexid) 
	{
		return incidentfaces[vertexid];
	}
	// two faces are adjacent if sharing at least one same vertex 
	unsigned int nAdjacentFaces(unsigned int faceid)  
	{ 
		return nneighborfaces[faceid];
	}
	unsigned int *getAdjacentFaces(unsigned int faceid)  
	{ 
		return neighborfaces[faceid];
	}
	bool areConnected(int faceid1, int faceid2) const
	{
		unsigned int iv1;
		for (iv1=0; iv1<3; iv1++)
		{
			unsigned int v1 = getFaceVertexId(faceid1, iv1);
			int iv2 = faceContainsVertex(faceid2, v1);
			if (iv2!=-1) return true;
		}
		return false;
	}
	bool areSharingEdge(unsigned int faceid1, unsigned int faceid2) const
	{
		unsigned int iv1;
		for (iv1=0; iv1<3; iv1++)
		{
			unsigned int v1 = getFaceVertexId(faceid1, iv1);
			int id = faceContainsVertex(faceid2, v1);
			unsigned int iv2 = (unsigned int)id;
			if (id!=-1)
			{
				if (getFaceVertexId(faceid1, iv1==0?2:iv1-1)==getFaceVertexId(faceid2, iv2==0?2:iv2-1)) return true;	
				if (getFaceVertexId(faceid1, iv1==2?0:iv1+1)==getFaceVertexId(faceid2, iv2==0?2:iv2-1)) return true;	
				if (getFaceVertexId(faceid1, iv1==2?0:iv1+1)==getFaceVertexId(faceid2, iv2==2?0:iv2+1)) return true;	
				if (getFaceVertexId(faceid1, iv1==0?2:iv1-1)==getFaceVertexId(faceid2, iv2==2?0:iv2+1)) return true;	
			}
		}
		return false;
	}
	bool getSharingEdge(unsigned int faceid1, unsigned int faceid2, unsigned int &vertexid1, unsigned int &vertexid2) const
	{
		unsigned int iv1;
		for (iv1=0; iv1<3; iv1++)
		{
			unsigned int v1 = getFaceVertexId(faceid1, iv1);
			int id = faceContainsVertex(faceid2, v1);
			unsigned int iv2 = (unsigned int)id;
			if (id!=-1)
			{
				if (getFaceVertexId(faceid1, iv1==0?2:iv1-1)==getFaceVertexId(faceid2, iv2==0?2:iv2-1)) {  vertexid1=(iv1==0?2:iv1-1); vertexid2=(iv2==0?2:iv2-1); return true; }	
				if (getFaceVertexId(faceid1, iv1==2?0:iv1+1)==getFaceVertexId(faceid2, iv2==0?2:iv2-1)) { vertexid1=iv1;  vertexid2=(iv2==0?2:iv2-1); return true; }	
				if (getFaceVertexId(faceid1, iv1==2?0:iv1+1)==getFaceVertexId(faceid2, iv2==2?0:iv2+1)) {  vertexid1=iv1; vertexid2=iv2; return true; }
				if (getFaceVertexId(faceid1, iv1==0?2:iv1-1)==getFaceVertexId(faceid2, iv2==2?0:iv2+1)) {  vertexid1=(iv1==0?2:iv1-1); vertexid2=iv2; return true; }	
			}
		}
		vertexid1=0; vertexid2=0;
		return false;
	}
	bool getVertexSharingEdge(unsigned int faceid1, unsigned int faceid2, unsigned int v1, unsigned int &v2) const
	{
		int id1 = faceContainsVertex(faceid1, v1);
		int id2 = faceContainsVertex(faceid2, v1);
		if (id1!=-1 && id2!=-1)
		{
			unsigned int iv1=(unsigned int)id1;
			unsigned int iv2=(unsigned int)id2;
			if (getFaceVertexId(faceid1, iv1==0?2:iv1-1)==getFaceVertexId(faceid2, iv2==0?2:iv2-1)) {  v2=getFaceVertexId(faceid1, iv1==0?2:iv1-1); return true; }	
			if (getFaceVertexId(faceid1, iv1==2?0:iv1+1)==getFaceVertexId(faceid2, iv2==0?2:iv2-1)) {  v2=getFaceVertexId(faceid1, iv1==2?0:iv1+1); return true;}	
			if (getFaceVertexId(faceid1, iv1==2?0:iv1+1)==getFaceVertexId(faceid2, iv2==2?0:iv2+1)) {  v2=getFaceVertexId(faceid1, iv1==2?0:iv1+1); return true; }
			if (getFaceVertexId(faceid1, iv1==0?2:iv1-1)==getFaceVertexId(faceid2, iv2==2?0:iv2+1)) {  v2=getFaceVertexId(faceid1, iv1==0?2:iv1-1); return true;}	
		}
		return false;
	}
	// orient face according to faceref, the faces must share an edge
	void reverseFaceOrientation(unsigned int faceid)
	{
		unsigned int swap=lface[faceid*3];
		lface[faceid*3]=lface[faceid*3+2];
		lface[faceid*3+2]=swap;	
	}

	bool orientFace(unsigned int faceref, unsigned int faceid, unsigned int iv1, unsigned int iv2)
	{
		if( getFaceVertexId(faceref,iv1)==getFaceVertexId(faceid, iv2) )
		{
			reverseFaceOrientation(faceid);
			return true;
		}
		return false;
	}
	bool areOrientedFaces(unsigned int faceref, unsigned int faceid) const
	{
		unsigned int iv1, iv2;
		if (getSharingEdge(faceref, faceid, iv1, iv2))
		{
			return getFaceVertexId(faceref,iv1)!=getFaceVertexId(faceid, iv2);
		}
		return false;
	}
	void _orientAllCells(bool *yet, unsigned int *queue, unsigned int &qlen)
	{
		int i;
		while (qlen>0) 
		{
			// pop front of queue
			unsigned int faceref = queue[0];
			for (i=0; i<(int)qlen-1; i++) queue[i]=queue[i+1];
			qlen--;
			if (neighborfaces[faceref]!=0)
			{
				for (i=0; i<(int)nneighborfaces[faceref]; i++)
				{
					unsigned int id = neighborfaces[faceref][i];
					if (!yet[id])
					{
						unsigned int iv1,iv2;
						if (getSharingEdge(faceref, id, iv1, iv2))
						{
							orientFace(faceref, id, iv1, iv2);
							yet[id]= true;
							// push at back of queue
							queue[qlen]=id; qlen++;
						}
					}
				}
			}
		}
	}
	int orientAllFaces(unsigned int faceref)
	{
		int i, n=1; unsigned int curr=faceref;
		bool *yet= (bool *)malloc(sizeof(bool)*nFaces());
		for (i=0; i<(int)nFaces(); i++) yet[i]=false;
		unsigned int *queue=(unsigned int *)malloc(sizeof(unsigned int)*nFaces());
		unsigned int qlen=0;
		bool cont=true;
		do {
			qlen=1;
			queue[0]=curr;
			yet[curr]=true;
			_orientAllCells(yet, queue, qlen);
			for (i=0; i<(int)nFaces(); i++) if (!yet[i]) break;
			if (i<(int)nFaces()) { n++; curr=(unsigned int)i; } else { cont=false; }
		} while(cont);
		return n;
	}
	bool computeTopology()
	{
		unsigned int i,j,k;
		//printf("Computing topology...\n");
		incidentfaces=(unsigned int **)malloc(sizeof(unsigned int *)*nVertices());
		nincidentfaces=(unsigned int *)malloc(sizeof(unsigned int)*nVertices());
		for (i=0; i<nVertices(); i++) { incidentfaces[i]=0; nincidentfaces[i]=0; }
		neighborfaces=(unsigned int **)malloc(sizeof(unsigned int *)*nFaces());
		nneighborfaces=(unsigned int *)malloc(sizeof(unsigned int)*nFaces());
		for (i=0; i<nFaces(); i++) { neighborfaces[i]=0; nneighborfaces[i]=0; }
		for (i=0; i<nFaces(); i++)
		{
			for (j=0; j<3; j++)
			{
				unsigned int id = getFaceVertexId(i,j);
				unsigned int *infa = incidentfaces[id];
				if (infa==0) { infa =(unsigned int *)malloc(sizeof(unsigned int)*MAX_INCIDENT_CELLS); incidentfaces[id]=infa; }
				if (nincidentfaces[id]==MAX_INCIDENT_CELLS) { printf("cannot compute topology!Too many incident faces on vertex %d\n",id); return false; }
				infa[nincidentfaces[id]]=i;
				nincidentfaces[id]++;
			}
		}
		for (i=0; i<nFaces(); i++)
		{
			for (j=0; j<3; j++)
			{
				unsigned int id = getFaceVertexId(i,j);
				unsigned int *infa = incidentfaces[id];
				for (k=0; k<nincidentfaces[id]; k++)
				{
					if (infa[k]!=i)
					{
						if (neighborfaces[i]==0) neighborfaces[i]=(unsigned int *)malloc(sizeof(unsigned int)*MAX_NEIGHBOR_CELLS);
						if (nneighborfaces[i]==MAX_NEIGHBOR_CELLS){ printf("cannot compute topology!Too many adjacent faces on face %d\n",i); return false; }
						unsigned int *nface = neighborfaces[i];
						bool deja=false;
						for (int ii=0; ii<(int)nneighborfaces[i] && !deja; ii++) if (nface[ii]==infa[k]) deja=true;
						if (!deja) { nface[nneighborfaces[i]]=infa[k]; nneighborfaces[i]++; }
					}
				}
			}
		}
		return true;
	}
	// Computes the normals on all vertices using the Phong principle (average normal for all incident cells)
	void computePhongNormals()
	{
		//printf("Computing phong normals...\n");
		normals = (T *)malloc(sizeof(T)*3*ns);
		for (unsigned int i=0; i<nVertices(); i++) 
		{
			vec3<T> nn = computePhongNormalVertex(i);
			normals[i*3]=nn.X(); normals[i*3+1]=nn.Y(); normals[i*3+2]=nn.Z();
		}
	}
	// Computes the average normal on a single vertex for all incident cells
	vec3<T> computePhongNormalVertex(unsigned int vid)
	{
		int i;
		bool yet[MAX_INCIDENT_CELLS];
		vec3<T> no[MAX_INCIDENT_CELLS], dir, bc1, bc2, dd1, dd2, no1, no2;

		unsigned int *lfa = getIncidentFaces(vid);
		unsigned int nlfa = nIncidentFaces(vid);
		if (nlfa==0) { printf("isolated vertex %d!\n", vid); return vec3<T>(); }
		for (i=0;  i<(int)nlfa; i++) yet[i]=false;
		unsigned int istart = lfa[0];
		no[0]=normalFace(istart);
		unsigned int pifa=istart;
		yet[0]=true;
		bool second=false, error=false, cont;
		unsigned int oldi=0, ifa=0, nis=0;
		if (nlfa>1)
		{
		do {
			cont=true;
			for (i=1; i<(int)nlfa; i++)
				{
				ifa = lfa[i];
				if (!yet[i] && getVertexSharingEdge(ifa, pifa, vid, nis) ) break;
				}
			if (i<(int)nlfa)
				{
				no[i]=normalFace(ifa);
				dir.PVec(getVertex(vid),getVertex(nis));
				bc1 = getFaceBary(ifa);
				bc2 = getFaceBary(pifa);
				dd1.PVec(getVertex(vid), bc1);
				dd2.PVec(getVertex(vid), bc2);
				no1.cross(dir, dd1);
				no1.normalize(no1.norm());
				no2.cross(dd2, dir);
				no2.normalize(no2.norm());
				if (no2.dot(no[oldi])>T(0))
					{
						if (no1.dot(no[i])<T(0)) { printf("reversing on vertex %d\n", vid); no[i].reverse(); }
					}
				else
					{
						if (no1.dot(no[i])>T(0)) { printf("reversing on vertex %d\n", vid); no[i].reverse(); }
					}
				pifa=ifa;
				oldi=i;
				yet[i]=true;
				cont=false;
				for (i=0; i<(int)nlfa; i++) if (!yet[i]) { cont=true; break; }
				}
			else
				{
				if (second)
					{
					cont=false;
					error=true;
					}
				else
					{
					second=true;
					pifa=istart;
					oldi=0;
					}
				}
			} while(cont);
		}
		if (error) 
			{
			printf("\n\nwarning: vertex %d is not locally orientable.\n", vid);
			}
		vec3<T> n;
		int nb=0;
		for (i=0; i<(int)nlfa; i++) if (yet[i]) { n+=no[i]; nb++; }
		n.scale(T(1)/T(nb));
		T norm=n.norm();
		if (norm==T(0)) { printf("warning : sum is null in Shaded normal!"); return no[0]; }
		n.normalize(norm);
		return n;
	}

	bool intersectRay(const line3<T> &dd, intersection3<T> &itmin, int &iface, triangle3<T> &trimin, 
		T err, bool withdmin, T distmin)
	{
		intersection3<T> it;
		int j;
		bool first=true;
		if (it.intersectRay(dd,getMin(), getMax())>0)
		{
			if (!withdmin || (withdmin && it.distance(0)<distmin))
			for (j=0; j<(int)nFaces(); j++)
			{
				vec3<T> a,b,c;
				a = getFaceVertex(j,0);
				b = getFaceVertex(j,1);
				c = getFaceVertex(j,2);
				triangle3<float> tri(a,b,c);
				if (tri.intersectRay(dd,err,it)>0)
				{
					if (!withdmin || (withdmin && it.distance(0)<distmin))
					{
						if (first) { itmin=it; iface=j; trimin=tri; first=false; }
						else
						{
							if (it.distance(0)<itmin.distance(0)) { itmin=it; iface=j; trimin=tri; }
						}
					}
				}
			}
		}	
		return !first;
	}
};


#endif // !defined(AFX_VEC3_H__A6721A04_8090_4DA7_A7BB_8AAB841F9C11__INCLUDED_)
