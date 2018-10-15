// BBH.h: interface for the bounding box-based space partitionning class (bounding box hierarchy).
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BBH_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_)
#define AFX_BBH_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "box3.h"

////////////////////////////////////////////
/*
bbh: a bounding box-based cluster consists of a space partition in dimension 3
using the principle of hierarchical bounding boxes.
It is an octree structure
The leafs of the tree contain unsigned int components indentifying on object index.
Each node contains a max of 8 child nodes
*/ 
////////////////////////////////////////////
////////////////////////////////////////////

////////////////////////////////////////////

////////////////////////////////////////////
template <class T> class bbh
{
protected:
	bbh<T>			*nodes[8];
	bbh<T>			*father;
	box3<T>			bbox;
	unsigned int	*leafs;
	unsigned int	nleafs;
public:
	bbh<T>() { father=0; nleafs =0; leafs=0; for(int i=0; i<8; i++) nodes[i]=0; } 
	bbh<T>(unsigned int nb, box3<T> *blist, unsigned int maxcomp)
	{
		father=0; nleafs =0; leafs=0; for(int i=0; i<8; i++) nodes[i]=0;
		unsigned int *index=new unsigned int [nb];
		for (unsigned int i=0; i<nb; i++) index[i]=i;
		this->build(nb,blist,index, 0,maxcomp);
		delete [] index;
	}

	 void build(unsigned int nb, box3<T> *blist, unsigned int *index, bbh<T> *ff, unsigned int maxcomp)
	 {
		 unsigned int i,j,lii;
		 box3<T> bf[8];
		 unsigned int nf[8];
		 unsigned int *lindex[8];
		 int lnindex[8];
		 unsigned int *_cl_iface= new unsigned int [nb];
		 for (i=0; i<nb; i++) _cl_iface[i]=-1;

		 father=ff;
		 bbox = blist[index[0]]; for (i=1; i<nb; i++) bbox.merge(bbox,blist[index[i]]);
		 if (nb<=maxcomp) 
		 {
			 nleafs=nb; leafs= new unsigned int [nb];
			 for (i=0; i<nb; i++) leafs[i]=index[i];
			 for (i=0; i<8; i++) nodes[i]=0;
			 return;
		 }
		nleafs = 0;
		leafs=0;
		bbox.splitt(bf);
		for (i=0; i<8; i++) nf[i]=0;
		for (i=0; i<nb; i++)
			{
			lii = index[i];
			vec3<T> bc = blist[lii].middle();
			for (j=0; j<8; j++)
				{
        			if (bf[j].isInside(bc,0)) { nf[j]++; _cl_iface[i]=(int)j; break; }
				}
			if (j==8)
				{
				printf("_build: into no box ??"); exit(0);
				}
			}
		int nrep=0;
		for (i=0; i<8; i++)
			{
				if (nf[i]>0) 
				{
					lindex[i]= new unsigned int [nf[i]];
					nrep++;
				}
				else lindex[i]=0;
				lnindex[i]=0;
			}
		if (nrep<2) 
		{
			for (i=0; i<8; i++) if (nf[i]>0) delete [] lindex[i];
			 nleafs=nb; leafs= new unsigned int [nb];
			 for (i=0; i<nb; i++) leafs[i]=index[i];
			 for (i=0; i<8; i++) nodes[i]=0;
			 return;
		}
		for (i=0; i<nb; i++)
			{
			j = _cl_iface[i];
			lii = index[i];
			lindex[j][lnindex[j]]=lii;
			lnindex[j]++;
			}
		for (i=0; i<8; i++)
			{
			if (lnindex[i]!=nf[i]) { printf("inconsistency error!\n"); exit(1); }
			if (lnindex[i]>0)
				{
					nodes[i]=new bbh<T>();
					nodes[i]->build(lnindex[i],blist,lindex[i],this,maxcomp);
					delete [] lindex[i];
				}
			else nodes[i]=0;
			}
		delete [] _cl_iface;
	 }

	 bbh<T> *getNode(int i) { return nodes[i]; }
	 unsigned int nLeafs() const { return nleafs; }
	 box3<T> getBox() const { return bbox; }
	 unsigned int getLeaf(int i) const { return leafs[i]; }
	 bbh<T> *getFather() { return father; }
};


#endif // !defined(AFX_CLUSTER_H__B6AC0A32_75EF_428E_BC10_6219F619FA29__INCLUDED_)
