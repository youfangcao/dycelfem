// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     delaunay.cpp        ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef DELAUNAY_CPP
#define DELAUNAY_CPP

#include <iostream>
#include <vector>
#include <cmath>
#include "delaunay.h"
#include "util.h"

using namespace std;

CXYZW::~CXYZW()
{
	/*S0 = NULL;
	p = NULL;*/
}

CEdge::~CEdge()
{
	/*A = NULL;
	B = NULL;
	Cell = NULL;
	Rev = NULL;
	Enx = NULL;
	Fnx = NULL;
	S1 = NULL;
	S2 = NULL;*/
}

void CEdge::set_center()
{
	double xi = A->X;
	double yi = A->Y;
	double xj = B->X;
	double yj = B->Y;
	double wi = A->W;
	double wj = B->W;

	double aux = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj);
	double i4 = wi*wi;
	double j4 = wj*wj;
	double mu = 0.5 - 0.5*(i4-j4)/aux;
	center[0] = mu*xi + (1-mu)*xj;
	center[1] = mu*yi + (1-mu)*yj;
}

CCell::~CCell()
{
	/*VV[0] = NULL;
	VV[1] = NULL;
	VV[2] = NULL;
	ab = NULL;
	bc = NULL;
	ca = NULL;
	CC[0] = NULL;
	CC[1] = NULL;
	CC[2] = NULL;*/
}

void CCell::set_center()
{
	double xi = VV[0]->X;
	double yi = VV[0]->Y;
	double xj = VV[1]->X;
	double yj = VV[1]->Y;
	double xk = VV[2]->X;
	double yk = VV[2]->Y;

	double i0 = 0.5*(VV[0]->W*VV[0]->W - xi*xi - yi*yi);
	double j0 = 0.5*(VV[1]->W*VV[1]->W - xj*xj - yj*yj);
	double k0 = 0.5*(VV[2]->W*VV[2]->W - xk*xk - yk*yk);

	double D0 = (xj*yk - yj*xk) - (xi*yk - yi*xk) + (xi*yj - yi*xj);
	double Dx = (j0*yk - yj*k0) - (i0*yk - yi*k0) + (i0*yj - yi*j0);
	double Dy = (xj*k0 - j0*xk) - (xi*k0 - i0*xk) + (xi*j0 - i0*xj);

	center[0] = -Dx/D0;
	center[1] = -Dy/D0;

	alphavalue = (center[0]-xi)*(center[0]-xi) + (center[1]-yi)*(center[1]-yi);
}

bool CCell::inside_complex()
{
	bool valid = false;
	if (VV[0]->P()->get_cell_boundary() &&
		VV[1]->P()->get_cell_boundary() &&
		VV[2]->P()->get_cell_boundary())
	{
		point3D *P1 = VV[0]->P();
		point3D *P2 = VV[1]->P();
		point3D *P3 = VV[2]->P();
		double x12 = 0.5*(P1->x() + P2->x());
		double y12 = 0.5*(P1->y() + P2->y());
		double x23 = 0.5*(P2->x() + P3->x());
		double y23 = 0.5*(P2->y() + P3->y());
		double x13 = 0.5*(P1->x() + P3->x());
		double y13 = 0.5*(P1->y() + P3->y());
		if (!P1->getab()->judge_vert_out_side(x23,y23) && 
			!P2->getab()->judge_vert_out_side(x13,y13) &&
			!P3->getab()->judge_vert_out_side(x12,y12))
		{
			valid = true;
		}
	}
	else
	{
		valid = true;
	}
	return valid;
}

CD3DW::CD3DW(double range[4])
{
	double xmin = range[0];
	double ymin = range[1];
	double xmax = range[2];
	double ymax = range[3];
	CXYZW *a = new CXYZW();
	CXYZW *b = new CXYZW();
	CXYZW *c = new CXYZW();
	V0 = a;
	V1 = b;
	V2 = c;
	V0->set(1.5*xmin-0.5*xmax,ymin,10);
	V1->set(0.5*xmin+0.5*xmax,2*ymax-ymin,10);
	V2->set(1.5*xmax-0.5*xmin,ymin,10);
	V0->N=-1;V1->N=-2;V2->N=-3;
	
	CEdge *ab,*bc,*ca,*ba,*cb,*ac;
	ab = new CEdge(V0,V1);Es.push_back(ab);//ab->~CEdge();
	ba = new CEdge(V1,V0);Es.push_back(ba);//ba->~CEdge();
	Es[0]->setA(V0);Es[0]->setB(V1);Es[0]->Rev = Es[1];
	Es[1]->Rev = Es[0];
	Es[0]->set_center();
	Es[0]->alphavalue=-10;
	Es[1]->set_center();
	Es[1]->alphavalue=-10;
	bc = new CEdge(V1,V2);Es.push_back(bc);//bc->~CEdge();
	cb = new CEdge(V2,V1);Es.push_back(cb);//cb->~CEdge();
	Es[2]->Rev = Es[3];
	Es[3]->Rev = Es[2];
	Es[2]->set_center();
	Es[2]->alphavalue=-10;
	Es[3]->set_center();
	Es[3]->alphavalue=-10;
	ca = new CEdge(V2,V0);Es.push_back(ca);//ca->~CEdge();
	ac = new CEdge(V0,V2);Es.push_back(ac);//ac->~CEdge();
	Es[4]->Rev = Es[5];
	Es[5]->Rev = Es[4];
	Es[4]->set_center();
	Es[4]->alphavalue=-10;
	Es[5]->set_center();
	Es[5]->alphavalue=-10;

	Es[0]->Enx=Es[2];Es[2]->Enx=Es[4];Es[4]->Enx=Es[0];Es[1]->Enx=Es[5];Es[5]->Enx=Es[3];Es[3]->Enx=Es[1];
	Es[0]->Fnx=Es[5];Es[5]->Fnx=Es[0];Es[2]->Fnx=Es[1];Es[1]->Fnx=Es[2];Es[4]->Fnx=Es[3];Es[3]->Fnx=Es[4];

	CCell *C = new CCell;
	Cs.push_back(C);
	Cs[0]->set_BN(1);
	Cs[0]->VV[0]=V0;
	Cs[0]->VV[1]=V1;
	Cs[0]->VV[2]=V2;
	Cs[0]->ab=Es[0];
	Cs[0]->bc=Es[2];
	Cs[0]->ca=Es[4];
	Cs[0]->set_center(); //C->alphavalue=-10;

	Es[0]->Cell=Cs[0];Es[2]->Cell=Cs[0];Es[4]->Cell=Cs[0];
	Es[1]->Cell=Cs[0];Es[3]->Cell=Cs[0];Es[5]->Cell=Cs[0];
	Es[0]->Flag=2;Es[2]->Flag=2;Es[4]->Flag=2;
	Es[1]->Flag=2;Es[3]->Flag=2;Es[5]->Flag=2;
}

CD3DW::~CD3DW()
{
	int i;
	for (i=0;i<(int)Vs.size();i++) {delete Vs[i];Vs[i] = NULL;}
	for (i=0;i<(int)Es.size();i++) {delete Es[i];Es[i] = NULL;}
	for (i=0;i<(int)Cs.size();i++) {delete Cs[i];Cs[i] = NULL;}
	Vs.clear();
	Es.clear();
	Cs.clear();
	delete V0;delete V1;delete V2;
	V0 = NULL;V1 = NULL;V2 = NULL;
}

void CD3DW::add(double x, double y, double w, int n, void *l)
{
	CXYZW *d = new CXYZW; 
	d->Label=l;
	d->set(x,y,w);
	d->N=n;
	Vs.push_back(d);

	CEdge *ab;
	ab = locate(d);
	if (!judge_on_edge(ab->Cell,d)) 
	{
		CXYZW *a,*b,*c;
		CEdge *bc,*ca;
		bc = ab->Enx;
		ca = bc->Enx;
		CEdge *ba,*cb,*ac;
		CCell *abc,*abd,*cad,*bcd;
		ba=ab->Rev;cb=bc->Rev;ac=ca->Rev;
		a=ab->A;b=ab->B;c=ab->Enx->B;

		CEdge *ad,*da,*bd,*db,*cd,*dc;
		ad = new CEdge(a,d);Es.push_back(ad);
		da = new CEdge(d,a);Es.push_back(da);
		ad->Rev=da;da->Rev=ad;
		bd = new CEdge(b,d);Es.push_back(bd);
		db = new CEdge(d,b);Es.push_back(db);
		bd->Rev=db;db->Rev=bd;
		cd = new CEdge(c,d);Es.push_back(cd);
		dc = new CEdge(d,c);Es.push_back(dc);
		cd->Rev=dc;dc->Rev=cd;

		ad->set_center();ad->alphavalue=-10;
		da->set_center();da->alphavalue=-10;
		bd->set_center();bd->alphavalue=-10;
		db->set_center();db->alphavalue=-10;
		cd->set_center();cd->alphavalue=-10;
		dc->set_center();dc->alphavalue=-10;

		ab->Enx=bd;bd->Enx=da;da->Enx=ab;
		bc->Enx=cd;cd->Enx=db;db->Enx=bc;
		ca->Enx=ad;ad->Enx=dc;dc->Enx=ca;

		ad->Fnx=ab->Fnx;ab->Fnx=ad;
		bd->Fnx=bc->Fnx;bc->Fnx=bd;
		cd->Fnx=ca->Fnx;ca->Fnx=cd;

		da->Fnx=db;db->Fnx=dc;dc->Fnx=da;

		abc=ab->Cell;
		abd = new CCell;
		Cs.push_back(abd);
		abd->ab=ab;abd->bc=bd;abd->ca=da;
		abd->VV[0]=a;abd->VV[1]=b;abd->VV[2]=d;
		abd->set_center();abd->alphavalue=-10;

		bcd = new CCell;
		Cs.push_back(bcd);
		bcd->ab=bc;bcd->bc=cd;bcd->ca=db;
		bcd->VV[0]=b;bcd->VV[1]=c;bcd->VV[2]=d;
		bcd->set_center();bcd->alphavalue=-10;

		cad = new CCell;
		Cs.push_back(cad);
		cad->ab=ca;cad->bc=ad;cad->ca=dc;
		cad->VV[0]=c;cad->VV[1]=a;cad->VV[2]=d;
		cad->set_center();cad->alphavalue=-10;

		abc->CC[0] = abd;abc->CC[1] = bcd;abc->CC[2] = cad;

		ab->Cell = abd;
		bc->Cell = bcd;
		ca->Cell = cad;

		ad->Cell = cad;da->Cell = abd;
		bd->Cell = abd;db->Cell = bcd;
		cd->Cell = bcd;dc->Cell = cad;

		if (ab->A->N>=0 || ab->B->N>=0) {St.push_back(ab);ab->Flag=1;}
		if (bc->A->N>=0 || bc->B->N>=0) {St.push_back(bc);bc->Flag=1;}
		if (ca->A->N>=0 || ca->B->N>=0) {St.push_back(ca);ca->Flag=1;}

		while (St.size())
		{
			ab = St[St.size()-1]; St.pop_back(); ab->Flag&=2;
			if (ab->Flag) {continue;}
			a = ab->A;b=ab->B;c=ab->Rev->Fnx->B;
			if (insphere(ab->Rev->Cell,d)) 
			{
				flip(ab);
				continue;
			}
		}
	}
	else
	{
		CXYZW *a,*b,*c,*e;
		ab = get_on_edge(ab->Cell,d);
		a = ab->A;
		b = ab->B;
		c=ab->Enx->B;
		e = ab->Rev->Enx->B;
		CEdge *ba,*bc,*cb,*ca,*ac,*ae,*ea,*eb,*be;
		CEdge *ad,*db,*dc,*da,*bd,*cd,*de,*ed;
		CCell *abc,*bae;
		CCell *adc,*dbc,*bde,*dae;

		ad = new CEdge(a,d);Es.push_back(ad);
		da = new CEdge(d,a);Es.push_back(da);
		ad->Rev=da;da->Rev=ad;
		db = new CEdge(d,b);Es.push_back(db);
		bd = new CEdge(b,d);Es.push_back(bd);
		db->Rev=bd;bd->Rev=db;
		dc = new CEdge(d,c);Es.push_back(dc);
		cd = new CEdge(c,d);Es.push_back(cd);
		dc->Rev=cd;cd->Rev=dc;
		de = new CEdge(d,e);Es.push_back(de);
		ed = new CEdge(e,d);Es.push_back(ed);
		de->Rev=ed;ed->Rev=de;

		ba=ab->Rev;
		bc=ab->Enx;cb=bc->Rev;
		ca=bc->Enx;ac=ca->Rev;
		ae=ba->Enx;ea=ae->Rev;
		eb=ae->Enx;be=eb->Rev;

		abc=ab->Cell;bae=ba->Cell;

		ad->set_center();ad->alphavalue=-10;
		da->set_center();da->alphavalue=-10;
		db->set_center();db->alphavalue=-10;
		bd->set_center();bd->alphavalue=-10;
		dc->set_center();dc->alphavalue=-10;
		cd->set_center();cd->alphavalue=-10;
		de->set_center();de->alphavalue=-10;
		ed->set_center();ed->alphavalue=-10;

		ad->Enx=dc;dc->Enx=ca;ca->Enx=ad;
		db->Enx=bc;bc->Enx=cd;cd->Enx=db;
		bd->Enx=de;de->Enx=eb;eb->Enx=bd;
		da->Enx=ae;ae->Enx=ed;ed->Enx=da;

		ad->Fnx=ac;ca->Fnx=cd;
		db->Fnx=dc;cd->Fnx=cb;
		bc->Fnx=bd;dc->Fnx=da;

		bd->Fnx=be;eb->Fnx=ed;
		da->Fnx=de;ae->Fnx=ad;
		de->Fnx=db;ed->Fnx=ea;

		adc = new CCell;
		Cs.push_back(adc);
		adc->ab=ad;adc->bc=dc;adc->ca=ca;
		adc->VV[0]=a;adc->VV[1]=d;adc->VV[2]=c;
		adc->set_center();adc->alphavalue=-10;

		dbc = new CCell;
		Cs.push_back(dbc);
		dbc->ab=db;dbc->bc=bc;dbc->ca=cd;
		dbc->VV[0]=d;dbc->VV[1]=b;dbc->VV[2]=c;
		dbc->set_center();dbc->alphavalue=-10;

		bde = new CCell;
		Cs.push_back(bde);
		bde->ab=bd;bde->bc=de;bde->ca=eb;
		bde->VV[0]=b;bde->VV[1]=d;bde->VV[2]=e;
		bde->set_center();bde->alphavalue=-10;

		dae = new CCell;
		Cs.push_back(dae);
		dae->ab=da;dae->bc=ae;dae->ca=ed;
		dae->VV[0]=d;dae->VV[1]=a;dae->VV[2]=e;
		dae->set_center();dae->alphavalue=-10;

		abc->CC[0]=adc;abc->CC[1]=dbc;
		bae->CC[0]=bde;bae->CC[1]=dae;

		ad->Cell=adc;dc->Cell=adc;ca->Cell=adc;
		db->Cell=dbc;bc->Cell=dbc;cd->Cell=dbc;
		bd->Cell=bde;de->Cell=bde;eb->Cell=bde;
		da->Cell=dae;ae->Cell=dae;ed->Cell=dae;

		if (bc->A->N>=0 || bc->B->N>=0) {St.push_back(bc);bc->Flag=1;}
		if (ca->A->N>=0 || ca->B->N>=0) {St.push_back(ca);ca->Flag=1;}
		if (ae->A->N>=0 || ae->B->N>=0) {St.push_back(ae);ae->Flag=1;}
		if (eb->A->N>=0 || eb->B->N>=0) {St.push_back(eb);eb->Flag=1;}

		while (St.size())
		{
			ab = St[St.size()-1]; St.pop_back(); ab->Flag&=2;
			if (ab->Flag) {continue;}
			a = ab->A;b=ab->B;c=ab->Rev->Fnx->B;
			if (insphere(ab->Rev->Cell,d)) 
			{
				flip(ab);
				continue;
			}
		}
	}
}

void CD3DW::add_resample(double x, double y, double w, int n, point3D *l)
{
	CXYZW *d;
	
	d = new CXYZW; 
	d->set_p(l);
	d->set(x,y,w);
	d->N=n;
	Vs.push_back(d);

	CEdge *ab;
	ab = locate(d);
	if (!judge_on_edge(ab->Cell,d)) 
	{
		CXYZW *a,*b,*c;
		CEdge *bc,*ca;
		bc = ab->Enx;ca = bc->Enx;
		CEdge *ba,*cb,*ac;
		CCell *abc,*abd,*cad,*bcd;
		ba=ab->Rev;cb=bc->Rev;ac=ca->Rev;
		a=ab->A;b=ab->B;c=ab->Enx->B;
		
		CEdge *ad,*da,*bd,*db,*cd,*dc;
		ad = new CEdge(a,d);Es.push_back(ad);
		da = new CEdge(d,a);Es.push_back(da);
		ad->Rev=da;da->Rev=ad;
		bd = new CEdge(b,d);Es.push_back(bd);
		db = new CEdge(d,b);Es.push_back(db);
		bd->Rev=db;db->Rev=bd;
		cd = new CEdge(c,d);Es.push_back(cd);
		dc = new CEdge(d,c);Es.push_back(dc);
		cd->Rev=dc;dc->Rev=cd;

		ad->set_center();ad->alphavalue=-10;
		da->set_center();da->alphavalue=-10;
		bd->set_center();bd->alphavalue=-10;
		db->set_center();db->alphavalue=-10;
		cd->set_center();cd->alphavalue=-10;
		dc->set_center();dc->alphavalue=-10;

		ab->Enx=bd;bd->Enx=da;da->Enx=ab;
		bc->Enx=cd;cd->Enx=db;db->Enx=bc;
		ca->Enx=ad;ad->Enx=dc;dc->Enx=ca;

		ad->Fnx=ab->Fnx;ab->Fnx=ad;
		bd->Fnx=bc->Fnx;bc->Fnx=bd;
		cd->Fnx=ca->Fnx;ca->Fnx=cd;

		da->Fnx=db;db->Fnx=dc;dc->Fnx=da;

		abc=ab->Cell;
		abd = new CCell;
		Cs.push_back(abd);
		abd->ab=ab;abd->bc=bd;abd->ca=da;
		abd->VV[0]=a;abd->VV[1]=b;abd->VV[2]=d;
		abd->set_center();//abd->alphavalue=-10;

		bcd = new CCell;
		Cs.push_back(bcd);
		bcd->ab=bc;bcd->bc=cd;bcd->ca=db;
		bcd->VV[0]=b;bcd->VV[1]=c;bcd->VV[2]=d;
		bcd->set_center();//bcd->alphavalue=-10;

		cad = new CCell;
		Cs.push_back(cad);
		cad->ab=ca;cad->bc=ad;cad->ca=dc;
		cad->VV[0]=c;cad->VV[1]=a;cad->VV[2]=d;
		cad->set_center();//cad->alphavalue=-10;

		abc->set_BN(0);
		abd->set_BN(1);
		bcd->set_BN(1);
		cad->set_BN(1);
		abc->CC[0] = abd;abc->CC[1] = bcd;abc->CC[2] = cad;

		ab->Cell = abd;
		bc->Cell = bcd;
		ca->Cell = cad;

		ad->Cell = cad;da->Cell = abd;
		bd->Cell = abd;db->Cell = bcd;
		cd->Cell = bcd;dc->Cell = cad;

		if (ab->A->N>=0 || ab->B->N>=0) {St.push_back(ab);ab->Flag=1;}
		if (bc->A->N>=0 || bc->B->N>=0) {St.push_back(bc);bc->Flag=1;}
		if (ca->A->N>=0 || ca->B->N>=0) {St.push_back(ca);ca->Flag=1;}

		while (St.size())
		{
			ab = St[St.size()-1]; St.pop_back(); ab->Flag&=2;
			if (ab->Flag) {continue;}
			a = ab->A;b=ab->B;c=ab->Rev->Fnx->B;
			if (insphere(ab->Rev->Cell,d)) 
			{
				flip_resample(ab);
				continue;
			}
		}
	}
	else
	{
		CXYZW *a,*b,*c,*e;
		ab = get_on_edge(ab->Cell,d);
		a = ab->A;
		b = ab->B;
		c=ab->Enx->B;
		e = ab->Rev->Enx->B;
		CEdge *ba,*bc,*cb,*ca,*ac,*ae,*ea,*eb,*be;
		CEdge *ad,*db,*dc,*da,*bd,*cd,*de,*ed;
		CCell *abc,*bae;
		CCell *adc,*dbc,*bde,*dae;

		ad = new CEdge(a,d);Es.push_back(ad);
		da = new CEdge(d,a);Es.push_back(da);
		ad->Rev=da;da->Rev=ad;
		db = new CEdge(d,b);Es.push_back(db);
		bd = new CEdge(b,d);Es.push_back(bd);
		db->Rev=bd;bd->Rev=db;
		dc = new CEdge(d,c);Es.push_back(dc);
		cd = new CEdge(c,d);Es.push_back(cd);
		dc->Rev=cd;cd->Rev=dc;
		de = new CEdge(d,e);Es.push_back(de);
		ed = new CEdge(e,d);Es.push_back(ed);
		de->Rev=ed;ed->Rev=de;

		ba=ab->Rev;
		bc=ab->Enx;cb=bc->Rev;
		ca=bc->Enx;ac=ca->Rev;
		ae=ba->Enx;ea=ae->Rev;
		eb=ae->Enx;be=eb->Rev;

		abc=ab->Cell;bae=ba->Cell;

		ad->set_center();ad->alphavalue=-10;
		da->set_center();da->alphavalue=-10;
		db->set_center();db->alphavalue=-10;
		bd->set_center();bd->alphavalue=-10;
		dc->set_center();dc->alphavalue=-10;
		cd->set_center();cd->alphavalue=-10;
		de->set_center();de->alphavalue=-10;
		ed->set_center();ed->alphavalue=-10;

		ad->Enx=dc;dc->Enx=ca;ca->Enx=ad;
		db->Enx=bc;bc->Enx=cd;cd->Enx=db;
		bd->Enx=de;de->Enx=eb;eb->Enx=bd;
		da->Enx=ae;ae->Enx=ed;ed->Enx=da;

		ad->Fnx=ac;ca->Fnx=cd;
		db->Fnx=dc;cd->Fnx=cb;
		bc->Fnx=bd;dc->Fnx=da;

		bd->Fnx=be;eb->Fnx=ed;
		da->Fnx=de;ae->Fnx=ad;
		de->Fnx=db;ed->Fnx=ea;

		adc = new CCell;
		Cs.push_back(adc);
		adc->ab=ad;adc->bc=dc;adc->ca=ca;
		adc->VV[0]=a;adc->VV[1]=d;adc->VV[2]=c;
		adc->set_center();//adc->alphavalue=-10;

		dbc = new CCell;
		Cs.push_back(dbc);
		dbc->ab=db;dbc->bc=bc;dbc->ca=cd;
		dbc->VV[0]=d;dbc->VV[1]=b;dbc->VV[2]=c;
		dbc->set_center();//dbc->alphavalue=-10;

		bde = new CCell;
		Cs.push_back(bde);
		bde->ab=bd;bde->bc=de;bde->ca=eb;
		bde->VV[0]=b;bde->VV[1]=d;bde->VV[2]=e;
		bde->set_center();//bde->alphavalue=-10;

		dae = new CCell;
		Cs.push_back(dae);
		dae->ab=da;dae->bc=ae;dae->ca=ed;
		dae->VV[0]=d;dae->VV[1]=a;dae->VV[2]=e;
		dae->set_center();//dae->alphavalue=-10;

		abc->set_BN(0);
		bae->set_BN(0);
		adc->set_BN(1);
		dbc->set_BN(1);
		bde->set_BN(1);
		dae->set_BN(1);
		abc->CC[0]=adc;abc->CC[1]=dbc;
		bae->CC[0]=bde;bae->CC[1]=dae;

		ad->Cell=adc;dc->Cell=adc;ca->Cell=adc;
		db->Cell=dbc;bc->Cell=dbc;cd->Cell=dbc;
		bd->Cell=bde;de->Cell=bde;eb->Cell=bde;
		da->Cell=dae;ae->Cell=dae;ed->Cell=dae;

		if (bc->A->N>=0 || bc->B->N>=0) {St.push_back(bc);bc->Flag=1;}
		if (ca->A->N>=0 || ca->B->N>=0) {St.push_back(ca);ca->Flag=1;}
		if (ae->A->N>=0 || ae->B->N>=0) {St.push_back(ae);ae->Flag=1;}
		if (eb->A->N>=0 || eb->B->N>=0) {St.push_back(eb);eb->Flag=1;}

		while (St.size())
		{
			ab = St[St.size()-1]; St.pop_back(); ab->Flag&=2;
			if (ab->Flag) {continue;}
			a = ab->A;b=ab->B;c=ab->Rev->Fnx->B;
			if (insphere(ab->Rev->Cell,d)) 
			{
				flip_resample(ab);
				continue;
			}
		}
	}
}

void CD3DW::add_interior(double l)
{
	int Csn = (int)Cs.size();
	bool loop=1;
	int flag = 0;
	double ra = 1.9;
	for (int i=0;i<Csn;i++)
	{
		if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
	}
	while (loop && flag<100)
	{
		loop = 0;
		CCell *G3 = NULL,*G2 = NULL,*G1 = NULL,*G0 = NULL;
		int types3 = -1; // 0: VV[0] and VV[1], 1: VV[0] and VV[2], 2: VV[1] and VV[2]
		int types2 = -1; // 0: VV[0] and VV[1], 1: VV[0] and VV[2], 2: VV[1] and VV[2]
		int types1 = -1; // 3: VV[0], 4: VV[1], 5: VV[2]
		int types0 = -1; // 3: VV[0], 4: VV[1], 5: VV[2]
		double alpha_max_3 = -10000.0;
		double alpha_max_2 = -10000.0;
		double alpha_max_1 = -10000.0;
		double alpha_max_0 = -10000.0;
		int Cs_n = (int)Cs.size();
		for (int i=0;i<Cs_n;i++)
		{
			if (Cs[i]->BN())
			{
				point3D *P0 = Cs[i]->get_VV(0)->P();
				point3D *P1 = Cs[i]->get_VV(1)->P();
				point3D *P2 = Cs[i]->get_VV(2)->P();
				double xc = Cs[i]->get_center(0);
				double yc = Cs[i]->get_center(1);
				if (P0->get_cell_boundary() &&
					P1->get_cell_boundary() &&
					P2->get_cell_boundary())
				{
					if      ((P0->getfp()==P1 && 
					  		 !P0->getab()->judge_vert_out_side(P2->x(),P2->y())) ||
							 (P0->getrp()==P1 && 
							 !P0->getba()->judge_vert_out_side(P2->x(),P2->y())))
					{
						Cs[i]->set_type(0);
						bool flagin = P0->getab()->Cell()->test_point_inside(xc,yc);
						if (flagin)
						{
							double alphav = Cs[i]->alphavalue;
							double multi = alphav/l;
							if (alphav>alpha_max_3 && multi>ra)
							{
								alpha_max_3 = alphav;
								G3 = Cs[i];
								types3 = 0;
								loop = 1;
							}
						}
					}
					else if ((P0->getfp()==P2 && 
							 !P0->getab()->judge_vert_out_side(P1->x(),P1->y())) || 
							 (P0->getrp()==P2 && 
							 !P0->getba()->judge_vert_out_side(P1->x(),P1->y())))
					{
						Cs[i]->set_type(0);
						bool flagin = P0->getab()->Cell()->test_point_inside(xc,yc);
						if (flagin)
						{
							double alphav = Cs[i]->alphavalue;
							double multi = alphav/l;
							if (alphav>alpha_max_3 && multi>ra)
							{
								alpha_max_3 = alphav;
								G3 = Cs[i];
								types3 = 1;
								loop = 1;
							}
						}
					}
					else if ((P1->getfp()==P2 && 
							 !P1->getab()->judge_vert_out_side(P0->x(),P0->y())) || 
							 (P1->getrp()==P2 && 
							 !P1->getba()->judge_vert_out_side(P0->x(),P0->y())))
					{
						Cs[i]->set_type(0);
						bool flagin = P0->getab()->Cell()->test_point_inside(xc,yc);
						if (flagin)
						{
							double alphav = Cs[i]->alphavalue;
							double multi = alphav/l;
							if (alphav>alpha_max_3 && multi>ra)
							{
								alpha_max_3 = alphav;
								G3 = Cs[i];
								types3 = 2;
								loop = 1;
							}
						}
					}
				}
				else if (P0->get_cell_boundary() &&
						 P1->get_cell_boundary() &&
						!P2->get_cell_boundary() &&
						(P0->getfp()==P1 && !P0->getab()->judge_vert_out_side(P2->x(),P2->y()) || 
						 P0->getrp()==P1 && !P0->getba()->judge_vert_out_side(P2->x(),P2->y())))
				{
					Cs[i]->set_type(0);
					bool flagin = P0->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_2 && multi>ra)
						{
							alpha_max_2 = alphav;
							G2 = Cs[i];
							types2 = 0;
							loop = 1;
						}
					}
				}
				else if (P0->get_cell_boundary() &&
						!P1->get_cell_boundary() &&
						 P2->get_cell_boundary() &&
						(P0->getfp()==P2 && !P0->getab()->judge_vert_out_side(P1->x(),P1->y()) || 
						 P0->getrp()==P2 && !P0->getba()->judge_vert_out_side(P1->x(),P1->y())))
				{
					Cs[i]->set_type(0);
					bool flagin = P0->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_2 && multi>ra)
						{
							alpha_max_2 = alphav;
							G2 = Cs[i];
							types2 = 1;
							loop = 1;
						}
					}
				}
				else if (!P0->get_cell_boundary() &&
						  P1->get_cell_boundary() &&
						  P2->get_cell_boundary() &&
						 (P1->getfp()==P2 && !P1->getab()->judge_vert_out_side(P0->x(),P0->y()) || 
						  P1->getrp()==P2 && !P1->getba()->judge_vert_out_side(P0->x(),P0->y())))
				{
					Cs[i]->set_type(0);
					bool flagin = P1->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_2 && multi>ra)
						{
							alpha_max_2 = alphav;
							G2 = Cs[i];
							types2 = 2;
							loop = 1;
						}
					}
				}
				else if (P0->get_cell_boundary() &&
						!P1->get_cell_boundary() &&
						!P2->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					bool flagin = P0->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_1 && multi>ra)
						{
							alpha_max_1 = alphav;
							G1 = Cs[i];
							types1 = 3;
							loop = 1;
						}
					}
				}
				else if (!P0->get_cell_boundary() &&
						  P1->get_cell_boundary() &&
						 !P2->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					bool flagin = P1->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_1 && multi>ra)
						{
							alpha_max_1 = alphav;
							G1 = Cs[i];
							types1 = 3;						
							loop = 1;
						}
					}
				}
				else if (!P0->get_cell_boundary() &&
						 !P1->get_cell_boundary() &&
						  P2->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					bool flagin = P2->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_1 && multi>ra)
						{
							alpha_max_1 = alphav;
							G1 = Cs[i];
							types1 = 3;
							loop = 1;
						}
					}
				}
				else if (!P0->get_cell_boundary() &&
						 !P1->get_cell_boundary() &&
						 !P2->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					bool flagin = Vs[3]->P()->getab()->Cell()->test_point_inside(xc,yc);
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_0 && multi>ra)
						{
							alpha_max_0 = alphav;
							G0 = Cs[i];
							types0 = 3;
							loop = 1;
						}
					}
				}
			}
		}
		///////////////////
		if      (alpha_max_3>0)
		{
			double multi = alpha_max_3/l;
			if (multi>=ra && types3>=0)
			{
				double cx = G3->get_center(0);
				double cy = G3->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				int Csnt = (int)Cs.size();
				for (int i=0;i<Csnt;i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		else if (alpha_max_3<0 && alpha_max_2>0)
		{
			double multi = alpha_max_2/l;
			if (multi>=ra && types2>=0)
			{
				double cx = G2->get_center(0);
				double cy = G2->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				int Csnt = (int)Cs.size();
				for (int i=0;i<Csnt;i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		else if (alpha_max_3<0 && alpha_max_2<0 && alpha_max_1>0)
		{
			double multi = alpha_max_1/l;
			if (multi>=ra && types1>=0)
			{
				double cx = G1->get_center(0);
				double cy = G1->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				int Csnt = (int)Cs.size();
				for (int i=0;i<Csnt;i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		else if (alpha_max_3<0 && alpha_max_2<0 && alpha_max_1<0 && alpha_max_0>0)
		{
			double multi = alpha_max_0/l;
			if (multi>=ra && types0>=0)
			{
				double cx = G0->get_center(0);
				double cy = G0->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				int Csnt = (int)Cs.size();
				for (int i=0;i<Csnt;i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
	}
	///////////////// for test purpose /////////////////
	int Cs_n1 = (int)Cs.size();
	for (int i=0;i<Cs_n1;i++) 
	{
		if (Cs[i]->BN())
		{
			if (Cs[i]->inside_complex())
			{
				Cs[i]->set_type(0);
			}
		}
	}
}

void CD3DW::add_interior_debug(double l)
{
	bool loop=1;
	int flag = 0;
	double ra = 1.6;
	for (int i=0;i<(int)Cs.size();i++)
	{
		if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
	}
	while (loop && flag<3)
	{
		loop = 0;
		CCell *G3,*G2,*G1,*G0;
		int types3 = -1; // 0: VV[0] and VV[1], 1: VV[0] and VV[2], 2: VV[1] and VV[2]
		int types2 = -1; // 0: VV[0] and VV[1], 1: VV[0] and VV[2], 2: VV[1] and VV[2]
		int types1 = -1; // 3: VV[0], 4: VV[1], 5: VV[2]
		int types0 = -1; // 3: VV[0], 4: VV[1], 5: VV[2]
		double alpha_max_3 = -10000.0;
		double alpha_max_2 = -10000.0;
		double alpha_max_1 = -10000.0;
		double alpha_max_0 = -10000.0;
		for (int i=0;i<(int)Cs.size();i++)
		{
			if (Cs[i]->BN())
			{
				if (Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
					Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
					Cs[i]->get_VV(2)->P()->get_cell_boundary())
				{
					if      ((Cs[i]->get_VV(0)->P()->getfp()==Cs[i]->get_VV(1)->P() && 
					  		 !Cs[i]->get_VV(0)->P()->getab()->judge_vert_out_side(Cs[i]->get_VV(2)->P()->x(),Cs[i]->get_VV(2)->P()->y())) ||
							 (Cs[i]->get_VV(0)->P()->getrp()==Cs[i]->get_VV(1)->P() && 
							 !Cs[i]->get_VV(0)->P()->getba()->judge_vert_out_side(Cs[i]->get_VV(2)->P()->x(),Cs[i]->get_VV(2)->P()->y())))
					{
						Cs[i]->set_type(0);
						///////////////////
						bool flagin = true;
						point3D *Ps = Cs[i]->get_VV(0)->P();
						point3D *Pt = Cs[i]->get_VV(0)->P();
						while (Pt!=Ps->getfp() && flagin)
						{
							if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
							{flagin = false;}
							Ps = Ps->getfp();
						}
						Ps = NULL;Pt = NULL;
						delete Ps;delete Pt;
						////////////////////
						if (flagin)
						{
							double alphav = Cs[i]->alphavalue;
							double multi = alphav/l;
							if (alphav>alpha_max_3 && multi>ra)
							{
								alpha_max_3 = alphav;
								G3 = Cs[i];
								types3 = 0;
								loop = 1;
							}
						}
					}
					else if ((Cs[i]->get_VV(0)->P()->getfp()==Cs[i]->get_VV(2)->P() && 
							 !Cs[i]->get_VV(0)->P()->getab()->judge_vert_out_side(Cs[i]->get_VV(1)->P()->x(),Cs[i]->get_VV(1)->P()->y())) || 
							 (Cs[i]->get_VV(0)->P()->getrp()==Cs[i]->get_VV(2)->P() && 
							 !Cs[i]->get_VV(0)->P()->getba()->judge_vert_out_side(Cs[i]->get_VV(1)->P()->x(),Cs[i]->get_VV(1)->P()->y())))
					{
						Cs[i]->set_type(0);
						///////////////////
						bool flagin = true;
						point3D *Ps = Cs[i]->get_VV(0)->P();
						point3D *Pt = Cs[i]->get_VV(0)->P();
						while (Pt!=Ps->getfp() && flagin)
						{
							if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
							{flagin = false;}
							Ps = Ps->getfp();
						}
						Ps = NULL;Pt = NULL;
						delete Ps;delete Pt;
						////////////////////
						if (flagin)
						{
							double alphav = Cs[i]->alphavalue;
							double multi = alphav/l;
							if (alphav>alpha_max_3 && multi>ra)
							{
								alpha_max_3 = alphav;
								G3 = Cs[i];
								types3 = 1;
								loop = 1;
							}
						}
					}
					else if ((Cs[i]->get_VV(1)->P()->getfp()==Cs[i]->get_VV(2)->P() && 
							 !Cs[i]->get_VV(1)->P()->getab()->judge_vert_out_side(Cs[i]->get_VV(0)->P()->x(),Cs[i]->get_VV(0)->P()->y())) || 
							 (Cs[i]->get_VV(1)->P()->getrp()==Cs[i]->get_VV(2)->P() && 
							 !Cs[i]->get_VV(1)->P()->getba()->judge_vert_out_side(Cs[i]->get_VV(0)->P()->x(),Cs[i]->get_VV(0)->P()->y())))
					{
						Cs[i]->set_type(0);
						///////////////////
						bool flagin = true;
						point3D *Ps = Cs[i]->get_VV(0)->P();
						point3D *Pt = Cs[i]->get_VV(0)->P();
						while (Pt!=Ps->getfp() && flagin)
						{
							if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
							{flagin = false;}
							Ps = Ps->getfp();
						}
						Ps = NULL;Pt = NULL;
						delete Ps;delete Pt;
						////////////////////
						if (flagin)
						{
							double alphav = Cs[i]->alphavalue;
							double multi = alphav/l;
							if (alphav>alpha_max_3 && multi>ra)
							{
								alpha_max_3 = alphav;
								G3 = Cs[i];
								types3 = 2;
								loop = 1;
							}
						}
					}
				}
				else if (Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						 Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						!Cs[i]->get_VV(2)->P()->get_cell_boundary() &&
						(Cs[i]->get_VV(0)->P()->getfp()==Cs[i]->get_VV(1)->P() && !Cs[i]->get_VV(0)->P()->getab()->judge_vert_out_side(Cs[i]->get_VV(2)->x(),Cs[i]->get_VV(2)->y()) || 
						 Cs[i]->get_VV(0)->P()->getrp()==Cs[i]->get_VV(1)->P() && !Cs[i]->get_VV(0)->P()->getba()->judge_vert_out_side(Cs[i]->get_VV(2)->x(),Cs[i]->get_VV(2)->y())))
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Cs[i]->get_VV(0)->P();
					point3D *Pt = Cs[i]->get_VV(0)->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_2 && multi>ra)
						{
							alpha_max_2 = alphav;
							G2 = Cs[i];
							types2 = 0;
							loop = 1;
						}
					}
				}
				else if (Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						!Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						 Cs[i]->get_VV(2)->P()->get_cell_boundary() &&
						(Cs[i]->get_VV(0)->P()->getfp()==Cs[i]->get_VV(2)->P() && !Cs[i]->get_VV(0)->P()->getab()->judge_vert_out_side(Cs[i]->get_VV(1)->x(),Cs[i]->get_VV(1)->y()) || 
						 Cs[i]->get_VV(0)->P()->getrp()==Cs[i]->get_VV(2)->P() && !Cs[i]->get_VV(0)->P()->getba()->judge_vert_out_side(Cs[i]->get_VV(1)->x(),Cs[i]->get_VV(1)->y())))
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Cs[i]->get_VV(0)->P();
					point3D *Pt = Cs[i]->get_VV(0)->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_2 && multi>ra)
						{
							alpha_max_2 = alphav;
							G2 = Cs[i];
							types2 = 1;
							loop = 1;
						}
					}
				}
				else if (!Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						  Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						  Cs[i]->get_VV(2)->P()->get_cell_boundary() &&
						 (Cs[i]->get_VV(1)->P()->getfp()==Cs[i]->get_VV(2)->P() && !Cs[i]->get_VV(1)->P()->getab()->judge_vert_out_side(Cs[i]->get_VV(0)->x(),Cs[i]->get_VV(0)->y()) || 
						  Cs[i]->get_VV(1)->P()->getrp()==Cs[i]->get_VV(2)->P() && !Cs[i]->get_VV(1)->P()->getba()->judge_vert_out_side(Cs[i]->get_VV(0)->x(),Cs[i]->get_VV(0)->y())))
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Cs[i]->get_VV(1)->P();
					point3D *Pt = Cs[i]->get_VV(1)->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_2 && multi>ra)
						{
							alpha_max_2 = alphav;
							G2 = Cs[i];
							types2 = 2;
							loop = 1;
						}
					}
				}
				else if (Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						!Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						!Cs[i]->get_VV(2)->P()->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Cs[i]->get_VV(0)->P();
					point3D *Pt = Cs[i]->get_VV(0)->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_1 && multi>ra)
						{
							alpha_max_1 = alphav;
							G1 = Cs[i];
							types1 = 3;
							loop = 1;
						}
					}
				}
				else if (!Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						  Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						 !Cs[i]->get_VV(2)->P()->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Cs[i]->get_VV(1)->P();
					point3D *Pt = Cs[i]->get_VV(1)->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_1 && multi>ra)
						{
							alpha_max_1 = alphav;
							G1 = Cs[i];
							types1 = 3;						
							loop = 1;
						}
					}
				}
				else if (!Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						 !Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						  Cs[i]->get_VV(2)->P()->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Cs[i]->get_VV(2)->P();
					point3D *Pt = Cs[i]->get_VV(2)->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_1 && multi>ra)
						{
							alpha_max_1 = alphav;
							G1 = Cs[i];
							types1 = 3;
							loop = 1;
						}
					}
				}
				else if (!Cs[i]->get_VV(0)->P()->get_cell_boundary() &&
						 !Cs[i]->get_VV(1)->P()->get_cell_boundary() &&
						 !Cs[i]->get_VV(2)->P()->get_cell_boundary())
				{
					Cs[i]->set_type(0);
					///////////////////
					bool flagin = true;
					point3D *Ps = Vs[3]->P();
					point3D *Pt = Vs[3]->P();
					while (Pt!=Ps->getfp() && flagin)
					{
						if (Ps->getab()->judge_vert_out_side(Cs[i]->get_center(0),Cs[i]->get_center(1))) 
						{flagin = false;}
						Ps = Ps->getfp();
					}
					Ps = NULL;Pt = NULL;
					delete Ps;delete Pt;
					////////////////////
					if (flagin)
					{
						double alphav = Cs[i]->alphavalue;
						double multi = alphav/l;
						if (alphav>alpha_max_0 && multi>ra)
						{
							alpha_max_0 = alphav;
							G0 = Cs[i];
							types0 = 3;
							loop = 1;
						}
					}
				}
			}
		}
		///////////////////
		if      (alpha_max_3>0)
		{
			double multi = alpha_max_3/l;
			if (multi>=ra && types3>=0)
			{
				/*double c1x = G3->get_center(0);
				double c1y = G3->get_center(1);
				double c0x = 0;
				double c0y = 0;
				if      (types3==0)
				{
					c0x = (G3->get_VV(0)->x() + G3->get_VV(1)->x())/2;
					c0y = (G3->get_VV(0)->y() + G3->get_VV(1)->y())/2;
				}
				else if (types3==1)
				{
					c0x = (G3->get_VV(0)->x() + G3->get_VV(2)->x())/2;
					c0y = (G3->get_VV(0)->y() + G3->get_VV(2)->y())/2;
				}
				else if (types3==2)
				{
					c0x = (G3->get_VV(1)->x() + G3->get_VV(2)->x())/2;
					c0y = (G3->get_VV(1)->y() + G3->get_VV(2)->y())/2;
				}
				double ta = 1 - 0.866*sqrt(l/((c0x-c1x)*(c0x-c1x)+(c0y-c1y)*(c0y-c1y)));
				double cx = (1-ta)*c1x+ta*c0x;
				double cy = (1-ta)*c1y+ta*c0y;
				//////////
				cx = G3->get_center(0);
				cy = G3->get_center(1);
				//////////*/
				double cx = G3->get_center(0);
				double cy = G3->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				for (int i=0;i<(int)Cs.size();i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		else if (alpha_max_3<0 && alpha_max_2>0)
		{
			double multi = alpha_max_2/l;
			if (multi>=ra && types2>=0)
			{
				/*double c1x = G2->get_center(0);
				double c1y = G2->get_center(1);
				double c0x = 0;
				double c0y = 0;
				if      (types2==0)
				{
					c0x = (G2->get_VV(0)->x() + G2->get_VV(1)->x())/2;
					c0y = (G2->get_VV(0)->y() + G2->get_VV(1)->y())/2;
				}
				else if (types2==1)
				{
					c0x = (G2->get_VV(0)->x() + G2->get_VV(2)->x())/2;
					c0y = (G2->get_VV(0)->y() + G2->get_VV(2)->y())/2;
				}
				else if (types2==2)
				{
					c0x = (G2->get_VV(1)->x() + G2->get_VV(2)->x())/2;
					c0y = (G2->get_VV(1)->y() + G2->get_VV(2)->y())/2;
				}
				double ta = 1 - 0.866*sqrt(l/((c0x-c1x)*(c0x-c1x)+(c0y-c1y)*(c0y-c1y)));
				double cx = (1-ta)*c1x+ta*c0x;
				double cy = (1-ta)*c1y+ta*c0y;
				//////////
				cx = G2->get_center(0);
				cy = G2->get_center(1);
				//////////*/
				double cx = G2->get_center(0);
				double cy = G2->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				//P0->~point3D();
				for (int i=0;i<(int)Cs.size();i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		else if (alpha_max_3<0 && alpha_max_2<0 && alpha_max_1>0)
		{
			double multi = alpha_max_1/l;
			if (multi>=ra && types1>=0)
			{
				/*double c1x = G1->get_center(0);
				double c1y = G1->get_center(1);
				double c0x = 0;
				double c0y = 0;
				if     (types1==3)
				{
					c0x = (G1->get_VV(1)->x() + G1->get_VV(2)->x())/2;
					c0y = (G1->get_VV(1)->y() + G1->get_VV(2)->y())/2;
				}
				else if (types1==4)
				{
					c0x = (G1->get_VV(0)->x() + G1->get_VV(2)->x())/2;
					c0y = (G1->get_VV(0)->y() + G1->get_VV(2)->y())/2;
				}
				else if (types1==5)
				{
					c0x = (G1->get_VV(0)->x() + G1->get_VV(1)->x())/2;
					c0y = (G1->get_VV(0)->y() + G1->get_VV(1)->y())/2;
				}
				double ta = 1 - 0.866*sqrt(l/((c0x-c1x)*(c0x-c1x)+(c0y-c1y)*(c0y-c1y)));
				double cx = (1-ta)*c1x+ta*c0x;
				double cy = (1-ta)*c1y+ta*c0y;
				//////////
				cx = G1->get_center(0);
				cy = G1->get_center(1);
				//////////*/
				double cx = G1->get_center(0);
				double cy = G1->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				//P0->~point3D();
				for (int i=0;i<(int)Cs.size();i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		else if (alpha_max_3<0 && alpha_max_2<0 && alpha_max_1<0 && alpha_max_0>0)
		{
			double multi = alpha_max_0/l;
			if (multi>=ra && types0>=0)
			{
				/*double c1x = G0->get_center(0);
				double c1y = G0->get_center(1);
				double c0x = 0;
				double c0y = 0;
				if      (types0==3)
				{
					c0x = (G0->get_VV(1)->x() + G0->get_VV(2)->x())/2;
					c0y = (G0->get_VV(1)->y() + G0->get_VV(2)->y())/2;
				}
				else if (types0==4)
				{
					c0x = (G0->get_VV(0)->x() + G0->get_VV(2)->x())/2;
					c0y = (G0->get_VV(0)->y() + G0->get_VV(2)->y())/2;
				}
				else if (types0==5)
				{
					c0x = (G0->get_VV(0)->x() + G0->get_VV(1)->x())/2;
					c0y = (G0->get_VV(0)->y() + G0->get_VV(1)->y())/2;
				}
				double ta = 1 - 0.866*sqrt(l/((c0x-c1x)*(c0x-c1x)+(c0y-c1y)*(c0y-c1y)));
				double cx = (1-ta)*c1x+ta*c0x;
				double cy = (1-ta)*c1y+ta*c0y;
				//////////
				cx = G0->get_center(0);
				cy = G0->get_center(1);
				//////////*/
				double cx = G0->get_center(0);
				double cy = G0->get_center(1);
				point3D *P0 = new point3D(flag, cx, cy);
				P0->set_cell_boundary(0);
				add_resample(cx,cy,l,flag,P0);
				//P0->~point3D();
				for (int i=0;i<(int)Cs.size();i++)
				{
					if (Cs[i]->get_VV(0)->N<0 || Cs[i]->get_VV(1)->N<0 || Cs[i]->get_VV(2)->N<0) {Cs[i]->set_BN(0);}
				}
				flag++;
			}
		}
		//////////
		G3 = NULL;
		G2 = NULL;
		G1 = NULL;
		G0 = NULL;
		delete G3;
		delete G2;
		delete G1;
		delete G0;
	}
	///////////////// for test purpose /////////////////
	for (int i=0;i<(int)Cs.size();i++) 
	{
		//cout<<i<<": "<<Cs[i]->get_center(0)<<" "<<Cs[i]->get_center(1)<<"     ";
		if (Cs[i]->BN())
		{
			//cout<<i<<": "<<Cs[i]->get_VV(0)->P()->id()<<" "<<Cs[i]->get_VV(1)->P()->id()<<" "<<Cs[i]->get_VV(2)->P()->id()<<"      ";
			if (Cs[i]->inside_complex())
			{
				Cs[i]->set_type(0);
			}
		}
	}
	/////////////////////////////////////////////////////
}

inline void CD3DW::flip(CEdge *ab)
{
	CXYZW *a,*b,*c,*d;
	CEdge *ba;
	CEdge *ac,*ca,*bc,*cb,*ad,*da,*bd,*db,*cd,*dc;

	/****************
	    d
	   /|\
	  /   \
	 /  |  \
    /       \
	a--->----b
    \   |   /
     \     /
      \ | /
	   \ /
		c	    
	****************/
	ba = ab->Rev;
	bd = ab->Enx;db = bd->Rev;
	ad = ab->Fnx;da = ad->Rev;
	ac = ba->Enx;ca = ac->Rev;
	bc = ba->Fnx;cb = bc->Rev;

	a=ab->A;b=ab->B;c=bc->B;d=bd->B;

	cd = new CEdge(c,d);Es.push_back(cd);
	dc = new CEdge(d,c);Es.push_back(dc);
	cd->Rev=dc;dc->Rev=cd;
	cd->set_center();cd->alphavalue=-10;
	dc->set_center();dc->alphavalue=-10;

	cd->Enx=da;da->Enx=ac;ac->Enx=cd;
	cd->Fnx=ca;da->Fnx=dc;ac->Fnx=ad;
	dc->Enx=cb;cb->Enx=bd;bd->Enx=dc;
	dc->Fnx=db;cb->Fnx=cd;bd->Fnx=bc;
	
	CCell *acd,*dcb,*abd,*acb;
	abd=ab->Cell;acb=ba->Cell;

	acd = new CCell;Cs.push_back(acd);
	acd->ab=ac;acd->bc=cd;acd->ca=da;
	acd->VV[0]=a;acd->VV[1]=c;acd->VV[2]=d;
	acd->set_center();acd->alphavalue=-10;

	dcb = new CCell;Cs.push_back(dcb);
	dcb->ab=dc;dcb->bc=cb;dcb->ca=bd;
	dcb->VV[0]=d;dcb->VV[1]=c;dcb->VV[2]=b;
	dcb->set_center();dcb->alphavalue=-10;

	abd->set_BN(0);
	acb->set_BN(0);
	acd->set_BN(1);
	dcb->set_BN(1);
	
	abd->CC[0]=acd;abd->CC[1]=dcb;
	acb->CC[0]=acd;acb->CC[1]=dcb;

	ac->Cell=acd;cd->Cell=acd;da->Cell=acd;
	dc->Cell=dcb;cb->Cell=dcb;bd->Cell=dcb;

	ab->Flag|=2;ba->Flag|=2;
	
	if (ac->A->N>=0 || ac->B->N>=0) {St.push_back(ac);ac->Flag=1;}
	if (cb->A->N>=0 || cb->B->N>=0) {St.push_back(cb);cb->Flag=1;}
}

inline void CD3DW::flip_resample(CEdge *ab)
{
	CXYZW *a,*b,*c,*d;
	CEdge *ba;
	CEdge *ac,*ca,*bc,*cb,*ad,*da,*bd,*db,*cd,*dc;

	/****************
	    d
	   /|\
	  /   \
	 /  |  \
    /       \
	a--->----b
    \   |   /
     \     /
      \ | /
	   \ /
		c	    
	****************/
	ba = ab->Rev;
	bd = ab->Enx;db = bd->Rev;
	ad = ab->Fnx;da = ad->Rev;
	ac = ba->Enx;ca = ac->Rev;
	bc = ba->Fnx;cb = bc->Rev;

	a=ab->A;b=ab->B;c=bc->B;d=bd->B;

	cd = new CEdge(c,d);Es.push_back(cd);
	dc = new CEdge(d,c);Es.push_back(dc);
	cd->Rev=dc;dc->Rev=cd;
	cd->set_center();cd->alphavalue=-10;
	dc->set_center();dc->alphavalue=-10;

	cd->Enx=da;da->Enx=ac;ac->Enx=cd;
	cd->Fnx=ca;da->Fnx=dc;ac->Fnx=ad;
	dc->Enx=cb;cb->Enx=bd;bd->Enx=dc;
	dc->Fnx=db;cb->Fnx=cd;bd->Fnx=bc;
	
	CCell *acd,*dcb,*abd,*acb;
	abd=ab->Cell;acb=ba->Cell;

	acd = new CCell;Cs.push_back(acd);
	acd->ab=ac;acd->bc=cd;acd->ca=da;
	acd->VV[0]=a;acd->VV[1]=c;acd->VV[2]=d;
	acd->set_center();//acd->alphavalue=-10;

	dcb = new CCell;Cs.push_back(dcb);
	dcb->ab=dc;dcb->bc=cb;dcb->ca=bd;
	dcb->VV[0]=d;dcb->VV[1]=c;dcb->VV[2]=b;
	dcb->set_center();//dcb->alphavalue=-10;

	abd->set_BN(0);
	acb->set_BN(0);
	acd->set_BN(1);
	dcb->set_BN(1);
	
	abd->CC[0]=acd;abd->CC[1]=dcb;
	acb->CC[0]=acd;acb->CC[1]=dcb;

	ac->Cell=acd;cd->Cell=acd;da->Cell=acd;
	dc->Cell=dcb;cb->Cell=dcb;bd->Cell=dcb;

	ab->Flag|=2;ba->Flag|=2;
	
	if (ac->A->N>=0 || ac->B->N>=0) {St.push_back(ac);ac->Flag=1;}
	if (cb->A->N>=0 || cb->B->N>=0) {St.push_back(cb);cb->Flag=1;}
}

inline double CD3DW::judge_inside(CCell *C,CXYZW *V)
{
	double inside_v = 0;

	double xa = C->VV[0]->X;
	double ya = C->VV[0]->Y;
	double xb = C->VV[1]->X;
	double yb = C->VV[1]->Y;
	double xc = C->VV[2]->X;
	double yc = C->VV[2]->Y;
	double xd = V->X;
	double yd = V->Y;

	double det_ab = (xb*yd - xd*yb) - (xa*yd - xd*ya) + (xa*yb - xb*ya);
	double det_bc = (xc*yd - xd*yc) - (xb*yd - xd*yb) + (xb*yc - xc*yb);
	double det_ca = (xa*yd - xd*ya) - (xc*yd - xd*yc) + (xc*ya - xa*yc);

	double c1 = abs(det_ab);
	double c2 = abs(det_bc);
	double c3 = abs(det_ca);

	double v1 = 0, v2 = 0, v3 = 0;
	if      (c1>c2 && c2>c3) {v1 = c1;v2 = c2;v3 = c3;}
	else if (c1>c3 && c3>c2) {v1 = c1;v2 = c3;v3 = c2;}
	else if (c2>c1 && c1>c3) {v1 = c2;v2 = c1;v3 = c3;}
	else if (c2>c3 && c3>c1) {v1 = c2;v2 = c3;v3 = c1;}
	else if (c3>c1 && c1>c2) {v1 = c3;v2 = c1;v3 = c2;}
	else if (c3>c2 && c2>c1) {v1 = c3;v2 = c2;v3 = c1;}

	inside_v = v3;

	return inside_v;
}

inline bool CD3DW::inside(CCell *C,CXYZW *V)
{
	bool valid = false;
	
	double xa = C->VV[0]->X;
	double ya = C->VV[0]->Y;
	double xb = C->VV[1]->X;
	double yb = C->VV[1]->Y;
	double xc = C->VV[2]->X;
	double yc = C->VV[2]->Y;
	double xd = V->X;
	double yd = V->Y;

	double det_ab = (xb*yd - xd*yb) - (xa*yd - xd*ya) + (xa*yb - xb*ya);
	double det_bc = (xc*yd - xd*yc) - (xb*yd - xd*yb) + (xb*yc - xc*yb);
	double det_ca = (xa*yd - xd*ya) - (xc*yd - xd*yc) + (xc*ya - xa*yc);

	if (abs(det_ab)<0.0000001) {det_ab=0;}
	if (abs(det_bc)<0.0000001) {det_bc=0;}
	if (abs(det_ca)<0.0000001) {det_ca=0;}

	double c_ab_bc = det_ab*det_bc;
	double c_bc_ca = det_bc*det_ca;
	double c_ca_ab = det_ca*det_ab;

	if (c_ab_bc>=0 && c_bc_ca>=0 && c_ca_ab>=0) {valid = true;}

	return valid;
}

inline bool CD3DW::insphere(CCell *C,CXYZW *V)
{
	bool valid = false;

	double xa = C->VV[0]->X;
	double ya = C->VV[0]->Y;
	double xb = C->VV[1]->X;
	double yb = C->VV[1]->Y;
	double xc = C->VV[2]->X;
	double yc = C->VV[2]->Y;
	double xd = V->X;
	double yd = V->Y;

	double D = 2*(xa*(yb-yc)+xb*(yc-ya)+xc*(ya-yb));
	double Ux = ((xa*xa+ya*ya)*(yb-yc)+(xb*xb+yb*yb)*(yc-ya)+(xc*xc+yc*yc)*(ya-yb))/D;
	double Uy = ((xa*xa+ya*ya)*(xc-xb)+(xb*xb+yb*yb)*(xa-xc)+(xc*xc+yc*yc)*(xb-xa))/D;

	double distUd = (Ux-xd)*(Ux-xd)+(Uy-yd)*(Uy-yd);
	double distUa = (Ux-xa)*(Ux-xa)+(Uy-ya)*(Uy-ya);

	double det = 0;

	if      (distUa>distUd) {det=1;}
	else if (distUa<=distUd) {det=-1;}

	if (det>0) {valid = true;}
	return valid;
}

inline CEdge* CD3DW::locate(CXYZW *V)
{
	CCell *C = Cs[0];
	while (C->CC[0])
	{
		if(!C->CC[1]) {C=C->CC[0];continue;}
		if(inside(C->CC[0],V)) {C=C->CC[0];continue;}
		if(!C->CC[2]) {C=C->CC[1];continue;}
		if(inside(C->CC[1],V)) {C=C->CC[1];continue;}
		C=C->CC[2];
	}
	return C->ab;
}

inline bool CD3DW::judge_on_edge(CCell *C,CXYZW *V)
{
	bool valid = false;
	double xa = C->VV[0]->X;
	double ya = C->VV[0]->Y;
	double xb = C->VV[1]->X;
	double yb = C->VV[1]->Y;
	double xc = C->VV[2]->X;
	double yc = C->VV[2]->Y;
	double xd = V->X;
	double yd = V->Y;

	double det_ab = (xb*yd - xd*yb) - (xa*yd - xd*ya) + (xa*yb - xb*ya);
	double det_bc = (xc*yd - xd*yc) - (xb*yd - xd*yb) + (xb*yc - xc*yb);
	double det_ca = (xa*yd - xd*ya) - (xc*yd - xd*yc) + (xc*ya - xa*yc);

	if (abs(det_ab)<0.0000001) {det_ab=0;}
	if (abs(det_bc)<0.0000001) {det_bc=0;}
	if (abs(det_ca)<0.0000001) {det_ca=0;}

	if (det_ab==0 || det_bc==0 || det_ca==0) {valid = true;}

	return valid;
}

inline CEdge* CD3DW::get_on_edge(CCell *C,CXYZW *V)
{
	CEdge *ab = NULL;
	double xa = C->VV[0]->X;
	double ya = C->VV[0]->Y;
	double xb = C->VV[1]->X;
	double yb = C->VV[1]->Y;
	double xc = C->VV[2]->X;
	double yc = C->VV[2]->Y;
	double xd = V->X;
	double yd = V->Y;

	double det_ab = (xb*yd - xd*yb) - (xa*yd - xd*ya) + (xa*yb - xb*ya);
	double det_bc = (xc*yd - xd*yc) - (xb*yd - xd*yb) + (xb*yc - xc*yb);
	double det_ca = (xa*yd - xd*ya) - (xc*yd - xd*yc) + (xc*ya - xa*yc);
	double abs_ab = abs(det_ab);
	double abs_bc = abs(det_bc);
	double abs_ca = abs(det_ca);

	if      (abs_ab<0.0000001) {ab = C->ab;}
	else if (abs_bc<0.0000001) {ab = C->bc;}
	else if (abs_ca<0.0000001) {ab = C->ca;}
	else
	{
		if      (abs_ab<abs_bc && abs_ab<abs_ca) {ab = C->ab;}
		else if (abs_bc<abs_ab && abs_ab<abs_ca) {ab = C->bc;}
		else if (abs_ca<abs_bc && abs_ab<abs_ab) {ab = C->ca;}
	}
	return ab;
}

void CComplex::sort(CSimplex *A)
{
	double alphav1 = A->alphavalue;
	int size = (int)P.size();
	int pos1 = 0;
	int pos2 = size-1;
	int pos3 = 0;
	int leng = 0,leng1 = 0;
	int flag = 0;

	if (size==0)
	{
		P.push_back(A);
		flag = 1;
	}
	else if (size==1)
	{
		if (P[0]->alphavalue<=alphav1)
		{
			P.push_back(A);
		}
		else {P.insert(P.begin(),A);}
		flag = 1;
	}
	else if (size==2)
	{
		flag = 2;
	}
	else if (size>2)
	{
		double nb = log(2.0);
		double ni = size;
		double nv = log(ni)/nb;
		int checkn = int(ceil(nv));

		while (checkn>=0 && flag==0)
		{
			if (P[pos1]->alphavalue>alphav1)
			{
				P.insert(P.begin()+pos1,A);
				flag = 1;
			}
			else if (P[pos2]->alphavalue<=alphav1)
			{
				if (pos2==P.size()-1) {P.push_back(A);}
				else {P.insert(P.begin()+pos2+1,A);}
				flag = 1;
			}
			else
			{
				leng = pos2-pos1+1;
				if (leng==2)
				{
					flag = 2;
				}
				else if (leng>2)
				{
					if (leng%2==0) {leng1=leng/2;}
					else {leng1=(leng-1)/2;}
					pos3 = pos1 + leng1;
					if (P[pos3]->alphavalue>alphav1)
					{
						pos2=pos3;
					}
					else if (P[pos3]->alphavalue<=alphav1)
					{
						pos1=pos3;
					}
				}
			}
			checkn--;
		}
	}

	if (flag==2)
	{
		if (P[pos1]->alphavalue>alphav1)
		{
			P.insert(P.begin()+pos1,A);
		}
		else if (P[pos2]->alphavalue<=alphav1)
		{
			if (pos2==P.size()-1) {P.push_back(A);}
			else {P.insert(P.begin()+pos2+1,A);}
		}
		else 
		{
			P.insert(P.begin()+pos2,A);
		}
	}
}

CComplex::CComplex(CD3DW *D)
{
	int i=0;
	CEdge *E,*F;
	for (i=0;i<(int)D->Es.size();i++) 
	{
		D->Es[i]->S1=0;D->Es[i]->S2=0;
	}
	for (i=0;i<(int)D->Es.size();i++)
	{
		E = D->Es[i];
		if (E->Flag&2) continue;
		if (E->A->N<0 || E->B->N<0 || E->Fnx->B->N<0) continue;
		if (E->S2) continue;
		int pn = (int)P.size();
		C2Simplex *S2 = new C2Simplex;
		P.push_back(S2);
		S2 = NULL;
		delete S2;
		P[pn]->E = E;
		
		((C2Simplex*)P[pn])->A = E->A;
		((C2Simplex*)P[pn])->a = E->A->N;
		((C2Simplex*)P[pn])->B = E->B;
		((C2Simplex*)P[pn])->b = E->B->N;
		((C2Simplex*)P[pn])->C = E->Fnx->B;
		((C2Simplex*)P[pn])->c = E->Fnx->B->N;
		E->S2 = ((C2Simplex*)P[pn]);E->Enx->S2 = ((C2Simplex*)P[pn]);E->Enx->Enx->S2 = ((C2Simplex*)P[pn]);
		F = E->Fnx->Rev;
		F->S2 = ((C2Simplex*)P[pn]);F->Enx->S2 = ((C2Simplex*)P[pn]);F->Enx->Enx->S2 = ((C2Simplex*)P[pn]);
		F = E->Enx;
		F->S2 = ((C2Simplex*)P[pn]);F->Enx->S2 = ((C2Simplex*)P[pn]);F->Enx->Enx->S2 = ((C2Simplex*)P[pn]);
		
		((C2Simplex*)P[pn])->center[0] = E->Cell->center[0];
		((C2Simplex*)P[pn])->center[1] = E->Cell->center[1];
		((C2Simplex*)P[pn])->alphavalue = E->Cell->alphavalue;
		((C2Simplex*)P[pn])->geometry();
	}
	for (i=0;i<(int)D->Es.size();i++)
	{
		E = D->Es[i];
		if (E->Flag&2) continue;
		if (E->A->N<0 || E->B->N<0) continue;
		if (E->S1) continue;
		int pn = (int)P.size();
		C1Simplex *S1 = new C1Simplex;
		P.push_back(S1);
		S1 = NULL;
		delete S1;
		((C1Simplex*)P[pn])->E = E;
		((C1Simplex*)P[pn])->A = E->A;
		((C1Simplex*)P[pn])->a = E->A->N;
		((C1Simplex*)P[pn])->B = E->B;
		((C1Simplex*)P[pn])->b = E->B->N;
		E->S1 = ((C1Simplex*)P[pn]);E->Rev->S1 = ((C1Simplex*)P[pn]);

		((C1Simplex*)P[pn])->center[0] = E->center[0];
		((C1Simplex*)P[pn])->center[1] = E->center[1];
		((C1Simplex*)P[pn])->alphavalue = E->alphavalue;
		((C1Simplex*)P[pn])->geometry();
	}
}

CComplex::~CComplex()
{
	int i;
	for (i=0;i<(int)P.size();i++) {delete P[i];P[i] = NULL;}
	P.clear();
}

CXYZW* C0Simplex::operator [](int i)
{
	switch(i)
	{
		case 0: return A;
	}
	return 0;
}

CXYZW* C1Simplex::operator [](int i)
{
	switch(i)
	{
		case 0: return A;
		case 1: return B;
	}
	return 0;
}

CXYZW* C2Simplex::operator [](int i)
{
	switch(i)
	{
		case 0: return A;
		case 1: return B;
		case 2: return C;
	}
	return 0;
}

void C1Simplex::geometry()
{
	double x1 = E->A->X;
	double y1 = E->A->Y;
	double w1 = E->A->W;
	int n1 = E->A->N;
	int n2 = E->B->N;
	
	double temp1 = (center[0]-x1)*(center[0]-x1)+(center[1]-y1)*(center[1]-y1)-w1*w1;
	if (temp1<=0)
	{
		alphavalue = max_val2(-sqrt(-temp1),alphavalue);
	}
	else
	{
		alphavalue = max_val2(sqrt(temp1),alphavalue);
	}
}

void C2Simplex::geometry()
{
	double temp1=-100000,temp2,temp3,temp4;
	
	double x1 = E->A->X;
	double y1 = E->A->Y;
	double w1 = E->A->W;
	double x2 = E->B->X;
	double y2 = E->B->Y;
	double w2 = E->B->W;
	double x3 = E->Fnx->B->X;
	double y3 = E->Fnx->B->Y;
	double w3 = E->Fnx->B->W;

	int n1 = E->A->N;
	int n2 = E->B->N;
	int n3 = E->Fnx->B->N;

	double dist_1_12 = (E->center[0]-x1)*(E->center[0]-x1)+(E->center[1]-y1)*(E->center[1]-y1);
	double dist_3_12 = (E->center[0]-x3)*(E->center[0]-x3)+(E->center[1]-y3)*(E->center[1]-y3);
	double dist_1_13 = (E->Fnx->center[0]-x1)*(E->Fnx->center[0]-x1)+(E->Fnx->center[1]-y1)*(E->Fnx->center[1]-y1);
	double dist_2_13 = (E->Fnx->center[0]-x2)*(E->Fnx->center[0]-x2)+(E->Fnx->center[1]-y2)*(E->Fnx->center[1]-y2);
	double dist_1_23 = (E->Enx->center[0]-x1)*(E->Enx->center[0]-x1)+(E->Enx->center[1]-y1)*(E->Enx->center[1]-y1);
	double dist_2_23 = (E->Enx->center[0]-x2)*(E->Enx->center[0]-x2)+(E->Enx->center[1]-y2)*(E->Enx->center[1]-y2);

	double temp11 = (center[0]-x1)*(center[0]-x1)+(center[1]-y1)*(center[1]-y1) - w1*w1; if (temp11>temp1) {temp1 = temp11;}
	double temp12 = (center[0]-x2)*(center[0]-x2)+(center[1]-y2)*(center[1]-y2) - w2*w2; if (temp12>temp1) {temp1 = temp12;}
	double temp13 = (center[0]-x3)*(center[0]-x3)+(center[1]-y3)*(center[1]-y3) - w3*w3; if (temp13>temp1) {temp1 = temp13;}
	temp2 = dist_1_12 - w1*w1 - dist_3_12 + w3*w3;
	temp3 = dist_1_13 - w1*w1 - dist_2_13 + w2*w2;
	temp4 = dist_2_23 - w2*w2 - dist_1_23 + w1*w1;

	if (temp1<=0)
	{
		alphavalue = max_val2(-sqrt(-temp1),alphavalue);
	}
	else if (temp1>0)
	{
		alphavalue = max_val2(sqrt(temp1),alphavalue);
	}
	E->Cell->alphavalue = alphavalue;
	if (temp2>0)
	{
		E->alphavalue = max_val3(alphavalue,E->alphavalue,E->Rev->alphavalue);
		E->Rev->alphavalue = E->alphavalue;
	}
	if (temp3>0)
	{
		E->Fnx->alphavalue = max_val3(alphavalue,E->Fnx->alphavalue,E->Fnx->Rev->alphavalue);
		E->Fnx->Rev->alphavalue = E->Fnx->alphavalue;
	}
	if (temp4>0)
	{
		E->Enx->alphavalue = max_val3(alphavalue,E->Enx->alphavalue,E->Enx->Rev->alphavalue);
		E->Enx->Rev->alphavalue = E->Enx->alphavalue;
	}
}

#endif

