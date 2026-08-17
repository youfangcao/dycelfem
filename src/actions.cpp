// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/***************************************
***    Project:  Cell Growth         ***
***    File:     actions.cpp         ***
***    Author:   Jieling Zhao        ***
***                                  ***
***    Created on March 19, 2012     ***
***************************************/

#ifndef ACTIONS_CPP
#define ACTIONS_CPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>
#include <algorithm>
#include "dbReader.h"
#include "biology.h"

using namespace std;

void dbReader::cell_division(cell *C, int life, int type_Co, int type_Cn, double div_angle)
{
	int cell_n = (int)cellList.size();
	cell *A = new cell(cell_n);
	C->set_cell_type(type_Co);
	A->set_cell_type(type_Cn);
	double ri = sqrt(C->get_area()/2/PI);
	A->set_r(ri);
	C->set_r(ri);
	cout<<"cell division: "<<C->id()<<", "<<A->id();
	double color_r_C = color_class_R[type_Co];
	double color_g_C = color_class_G[type_Co];
	double color_b_C = color_class_B[type_Co];
	double color_r_A = color_class_R[type_Cn];
	double color_g_A = color_class_G[type_Cn];
	double color_b_A = color_class_B[type_Cn];
	A->set_life(life);
	C->set_life(life);
	A->set_initial_area(C->get_initial_area());
	A->setup_tcoef(tension_edge_class[type_Cn]);
	A->setup_pcoef(pressure_edge_class[type_Cn]);
	A->setup_mcoef(mass_edge_class[type_Cn]);
	A->set_initial_area(AR/2);
	A->setup_tcoef(tension_edge_class[type_Cn]);
	A->setup_pcoef(pressure_edge_class[type_Cn]);
	A->setup_mcoef(mass_edge_class[type_Cn]);
	A->set_lame(0,elastic_MU[type_Cn]);
	A->set_lame(1,elastic_LAMBDA[type_Cn]);
	A->set_GR(growth_rate_class[type_Cn]);
	C->set_GR(growth_rate_class[type_Co]);
	A->set_set_dead(C->Set_Dead());
	C->set_cell_color(color_r_C, color_g_C, color_b_C);
	A->set_cell_color(color_r_A, color_g_A, color_b_A);
	cellList.push_back(A);
	for (int i=0;i<biology->getNumSpecies();i++)
	{
		C->set_cell_state(i,C->get_cell_state(i)/2);
		A->set_cell_state(i,C->get_cell_state(i));
	}
	C->setup_Circle();
	double angle = 0;
	if (div_angle>=0 && div_angle<=360)
	{
		angle = div_angle;
	}
	else
	{
		angle = C->shape_based_division_angle();
	}
	C->center_refresh();
	double xc = C->get_center(0);
	double yc = C->get_center(1);
	point3D *s = NULL,*t = NULL,*e = NULL,*f = NULL,*g = NULL,*h = NULL,*gt = NULL,*ht = NULL,*p = NULL;
	/************************
	           \
	          o-\*----o
		      |  \     \
	    o-----o   \     o
		 \         O     \
		  o-----o   \     o
		         \   \    |
		          o---\*--o
		               \
     to divide the cell into 
	 two daughtors
	************************/
	if      (angle==90 || angle==270)
	{
		vector<point3D*> pp;
		int sn = (int)C->get_sides()->size();
		for (int i=0;i<sn;i++)
		{
			double x1 = C->get_sides()->at(i)->p1()->x();
			double x2 = C->get_sides()->at(i)->p2()->x();
			double t1 = x1 - xc;
			double t2 = x2 - xc;
			if ((t1>0 && t2<0) || (t1<0 && t2>0))
			{
				if (abs(t1)<=abs(t2)) {pp.push_back(C->get_sides()->at(i)->p1());}
				else                  {pp.push_back(C->get_sides()->at(i)->p2());}
			}
			else if (t1==0) {pp.push_back(C->get_sides()->at(i)->p1());}
		}
		if (pp.size()==2)
		{
			s = pp[0];
			t = pp[1];
		}
		pp.clear();
	}
	else
	{
		double slope = tan(angle*PI/180);
		double intercept = yc - slope*xc;
		vector<point3D*> pp; 
		for (int i=0;i<(int)C->get_sides()->size();i++)
		{
			double x1 = C->get_sides()->at(i)->p1()->x();
			double y1 = C->get_sides()->at(i)->p1()->y();
			double x2 = C->get_sides()->at(i)->p2()->x();
			double y2 = C->get_sides()->at(i)->p2()->y();
			double t1 = y1 - slope*x1 - intercept;
			double t2 = y2 - slope*x2 - intercept;
			if ((t1>0 && t2<0) || (t1<0 && t2>0))
			{
				if (abs(t1)<=abs(t2)) {pp.push_back(C->get_sides()->at(i)->p1());}
				else				  {pp.push_back(C->get_sides()->at(i)->p2());}
			}
		}
		if (pp.size()==2)
		{
			s = pp[0];
			t = pp[1];
		}
		pp.clear();
	}
	/***************************************
	 points on division wall should be 
	 avoided 
	***************************************/
	int s_nei_n = (int)s->getmp()->size();
	int t_nei_n = (int)t->getmp()->size();
	while (s_nei_n>=2 || t_nei_n>=2)
	{
		s = s->getfp();
		t = t->getfp();
		s_nei_n = (int)s->getmp()->size();
		t_nei_n = (int)t->getmp()->size();
	}
	/***************************************
	  update vertex direction relationship 
	***************************************/
	p = t;
	while (p!=s)
	{
		C->removeside(p->getab());
		p->getab()->reset_Cell(cellList[cell_n]);
		cellList[cell_n]->pushside(p->getab());
		p = p->getfp();
	}
	/*******************************
	           s  V1=e
	       o-<--o*--<-o
	      /     ||     \
         o      **      o
		 |  C   ||   A  |
		 o   V3 ** V4   o
		  \     ||     /
		   o->--o*-->-o
		       t  V2=f
	*******************************/
	int vertex_n = (int)vertexList.size();
	int edge_n = (int)edgeList.size();
	point3D *V1 = new point3D(vertex_n,s->x(),s->y());
	point3D *V2 = new point3D(vertex_n+1,t->x(),t->y());
	vertexList.push_back(V1);
	vertexList.push_back(V2);
	e = vertexList[vertex_n];
	f = vertexList[vertex_n+1];
	/////////////////////////////////
	if ((int)s->getmp()->size()==0)
	{
		s->pushmp(e);
		e->pushmp(s);
		s->get_node_id()->pushvertex(e);
		e->set_node_id(s->get_node_id());
	}
	else
	{
		int smn = (int)s->getmp()->size();
		for (int i=0;i<smn;i++)
		{
			s->getmp()->at(i)->pushmp(e);
			e->pushmp(s->getmp()->at(i));
		}
		s->pushmp(e);
		e->pushmp(s);
		s->get_node_id()->pushvertex(e);
		e->set_node_id(s->get_node_id());
	}
	if (t->getmp()->size()==0)
	{
		t->pushmp(f);
		f->pushmp(t);
		t->get_node_id()->pushvertex(f);
		f->set_node_id(t->get_node_id());
	}
	else
	{
		int tmn = (int)t->getmp()->size();
		for (int i=0;i<tmn;i++)
		{
			t->getmp()->at(i)->pushmp(f);
			f->pushmp(t->getmp()->at(i));
		}
		t->pushmp(f);
		f->pushmp(t);
		t->get_node_id()->pushvertex(f);
		f->set_node_id(t->get_node_id());
	}
	/*****************************
	 update neighbor relationship
	*****************************/
	/** just for use once **/
	e->setba(s->getba());
	f->setab(t->getab());
	//////////////////////////////
	// add pairlist index
	//////////////////////////////
	if (C->get_pairs()->size()>0)
	{
		int cpn = (int)C->get_pairs()->size();
		for (int i=0;i<cpn;i++)
		{
			interpair *Cp = C->get_pairs()->at(i);
			if (Cp->get_redundant()) continue;
			cell *Cp1 = Cp->get_I1();
			cell *Cp2 = Cp->get_I2();
			if      (Cp1==C)
			{
				int cpmn = (int)Cp->MP()->size();
				if (cpmn==0)
				{
				}
				else
				{
					for (int j=0;j<cpmn;j++)
					{
						mergepair *CpMp = Cp->MP()->at(j);
						point3D *s1 = CpMp->S1();
						point3D *t1 = CpMp->T1();
						point3D *s2 = CpMp->S2();
						point3D *t2 = CpMp->T2();
						point3D *s1t = CpMp->S1();
						point3D *t1t = CpMp->T1();
						if (s1==t) {s1 = f;}
						if (t1==s) {t1 = e;}
						/**************************************
						          |         |
								  o         o
								  |  C      |
							   s1 o---->----o t1 
						o---<----oo----<----oo----<-----o 
						|	     || t2    s2||          |
						o		 oo         oo          o
						**************************************/
						if (s2==t2)
						{
							if      (s1==f)
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp2);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,f,f,s2,s2);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp2->pushpair(collisionpairList[pair_n1]);
							}
							else if (t1==e)
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp2);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,e,e,t2,t2);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp2->pushpair(collisionpairList[pair_n1]);
							}
							else
							{
								if (s1->getab()->Cell()==cellList[cell_n] && t1->getba()->Cell()==cellList[cell_n] && (s2->getba()->get_attach() && t2->getab()->get_attach()))
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s1,t1,s2,t2);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
									Cp->remove_MP(CpMp);
									if ((int)Cp->MP()->size()==0) {Cp->set_redundant(1);}
								}
								else if (s1->getab()->Cell()==cellList[cell_n] && t1->getba()->Cell()==cellList[cell_n] && (!t2->getab()->get_attach() || !s2->getba()->get_attach()))
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(0,s1,t1,s2,t2);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
									Cp->remove_MP(CpMp);
								}
							}
						}
						else
						{
							if      (s1->getab()->Cell()==C && t1->getba()->Cell()==C && (s2->getba()->get_attach() && t2->getab()->get_attach()))
							{
								if      (s1t==s || t1t==t)
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									if      (s1t==s) 
									{
										mergepair *Mn = new mergepair(1,e,e,t2,t2);
										collisionpairList[pair_n1]->push_MP(Mn);
									}
									else if (t1t==t) 
									{
										mergepair *Mn = new mergepair(1,f,f,s2,s2);
										collisionpairList[pair_n1]->push_MP(Mn);
									}
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
								}
							}
							else if (s1->getab()->Cell()==C && t1->getba()->Cell()==C && (!t2->getab()->get_attach() || !s2->getba()->get_attach()))
							{
								if      (s1t==s || t1t==t)
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									if (s1t==s)
									{
										mergepair *Mn = new mergepair(1,e,e,t2,t2);
										collisionpairList[pair_n1]->push_MP(Mn);
									}
									else if (t1t==t)
									{
										mergepair *Mn = new mergepair(1,f,f,s2,s2);
										collisionpairList[pair_n1]->push_MP(Mn);
									}
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
								}
							}
							else if (s1->getab()->Cell()==cellList[cell_n] && t1->getba()->Cell()==cellList[cell_n] && (s2->getba()->get_attach() && t2->getab()->get_attach()))
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp2);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,s1,t1,s2,t2);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp2->pushpair(collisionpairList[pair_n1]);
								if (s1t->getfp()==t1t) 
								{
									s1t->getab()->set_Pair(collisionpairList[pair_n1]);
									s1t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
								}
								else
								{
									p = s1t;
									while (p!=t1t)
									{
										p->getab()->set_Pair(collisionpairList[pair_n1]);
										p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
										p = p->getfp();
									}
								}
								if (s1t==t)
								{
									CpMp->setS1(t);
									CpMp->setT1(t);
									CpMp->setS2(t2);
								}
								else if (t1t==s)
								{
									CpMp->setS1(s);
									CpMp->setT1(s);
									CpMp->setT2(s2);
								}
								else
								{
									Cp->remove_MP(CpMp);
								}
								if ((int)Cp->MP()->size()==0) {Cp->set_redundant(1);}
							}
							else if (s1->getab()->Cell()==cellList[cell_n] && t1->getba()->Cell()==cellList[cell_n] && (!t2->getab()->get_attach() || !s2->getba()->get_attach()))
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp2);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(0,s1,t1,s2,t2);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp2->pushpair(collisionpairList[pair_n1]);
								if (s1t->getfp()==t1t) 
								{
									s1t->getab()->set_Pair(collisionpairList[pair_n1]);
									s1t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
								}
								else
								{
									p = s1t;
									while (p!=t1t)
									{
										p->getab()->set_Pair(collisionpairList[pair_n1]);
										p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
										p = p->getfp();
									}
								}
								if (s1t==t)
								{
									CpMp->setS1(t);
									CpMp->setT1(t);
									CpMp->setS2(t2);
								}
								else if (t1t==s)
								{
									CpMp->setS1(s);
									CpMp->setT1(s);
									CpMp->setT2(s2);
								}
								else
								{
									Cp->remove_MP(CpMp);
								}
							}
							else if (s1->getab()->Cell()==C && t1->getba()->Cell()==cellList[cell_n])
							{
								int flagst = 0; // s = 0;t = 1;
								p = s1t;
								while (p!=t1t)
								{
									if (p==t) {flagst = 1;break;}
									p = p->getfp();
								}
								if (flagst==0) // s
								{
									CpMp->setT1(e);
									for (int k=0;k<(int)s->getmp()->size();k++)
									{
										if (s->getmp()->at(k)->getba()->Cell()==Cp2)
										{
											p = s->getmp()->at(k);break;
										}
									}
									CpMp->setS2(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(0,s,t1,s2,p);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
									if (s->getfp()==t1t)
									{
										s->getab()->set_Pair(collisionpairList[pair_n1]);
										s->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = s;
										while (p!=t1t)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
								else // t
								{
									CpMp->setT1(t);
									for (int k=0;k<(int)t->getmp()->size();k++)
									{
										if (t->getmp()->at(k)->getab()->Cell()==Cp2)
										{
											p = t->getmp()->at(k);break;
										}
									}
									CpMp->setS2(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(0,f,t1,s2,p);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
									if (t->getfp()==t1t)
									{
										t->getab()->set_Pair(collisionpairList[pair_n1]);
										t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = t;
										while (p!=t1t)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
							}
							else if (s1->getab()->Cell()==cellList[cell_n] && t1->getba()->Cell()==C)
							{
								int flagst = 0; // s = 0; t = 1; 
								p = s1t;
								while (p!=t1t)
								{
									if (p==t) {flagst = 1;break;}
									p = p->getfp();
								}
								if (flagst==0) // s
								{
									CpMp->setS1(s);
									int smn1 = (int)s->getmp()->size();
									for (int k=0;k<smn1;k++)
									{
										if (s->getmp()->at(k)->getba()->Cell()==Cp2)
										{
											p = s->getmp()->at(k);break;
										}
									}
									CpMp->setT2(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s1,e,p,t2);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
									if (s1t->getfp()==s)
									{
										s1t->getab()->set_Pair(collisionpairList[pair_n1]);
										s1t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = s1t;
										while (p!=s)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
								else // t
								{
									CpMp->setS1(f);
									int tmn1 = (int)t->getmp()->size();
									for (int k=0;k<tmn1;k++)
									{
										if (t->getmp()->at(k)->getab()->Cell()==Cp2)
										{
											p = t->getmp()->at(k);break;
										}
									}
									CpMp->setT2(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp2);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp2);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp2->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s1,t,p,t2);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp2->pushpair(collisionpairList[pair_n1]);
									if (s1t->getfp()==t)
									{
										s1t->getab()->set_Pair(collisionpairList[pair_n1]);
										s1t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = s1t;
										while (p!=t)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
							}
						}
					}
				}
			}
			else if (Cp2==C)
			{
				int cpmn = (int)C->get_pairs()->at(i)->MP()->size();
				if (cpmn==0) {}
				else
				{
					for (int j=0;j<cpmn;j++)
					{
						mergepair *CpMp = Cp->MP()->at(j);
						point3D *s1 = CpMp->S1();
						point3D *t1 = CpMp->T1();
						point3D *s2 = CpMp->S2();
						point3D *t2 = CpMp->T2();
						point3D *s2t = CpMp->S2();
						point3D *t2t = CpMp->T2();
						if (t2==s) {t2 = e;}
						if (s2==t) {s2 = f;}
						/*********************************************/
						if (s1==t1)
						{
							if      (s2==f)
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp1);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,f,f,s1,s1);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp1->pushpair(collisionpairList[pair_n1]);
							}
							else if (t2==e)
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp1);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,e,e,t1,t1);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp1->pushpair(collisionpairList[pair_n1]);
							}
							else
							{
								if (s2->getab()->Cell()==cellList[cell_n] && t2->getba()->Cell()==cellList[cell_n] && (s1->getba()->get_attach() && t1->getab()->get_attach()))
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s2,t2,s1,t1);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
									Cp->remove_MP(CpMp);
									if ((int)Cp->MP()->size()==0) 
									{
										Cp->set_redundant(1);
									}
								}
								else if (s2->getab()->Cell()==cellList[cell_n] && t2->getba()->Cell()==cellList[cell_n] && (!t1->getab()->get_attach() || !s1->getba()->get_attach()))
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s2,t2,s1,t1);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
									Cp->remove_MP(CpMp);
								}
							}
						}
						else
						{
							if      (t2->getba()->Cell()==C && s2->getab()->Cell()==C && (t1->getab()->get_attach() && s1->getba()->get_attach()))
							{
								if (s2t==s || t2t==t)
								{
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									if      (s2t==s) 
									{
										mergepair *Mn = new mergepair(1,e,e,t1,t1);
										collisionpairList[pair_n1]->push_MP(Mn);
									}
									else if (t2t==t)
									{
										mergepair *Mn = new mergepair(1,f,f,s1,s1);
										collisionpairList[pair_n1]->push_MP(Mn);
									}
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
								}
							}
							else if (t2->getba()->Cell()==C && s2->getab()->Cell()==C && (!t1->getab()->get_attach() || !s1->getba()->get_attach()))
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp1);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								if (s2t==s)
								{
									mergepair *Mn = new mergepair(1,e,e,t1,t1);
									collisionpairList[pair_n1]->push_MP(Mn);
								}
								else if (t2t==t)
								{
									mergepair *Mn = new mergepair(1,f,f,s1,s1);
									collisionpairList[pair_n1]->push_MP(Mn);
								}
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp1->pushpair(collisionpairList[pair_n1]);
							}
							else if (s2->getab()->Cell()==cellList[cell_n] && t2->getba()->Cell()==cellList[cell_n] && (s1->getba()->get_attach() && t1->getab()->get_attach()))
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp1);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,s2,t2,s1,t1);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp1->pushpair(collisionpairList[pair_n1]);
								if (s2t->getfp()==t2t)
								{
									s2t->getab()->set_Pair(collisionpairList[pair_n1]);
									s2t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
								}
								else
								{
									p = s2t;
									while (p!=t2t)
									{
										p->getab()->set_Pair(collisionpairList[pair_n1]);
										p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
										p = p->getfp();
									}
								}
								if      (s2t==t)
								{
									CpMp->setS1(t1);
									CpMp->setS2(t);
									CpMp->setT2(t);
								}
								else if (t2t==s)
								{
									CpMp->setT1(s1);
									CpMp->setS2(s);
									CpMp->setT2(s);
								}
								else
								{
									Cp->remove_MP(CpMp);
								}
								if ((int)Cp->MP()->size()==0) 
								{
									Cp->set_redundant(1);
								}
							}
							else if (s2->getab()->Cell()==cellList[cell_n] && t2->getba()->Cell()==cellList[cell_n] && (!t1->getab()->get_attach() || !s1->getba()->get_attach()))
							{
								int pair_n1 = (int)collisionpairList.size();
								interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
								collisionpairList.push_back(Pn);
								//int iptype = interpair_type(cellList[cell_n],Cp1);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
								double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
								collisionpairList[pair_n1]->set_adhesion(adhesion_force);
								collisionpairList[pair_n1]->set_friction(adhesion_force);
								mergepair *Mn = new mergepair(1,s2,t2,s1,t1);
								collisionpairList[pair_n1]->push_MP(Mn);
								//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
								cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
								Cp1->pushpair(collisionpairList[pair_n1]);
								if (s2t->getfp()==t2t)
								{
									s2t->getab()->set_Pair(collisionpairList[pair_n1]);
									s2t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
								}
								else
								{
									p = s2t;
									while (p!=t2t)
									{
										p->getab()->set_Pair(collisionpairList[pair_n1]);
										p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
										p = p->getfp();
									}
								}
								if      (s2t==t)
								{
									CpMp->setS1(t1);
									CpMp->setS2(t);
									CpMp->setT2(t);
								}
								else if (t2t==s)
								{
									CpMp->setT1(s1);
									CpMp->setS2(s);
									CpMp->setT2(s);
								}
								else 
								{
									Cp->remove_MP(Cp->MP()->at(j));
								}
							}
							else if (s2->getab()->Cell()==C && t2->getba()->Cell()==cellList[cell_n])
							{
								int flagst = 0; // s = 0;t = 1;
								p = s2t;
								while (p!=t2t)
								{
									if (p==t) {flagst = 1;break;}
									p = p->getfp();
								}
								if (flagst==0) // s
								{
									CpMp->setT2(e);
									int smn2 = (int)s->getmp()->size();
									for (int k=0;k<smn2;k++)
									{
										if (s->getmp()->at(k)->getba()->Cell()==Cp1)
										{
											p = s->getmp()->at(k);break;
										}
									}
									Cp->MP()->at(j)->setS1(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s,t2,s1,p);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
									if (s->getfp()==t2t)
									{
										s->getab()->set_Pair(collisionpairList[pair_n1]);
										s->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = s;
										while (p!=t2t)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
								else // t
								{
									CpMp->setT2(t);
									int tmn2 = (int)t->getmp()->size();
									for (int k=0;k<tmn2;k++)
									{
										if (t->getmp()->at(k)->getab()->Cell()==Cp1)
										{
											p = t->getmp()->at(k);break;
										}
									}
									CpMp->setS1(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,f,t2,s1,p);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
									if (t->getfp()==t2t)
									{
										t->getab()->set_Pair(collisionpairList[pair_n1]);
										t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = t;
										while (p!=t2t)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
							}
							else if (s2->getab()->Cell()==cellList[cell_n] && t2->getba()->Cell()==C)
							{
								int flagst = 0; // s = 0; t = 1; 
								p = s2t;
								while (p!=t2t)
								{
									if (p==t) {flagst = 1;break;}
									p = p->getfp();
								}
								if (flagst==0) // s
								{
									CpMp->setS2(s);
									int smn2 = (int)s->getmp()->size();
									for (int k=0;k<smn2;k++)
									{
										if (s->getmp()->at(k)->getba()->Cell()==Cp1)
										{
											p = s->getmp()->at(k);break;
										}
									}
									CpMp->setT1(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s2,e,p,t1);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
									if (s2t->getfp()==s)
									{
										s2t->getab()->set_Pair(collisionpairList[pair_n1]);
										s2t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = s2t;
										while (p!=s)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
								else // t
								{
									CpMp->setS2(f);
									int tmn2 = (int)t->getmp()->size();
									for (int k=0;k<tmn2;k++)
									{
										if (t->getmp()->at(k)->getab()->Cell()==Cp1)
										{
											p = t->getmp()->at(k);break;
										}
									}
									CpMp->setT1(p);
									int pair_n1 = (int)collisionpairList.size();
									interpair *Pn = new interpair(pair_n1,cellList[cell_n],Cp1);
									collisionpairList.push_back(Pn);
									//int iptype = interpair_type(cellList[cell_n],Cp1);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[iptype]);
									double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][Cp1->get_cell_type()];
									collisionpairList[pair_n1]->set_adhesion(adhesion_force);
									collisionpairList[pair_n1]->set_friction(adhesion_force);
									mergepair *Mn = new mergepair(1,s2,t,p,t1);
									collisionpairList[pair_n1]->push_MP(Mn);
									//collisionpairList[pair_n1]->set_adhesion(adhesion_pair_class[0]);
									cellList[cell_n]->pushpair(collisionpairList[pair_n1]);
									Cp1->pushpair(collisionpairList[pair_n1]);
									if (s2t->getfp()==t)
									{
										s2t->getab()->set_Pair(collisionpairList[pair_n1]);
										s2t->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
									}
									else
									{
										p = p->getfp();
										while (p!=t)
										{
											p->getab()->set_Pair(collisionpairList[pair_n1]);
											p->getab()->get_Neighbor()->set_Pair(collisionpairList[pair_n1]);
											p = p->getfp();
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	/********* add one connection between the two cells **********/
	int pair_n = (int)collisionpairList.size();
	interpair *Pl = new interpair(pair_n,C,cellList[cell_n]);
	collisionpairList.push_back(Pl);
	//int iptype = interpair_type(C,cellList[cell_n]);
	//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
	double adhesion_force = adhesion_pair_array[cellList[cell_n]->get_cell_type()][C->get_cell_type()];
	collisionpairList[pair_n]->set_adhesion(adhesion_force);
	collisionpairList[pair_n]->set_friction(adhesion_force);
	/*****************************
	  update the partition wall
	*****************************/
	t->getfp()->setrp(f);
	f->setfp(t->getfp());
	s->getrp()->setfp(e);
	e->setrp(s->getrp());
	f->setab(f->getfp()->getba());
	f->getab()->reset_p1(f);
	e->setba(e->getrp()->getab());
	e->getba()->reset_p2(e);
	double BD = 2*PI*CR/SN;
	double dist = sqrt((s->x()-t->x())*(s->x()-t->x())+(s->y()-t->y())*(s->y()-t->y()));
	int num = (int)(dist/BD);
	if      (num<=1)
	{
		edge *E1 = new edge(edge_n,t,s,C);
		edge *E2 = new edge(edge_n+1,e,f,cellList[cell_n]);
		edgeList.push_back(E1);
		edgeList.push_back(E2);
		E1->set_Neighbor(E2);
		E2->set_Neighbor(E1);
		E1->set_attach(1);
		E2->set_attach(1);
		E1->set_Pair(collisionpairList[pair_n]);
		E2->set_Pair(collisionpairList[pair_n]);
		t->setab(E1);
		s->setba(E1);
		e->setab(E2);
		f->setba(E2);
		t->getba()->Cell()->pushside(E1);
		f->getab()->Cell()->pushside(E2);
		C->pushside(E1);
		cellList[cell_n]->pushside(E2);
		t->setfp(s);
		s->setrp(t);
		e->setfp(f);
		f->setrp(e);
	}
	else if (num>1)
	{
		double lambda = 1.0/(num);
		vertex_n = (int)vertexList.size();
		vertex_n++;
		double xk = lambda*s->x() + (1 - lambda)*t->x();
		double yk = lambda*s->y() + (1 - lambda)*t->y();
		point3D *V3 = new point3D(vertex_n-1,xk,yk);
		vertexList.push_back(V3);

		vertex_n++;
		point3D *V4 = new point3D(vertex_n-1,xk,yk);
		vertexList.push_back(V4);
		g = vertexList[vertex_n-2];
		h = vertexList[vertex_n-1];
		g->pushmp(h);
		h->pushmp(g);
		int node_n = (int)nodeList.size();
		node *N3 = new node(node_n);
		nodeList.push_back(N3);
		nodeList[node_n]->pushvertex(g);
		nodeList[node_n]->pushvertex(h);
		g->set_node_id(nodeList[node_n]);
		h->set_node_id(nodeList[node_n]);

		edge_n++;
		edge *E1 = new edge(edge_n-1,t,g,C);
		edge_n++;
		edge *E2 = new edge(edge_n-1,h,f,cellList[cell_n]);
		edgeList.push_back(E1);
		edgeList.push_back(E2);
		E1->set_Neighbor(E2);
		E2->set_Neighbor(E1);
		E1->set_attach(1);
		E2->set_attach(1);
		E1->set_Pair(collisionpairList[pair_n]);
		E2->set_Pair(collisionpairList[pair_n]);
		t->setfp(g);
		g->setrp(t);
		h->setfp(f);
		f->setrp(h);
		t->setab(E1);
		g->setba(E1);
		f->setba(E2);
		h->setab(E2);
		C->pushside(E1);
		cellList[cell_n]->pushside(E2);
	
		if (num>2)
		{
			for (int i=2;i<num;i++)
			{
				gt = g;
				ht = h;
				vertex_n++;
				xk = lambda*i*s->x() + (1 - lambda*i)*t->x();
				yk = lambda*i*s->y() + (1 - lambda*i)*t->y();
				point3D *V5 = new point3D(vertex_n-1,xk,yk);
				vertexList.push_back(V5);

				vertex_n++;
				point3D *V6 = new point3D(vertex_n-1,xk,yk);
				vertexList.push_back(V6);
				g = vertexList[vertex_n-2];
				h = vertexList[vertex_n-1];
				g->pushmp(h);
				h->pushmp(g);
				node_n++;
				node *N4 = new node(node_n);
				nodeList.push_back(N4);
				nodeList[node_n]->pushvertex(g);
				nodeList[node_n]->pushvertex(h);
				g->set_node_id(nodeList[node_n]);
				h->set_node_id(nodeList[node_n]);
				gt->setfp(g);
				g->setrp(gt);
				ht->setrp(h);
				h->setfp(ht);
				edge_n++;
				edge *E3 = new edge(edge_n-1,gt,g,C);
				edge_n++;
				edge *E4 = new edge(edge_n-1,h,ht,cellList[cell_n]);
				edgeList.push_back(E3);
				edgeList.push_back(E4);
				E3->set_Neighbor(E4);
				E4->set_Neighbor(E3);
				E3->set_attach(1);
				E4->set_attach(1);
				E3->set_Pair(collisionpairList[pair_n]);
				E4->set_Pair(collisionpairList[pair_n]);
				gt->setab(E3);
				g->setba(E3);
				h->setab(E4);
				ht->setba(E4);
				C->pushside(E3);
				cellList[cell_n]->pushside(E4);
			}
			edge_n++;
			edge *E5 = new edge(edge_n-1,g,s,C);
			edge_n++;
			edge *E6 = new edge(edge_n-1,e,h,cellList[cell_n]);
			edgeList.push_back(E5);
			edgeList.push_back(E6);
			E5->set_Neighbor(E6);
			E6->set_Neighbor(E5);
			E5->set_attach(1);
			E6->set_attach(1);
			E5->set_Pair(collisionpairList[pair_n]);
			E6->set_Pair(collisionpairList[pair_n]);
			g->setfp(s);
			s->setrp(g);
			h->setrp(e);
			e->setfp(h);
			g->setab(E5);
			s->setba(E5);
			h->setba(E6);
			e->setab(E6);
			C->pushside(E5);
			cellList[cell_n]->pushside(E6);
		}
		else
		{
			edge_n++;
			edge *E3 = new edge(edge_n-1,g,s,C);
			edge_n++;
			edge *E4 = new edge(edge_n-1,e,h,cellList[cell_n]);
			edgeList.push_back(E3);
			edgeList.push_back(E4);
			E3->set_Neighbor(E4);
			E4->set_Neighbor(E3);
			E3->set_attach(1);
			E4->set_attach(1);
			E3->set_Pair(collisionpairList[pair_n]);
			E4->set_Pair(collisionpairList[pair_n]);
			g->setfp(s);
			s->setrp(g);
			h->setrp(e);
			e->setfp(h);
			g->setab(E3);
			s->setba(E3);
			h->setba(E4);
			e->setab(E4);
			C->pushside(E3);
			cellList[cell_n]->pushside(E4);
		}
	}
	mergepair *Ml = new mergepair(1,t,s,e,f);
	collisionpairList[pair_n]->push_MP(Ml);
	collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[0]);
	C->pushpair(collisionpairList[pair_n]);
	cellList[cell_n]->pushpair(collisionpairList[pair_n]);
	cellList[cell_n]->refresh_AABB();
	C->refresh_AABB();
	cell_interior_resample(C,life);
	cell_interior_resample(cellList[cell_n],life);
	vector<cell*> neighbor_c;
	int cpn = (int)C->get_pairs()->size();
	for (int i=0;i<cpn;i++)
	{
		interpair *Cp = C->get_pairs()->at(i);
		if (Cp->get_redundant()) continue;
		if (C->get_pairs()->at(i)->get_I1()==C)
		{
			neighbor_c.push_back(Cp->get_I2());
		}
		else
		{
			neighbor_c.push_back(Cp->get_I1());
		}
	}
	double range[4]; // 0:xmin 1:ymin 2:xmax 3:ymax
	range[0] = 1000000000000;
	range[1] = 1000000000000;
	range[2] =-1000000000000;
	range[3] =-1000000000000;

	int ncn = (int)neighbor_c.size();
	for (int i=0;i<ncn;i++)
	{
		int ncpn = (int)neighbor_c[i]->get_inner_p()->size();
		for (int j=0;j<ncpn;j++)
		{
			double vx = neighbor_c[i]->get_inner_p()->at(j)->x();
			double vy = neighbor_c[i]->get_inner_p()->at(j)->y();
			if (range[0]>vx) {range[0] = vx;}
			if (range[1]>vy) {range[1] = vy;}
			if (range[2]<vx) {range[2] = vx;}
			if (range[3]<vy) {range[3] = vy;}
		}
	}
	int Cpn1 = (int)C->get_inner_p()->size();
	for (int i=0;i<Cpn1;i++)
	{
		double vx = C->get_inner_p()->at(i)->x();
		double vy = C->get_inner_p()->at(i)->y();
		if (range[0]>vx) {range[0] = vx;}
		if (range[1]>vy) {range[1] = vy;}
		if (range[2]<vx) {range[2] = vx;}
		if (range[3]<vy) {range[3] = vy;}
	}
	range[0] -= 1000;
	range[1] -= 1000;
	range[2] += 1000;
	range[3] += 1000;

	CD3DW *Del = new CD3DW(range);
	int c_in = 0;
	for (int i=0;i<ncn;i++)
	{
		int ncpn = (int)neighbor_c[i]->get_inner_p()->size();
		if (ncpn>0)
		{
			for (int j=0;j<ncpn;j++)
			{
				Del->add_resample(neighbor_c[i]->get_inner_p()->at(j)->x(),
					neighbor_c[i]->get_inner_p()->at(j)->y(),
					BD*BD,
					c_in++,
					neighbor_c[i]->get_inner_p()->at(j));
				c_in++;
			}
		}
		else
		{
			int ncpsn = (int)neighbor_c[i]->get_inner_p_slip()->size();
			for (int j=0;j<ncpsn;j++)
			{
				Del->add_resample(neighbor_c[i]->get_inner_p_slip()->at(j)->x(),
					neighbor_c[i]->get_inner_p_slip()->at(j)->y(),
					BD*BD,
					c_in++,
					neighbor_c[i]->get_inner_p_slip()->at(j));
				c_in++;
			}
		}
	}
	int cpn2 = (int)C->get_inner_p()->size();
	if (cpn2>0)
	{
		for (int i=0;i<cpn2;i++)
		{
			if (C->get_inner_p()->at(i)->get_burry()) continue;
			Del->add_resample(C->get_inner_p()->at(i)->x(),
				C->get_inner_p()->at(i)->y(),
				BD*BD,
				c_in++,
				C->get_inner_p()->at(i));
			c_in++;
		}
	}
	else
	{
		int cpsn2 = (int)C->get_inner_p_slip()->size();
		for (int i=0;i<cpsn2;i++)
		{
			Del->add_resample(C->get_inner_p_slip()->at(i)->x(),
				C->get_inner_p_slip()->at(i)->y(),
				BD*BD,
				c_in++,
				C->get_inner_p_slip()->at(i));
			c_in++;
		}
	}
	int cs_n1 = (int)Del->get_Cs()->size();
	vector<cell*> nPairi; // the new cell: cellList[cell_n]
	vector<cell*> nPairj; // the new cell: cellList[cell_n]
	vector<cell*> oPairi; // the old cell: C
	vector<cell*> oPairj; // the old cell: C
	for (int i=0;i<cs_n1;i++)
	{
		CCell *Ci = Del->get_Cs()->at(i);
		CXYZW *V0 = Ci->get_VV(0);
		CXYZW *V1 = Ci->get_VV(1);
		CXYZW *V2 = Ci->get_VV(2);
		if (Ci->BN() &&
			V0->n()>=0 &&
			V1->n()>=0 &&
			V2->n()>=0)
		{
			cell *CC0 = V0->P()->Cell();
			cell *CC1 = V1->P()->Cell();
			cell *CC2 = V2->P()->Cell();
			if      (CC0==cellList[cell_n] &&
					 CC1!=cellList[cell_n] &&
					 CC1==CC2)
			{
				bool exist12 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC0 && nPairj[j]==CC1) ||
						(nPairi[j]==CC1 && nPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC0);
					nPairj.push_back(CC1);
				}
			}
			else if (CC1==cellList[cell_n] &&
					 CC0!=cellList[cell_n] &&
					 CC0==CC2)
			{
				bool exist12 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC0 && nPairj[j]==CC1) ||
						(nPairi[j]==CC1 && nPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC1);
					nPairj.push_back(CC0);
				}
			}
			else if (CC2==cellList[cell_n] &&
					 CC0!=cellList[cell_n] &&
					 CC0==CC1)
			{
				bool exist12 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC2 && nPairj[j]==CC0) ||
						(nPairi[j]==CC0 && nPairj[j]==CC2)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC2);
					nPairj.push_back(CC0);
				}
			}
			else if (CC0==cellList[cell_n] &&
					 CC1!=cellList[cell_n] &&
					 CC2!=cellList[cell_n] &&
					 CC1!=CC2)
			{
				bool exist12 = false;
				bool exist13 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC0 && nPairj[j]==CC1) ||
						(nPairi[j]==CC1 && nPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC0);
					nPairj.push_back(CC1);
				}
				nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC0 && nPairj[j]==CC2) ||
						(nPairi[j]==CC2 && nPairj[j]==CC0)) {exist13 = true;break;}
				}
				if (!exist13)
				{
					nPairi.push_back(CC0);
					nPairj.push_back(CC2);
				}
			}
			else if (CC1==cellList[cell_n] &&
					 CC0!=cellList[cell_n] &&
					 CC2!=cellList[cell_n] &&
					 CC0!=CC2)
			{
				bool exist12 = false;
				bool exist13 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC1 && nPairj[j]==CC0) ||
						(nPairi[j]==CC0 && nPairj[j]==CC1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC1);
					nPairj.push_back(CC0);
				}
				nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC1 && nPairj[j]==CC2) ||
						(nPairi[j]==CC2 && nPairj[j]==CC1)) {exist13 = true;break;}
				}
				if (!exist13)
				{
					nPairi.push_back(CC1);
					nPairj.push_back(CC2);
				}
			}
			else if (CC2==cellList[cell_n] &&
					 CC0!=cellList[cell_n] &&
					 CC1!=cellList[cell_n] &&
					 CC0!=CC1)
			{
				bool exist12 = false;
				bool exist13 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC2 && nPairj[j]==CC0) ||
						(nPairi[j]==CC0 && nPairj[j]==CC2)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC2);
					nPairj.push_back(CC0);
				}
				nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC2 && nPairj[j]==CC1) ||
						(nPairi[j]==CC1 && nPairj[j]==CC2)) {exist13 = true;break;}
				}
				if (!exist13)
				{
					nPairi.push_back(CC2);
					nPairj.push_back(CC1);
				}
			}
			else if (CC0==cellList[cell_n] &&
					 CC1==cellList[cell_n] &&
					 CC2!=cellList[cell_n])
			{
				bool exist12 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC0 && nPairj[j]==CC2) ||
						(nPairi[j]==CC2 && nPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC0);
					nPairj.push_back(CC2);
				}
			}
			else if (CC0==cellList[cell_n] &&
					 CC2==cellList[cell_n] &&
					 CC1!=cellList[cell_n])
			{
				bool exist12 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC0 && nPairj[j]==CC1) ||
						(nPairi[j]==CC1 && nPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC0);
					nPairj.push_back(CC1);
				}
			}
			else if (CC1==cellList[cell_n] &&
					 CC2==cellList[cell_n] &&
					 CC0!=cellList[cell_n])
			{
				bool exist12 = false;
				int nii = (int)nPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((nPairi[j]==CC1 && nPairj[j]==CC0) ||
						(nPairi[j]==CC0 && nPairj[j]==CC1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					nPairi.push_back(CC1);
					nPairj.push_back(CC0);
				}
			}
			if      (CC0==C &&
					 CC1!=C &&
					 CC1==CC2)
			{
				bool exist12 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC0 && oPairj[j]==CC1) || 
						(oPairi[j]==CC1 && oPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC0);
					oPairj.push_back(CC1);
				}
			}
			else if (CC1==C &&
					 CC0!=C &&
					 CC0==CC2)
			{
				bool exist12 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC1 && oPairj[j]==CC0) || 
						(oPairi[j]==CC0 && oPairj[j]==CC1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC1);
					oPairj.push_back(CC0);
				}
			}
			else if (CC2==C &&
					 CC0!=C &&
					 CC0==CC1)
			{
				bool exist12 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC2 && oPairj[j]==CC0) || 
						(oPairi[j]==CC0 && oPairj[j]==CC2)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC2);
					oPairj.push_back(CC0);
				}
			}
			else if (CC0==C &&
					 CC1!=C &&
					 CC2!=C &&
					 CC1!=CC2)
			{
				bool exist12 = false;
				bool exist13 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC0 && oPairj[j]==CC1) || 
						(oPairi[j]==CC1 && oPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC0);
					oPairj.push_back(CC1);
				}
				nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC0 && oPairj[j]==CC2) ||
						(oPairi[j]==CC2 && oPairj[j]==CC0)) {exist13 = true;break;}
				}
				if (!exist13)
				{
					oPairi.push_back(CC0);
					oPairj.push_back(CC2);
				}
			}
			else if (CC1==C &&
					 CC0!=C &&
					 CC2!=C &&
					 CC0!=CC2)
			{
				bool exist12 = false;
				bool exist13 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC1 && oPairj[j]==CC0) || 
						(oPairi[j]==CC0 && oPairj[j]==CC1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC1);
					oPairj.push_back(CC0);
				}
				nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC1 && oPairj[j]==CC2) ||
						(oPairi[j]==CC2 && oPairj[j]==CC1)) {exist13 = true;break;}
				}
				if (!exist13)
				{
					oPairi.push_back(CC1);
					oPairj.push_back(CC2);
				}
			}
			else if (CC2==C &&
					 CC0!=C &&
					 CC1!=C &&
					 CC0!=CC1)
			{
				bool exist12 = false;
				bool exist13 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC2 && oPairj[j]==CC0) || 
						(oPairi[j]==CC0 && oPairj[j]==CC2)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC2);
					oPairj.push_back(CC0);
				}
				nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC2 && oPairj[j]==CC1) ||
						(oPairi[j]==CC1 && oPairj[j]==CC2)) {exist13 = true;break;}
				}
				if (!exist13)
				{
					oPairi.push_back(CC2);
					oPairj.push_back(CC1);
				}
			}
			else if (CC0==C &&
					 CC1==C &&
					 CC2!=C)
			{
				bool exist12 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC0 && oPairj[j]==CC2) || 
						(oPairi[j]==CC2 && oPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC0);
					oPairj.push_back(CC2);
				}
			}
			else if (CC0==C &&
					 CC2==C &&
					 CC1!=C)
			{
				bool exist12 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC0 && oPairj[j]==CC1) || 
						(oPairi[j]==CC1 && oPairj[j]==CC0)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC0);
					oPairj.push_back(CC1);
				}
			}
			else if (CC1==C &&
					 CC2==C &&
					 CC0!=C)
			{
				bool exist12 = false;
				int nii = (int)oPairi.size();
				for (int j=0;j<nii;j++)
				{
					if ((oPairi[j]==CC1 && oPairj[j]==CC0) || 
						(oPairi[j]==CC0 && oPairj[j]==CC1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					oPairi.push_back(CC1);
					oPairj.push_back(CC0);
				}
			}
		}
	}
	int npn1 = (int)nPairi.size(); // cellList[cell_n]
	int opn2 = (int)oPairi.size(); // C
	for (int i=0;i<npn1;i++)
	{
		int cnpn5 = (int)cellList[cell_n]->get_pairs()->size();
		bool exist = false;
		for (int j=0;j<cnpn5;j++)
		{
			cell *Cpn1 = cellList[cell_n]->get_pairs()->at(j)->get_I1();
			cell *Cpn2 = cellList[cell_n]->get_pairs()->at(j)->get_I2();
			if ((Cpn1==nPairi[i] && Cpn2==nPairj[i]) ||
				(Cpn2==nPairi[i] && Cpn1==nPairj[i])) {exist = true;break;}
		}
		if (!exist)
		{
			int pair_n = (int)collisionpairList.size();
			interpair *Pl = new interpair(pair_n,nPairi[i],nPairj[i]);
			collisionpairList.push_back(Pl);
			//int iptype = interpair_type(nPairi[i],nPairj[i]);
			//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
			//collisionpairList[pair_n]->set_friction(friction_pair_class[iptype]);
			double adhesion_force = adhesion_pair_array[nPairi[i]->get_cell_type()][nPairj[i]->get_cell_type()];
			collisionpairList[pair_n]->set_adhesion(adhesion_force);
			collisionpairList[pair_n]->set_friction(adhesion_force);
			nPairi[i]->pushpair(collisionpairList[pair_n]);
			nPairj[i]->pushpair(collisionpairList[pair_n]);
		}
	}
	for (int i=0;i<opn2;i++)
	{
		int copn5 = (int)C->get_pairs()->size();
		bool exist = false;
		for (int j=0;j<copn5;j++)
		{
			cell *Cpo1 = C->get_pairs()->at(j)->get_I1();
			cell *Cpo2 = C->get_pairs()->at(j)->get_I2();
			if ((oPairi[i]==Cpo1 && oPairj[i]==Cpo2) ||
				(oPairi[i]==Cpo2 && oPairj[i]==Cpo1)) {exist = true;break;}
		}
		if (!exist)
		{
			int pair_n = (int)collisionpairList.size();
			interpair *Pl = new interpair(pair_n,oPairi[i],oPairj[i]);
			collisionpairList.push_back(Pl);
			//int iptype = interpair_type(nPairi[i],nPairj[i]);
			//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
			//collisionpairList[pair_n]->set_friction(friction_pair_class[iptype]);
			double adhesion_force = adhesion_pair_array[nPairi[i]->get_cell_type()][nPairj[i]->get_cell_type()];
			collisionpairList[pair_n]->set_adhesion(adhesion_force);
			collisionpairList[pair_n]->set_friction(adhesion_force);
			oPairi[i]->pushpair(collisionpairList[pair_n]);
			oPairj[i]->pushpair(collisionpairList[pair_n]);
		}
	}
	nPairi.clear();
	nPairj.clear();
	oPairi.clear();
	oPairj.clear();
	Del->~CD3DW();
	vector<interpair*> C_pair;
	vector<interpair*> Cn_pair;
	int cpn5 = (int)C->get_pairs()->size();
	int cpn6 = (int)cellList[cell_n]->get_pairs()->size();
	////// Pair refresh for burried cell ///////
	C->setup_Burry();
	cellList[cell_n]->setup_Burry();
	for (int i=cpn5;i>0;i--)
	{
		if (C->get_pairs()->at(i-1)->get_redundant()) 
		{
			C_pair.push_back(C->get_pairs()->at(i-1));
		}
		if (C->Burry() && (int)C->get_pairs()->at(i-1)->MP()->size()==0)
		{
			C->get_pairs()->at(i-1)->set_redundant(1);
			C_pair.push_back(C->get_pairs()->at(i-1));
		}
	}
	for (int i=0;i<(int)C_pair.size();i++)
	{
		cell *C11 = C_pair[i]->get_I1();
		cell *C12 = C_pair[i]->get_I2();
		C11->removepair(C_pair[i]);
		C12->removepair(C_pair[i]);
	}
	for (int i=cpn6;i>0;i--)
	{
		if (cellList[cell_n]->get_pairs()->at(i-1)->get_redundant()) 
		{
			Cn_pair.push_back(cellList[cell_n]->get_pairs()->at(i-1));
		}
		if (cellList[cell_n]->Burry() && (int)cellList[cell_n]->get_pairs()->at(i-1)->MP()->size()==0)
		{
			cellList[cell_n]->get_pairs()->at(i-1)->set_redundant(1);
			Cn_pair.push_back(cellList[cell_n]->get_pairs()->at(i-1));
		}
	}
	for (int i=0;i<(int)Cn_pair.size();i++)
	{
		cell *C11 = Cn_pair[i]->get_I1();
		cell *C12 = Cn_pair[i]->get_I2();
		C11->removepair(Cn_pair[i]);
		C12->removepair(Cn_pair[i]);
	}
	neighbor_c.clear();
	C_pair.clear();
	Cn_pair.clear();

	/***************************************
	       update migrating vertices
	***************************************/
	if (C->Migrate())
	{
		cellList[cell_n]->set_migrate_relax(1);
		point3D *m_s = C->get_migrate_p(0);
		point3D *m_t = C->get_migrate_p(1);
		bool one_side = false;
		if (m_s!=NULL &&
			m_t!=NULL &&
			m_s->getab()->Cell()==C &&
			m_t->getab()->Cell()==C)
		{one_side = true;}
		if (!one_side)
		{
			cellList[cell_n]->set_migrate(0);
			C->set_migrate(0);
			cellList[cell_n]->set_migrate_p(NULL,NULL);
			C->set_migrate_p(NULL,NULL);
		}
		cout<<", migrating status: "<<C->Migrate()<<", "<<cellList[cell_n]->Migrate()<<endl;
	}
	else
	{
		cout<<endl;
	}
}

void dbReader::cell_death_correction()
{
	map< int, vector<int> > dead_cell_index; // to store dying cells
	map< int, vector<int> > dead_cell_neighbor_index; // to store neighbor cells of dying cells
	int cn = (int)cellList.size();
	int dc_index = -1;
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (cellList[i]->Set_Dead())
		{
			cellList[i]->setup_Dead(1);
			dc_index++;
			if (dead_cell_index.size()==0)
			{
				vector<int> dead_cell_0;
				vector<int> dead_cell_neighbor_0;
				dead_cell_0.push_back(i);
				dead_cell_index[0] = dead_cell_0;
				int cpn = (int)cellList[i]->get_pairs()->size();
				for (int j=0;j<cpn;j++)
				{
					interpair *Ipn = cellList[i]->get_pairs()->at(j);
					if (Ipn->get_redundant()) continue;
					if      (Ipn->get_I1()==cellList[i] && !Ipn->get_I2()->Set_Dead())
					{dead_cell_neighbor_0.push_back(Ipn->get_I2()->id());}
					else if (Ipn->get_I2()==cellList[i] && !Ipn->get_I1()->Set_Dead())
					{dead_cell_neighbor_0.push_back(Ipn->get_I1()->id());}
				}
				dead_cell_neighbor_index[0] = dead_cell_neighbor_0;
			}
			else
			{
				bool flagc = false;
				vector<int> map_index;
				map<int, vector<int> >::iterator its;
				for (its=dead_cell_index.begin();its!=dead_cell_index.end();its++)
				{
					int vsn = (int)its->second.size();
					bool flagcid = false;
					for (int j=0;j<vsn;j++)
					{
						int cid = its->second[j];
						int cpn = cellList[cid]->get_pairs()->size();
						for (int k=0;k<cpn;k++)
						{
							if (cellList[cid]->get_pairs()->at(k)->get_redundant()) continue;
							if ((int)cellList[cid]->get_pairs()->at(k)->MP()->size()==0) continue;
							if (cellList[cid]->get_pairs()->at(k)->get_I1()==cellList[i] ||
								cellList[cid]->get_pairs()->at(k)->get_I2()==cellList[i])
							{
								flagcid = true;break;
							}
						}
						if (flagcid) break;
					}
					if (flagcid)
					{
						map_index.push_back(its->first);
					}
				}
				int map_index_n = (int)map_index.size();
				if (map_index_n==0)
				{
					vector<int> dead_cell_0;
					vector<int> dead_cell_neighbor_0;
					dead_cell_0.push_back(i);
					dead_cell_index[dc_index] = dead_cell_0;
					int cpn = (int)cellList[i]->get_pairs()->size();
					for (int j=0;j<cpn;j++)
					{
						interpair *Ipn = cellList[i]->get_pairs()->at(j);
						if (Ipn->get_redundant()) continue;
						if      (Ipn->get_I1()==cellList[i] && !Ipn->get_I2()->Set_Dead()) 
						{dead_cell_neighbor_0.push_back(Ipn->get_I2()->id());}
						else if (Ipn->get_I2()==cellList[i] && !Ipn->get_I1()->Set_Dead()) 
						{dead_cell_neighbor_0.push_back(Ipn->get_I1()->id());}
					}
					dead_cell_neighbor_index[dc_index] = dead_cell_neighbor_0;
				}
				else if (map_index_n==1)
				{
					map<int, vector<int> >::iterator itt;
					itt = dead_cell_index.find(map_index[0]);
					itt->second.push_back(i);
					int cpn = cellList[i]->get_pairs()->size();
					for (int j=0;j<cpn;j++)
					{
						interpair *Ipn = cellList[i]->get_pairs()->at(j);
						if (Ipn->get_redundant()) continue;
						int cnid = 0;
						if (Ipn->get_I1()==cellList[i]) {cnid = Ipn->get_I2()->id();}
						else                            {cnid = Ipn->get_I1()->id();}
						if (cellList[cnid]->Set_Dead()) continue;
						map<int, vector<int> >::iterator ittn;
						ittn = dead_cell_neighbor_index.find(map_index[0]);
						bool flagk = false;
						for (int k=0;k<(int)ittn->second.size();k++)
						{
							int cid = ittn->second[k];
							int cpn = (int)cellList[i]->get_pairs()->size();
							if (cellList[cnid]==cellList[cid])
							{
								flagk = true;break;
							}
						}
						if (!flagk)
						{
							ittn->second.push_back(cnid);
						}
					}
				}
				else if (map_index_n>1)
				{
					for (int j=1;j<map_index_n;j++)
					{
						map<int, vector<int> >::iterator itto;
						itto = dead_cell_index.find(map_index[0]);
						map<int, vector<int> >::iterator ittt;
						ittt = dead_cell_index.find(map_index[j]);
						for (int k=0;k<(int)ittt->second.size();k++)
						{
							itto->second.push_back(ittt->second[k]);
						}
						dead_cell_index.erase(map_index[j]);
					}
					map<int, vector<int> >::iterator itt;
					itt = dead_cell_index.find(map_index[0]);
					itt->second.push_back(i);
					for (int j=1;j<map_index_n;j++)
					{
						map<int, vector<int> >::iterator itto;
						itto = dead_cell_neighbor_index.find(map_index[0]);
						map<int, vector<int> >::iterator ittt;
						ittt = dead_cell_neighbor_index.find(map_index[j]);
						for (int k=0;k<(int)ittt->second.size();k++)
						{
							bool flagl = false;
							for (int l=0;l<(int)itto->second.size();l++)
							{
								if (ittt->second[k]==itto->second[l])
								{flagl = true;break;}
							}
							if (!flagl) {itto->second.push_back(ittt->second[k]);}
						}
						dead_cell_neighbor_index.erase(map_index[j]);
					}
					map<int, vector<int> >::iterator ittt;
					ittt = dead_cell_neighbor_index.find(map_index[0]);
					int cpn = (int)cellList[i]->get_pairs()->size();
					for (int j=0;j<cpn;j++)
					{
						interpair *Ipn = cellList[i]->get_pairs()->at(j);
						if (Ipn->get_redundant()) continue;
						int cid = 0;
						if      (Ipn->get_I1()==cellList[i])
						{cid = Ipn->get_I2()->id();}
						else if (Ipn->get_I2()==cellList[i])
						{cid = Ipn->get_I1()->id();}
						if (cellList[cid]->Set_Dead()) continue;
						bool flagk = false;
						for (int k=0;k<(int)ittt->second.size();k++)
						{
							if (ittt->second[k]==cid)
							{flagk = true;break;}
						}
						if (!flagk)
						{ittt->second.push_back(cid);}
					}
				}
			}
		}
	}

	if ((int)dead_cell_index.size()>0)
	{
		map<int, vector<int> >::iterator its;
		for (its = dead_cell_index.begin();its!=dead_cell_index.end();its++)
		{
			int dc_n = (int)its->second.size();
			for (int i=0;i<dc_n;i++)
			{
				int c_n = its->second[i];
				//cout<<"cell "<<c_n<<" is dead"<<endl;
				cellList[c_n]->set_set_dead(0);
				////////// update side relationship ///////////
				for (int j=0;j<(int)cellList[c_n]->get_sides()->size();j++)
				{
					if (cellList[c_n]->get_sides()->at(j)->p1()->getmp()->size()>0) // attached with other vertices
					{
						cellList[c_n]->get_sides()->at(j)->p1()->get_node_id()->removevertex(cellList[c_n]->get_sides()->at(j)->p1());
						cellList[c_n]->get_sides()->at(j)->p1()->set_node_id(NULL);
						cellList[c_n]->get_sides()->at(j)->p1()->set_in_cell(0);
						for (int k=0;k<(int)cellList[c_n]->get_sides()->at(j)->p1()->getmp()->size();k++)
						{
							cellList[c_n]->get_sides()->at(j)->p1()->getmp()->at(k)->removemp(cellList[c_n]->get_sides()->at(j)->p1());
						}
						if (cellList[c_n]->get_sides()->at(j)->get_attach())
						{
							cellList[c_n]->get_sides()->at(j)->get_Neighbor()->set_attach(0);
							cellList[c_n]->get_sides()->at(j)->get_Neighbor()->remove_Neighbor();
						}
					}
					else // single vertex
					{
						cellList[c_n]->get_sides()->at(j)->p1()->set_in_cell(0);
						cellList[c_n]->get_sides()->at(j)->p1()->get_node_id()->set_in_global(0);
					}
				}
				/////////////// update neighbor relationship //////////////
				for (int j=0;j<(int)cellList[c_n]->get_pairs()->size();j++)
				{
					cellList[c_n]->get_pairs()->at(j)->set_redundant(1);
					if      (cellList[c_n]->get_pairs()->at(j)->get_I1()==cellList[c_n] &&
							!cellList[c_n]->get_pairs()->at(j)->get_I2()->Dead())
					{
						cellList[c_n]->get_pairs()->at(j)->get_I2()->removepair(cellList[c_n]->get_pairs()->at(j));
					}
					else if (cellList[c_n]->get_pairs()->at(j)->get_I2()==cellList[c_n] &&
							!cellList[c_n]->get_pairs()->at(j)->get_I1()->Dead())
					{
						cellList[c_n]->get_pairs()->at(j)->get_I1()->removepair(cellList[c_n]->get_pairs()->at(j));
					}
				}
			}
		}
		/**********************************************
		  for migrating cell, renew all possible 
		  collision pair for migration which might 
		  break the regular interrelationship
		**********************************************/
		map<int, vector<int> >::iterator itt;
		for (itt = dead_cell_neighbor_index.begin();itt!=dead_cell_neighbor_index.end();itt++)
		{
			vector<cell*> migrate_cell;
			vector<cell*> static_cell;
			int dcn_n = (int)itt->second.size();
			for (int i=0;i<dcn_n;i++)
			{
				int cid = itt->second[i];
				if (cellList[cid]->Migrate()) {migrate_cell.push_back(cellList[cid]);}
				else                          {static_cell.push_back(cellList[cid]);}
			}
			int mcn = (int)migrate_cell.size();
			int scn = (int)static_cell.size();
			for (int i=0;i<mcn;i++)
			{
				for (int j=0;j<scn;j++)
				{
					cell *Cm = migrate_cell[i];
					cell *Cs = static_cell[j];
					int cpn = (int)Cm->get_pairs()->size();
					bool flagc = false;
					for (int k=0;k<cpn;k++)
					{
						if ((Cm->get_pairs()->at(k)->get_I1()==Cm && Cm->get_pairs()->at(k)->get_I2()==Cs) || 
							(Cm->get_pairs()->at(k)->get_I2()==Cm && Cm->get_pairs()->at(k)->get_I1()==Cs))
						{flagc = true;break;}
					}
					if (!flagc)
					{
						int pair_n = (int)collisionpairList.size();
						interpair *Pl = new interpair(pair_n,Cm,Cs);
						collisionpairList.push_back(Pl);
						//int iptype = interpair_type(Cm,Cs);
						//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
						double adhesion_force = adhesion_pair_array[Cm->get_cell_type()][Cs->get_cell_type()];
						collisionpairList[pair_n]->set_adhesion(adhesion_force);
						collisionpairList[pair_n]->set_friction(adhesion_force);
						Cm->pushpair(collisionpairList[pair_n]);
						Cs->pushpair(collisionpairList[pair_n]);
					}
				}
			}
			if (mcn>1)
			{
				for (int i=0;i<mcn-1;i++)
				{
					for (int j=i+1;j<mcn;j++)
					{
						cell *Cmi = migrate_cell[i];
						cell *Cmj = migrate_cell[j];
						int cpn = (int)Cmi->get_pairs()->size();
						bool flagc = false;
						for (int k=0;k<cpn;k++)
						{
							if ((Cmi->get_pairs()->at(k)->get_I1()==Cmi && Cmi->get_pairs()->at(k)->get_I2()==Cmj) ||
								(Cmi->get_pairs()->at(k)->get_I2()==Cmi && Cmi->get_pairs()->at(k)->get_I1()==Cmj))
							{flagc = true;break;}
						}
						if (!flagc)
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pl = new interpair(pair_n,Cmi,Cmj);
							collisionpairList.push_back(Pl);
							//int iptype = interpair_type(Cmi,Cmj);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[Cmi->get_cell_type()][Cmj->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							Cmi->pushpair(collisionpairList[pair_n]);
							Cmj->pushpair(collisionpairList[pair_n]);
						}
					}
				}
			}
			migrate_cell.clear();
			static_cell.clear();
		}
		//////////////// renew neighbor relationship ///////////////
		for (itt = dead_cell_neighbor_index.begin();itt!=dead_cell_neighbor_index.end();itt++) 
		{
			double range[4]; // 0:xmin 1:ymin 2:xmax 3:ymax
			range[0] = 1000000000000;
			range[1] = 1000000000000;
			range[2] =-1000000000000;
			range[3] =-1000000000000;
			for (int i=0;i<(int)itt->second.size();i++)
			{
				int c_n = itt->second[i];
				for (int j=0;j<(int)cellList[c_n]->get_inner_p()->size();j++)
				{
					double vx = cellList[c_n]->get_inner_p()->at(j)->x();
					double vy = cellList[c_n]->get_inner_p()->at(j)->y();
					if (range[0]>vx) {range[0] = vx;}
					if (range[1]>vy) {range[1] = vy;}
					if (range[2]<vx) {range[2] = vx;}
					if (range[3]<vy) {range[3] = vy;}
				}
			}
			range[0] -= 1000;
			range[1] -= 1000;
			range[2] += 1000;
			range[3] += 1000;

			CD3DW *Del = new CD3DW(range);
			int flag = 0;
			for (int i=0;i<(int)itt->second.size();i++)
			{
				int c_n = itt->second[i];
				if ((int)cellList[c_n]->get_inner_p()->size()>0)
				{
					for (int j=0;j<(int)cellList[c_n]->get_inner_p()->size();j++)
					{
						if (cellList[c_n]->get_inner_p()->at(j)->get_burry()) continue;
						Del->add_resample(cellList[c_n]->get_inner_p()->at(j)->x(),
							cellList[c_n]->get_inner_p()->at(j)->y(),
							BD*BD,
							flag,
							cellList[c_n]->get_inner_p()->at(j));
						flag++;
					}
				}
				else
				{
					for (int j=0;j<(int)cellList[c_n]->get_inner_p_slip()->size();j++)
					{
						Del->add_resample(cellList[c_n]->get_inner_p_slip()->at(j)->x(),
							cellList[c_n]->get_inner_p_slip()->at(j)->y(),
							BD*BD,
							flag,
							cellList[c_n]->get_inner_p_slip()->at(j));
						flag++;
					}
				}
			}
			for (int i=0;i<(int)Del->get_Cs()->size();i++)
			{
				if (Del->get_Cs()->at(i)->BN() &&
					Del->get_Cs()->at(i)->get_VV(0)->n()>=0 &&
					Del->get_Cs()->at(i)->get_VV(1)->n()>=0 &&
					Del->get_Cs()->at(i)->get_VV(2)->n()>=0)
				{
					if		(Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()) 
					{
						cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
						cell* C2 = Del->get_Cs()->at(i)->get_VV(1)->P()->Cell();
						cell* C3 = Del->get_Cs()->at(i)->get_VV(2)->P()->Cell();
						bool exist12 = false;
						bool exist13 = false;
						bool exist23 = false;
						for (int j=0;j<(int)C1->get_pairs()->size();j++) 
						{
							if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
								(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;}
							if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C3) ||
								(C1->get_pairs()->at(j)->get_I1()==C3 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist13 = true;}
						}
						for (int j=0;j<(int)C2->get_pairs()->size();j++) 
						{
							if ((C2->get_pairs()->at(j)->get_I1()==C2 && C2->get_pairs()->at(j)->get_I2()==C3) ||
								(C2->get_pairs()->at(j)->get_I1()==C3 && C2->get_pairs()->at(j)->get_I2()==C2)) {exist23 = true;break;}
						}
						if (!exist12) 
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pn = new interpair(pair_n,C1,C2);
							collisionpairList.push_back(Pn);
							//int iptype = interpair_type(C1,C2);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							C1->pushpair(collisionpairList[pair_n]);
							C2->pushpair(collisionpairList[pair_n]);
						}
						if (!exist13) 
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pn = new interpair(pair_n,C1,C3);
							collisionpairList.push_back(Pn);
							//int iptype = interpair_type(C1,C3);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C3->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							C1->pushpair(collisionpairList[pair_n]);
							C3->pushpair(collisionpairList[pair_n]);
						}
						if (!exist23) 
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pn = new interpair(pair_n,C2,C3);
							collisionpairList.push_back(Pn);
							//int iptype = interpair_type(C2,C3);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[C2->get_cell_type()][C3->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							C2->pushpair(collisionpairList[pair_n]);
							C3->pushpair(collisionpairList[pair_n]);
						}
					}
					else if (Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()==Del->get_Cs()->at(i)->get_VV(1)->P()->Cell())
					{
						cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
						cell* C2 = Del->get_Cs()->at(i)->get_VV(2)->P()->Cell();
						bool exist12 = false;
						for (int j=0;j<(int)C1->get_pairs()->size();j++) 
						{
							if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
								(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;break;}
						}
						if (!exist12)
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pn = new interpair(pair_n,C1,C2);
							collisionpairList.push_back(Pn);
							//int iptype = interpair_type(C1,C2);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							C1->pushpair(collisionpairList[pair_n]);
							C2->pushpair(collisionpairList[pair_n]);
						}
					}
					else if (Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(1)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(2)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(1)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()==Del->get_Cs()->at(i)->get_VV(2)->P()->Cell())
					{
						cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
						cell* C2 = Del->get_Cs()->at(i)->get_VV(1)->P()->Cell();
						bool exist12 = false;
						for (int j=0;j<(int)C1->get_pairs()->size();j++) 
						{
							if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
								(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;break;}
						}
						if (!exist12)
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pn = new interpair(pair_n,C1,C2);
							collisionpairList.push_back(Pn);
							//int iptype = interpair_type(C1,C2);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							C1->pushpair(collisionpairList[pair_n]);
							C2->pushpair(collisionpairList[pair_n]);
						}
					}
					else if (Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(0)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(2)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(0)->P()->Cell() &&
							 Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()==Del->get_Cs()->at(i)->get_VV(2)->P()->Cell())
					{
						cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
						cell* C2 = Del->get_Cs()->at(i)->get_VV(2)->P()->Cell();
						bool exist12 = false;
						for (int j=0;j<(int)C1->get_pairs()->size();j++) 
						{
							if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
								(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;break;}
						}
						if (!exist12)
						{
							int pair_n = (int)collisionpairList.size();
							interpair *Pn = new interpair(pair_n,C1,C2);
							collisionpairList.push_back(Pn);
							//int iptype = interpair_type(C1,C2);
							//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
							double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
							collisionpairList[pair_n]->set_adhesion(adhesion_force);
							collisionpairList[pair_n]->set_friction(adhesion_force);
							C1->pushpair(collisionpairList[pair_n]);
							C2->pushpair(collisionpairList[pair_n]);
						}
					}
				}
			}
			Del->~CD3DW();
		}
	}
	dead_cell_index.clear();
	dead_cell_neighbor_index.clear();
}

void dbReader::cell_death_response(cell *C)
{
	C->set_set_dead(0);
	C->setup_Dead(1);
	int sn = C->get_sides()->size();
	int pn = C->get_pairs()->size();
	vector<cell*> neighbor_cells;
	////////// update side relationship ///////////
	for (int i=0;i<sn;i++)
	{
		point3D *p1 = C->get_sides()->at(i)->p1();
		if (p1->getmp()->size()>0)
		{
			p1->get_node_id()->removevertex(p1);
			p1->set_node_id(NULL);
			p1->set_in_cell(0);
			for (int j=0;j<(int)p1->getmp()->size();j++)
			{
				p1->getmp()->at(j)->removemp(p1);
			}
			if (p1->getab()->get_attach())
			{
				p1->getab()->get_Neighbor()->set_attach(0);
				p1->getab()->get_Neighbor()->remove_Neighbor();
			}
		}
		else
		{
			p1->set_in_cell(0);
			p1->get_node_id()->set_in_global(0);
		}
	}
	////////// update neighbor relationship ///////////
	for (int i=0;i<pn;i++)
	{
		interpair *Ip = C->get_pairs()->at(i);
		Ip->set_redundant(1);
		if      (Ip->get_I1()==C &&
				!Ip->get_I2()->Dead())
		{
			Ip->get_I2()->removepair(Ip);
			neighbor_cells.push_back(Ip->get_I2());
		}
		else if (Ip->get_I2()==C &&
				!Ip->get_I1()->Dead())
		{
			Ip->get_I1()->removepair(Ip);
			neighbor_cells.push_back(Ip->get_I1());
		}
	}
	/*********************************************
	   for migrating cell, renew all possible 
	   collision pair for migration which might 
	   break the regular interrelationship
	*********************************************/
	vector<cell*> migrate_cell;
	vector<cell*> static_cell;
	int dcn_n = (int)neighbor_cells.size();
	for (int i=0;i<dcn_n;i++)
	{
		if (neighbor_cells[i]->Migrate()) {migrate_cell.push_back(neighbor_cells[i]);}
		else                              {static_cell.push_back(neighbor_cells[i]);}
	}
	int mcn = (int)migrate_cell.size();
	int scn = (int)static_cell.size();
	for (int i=0;i<mcn;i++)
	{
		for (int j=0;j<scn;j++)
		{
			cell *Cm = migrate_cell[i];
			cell *Cs = static_cell[j];
			int cpn = (int)Cm->get_pairs()->size();
			bool flagc = false;
			for (int k=0;k<cpn;k++)
			{
				if ((Cm->get_pairs()->at(k)->get_I1()==Cm && Cm->get_pairs()->at(k)->get_I2()==Cs) || 
					(Cm->get_pairs()->at(k)->get_I2()==Cm && Cm->get_pairs()->at(k)->get_I1()==Cs))
				{flagc = true;break;}
			}
			if (!flagc)
			{
				int pair_n = (int)collisionpairList.size();
				interpair *Pl = new interpair(pair_n,Cm,Cs);
				collisionpairList.push_back(Pl);
				//int iptype = interpair_type(Cm,Cs);
				//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
				double adhesion_force = adhesion_pair_array[Cm->get_cell_type()][Cs->get_cell_type()];
				collisionpairList[pair_n]->set_adhesion(adhesion_force);
				collisionpairList[pair_n]->set_friction(adhesion_force);
				Cm->pushpair(collisionpairList[pair_n]);
				Cs->pushpair(collisionpairList[pair_n]);
			}
		}
	}
	if (mcn>1)
	{
		for (int i=0;i<mcn-1;i++)
		{
			for (int j=i+1;j<mcn;j++)
			{
				cell *Cmi = migrate_cell[i];
				cell *Cmj = migrate_cell[j];
				int cpn = (int)Cmi->get_pairs()->size();
				bool flagc = false;
				for (int k=0;k<cpn;k++)
				{
					if ((Cmi->get_pairs()->at(k)->get_I1()==Cmi && Cmi->get_pairs()->at(k)->get_I2()==Cmj) ||
						(Cmi->get_pairs()->at(k)->get_I2()==Cmi && Cmi->get_pairs()->at(k)->get_I1()==Cmj))
					{flagc = true;break;}
				}
				if (!flagc)
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pl = new interpair(pair_n,Cmi,Cmj);
					collisionpairList.push_back(Pl);
					//int iptype = interpair_type(Cmi,Cmj);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[Cmi->get_cell_type()][Cmj->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					Cmi->pushpair(collisionpairList[pair_n]);
					Cmj->pushpair(collisionpairList[pair_n]);
				}
			}
		}
	}
	migrate_cell.clear();
	static_cell.clear();
	////////////// renew neighbor relationship ////////////////
	double range[4]; // 0:xmin 1:ymin 2:xmax 3:ymax
	range[0] = 1000000000000;
	range[1] = 1000000000000;
	range[2] =-1000000000000;
	range[3] =-1000000000000;
	for (int i=0;i<dcn_n;i++)
	{
		for (int j=0;j<(int)neighbor_cells[i]->get_inner_p()->size();j++)
		{
			double vx = neighbor_cells[i]->get_inner_p()->at(j)->x();
			double vy = neighbor_cells[i]->get_inner_p()->at(j)->y();
			if (range[0]>vx) {range[0] = vx;}
			if (range[1]>vy) {range[1] = vy;}
			if (range[2]<vx) {range[2] = vx;}
			if (range[3]<vy) {range[3] = vy;}
		}
	}
	range[0] -= 1000;
	range[1] -= 1000;
	range[2] += 1000;
	range[3] += 1000;

	CD3DW *Del = new CD3DW(range);
	int flag = 0;
	for (int i=0;i<(int)neighbor_cells.size();i++)
	{
		if ((int)neighbor_cells[i]->get_inner_p()->size()>0)
		{
			for (int j=0;j<(int)neighbor_cells[i]->get_inner_p()->size();j++)
			{
				if (neighbor_cells[i]->get_inner_p()->at(j)->get_burry()) continue;
				Del->add_resample(neighbor_cells[i]->get_inner_p()->at(j)->x(),
					neighbor_cells[i]->get_inner_p()->at(j)->y(),
					BD*BD,
					flag,
					neighbor_cells[i]->get_inner_p()->at(j));
				flag++;
			}
		}
		else
		{
			for (int j=0;j<(int)neighbor_cells[i]->get_inner_p_slip()->size();j++)
			{
				Del->add_resample(neighbor_cells[i]->get_inner_p_slip()->at(j)->x(),
					neighbor_cells[i]->get_inner_p_slip()->at(j)->y(),
					BD*BD,
					flag,
					neighbor_cells[i]->get_inner_p_slip()->at(j));
				flag++;
			}
		}
	}
	for (int i=0;i<(int)Del->get_Cs()->size();i++)
	{
		if (Del->get_Cs()->at(i)->BN() &&
			Del->get_Cs()->at(i)->get_VV(0)->n()>=0 &&
			Del->get_Cs()->at(i)->get_VV(1)->n()>=0 &&
			Del->get_Cs()->at(i)->get_VV(2)->n()>=0)
		{
			if		(Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()) 
			{
				cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
				cell* C2 = Del->get_Cs()->at(i)->get_VV(1)->P()->Cell();
				cell* C3 = Del->get_Cs()->at(i)->get_VV(2)->P()->Cell();
				bool exist12 = false;
				bool exist13 = false;
				bool exist23 = false;
				for (int j=0;j<(int)C1->get_pairs()->size();j++) 
				{
					if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
						(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;}
					if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C3) ||
						(C1->get_pairs()->at(j)->get_I1()==C3 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist13 = true;}
				}
				for (int j=0;j<(int)C2->get_pairs()->size();j++) 
				{
					if ((C2->get_pairs()->at(j)->get_I1()==C2 && C2->get_pairs()->at(j)->get_I2()==C3) ||
						(C2->get_pairs()->at(j)->get_I1()==C3 && C2->get_pairs()->at(j)->get_I2()==C2)) {exist23 = true;break;}
				}
				if (!exist12) 
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pn = new interpair(pair_n,C1,C2);
					collisionpairList.push_back(Pn);
					//int iptype = interpair_type(C1,C2);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					C1->pushpair(collisionpairList[pair_n]);
					C2->pushpair(collisionpairList[pair_n]);
				}
				if (!exist13) 
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pn = new interpair(pair_n,C1,C3);
					collisionpairList.push_back(Pn);
					//int iptype = interpair_type(C1,C3);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C3->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					C1->pushpair(collisionpairList[pair_n]);
					C3->pushpair(collisionpairList[pair_n]);
				}
				if (!exist23) 
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pn = new interpair(pair_n,C2,C3);
					collisionpairList.push_back(Pn);
					//int iptype = interpair_type(C2,C3);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[C2->get_cell_type()][C3->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					C2->pushpair(collisionpairList[pair_n]);
					C3->pushpair(collisionpairList[pair_n]);
				}
			}
			else if (Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(2)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()==Del->get_Cs()->at(i)->get_VV(1)->P()->Cell())
			{
				cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
				cell* C2 = Del->get_Cs()->at(i)->get_VV(2)->P()->Cell();
				bool exist12 = false;
				for (int j=0;j<(int)C1->get_pairs()->size();j++) 
				{
					if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
						(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pn = new interpair(pair_n,C1,C2);
					collisionpairList.push_back(Pn);
					//int iptype = interpair_type(C1,C2);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					C1->pushpair(collisionpairList[pair_n]);
					C2->pushpair(collisionpairList[pair_n]);
				}
			}
			else if (Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(1)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(2)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(1)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(0)->P()->Cell()==Del->get_Cs()->at(i)->get_VV(2)->P()->Cell())
			{
				cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
				cell* C2 = Del->get_Cs()->at(i)->get_VV(1)->P()->Cell();
				bool exist12 = false;
				for (int j=0;j<(int)C1->get_pairs()->size();j++) 
				{
					if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
						(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pn = new interpair(pair_n,C1,C2);
					collisionpairList.push_back(Pn);
					//int iptype = interpair_type(C1,C2);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					C1->pushpair(collisionpairList[pair_n]);
					C2->pushpair(collisionpairList[pair_n]);
				}
			}
			else if (Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(0)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(2)->P()->Cell()!=Del->get_Cs()->at(i)->get_VV(0)->P()->Cell() &&
					 Del->get_Cs()->at(i)->get_VV(1)->P()->Cell()==Del->get_Cs()->at(i)->get_VV(2)->P()->Cell())
			{
				cell* C1 = Del->get_Cs()->at(i)->get_VV(0)->P()->Cell();
				cell* C2 = Del->get_Cs()->at(i)->get_VV(2)->P()->Cell();
				bool exist12 = false;
				for (int j=0;j<(int)C1->get_pairs()->size();j++) 
				{
					if ((C1->get_pairs()->at(j)->get_I1()==C1 && C1->get_pairs()->at(j)->get_I2()==C2) ||
						(C1->get_pairs()->at(j)->get_I1()==C2 && C1->get_pairs()->at(j)->get_I2()==C1)) {exist12 = true;break;}
				}
				if (!exist12)
				{
					int pair_n = (int)collisionpairList.size();
					interpair *Pn = new interpair(pair_n,C1,C2);
					collisionpairList.push_back(Pn);
					//int iptype = interpair_type(C1,C2);
					//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
					double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
					collisionpairList[pair_n]->set_adhesion(adhesion_force);
					collisionpairList[pair_n]->set_friction(adhesion_force);
					C1->pushpair(collisionpairList[pair_n]);
					C2->pushpair(collisionpairList[pair_n]);
				}
			}
		}
	}
	Del->~CD3DW();
	neighbor_cells.clear();
}

void dbReader::cell_migrate_response_0(cell *C, double angle, vector<int>* d, double F)
{
	/************************************
	                 A
	                 |
	           o---o---o---o
	          /      E      \
	         o               o
	      - -| - - - o- - - -|- -
	    o---oo               oo---o
	   /     \\             //     \
	  o       oo---o---o---oo       o
	  |       ||           ||       |

	************************************/
	C->center_refresh();
	C->setup_Mark_number(1); // 1: migrate
	C->set_interior_refresh(1);
	double xc = C->get_center(0);
	double yc = C->get_center(1);
	edge *E = C->get_sides()->at(0);
	bool force_control = false;
	double dx = 0;
	double dy = 1;
	if (angle==90 || angle==270)
	{
		if (angle==270) {dy = -1;}
		bool flagE = false;
		int sn = (int)C->get_sides()->size();
		for (int i=0;i<sn;i++)
		{
			double x1 = C->get_sides()->at(i)->p1()->x();
			double x2 = C->get_sides()->at(i)->p2()->x();
			if (x1<=xc && x2>xc)
			{
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = dy*ny;
				if (plot_dn>0)
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	else
	{
		dx = cos(angle*PI/180);
		dy = sin(angle*PI/180);
		double slope = tan(angle*PI/180);
		double intercept = yc - slope*xc;
		int sn = (int)C->get_sides()->size();
		bool flagE = false;
		for (int i=0;i<sn;i++)
		{
			point3D *p1 = C->get_sides()->at(i)->p1();
			point3D *p2 = C->get_sides()->at(i)->p2();
			double x1 = p1->x();
			double y1 = p1->y();
			double x2 = p2->x();
			double y2 = p2->y();
			double t1 = y1 - slope*x1 - intercept;
			double t2 = y2 - slope*x2 - intercept;
			if ((t1>=0 && t2<0) || (t1<0 && t2>=0))
			{
				double nx = C->get_sides()->at(i)->getnormal(0);
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = nx*dx + ny*dy;
				if (plot_dn>0) 
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	point3D *sp = E->p1();
	point3D *tp = E->p2();
	point3D *tp_ex = NULL;
	point3D *sp_ex = NULL;
	if (!E->get_attach()) // single
	{
		force_control = true;
		bool flagst = false;
		while (!flagst)
		{
			if      (!sp->getab()->get_attach() && !sp->getba()->get_attach())
			{
				double nabx = sp->getab()->getnormal(0);
				double naby = sp->getab()->getnormal(1);
				double nbax = sp->getba()->getnormal(0);
				double nbay = sp->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double sp_x = nab_nx + nba_nx;
				double sp_y = nab_ny + nba_ny;
				double sp_d = sqrt(sp_x*sp_x + sp_y*sp_y);
				double sp_nx = sp_x/sp_d;
				double sp_ny = sp_y/sp_d;
				double plot_spd = sp_nx*dx + sp_ny*dy;
				if (plot_spd<0.17) {flagst = true;break;} // the angle is less than 80 degrees
				else               {sp = sp->getrp();}
			}
			else if (!sp->getab()->get_attach() && sp->getba()->get_attach())
			{flagst = true;break;}
		}
		bool flagtt = false;
		while (!flagtt)
		{
			if      (!tp->getab()->get_attach() && !tp->getba()->get_attach())
			{
				double nabx = tp->getab()->getnormal(0);
				double naby = tp->getab()->getnormal(1);
				double nbax = tp->getba()->getnormal(0);
				double nbay = tp->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double tp_x = nab_nx + nba_nx;
				double tp_y = nab_ny + nba_ny;
				double tp_d = sqrt(tp_x*tp_x + tp_y*tp_y);
				double tp_nx = tp_x/tp_d;
				double tp_ny = tp_y/tp_d;
				double plot_tpd = tp_nx*dx + tp_ny*dy;
				if (plot_tpd<0.17) {flagtt = true;break;} // the angle is less than 80 degrees
				else              {tp = tp->getfp();}
			}
			else if (tp->getab()->get_attach() && !tp->getba()->get_attach())
			{flagtt = true;break;}
		}
	}
	else // attached
	{
		interpair *Ip = E->get_Pair();
		cell *C1 = Ip->get_I1();
		cell *C2 = Ip->get_I2();
		cell *Cn = NULL; // the neighbor of C: C towards Cn
		bool Cn_death = false;
		bool Cn_move = false; // if the cell can or can't degrade the cell on the way
		if (Ip->get_I1()==C)
		{
			if (find(d->begin(),d->end(),C2->get_cell_type())!=d->end())
			{
				Cn_move = true;
				Cn = C2;
				/// test the soften cell thickness ///
				point3D *sc = Ip->MP()->at(0)->S2();
				point3D *tc = Ip->MP()->at(0)->T2();
				if (sc!=tc)
				{
					vector<point3D*> scp;
					vector<double> scp_plot;
					point3D *max_plot_index = NULL;
					double max_scp_plot = -1000000;
					point3D *pmc = sc;
					while (pmc!=tc)
					{
						double plot_c = pmc->x()*dx + pmc->y()*dy;
						if (max_scp_plot<plot_c)
						{
							max_scp_plot = plot_c;
							max_plot_index = pmc;
						}
						pmc = pmc->getfp();
					}
					if (max_plot_index!=NULL)
					{
						for (int i=0;i<(int)max_plot_index->get_inner_at()->size();i++)
						{
							triangle *Tt = max_plot_index->get_inner_at()->at(i);
							point3D *Ta = Tt->getA();
							point3D *Tb = Tt->getB();
							point3D *Tc = Tt->getC();
							double alphav = Tt->alpha();
							if (Ta->get_cell_boundary() && 
								Tb->get_cell_boundary() && 
								Tc->get_cell_boundary() &&
								alphav<BD*BD)
							{
								point3D *p1 = NULL;
								point3D *p2 = NULL;
								point3D *p3 = NULL;
								if      (Ta==max_plot_index) {p1 = Ta;p2 = Tb;p3 = Tc;}
								else if (Tb==max_plot_index) {p1 = Tb;p2 = Ta;p3 = Tc;}
								else if (Tc==max_plot_index) {p1 = Tc;p2 = Ta;p3 = Tb;}
								double n2abx = p2->getab()->getnormal(0);
								double n2aby = p2->getab()->getnormal(1);
								double n2bax = p2->getba()->getnormal(0);
								double n2bay = p2->getba()->getnormal(1);
								double n3abx = p3->getab()->getnormal(0);
								double n3aby = p3->getab()->getnormal(1);
								double n3bax = p3->getba()->getnormal(0);
								double n3bay = p3->getba()->getnormal(1);
								double n2abd = sqrt(n2abx*n2abx + n2aby*n2aby);
								double n2bad = sqrt(n2bax*n2bax + n2bay*n2bay);
								double n3abd = sqrt(n3abx*n3abx + n3aby*n3aby);
								double n3bad = sqrt(n3bax*n3bax + n3bay*n3bay);
								double n2x = n2abx/n2abd + n2bax/n2bad;
								double n2y = n2aby/n2abd + n2bay/n2bad;
								double n3x = n3abx/n3abd + n3bax/n3bad;
								double n3y = n3aby/n3abd + n3bay/n3bad;
								double n2d = sqrt(n2x*n2x + n2y*n2y);
								double n3d = sqrt(n3x*n3x + n3y*n3y);
								double n2nx = n2x/n2d;
								double n2ny = n2y/n2d;
								double n3nx = n3x/n3d;
								double n3ny = n3y/n3d;
								double plot_n2 = n2nx*dx + n2ny*dy;
								double plot_n3 = n3nx*dx + n3ny*dy;
								if (plot_n2>0 || plot_n3>0)
								{
									Cn_death = true;break;
								}
							}
						}
					}
				}
			}
			sp = Ip->MP()->at(0)->S1();
			tp = Ip->MP()->at(0)->T1();
		}
		else
		{
			if (find(d->begin(),d->end(),C1->get_cell_type())!=d->end())
			{
				Cn = C1;
				Cn_move = true;
				/// test the soften cell thickness ///
				point3D *sc = Ip->MP()->at(0)->S1();
				point3D *tc = Ip->MP()->at(0)->T1();
				if (sc!=tc)
				{
					vector<point3D*> scp;
					vector<double> scp_plot;
					point3D *max_plot_index = NULL;
					double max_scp_plot = -1000000;
					point3D *pmc = sc;
					while (pmc!=tc)
					{
						double plot_c = pmc->x()*dx + pmc->y()*dy;
						if (max_scp_plot<plot_c)
						{
							max_scp_plot = plot_c;
							max_plot_index = pmc;
						}
						pmc = pmc->getfp();
					}
					if (max_plot_index!=NULL)
					{
						for (int i=0;i<(int)max_plot_index->get_inner_at()->size();i++)
						{
							triangle *Tt = max_plot_index->get_inner_at()->at(i);
							point3D *Ta = Tt->getA();
							point3D *Tb = Tt->getB();
							point3D *Tc = Tt->getC();
							double alphav = Tt->alpha();
							if (Ta->get_cell_boundary() && 
								Tb->get_cell_boundary() && 
								Tc->get_cell_boundary() &&
								alphav<BD*BD)
							{
								point3D *p1 = NULL;
								point3D *p2 = NULL;
								point3D *p3 = NULL;
								if      (Ta==max_plot_index) {p1 = Ta;p2 = Tb;p3 = Tc;}
								else if (Tb==max_plot_index) {p1 = Tb;p2 = Ta;p3 = Tc;}
								else if (Tc==max_plot_index) {p1 = Tc;p2 = Ta;p3 = Tb;}
								double n2abx = p2->getab()->getnormal(0);
								double n2aby = p2->getab()->getnormal(1);
								double n2bax = p2->getba()->getnormal(0);
								double n2bay = p2->getba()->getnormal(1);
								double n3abx = p3->getab()->getnormal(0);
								double n3aby = p3->getab()->getnormal(1);
								double n3bax = p3->getba()->getnormal(0);
								double n3bay = p3->getba()->getnormal(1);
								double n2abd = sqrt(n2abx*n2abx + n2aby*n2aby);
								double n2bad = sqrt(n2bax*n2bax + n2bay*n2bay);
								double n3abd = sqrt(n3abx*n3abx + n3aby*n3aby);
								double n3bad = sqrt(n3bax*n3bax + n3bay*n3bay);
								double n2x = n2abx/n2abd + n2bax/n2bad;
								double n2y = n2aby/n2abd + n2bay/n2bad;
								double n3x = n3abx/n3abd + n3bax/n3bad;
								double n3y = n3aby/n3abd + n3bay/n3bad;
								double n2d = sqrt(n2x*n2x + n2y*n2y);
								double n3d = sqrt(n3x*n3x + n3y*n3y);
								double n2nx = n2x/n2d;
								double n2ny = n2y/n2d;
								double n3nx = n3x/n3d;
								double n3ny = n3y/n3d;
								double plot_n2 = n2nx*dx + n2ny*dy;
								double plot_n3 = n3nx*dx + n3ny*dy;
								if (plot_n2>0 || plot_n3>0)
								{
									Cn_death = true;break;
								}
							}
						}
					}
				}
			}
			sp = Ip->MP()->at(0)->S2();
			tp = Ip->MP()->at(0)->T2();
		}
		sp_ex = sp->getrp();
		tp_ex = tp->getfp();
		bool flagst = false;
		while (!flagst)
		{
			if      (!sp_ex->getab()->get_attach())
			{
				double nabx = sp_ex->getab()->getnormal(0);
				double naby = sp_ex->getab()->getnormal(1);
				double nbax = sp_ex->getba()->getnormal(0);
				double nbay = sp_ex->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double sp_x = nab_nx + nba_nx;
				double sp_y = nab_ny + nba_ny;
				double sp_d = sqrt(sp_x*sp_x + sp_y*sp_y);
				double sp_nx = sp_x/sp_d;
				double sp_ny = sp_y/sp_d;
				double plot_sp_abd = nab_nx*dx + nab_ny*dy;
				double plot_sp_bad = nba_nx*dx + nba_ny*dy;
				double plot_spd = sp_nx*dx + sp_ny*dy;
				if (plot_sp_abd<0.17 && plot_sp_bad<0.17) {flagst = true;break;}
				else
				{
					sp_ex->set_migration_rate(6*GAR);
					sp_ex = sp_ex->getrp();
				}
			}
			else if (sp_ex->getab()->get_attach())
			{flagst = true;break;}
		}
		bool flagtt = false;
		while (!flagtt)
		{
			if      (!tp_ex->getba()->get_attach())
			{
				double nabx = tp_ex->getab()->getnormal(0);
				double naby = tp_ex->getab()->getnormal(1);
				double nbax = tp_ex->getba()->getnormal(0);
				double nbay = tp_ex->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double tp_x = nab_nx + nba_nx;
				double tp_y = nab_ny + nba_ny;
				double tp_d = sqrt(tp_x*tp_x + tp_y*tp_y);
				double tp_nx = tp_x/tp_d;
				double tp_ny = tp_y/tp_d;
				double plot_tp_abd = nab_nx*dx + nab_ny*dy;
				double plot_tp_bad = nba_nx*dx + nba_ny*dy;
				double plot_tpd = tp_nx*dx + tp_ny*dy;
				if (plot_tp_abd<0.17 && plot_tp_bad<0.17) {flagtt = true;break;}
				else
				{
					tp_ex->set_migration_rate(6*GAR);
					tp_ex = tp_ex->getfp();
				}
			}
			else if (tp_ex->getba()->get_attach())
			{flagtt = true;break;}
		}
		///// push contact cell into migrating dead cell list of C /////
		if (Cn_move) // only for cell which can be degraded
		{
			int mdcn = C->get_migrate_filling_cells()->size();
			if (mdcn==0)
			{
				rear_fc *FC = new rear_fc(mdcn);
				FC->set_angle(angle);
				FC->set_C(C);
				FC->set_center(xc,yc);
				FC->set_initial(1);
				int Cpn = (int)C->get_pairs()->size();
				for (int i=0;i<Cpn;i++)
				{
					interpair *Ipi = C->get_pairs()->at(i);
					if (Ipi->get_redundant()) continue;
					if      (Ipi->get_I1()==C && Ipi->get_I2()!=Cn) {FC->push_NCs(Ipi->get_I2());}
					else if (Ipi->get_I2()==C && Ipi->get_I1()!=Cn) {FC->push_NCs(Ipi->get_I1());}
				}
				C->push_filling_cell(FC);
			}
			mdcn = C->get_migrate_filling_cells()->size();
			if (mdcn>0)
			{
				bool addC = true;
				if (mdcn>1)
				{
					for (int i=1;i<mdcn;i++)
					{
						double ccx = C->get_migrate_filling_cells()->at(i)->get_center(0);
						double ccy = C->get_migrate_filling_cells()->at(i)->get_center(1);
						double distcc = sqrt((ccx - xc)*(ccx - xc) + (ccy - yc)*(ccy - yc));
						if (C->get_migrate_filling_cells()->at(i)->getC()==Cn || distcc<CR)
						{
							addC = false;break;
						}
					}
				}
				else if (mdcn==1)
				{
					if (C->get_migrate_filling_cells()->at(0)->getC()==Cn)
					{addC = false;}
				}
				if (addC)
				{
					rear_fc *FC = new rear_fc(mdcn);
					FC->set_angle(angle);
					FC->set_C(Cn);
					FC->set_center(xc,yc);
					int Cnpn = (int)Cn->get_pairs()->size();
					for (int i=0;i<Cnpn;i++)
					{
						interpair *Ipi = Cn->get_pairs()->at(i);
						if (Ipi->get_redundant()) continue;
						if      (Ipi->get_I1()==Cn && Ipi->get_I2()!=C) {FC->push_NCs(Ipi->get_I2());}
						else if (Ipi->get_I2()==Cn && Ipi->get_I1()!=C) {FC->push_NCs(Ipi->get_I1());}
					}
					C->push_filling_cell(FC);
				}
			}
		}
		if (Cn_death)
		{
			force_control = true;
			Cn->set_set_dead(1);
		}
		else if (!Cn_death && Cn!=NULL)
		{
			if (!Cn->Set_Dead()) {cout<<"    -> cell "<<C->id()<<" migrates towards soften cell "<<Cn->id()<<endl;}
			if (find(d->begin(),d->end(),Cn->get_cell_type())!=d->end())
			{
				Cn->set_soften(1);
				Cn->setup_pcoef(pressure_edge_soften);
				Cn->setup_tcoef(tension_edge_soften);
				Cn->set_lame(0,MU_soft);
				Cn->set_lame(1,LAMBDA_soft);
				int cnpn = (int)Cn->get_pairs()->size();
				for (int i=0;i<cnpn;i++)
				{
					Cn->get_pairs()->at(i)->set_adhesion(adhesion_pair_class[15]);
					Cn->get_pairs()->at(i)->set_friction(friction_pair_class[15]);
				}
			}
		}
	}
	C->set_migrate_p(sp,tp);
	point3D *pm = sp;
	double vol = F;
	double vol_t = 6*GAR;
	if (force_control)
	{
		if (vol>vol_t) {vol = vol_t;}
	}
	while (pm!=tp)
	{
		pm->set_migrate_mark(1);
		pm->set_migration_rate(vol);
		pm = pm->getfp();
	}
	pm->set_migration_rate(vol);
	pm->set_migrate_mark(1);
}

void dbReader::cell_migrate_response_1(cell *C, double angle, vector<int>* d, double F)
{
	/************************************
	                 A
	                 |
	           o---o---o---o
	          /      E      \
	         o               o
	      - -| - - - o- - - -|- -
	    o---oo               oo---o
	   /     \\             //     \
	  o       oo---o---o---oo       o
	  |       ||           ||       |

	************************************/
	C->center_refresh();
	C->setup_Mark_number(1); // 1: migrate
	C->set_interior_refresh(1);
	double xc = C->get_center(0);
	double yc = C->get_center(1);
	edge *E = C->get_sides()->at(0);
	bool force_control = false;
	double dx = 0;
	double dy = 1;
	if (angle==90 || angle==270)
	{
		if (angle==270) {dy = -1;}
		bool flagE = false;
		int sn = (int)C->get_sides()->size();
		for (int i=0;i<sn;i++)
		{
			double x1 = C->get_sides()->at(i)->p1()->x();
			double x2 = C->get_sides()->at(i)->p2()->x();
			if (x1<=xc && x2>xc)
			{
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = dy*ny;
				if (plot_dn>0)
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	else
	{
		dx = cos(angle*PI/180);
		dy = sin(angle*PI/180);
		double slope = tan(angle*PI/180);
		double intercept = yc - slope*xc;
		int sn = (int)C->get_sides()->size();
		bool flagE = false;
		for (int i=0;i<sn;i++)
		{
			point3D *p1 = C->get_sides()->at(i)->p1();
			point3D *p2 = C->get_sides()->at(i)->p2();
			double x1 = p1->x();
			double y1 = p1->y();
			double x2 = p2->x();
			double y2 = p2->y();
			double t1 = y1 - slope*x1 - intercept;
			double t2 = y2 - slope*x2 - intercept;
			if ((t1>=0 && t2<0) || (t1<0 && t2>=0))
			{
				double nx = C->get_sides()->at(i)->getnormal(0);
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = nx*dx + ny*dy;
				if (plot_dn>0) 
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	point3D *sp = E->p1();
	point3D *tp = E->p2();
	point3D *tp_ex = NULL;
	point3D *sp_ex = NULL;
	if (!E->get_attach()) // single
	{
		force_control = true;
		bool flagst = false;
		while (!flagst)
		{
			if      (!sp->getab()->get_attach() && !sp->getba()->get_attach())
			{
				double nabx = sp->getab()->getnormal(0);
				double naby = sp->getab()->getnormal(1);
				double nbax = sp->getba()->getnormal(0);
				double nbay = sp->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double sp_x = nab_nx + nba_nx;
				double sp_y = nab_ny + nba_ny;
				double sp_d = sqrt(sp_x*sp_x + sp_y*sp_y);
				double sp_nx = sp_x/sp_d;
				double sp_ny = sp_y/sp_d;
				double plot_spd = sp_nx*dx + sp_ny*dy;
				if (plot_spd<0.17) {flagst = true;break;} // the angle is less than 80 degrees
				else               {sp = sp->getrp();}
			}
			else if (!sp->getab()->get_attach() && sp->getba()->get_attach())
			{flagst = true;break;}
		}
		bool flagtt = false;
		while (!flagtt)
		{
			if      (!tp->getab()->get_attach() && !tp->getba()->get_attach())
			{
				double nabx = tp->getab()->getnormal(0);
				double naby = tp->getab()->getnormal(1);
				double nbax = tp->getba()->getnormal(0);
				double nbay = tp->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double tp_x = nab_nx + nba_nx;
				double tp_y = nab_ny + nba_ny;
				double tp_d = sqrt(tp_x*tp_x + tp_y*tp_y);
				double tp_nx = tp_x/tp_d;
				double tp_ny = tp_y/tp_d;
				double plot_tpd = tp_nx*dx + tp_ny*dy;
				if (plot_tpd<0.17) {flagtt = true;break;} // the angle is less than 80 degrees
				else              {tp = tp->getfp();}
			}
			else if (tp->getab()->get_attach() && !tp->getba()->get_attach())
			{flagtt = true;break;}
		}
	}
	else // attached
	{
		interpair *Ip = E->get_Pair();
		cell *C1 = Ip->get_I1();
		cell *C2 = Ip->get_I2();
		cell *Cn = NULL; // the neighbor of C: C towards Cn
		bool Cn_death = false;
		bool Cn_move = false; // if the cell can or can't degrade the cell on the way
		if (Ip->get_I1()==C)
		{
			if (find(d->begin(),d->end(),C2->get_cell_type())!=d->end())
			{
				Cn_move = true;
				Cn = C2;
			}
			sp = Ip->MP()->at(0)->S1();
			tp = Ip->MP()->at(0)->T1();
		}
		else
		{
			if (find(d->begin(),d->end(),C1->get_cell_type())!=d->end())
			{
				Cn = C1;
				Cn_move = true;
			}
			sp = Ip->MP()->at(0)->S2();
			tp = Ip->MP()->at(0)->T2();
		}
		sp_ex = sp->getrp();
		tp_ex = tp->getfp();
		bool flagst = false;
		while (!flagst)
		{
			if      (!sp_ex->getab()->get_attach())
			{
				double nabx = sp_ex->getab()->getnormal(0);
				double naby = sp_ex->getab()->getnormal(1);
				double nbax = sp_ex->getba()->getnormal(0);
				double nbay = sp_ex->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double sp_x = nab_nx + nba_nx;
				double sp_y = nab_ny + nba_ny;
				double sp_d = sqrt(sp_x*sp_x + sp_y*sp_y);
				double sp_nx = sp_x/sp_d;
				double sp_ny = sp_y/sp_d;
				double plot_sp_abd = nab_nx*dx + nab_ny*dy;
				double plot_sp_bad = nba_nx*dx + nba_ny*dy;
				double plot_spd = sp_nx*dx + sp_ny*dy;
				if (plot_sp_abd<0.17 && plot_sp_bad<0.17) {flagst = true;break;}
				else
				{
					sp_ex->set_migration_rate(6*GAR);
					sp_ex = sp_ex->getrp();
				}
			}
			else if (sp_ex->getab()->get_attach())
			{flagst = true;break;}
		}
		bool flagtt = false;
		while (!flagtt)
		{
			if      (!tp_ex->getba()->get_attach())
			{
				double nabx = tp_ex->getab()->getnormal(0);
				double naby = tp_ex->getab()->getnormal(1);
				double nbax = tp_ex->getba()->getnormal(0);
				double nbay = tp_ex->getba()->getnormal(1);
				double nabd = sqrt(nabx*nabx + naby*naby);
				double nbad = sqrt(nbax*nbax + nbay*nbay);
				double nab_nx = nabx/nabd;
				double nab_ny = naby/nabd;
				double nba_nx = nbax/nbad;
				double nba_ny = nbay/nbad;
				double tp_x = nab_nx + nba_nx;
				double tp_y = nab_ny + nba_ny;
				double tp_d = sqrt(tp_x*tp_x + tp_y*tp_y);
				double tp_nx = tp_x/tp_d;
				double tp_ny = tp_y/tp_d;
				double plot_tp_abd = nab_nx*dx + nab_ny*dy;
				double plot_tp_bad = nba_nx*dx + nba_ny*dy;
				double plot_tpd = tp_nx*dx + tp_ny*dy;
				if (plot_tp_abd<0.17 && plot_tp_bad<0.17) {flagtt = true;break;}
				else
				{
					tp_ex->set_migration_rate(6*GAR);
					tp_ex = tp_ex->getfp();
				}
			}
			else if (tp_ex->getba()->get_attach())
			{flagtt = true;break;}
		}
		///// push contact cell into migrating dead cell list of C /////
		if (Cn_move) // only for cell which can be degraded
		{
			int mdcn = C->get_migrate_filling_cells()->size();
			if (mdcn==0)
			{
				rear_fc *FC = new rear_fc(mdcn);
				FC->set_angle(angle);
				FC->set_C(C);
				FC->set_center(xc,yc);
				FC->set_initial(1);
				int Cpn = (int)C->get_pairs()->size();
				for (int i=0;i<Cpn;i++)
				{
					interpair *Ipi = C->get_pairs()->at(i);
					if (Ipi->get_redundant()) continue;
					if      (Ipi->get_I1()==C && Ipi->get_I2()!=Cn) {FC->push_NCs(Ipi->get_I2());}
					else if (Ipi->get_I2()==C && Ipi->get_I1()!=Cn) {FC->push_NCs(Ipi->get_I1());}
				}
				C->push_filling_cell(FC);
			}
			mdcn = C->get_migrate_filling_cells()->size();
			if (mdcn>0)
			{
				bool addC = true;
				if (mdcn>1)
				{
					for (int i=1;i<mdcn;i++)
					{
						double ccx = C->get_migrate_filling_cells()->at(i)->get_center(0);
						double ccy = C->get_migrate_filling_cells()->at(i)->get_center(1);
						double distcc = sqrt((ccx - xc)*(ccx - xc) + (ccy - yc)*(ccy - yc));
						if (C->get_migrate_filling_cells()->at(i)->getC()==Cn || distcc<CR)
						{
							addC = false;break;
						}
					}
				}
				else if (mdcn==1)
				{
					if (C->get_migrate_filling_cells()->at(0)->getC()==Cn)
					{addC = false;}
				}
				if (addC)
				{
					rear_fc *FC = new rear_fc(mdcn);
					FC->set_angle(angle);
					FC->set_C(Cn);
					FC->set_center(xc,yc);
					int Cnpn = (int)Cn->get_pairs()->size();
					for (int i=0;i<Cnpn;i++)
					{
						interpair *Ipi = Cn->get_pairs()->at(i);
						if (Ipi->get_redundant()) continue;
						if      (Ipi->get_I1()==Cn && Ipi->get_I2()!=C) {FC->push_NCs(Ipi->get_I2());}
						else if (Ipi->get_I2()==Cn && Ipi->get_I1()!=C) {FC->push_NCs(Ipi->get_I1());}
					}
					C->push_filling_cell(FC);
				}
			}
		}
		if (Cn_death)
		{
			force_control = true;
			Cn->set_set_dead(1);
		}
		else if (!Cn_death && Cn!=NULL)
		{
			//if       (!Cn->Set_Dead() && Cn_move) {cout<<"    -> cell "<<C->id()<<" migrates towards soften cell "<<Cn->id()<<endl;}
			//else if (!Cn->Set_Dead() && !Cn_move) {cout<<"    -> cell "<<C->id()<<" migrates towards solid cell "<<Cn->id()<<endl;}
			if (find(d->begin(),d->end(),Cn->get_cell_type())!=d->end())
			{
				/*Cn->set_soften(1);
				Cn->setup_pcoef(pressure_edge_soften);
				Cn->setup_tcoef(tension_edge_soften);
				Cn->set_lame(0,MU_soft);
				Cn->set_lame(1,LAMBDA_soft);
				int cnpn = (int)Cn->get_pairs()->size();
				for (int i=0;i<cnpn;i++)
				{
					Cn->get_pairs()->at(i)->set_adhesion(adhesion_pair_class[15]);
					Cn->get_pairs()->at(i)->set_friction(friction_pair_class[15]);
				}*/
			}
		}
	}
	C->set_migrate_p(sp,tp);
	point3D *pm = sp;
	double vol = F;
	double vol_t = 6*GAR;
	if (force_control)
	{
		if (vol>vol_t) {vol = vol_t;}
	}
	while (pm!=tp)
	{
		pm->set_migrate_mark(1);
		pm->set_migration_rate(vol);
		pm = pm->getfp();
	}
	pm->set_migration_rate(vol);
	pm->set_migrate_mark(1);
}

cell* dbReader::cell_migrate_angle(cell *C, double angle)
{
	cell* Cr = NULL;
	double xc = C->get_center(0);
	double yc = C->get_center(1);
	edge *E = C->get_sides()->at(0);
	double dx = 0;
	double dy = 1;
	if (angle==90 || angle==270)
	{
		if (angle==270) {dy = -1;}
		bool flagE = false;
		int sn = (int)C->get_sides()->size();
		for (int i=0;i<sn;i++)
		{
			double x1 = C->get_sides()->at(i)->p1()->x();
			double x2 = C->get_sides()->at(i)->p2()->x();
			if (x1<=xc && x2>xc)
			{
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = dy*ny;
				if (plot_dn>0)
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	else
	{
		dx = cos(angle*PI/180);
		dy = sin(angle*PI/180);
		double slope = tan(angle*PI/180);
		double intercept = yc - slope*xc;
		int sn = (int)C->get_sides()->size();
		bool flagE = false;
		for (int i=0;i<sn;i++)
		{
			point3D *p1 = C->get_sides()->at(i)->p1();
			point3D *p2 = C->get_sides()->at(i)->p2();
			double x1 = p1->x();
			double y1 = p1->y();
			double x2 = p2->x();
			double y2 = p2->y();
			double t1 = y1 - slope*x1 - intercept;
			double t2 = y2 - slope*x2 - intercept;
			if ((t1>=0 && t2<0) || (t1<0 && t2>=0))
			{
				double nx = C->get_sides()->at(i)->getnormal(0);
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = nx*dx + ny*dy;
				if (plot_dn>0) 
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	if (E->get_attach())
	{
		if (E->get_Pair()->get_I1()==C) {Cr = E->get_Pair()->get_I2();}
		else                            {Cr = E->get_Pair()->get_I1();}
	}
	return Cr;
}

void dbReader::edge_response_b0(edge *E)
{
	/**************************
	 E breaks at p1 and p4
	---------------------------
			 o
			 |
		 o->-o
		 o-<-oo-<-o
		  p2 || p3
			 ||
		C1 E AV E1 C2
			 ||
		  p1 || p4
		 o->-oo->-o
			  o-<-o
			  |
			  o
	**************************/
	if      (!E->get_attach()) {cout<<"NULL"<<endl;}
	else if (E->get_attach())
	{
		edge *E1 = E->get_Neighbor();
		point3D *p1 = E->p1();
		point3D *p2 = E->p2();
		point3D *p3 = E1->p1();
		point3D *p4 = E1->p2();
		cell *C1 = E->Cell();
		cell *C2 = E1->Cell();
		//cout<<"edge break between Cell "<<C1->id()<<" and Cell "<<C2->id()<<" at edge "<<p1->id()<<"-"<<p2->id()<<endl;
		interpair *Pl1 = E->get_Pair();
		double break_force = Pl1->get_adhesion()/20;
		bool p23test = false;
		if (!p2->getab()->get_attach() || !p3->getba()->get_attach()) {p23test = true;}
		if ((Pl1->MP()->at(0)->S1()==p1 && Pl1->MP()->at(0)->T1()==p2 && p23test) ||
			(Pl1->MP()->at(0)->S2()==p1 && Pl1->MP()->at(0)->T2()==p2 && p23test))
		{
			if (!p1->getba()->get_attach() && !p4->getab()->get_attach())
			{
				p1->removemp(p4);
				p4->removemp(p1);
				node_remove_0(p4);
			}
			else if (!p1->getba()->get_attach() && p4->getab()->get_attach())
			{
				/***************************
				C1      C2
				   oo
				p1 || p4
			   o->-oo->-o
			        o-<-o
			        | p5
				    o    C3
				***************************/
				point3D *p5 = p4->getab()->get_Neighbor()->p2();
				cell *C3 = p5->getab()->Cell();
				p1->removemp(p4);p1->removemp(p5);
				p4->removemp(p1);
				p5->removemp(p1);
				for (int i=0;i<(int)C1->get_pairs()->size();i++)
				{
					interpair *C1pa = C1->get_pairs()->at(i);
					if ((C1pa->get_I1()==C1 && C1pa->get_I2()==C3) ||
						(C1pa->get_I2()==C1 && C1pa->get_I1()==C3)) {C1pa->remove_MP(C1pa->MP()->at(0));break;}
				}
				node_remove_0(p1);
			}
			else if (p1->getba()->get_attach() && !p4->getab()->get_attach())
			{
				/***************************
				C1      C2
				    oo
				 p1 || p4
				o->-oo->-o
				o-<-o
				p5	|
				C3	o
				***************************/
				point3D *p5 = p1->getba()->get_Neighbor()->p1();
				cell *C3 = p5->getab()->Cell();
				p1->removemp(p4);
				p5->removemp(p4);
				p4->removemp(p1);p4->removemp(p5);
				for (int i=0;i<(int)C2->get_pairs()->size();i++)
				{
					interpair *C2pa = C2->get_pairs()->at(i);
					if ((C2pa->get_I1()==C2 && C2pa->get_I2()==C3) ||
						(C2pa->get_I2()==C2 && C2pa->get_I1()==C3)) {C2pa->remove_MP(C2pa->MP()->at(0));break;}
				}
				node_remove_0(p4);
			}
			if (!p2->getab()->get_attach() && !p3->getba()->get_attach())
			{
				p3->removemp(p2);p2->removemp(p3);
				node_remove_0(p3);
			}
			else if (!p2->getab()->get_attach() && p3->getba()->get_attach())
			{
				/*********************************
				     o   C3
				     | p5
				     o->-o
				o-<-oo-<-o
				 p2 || p3
				C1	oo  C2
				*********************************/
				point3D *p5 = p3->getba()->get_Neighbor()->p1();
				cell *C3 = p5->getab()->Cell();
				p2->removemp(p3);p2->removemp(p5);
				p3->removemp(p2);
				p5->removemp(p2);
				for (int i=0;i<(int)C1->get_pairs()->size();i++)
				{
					interpair *C1pa = C1->get_pairs()->at(i);
					if ((C1pa->get_I1()==C1 && C1pa->get_I2()==C3) ||
						(C1pa->get_I2()==C1 && C1pa->get_I1()==C3)) {C1pa->remove_MP(C1pa->MP()->at(0));break;}
				}
				node_remove_0(p2);
			}
			else if (p2->getab()->get_attach() && !p3->getba()->get_attach())
			{
				/*********************************
			  C3    o
				 p5	|
				o->-o
				o-<-oo-<-o
				 p2 || p3
				C1	oo   C2
				*********************************/
				point3D *p5 = p2->getab()->get_Neighbor()->p2();
				cell *C3 = p5->getab()->Cell();
				p2->removemp(p3);
				p3->removemp(p2);p3->removemp(p5);
				p5->removemp(p3);
				for (int i=0;i<(int)C2->get_pairs()->size();i++)
				{
					interpair *C2pa = C2->get_pairs()->at(i);
					if ((C2pa->get_I1()==C2 && C2pa->get_I2()==C3) ||
						(C2pa->get_I2()==C2 && C2pa->get_I1()==C3)) {C2pa->remove_MP(C2pa->MP()->at(0));break;}
				}
				node_remove_0(p3);
			}
			E->set_attach(0);E1->set_attach(0);
			E->set_Neighbor(NULL);E1->set_Neighbor(NULL);
			Pl1->remove_MP(Pl1->MP()->at(0));
			
			double p1_tension_x = p1->getTension_ab(0) + p1->getTension_ba(0);
			double p1_tension_y = p1->getTension_ab(1) + p1->getTension_ba(1);
			double p2_tension_x = p2->getTension_ab(0) + p2->getTension_ba(0);
			double p2_tension_y = p2->getTension_ab(1) + p2->getTension_ba(1);
			double p3_tension_x = p3->getTension_ab(0) + p3->getTension_ba(0);
			double p3_tension_y = p3->getTension_ab(1) + p3->getTension_ba(1);
			double p4_tension_x = p4->getTension_ab(0) + p4->getTension_ba(0);
			double p4_tension_y = p4->getTension_ab(1) + p4->getTension_ba(1);
			double ld = 0.1;
			double p1_dx = ld*p1_tension_x + p1->x();
			double p1_dy = ld*p1_tension_y + p1->y();
			double p2_dx = ld*p2_tension_x + p2->x();
			double p2_dy = ld*p2_tension_y + p2->y();
			double p3_dx = ld*p3_tension_x + p3->x();
			double p3_dy = ld*p3_tension_y + p3->y();
			double p4_dx = ld*p4_tension_x + p4->x();
			double p4_dy = ld*p4_tension_y + p4->y();
			p1->setx(p1_dx);
			p1->sety(p1_dy);
			for (int i=0;i<(int)p1->getmp()->size();i++) {p1->getmp()->at(i)->setx(p1_dx);p1->getmp()->at(i)->sety(p1_dy);}
			p2->setx(p2_dx);
			p2->sety(p2_dy);
			for (int i=0;i<(int)p2->getmp()->size();i++) {p2->getmp()->at(i)->setx(p2_dx);p2->getmp()->at(i)->sety(p2_dy);}
			p3->setx(p3_dx);
			p3->sety(p3_dy);
			for (int i=0;i<(int)p3->getmp()->size();i++) {p3->getmp()->at(i)->setx(p3_dx);p3->getmp()->at(i)->sety(p3_dy);}
			p4->setx(p4_dx);
			p4->sety(p4_dy);
			for (int i=0;i<(int)p4->getmp()->size();i++) {p4->getmp()->at(i)->setx(p4_dx);p4->getmp()->at(i)->sety(p4_dy);}
			
			double E_nx = -E->getnormal(0);
			double E_ny = -E->getnormal(1);
			double E_nd = sqrt(E_nx*E_nx + E_ny*E_ny);
			double E_normal_x = E_nx/E_nd;
			double E_normal_y = E_ny/E_nd;
			double bf_x = break_force*E_normal_x;
			double bf_y = break_force*E_normal_y;
			double bf_x1 = bf_x + p1->get_adhesion_break(0);
			double bf_y1 = bf_y + p1->get_adhesion_break(1);
			double bf_x2 = bf_x + p2->get_adhesion_break(0);
			double bf_y2 = bf_y + p2->get_adhesion_break(1);
			double bf_x3 = -bf_x + p3->get_adhesion_break(0);
			double bf_y3 = -bf_y + p3->get_adhesion_break(1);
			double bf_x4 = -bf_x + p4->get_adhesion_break(0);
			double bf_y4 = -bf_y + p4->get_adhesion_break(1);
			p1->set_adhesion_break(bf_x1, bf_y1);
			p2->set_adhesion_break(bf_x2, bf_y2);
			p3->set_adhesion_break(bf_x3, bf_y3);
			p4->set_adhesion_break(bf_x4, bf_y4);
		}
		else
		{
			if (!p1->getba()->get_attach() && !p4->getab()->get_attach())
			{
				p1->removemp(p4);
				p4->removemp(p1);
				node_remove_0(p4);
			}
			else if (!p1->getba()->get_attach() && p4->getab()->get_attach())
			{
				/********************************
				    oo
				C1`	||    C2
				p1	|| p4
				o->-oo->-o
				     o-<-o
				     | p5
				     o   C3
				********************************/
				point3D *p5 = p4->getab()->get_Neighbor()->p2();
				cell *C3 = p5->getab()->Cell();
				p1->removemp(p4);p1->removemp(p5);
				p4->removemp(p1);
				p5->removemp(p1);
				for (int i=0;i<(int)C1->get_pairs()->size();i++)
				{
					interpair *C1pa = C1->get_pairs()->at(i);
					if ((C1pa->get_I1()==C1 && C1pa->get_I2()==C3) ||
						(C1pa->get_I2()==C1 && C1pa->get_I1()==C3))
					{
						C1pa->remove_MP(C1pa->MP()->at(0));
						break;
					}
				}
				node_remove_0(p1);
			}
			else if (p1->getba()->get_attach() && !p4->getab()->get_attach())
			{
				/********************************
				    oo
				C1`	||    C2
				p1	|| p4
				o->-oo->-o
				o-<-o
				 p5 |
				C3	o
				********************************/
				point3D *p5 = p1->getba()->get_Neighbor()->p1();
				cell *C3 = p5->getab()->Cell();
				p1->removemp(p4);
				p5->removemp(p4);
				p4->removemp(p1);p4->removemp(p5);
				for (int i=0;i<(int)C2->get_pairs()->size();i++)
				{
					interpair *C2pa = C2->get_pairs()->at(i);
					if ((C2pa->get_I1()==C2 && C2pa->get_I2()==C3) ||
						(C2pa->get_I2()==C2 && C2pa->get_I1()==C3)) {C2pa->remove_MP(C2pa->MP()->at(0));break;}
				}
				node_remove_0(p4);
			}
			if (Pl1->MP()->at(0)->S1()==p1)
			{
				Pl1->MP()->at(0)->setS1(p2);
				Pl1->MP()->at(0)->setT2(p3);
				double sF_p2 = p2->get_stored_apart_force();
				double sF_p3 = p3->get_stored_apart_force();
				Pl1->MP()->at(0)->set_F_s1(sF_p2);
				Pl1->MP()->at(0)->set_F_t2(sF_p3);
			}
			if (Pl1->MP()->at(0)->S2()==p1)
			{
				Pl1->MP()->at(0)->setS2(p2);
				Pl1->MP()->at(0)->setT1(p3);
				double sF_p2 = p2->get_stored_apart_force();
				double sF_p3 = p3->get_stored_apart_force();
				Pl1->MP()->at(0)->set_F_s2(sF_p2);
				Pl1->MP()->at(0)->set_F_t1(sF_p3);
			}
			E->set_attach(0);
			E1->set_attach(0);
			E->set_Neighbor(NULL);
			E1->set_Neighbor(NULL);
			double p1_tension_x = p1->getTension_ab(0) + p1->getTension_ba(0);
			double p1_tension_y = p1->getTension_ab(1) + p1->getTension_ba(1);
			double p4_tension_x = p4->getTension_ab(0) + p4->getTension_ba(0);
			double p4_tension_y = p4->getTension_ab(1) + p4->getTension_ba(1);
			double ld = 0.1;
			double p1_dx = ld*p1_tension_x + p1->x();
			double p1_dy = ld*p1_tension_y + p1->y();
			double p4_dx = ld*p4_tension_x + p4->x();
			double p4_dy = ld*p4_tension_y + p4->y();
			p1->setx(p1_dx);
			p1->sety(p1_dy);
			for (int i=0;i<(int)p1->getmp()->size();i++) {p1->getmp()->at(i)->setx(p1_dx);p1->getmp()->at(i)->sety(p1_dy);}
			p4->setx(p4_dx);
			p4->sety(p4_dy);
			for (int i=0;i<(int)p4->getmp()->size();i++) {p4->getmp()->at(i)->setx(p4_dx);p4->getmp()->at(i)->sety(p4_dy);}
			double E_nx = -E->getnormal(0);
			double E_ny = -E->getnormal(1);
			double E_nd = sqrt(E_nx*E_nx + E_ny*E_ny);
			double E_normal_x = E_nx/E_nd;
			double E_normal_y = E_ny/E_nd;
			double bf_x = break_force*E_normal_x;
			double bf_y = break_force*E_normal_y;
			double bf_x1 = bf_x + p1->get_adhesion_break(0);
			double bf_y1 = bf_y + p1->get_adhesion_break(1);
			double bf_x4 = -bf_x + p4->get_adhesion_break(0);
			double bf_y4 = -bf_y + p4->get_adhesion_break(1);
			p1->set_adhesion_break(bf_x1, bf_y1);
			p4->set_adhesion_break(bf_x4, bf_y4);
		}
	}
}

void dbReader::edge_response_b1(edge *E)
{
	/**************************
	 E breaks at p1 and p4
	 where C3 protrudes into
	 C1 and C2;
	---------------------------
	         C4
	      \\ p6 //
           \\  //
	        \\//
	      p2 oo p3  
	         || 
	    C1 E AV E1 C2
	         || 
	      p1 oo p4    
	        //\\ 
           //  \\ 
          // p5 \\ 
             C3
	**************************/
	if      (!E->get_attach()) {cout<<"  -> single edge no action is needed!"<<endl;}
	else if (E->get_attach())
	{
		edge *E1 = E->get_Neighbor();
		point3D *p1 = E->p1();
		point3D *p2 = E->p2();
		if (p1->getfp()==p2 && p2->getrp()==p1)
		{
			if (p1->get_in_cell() && p2->get_in_cell())
			{
				point3D *p3 = E1->p1();
				point3D *p4 = E1->p2();
				if      (p1->getba()->get_attach() && !p4->getab()->get_attach())
				{
					edge_response_b0(E);
				}
				else if (p1->getba()->get_attach() && p4->getab()->get_attach())
				{
					point3D *p5 = p1->getba()->get_Neighbor()->p1();
					cell *C1 = E->Cell();
					cell *C2 = E1->Cell();
					cell *C3 = p1->getba()->get_Neighbor()->Cell();
					cell *C3c = p4->getab()->get_Neighbor()->Cell();
					//cout<<"edge breaks between cell "<<C1->id()<<" and cell "<<C2->id()<<" at point "<<p1->id()<<endl;
					if (C3!=C3c)
					{
						interpair *Pl1 = E->get_Pair();
						if (Pl1->MP()->at(0)->S1()==p1)
						{
							double fs1b = Pl1->MP()->at(0)->get_F_s1()*0.5;
							double ft2b = Pl1->MP()->at(0)->get_F_t2()*0.5;
							Pl1->MP()->at(0)->set_F_s1(fs1b);
							Pl1->MP()->at(0)->set_F_t2(ft2b);
						}
					}
					else
					{
						interpair *Pl1 = E->get_Pair();
						interpair *Pl2 = p1->getba()->get_Pair();
						interpair *Pl3 = p4->getab()->get_Pair();
						double xm = (p1->x() + p2->x())*0.5;
						double ym = (p1->y() + p2->y())*0.5;
						double x1r = p1->getrp()->x();
						double y1r = p1->getrp()->y();
						double x4f = p4->getfp()->x();
						double y4f = p4->getfp()->y();
						double x2f = p2->getfp()->x();
						double y2f = p2->getfp()->y();
						double x3r = p3->getrp()->x();
						double y3r = p3->getrp()->y();
						double el1r = sqrt((xm - x1r)*(xm - x1r) + (ym - y1r)*(ym - y1r));
						double el2f = sqrt((xm - x2f)*(xm - x2f) + (ym - y2f)*(ym - y2f));
						double el4f = sqrt((xm - x4f)*(xm - x4f) + (ym - y4f)*(ym - y4f));
						double el3r = sqrt((xm - x3r)*(xm - x3r) + (ym - y3r)*(ym - y3r));
						int p2nei = p2->getmp()->size();
						bool burried = false;
						if (p2nei>=2)
						{
							point3D *p6 = NULL;
							for (int i6=0;i6<(int)p2->getmp()->size();i6++)
							{
								if (p2->getmp()->at(i6)!=p3 &&
									p2->getmp()->at(i6)!=p2) 
								{p6 = p2->getmp()->at(i6);break;}
							}
							cell *C4 = p6->getab()->Cell();
							int C3pn = (int)C3->get_pairs()->size();
							for (int i6=0;i6<C3pn;i6++)
							{
								if (C3->get_pairs()->at(i6)->get_redundant()) continue;
								if ((C3->get_pairs()->at(i6)->get_I1()==C3 && C3->get_pairs()->at(i6)->get_I2()==C4) ||
									(C3->get_pairs()->at(i6)->get_I1()==C4 && C3->get_pairs()->at(i6)->get_I2()==C3))
								{
									if (C3->get_pairs()->at(i6)->MP()->size()>0) {burried = true;break;}
								}
							}
						}
						if (!burried)
						{
							if ((p2nei>=2 && el1r<1.8*BD && el2f<1.8*BD && el4f<1.8*BD && el3r<1.8*BD) ||
								(p2nei==1))
							{
								p5->removemp(p1);
								p5->removemp(p4);
								node_remove_1(p1);
								node_remove_1(p4);
								p1->getrp()->setfp(p2);
								p2->setrp(p1->getrp());
								E->reset_p1(p1->getrp());
								C1->removeside(p1->getba());
								p1->getrp()->setab(E);
								p4->getfp()->setrp(p3);
								p3->setfp(p4->getfp());
								E1->reset_p2(p4->getfp());
								C2->removeside(p4->getab());
								p4->getfp()->setba(E1);
								E->set_Neighbor(p5->getab());
								p5->getab()->set_Neighbor(E);
								E1->set_Neighbor(p5->getba());
								p5->getba()->set_Neighbor(E1);
								E->set_Pair(Pl2);
								E1->set_Pair(Pl3);
								p5->pushmp(p2);p2->pushmp(p5);
								p5->pushmp(p3);p3->pushmp(p5);
								node_merge(p5,p2);
								if (Pl1->MP()->at(0)->S1()==p1) 
								{
									Pl1->MP()->at(0)->setS1(p2);
									Pl1->MP()->at(0)->setT2(p3);
									if (Pl1->MP()->at(0)->T1()!=p2)
									{
										double sF_s1 = p2->get_stored_apart_force();
										double sF_t2 = p3->get_stored_apart_force();
										Pl1->MP()->at(0)->set_F_s1(sF_s1);
										Pl1->MP()->at(0)->set_F_t2(sF_t2);
									}
									else
									{
										double sF_s1 = Pl1->MP()->at(0)->get_F_t1();
										double sF_t2 = Pl1->MP()->at(0)->get_F_s2();
										Pl1->MP()->at(0)->set_F_s1(sF_s1);
										Pl1->MP()->at(0)->set_F_t2(sF_t2);
									}
								}
								if (Pl1->MP()->at(0)->S2()==p1)
								{
									Pl1->MP()->at(0)->setT1(p3);
									Pl1->MP()->at(0)->setS2(p2);
									if (Pl1->MP()->at(0)->T2()!=p2)
									{
										double sF_t1 = p3->get_stored_apart_force();
										double sF_s2 = p2->get_stored_apart_force();
										Pl1->MP()->at(0)->set_F_t1(sF_t1);
										Pl1->MP()->at(0)->set_F_s2(sF_s2);
									}
									else
									{
										double sF_s2 = Pl1->MP()->at(0)->get_F_t2();
										double sF_t1 = Pl1->MP()->at(0)->get_F_s1();
										Pl1->MP()->at(0)->set_F_s2(sF_s2);
										Pl1->MP()->at(0)->set_F_t1(sF_t1);
									}
								}
								if (Pl2->MP()->at(0)->T1()==p1) 
								{
									Pl2->MP()->at(0)->setT1(p2);
									Pl2->MP()->at(0)->set_F_t1(0);
									Pl2->MP()->at(0)->set_F_s2(0);
								}
								if (Pl2->MP()->at(0)->T2()==p1) 
								{
									Pl2->MP()->at(0)->setT2(p2);
									Pl2->MP()->at(0)->set_F_t2(0);
									Pl2->MP()->at(0)->set_F_s1(0);
								}
								if (Pl3->MP()->at(0)->S1()==p4) 
								{
									Pl3->MP()->at(0)->setS1(p3);
									Pl3->MP()->at(0)->set_F_s1(0);
									Pl3->MP()->at(0)->set_F_t2(0);
								}
								if (Pl3->MP()->at(0)->S2()==p4) 
								{
									Pl3->MP()->at(0)->setS2(p3);
									Pl3->MP()->at(0)->set_F_s2(0);
									Pl3->MP()->at(0)->set_F_t1(0);
								}
								p2->setx(xm);
								p2->sety(ym);
								int p2mn = (int)p2->getmp()->size();
								if (p2mn>2)
								{
									for (int i=0;i<p2mn;i++)
									{
										point3D *Pt = p2->getmp()->at(i);
										Pt->setx(xm);
										Pt->sety(ym);
										if (Pt!=p5 && Pt!=p3)
										{
											p5->pushmp(Pt);
											Pt->pushmp(p5);
											cell *Ct = Pt->getba()->Cell();
											int cp3n = (int)C3->get_pairs()->size();
											int flag3n = false;
											for (int j=0;j<cp3n;j++)
											{
												if ((C3->get_pairs()->at(j)->get_I1()==C3 && C3->get_pairs()->at(j)->get_I2()==Ct) ||
													(C3->get_pairs()->at(j)->get_I2()==C3 && C3->get_pairs()->at(j)->get_I1()==Ct))
												{
													if (C3->get_pairs()->at(j)->MP()->size()==0)
													{
														if (C3->get_pairs()->at(j)->get_I1()==C3)
														{
															mergepair *Mp = new mergepair(0,p5,p5,Pt,Pt);
															C3->get_pairs()->at(j)->push_MP(Mp);
														}
														else
														{
															mergepair *Mp = new mergepair(0,Pt,Pt,p5,p5);
															C3->get_pairs()->at(j)->push_MP(Mp);
														}
													}
													else
													{
														if (C3->get_pairs()->at(j)->get_I1()==C3)
														{
															C3->get_pairs()->at(j)->MP()->at(0)->setS1(p5);
															C3->get_pairs()->at(j)->MP()->at(0)->setT1(p5);
															C3->get_pairs()->at(j)->MP()->at(0)->setS2(Pt);
															C3->get_pairs()->at(j)->MP()->at(0)->setT2(Pt);
														}
														else
														{
															C3->get_pairs()->at(j)->MP()->at(0)->setS1(Pt);
															C3->get_pairs()->at(j)->MP()->at(0)->setT1(Pt);
															C3->get_pairs()->at(j)->MP()->at(0)->setS2(p5);
															C3->get_pairs()->at(j)->MP()->at(0)->setT2(p5);
														}
													}
													flag3n = true;
													break;
												}
											}
											if (!flag3n)
											{
												int pair_n = (int)collisionpairList.size();
												interpair *PL = new interpair(pair_n,C3,Ct);
												collisionpairList.push_back(PL);
												//int iptype = interpair_type(C3,Ct);
												//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[iptype]);
												double adhesion_force = adhesion_pair_array[C3->get_cell_type()][Ct->get_cell_type()];
												collisionpairList[pair_n]->set_adhesion(adhesion_force);
												collisionpairList[pair_n]->set_friction(adhesion_force);
												mergepair *Mp = new mergepair(0,p5,p5,Pt,Pt);
												collisionpairList[pair_n]->push_MP(Mp);
												//collisionpairList[pair_n]->set_adhesion(adhesion_pair_class[0]);
												C3->pushpair(collisionpairList[pair_n]);
												Ct->pushpair(collisionpairList[pair_n]);
											}
										}
									}
									double p2ab = p2->getab()->getlength();
									if (p2ab>1.8*BD) {edge_response_d(p2->getab());}
									double p2ba = p2->getba()->getlength();
									if (p2ba>1.8*BD) {edge_response_d(p2->getba());}
									double p3ab = p3->getab()->getlength();
									if (p3ab>1.8*BD) {edge_response_d(p3->getab());}
									double p3ba = p3->getba()->getlength();
									if (p3ba>1.8*BD) {edge_response_d(p3->getba());}
									node *node_v = p2->get_node_id();
									int node_vn = (int)node_v->getvertex()->size();
									point3D *np0 = node_v->getvertex()->at(0);
									double angle0 = np0->point_angle();
									for (int jj=1;jj<node_vn;jj++)
									{
										point3D *npjj = node_v->getvertex()->at(jj);
										double anglejj = npjj->point_angle();
										if (angle0>anglejj) {np0 = npjj;angle0 = anglejj;}
									}
									cout<<"    -> additional ";
									int np0_p = pick_correction(np0);
									if (np0_p>=0) {cout<<"  -> error: pick correction collapsed!"<<endl;}
								}
								else
								{
									p3->setx(xm);
									p3->sety(ym);
									p5->setx(xm);
									p5->sety(ym);
									double p2ab = p2->getab()->getlength();
									if (p2ab>1.8*BD) {edge_response_d(p2->getab());}
									double p2ba = p2->getba()->getlength();
									if (p2ba>1.8*BD) {edge_response_d(p2->getba());}
									double p3ab = p3->getab()->getlength();
									if (p3ab>1.8*BD) {edge_response_d(p3->getab());}
									double p3ba = p3->getba()->getlength();
									if (p3ba>1.8*BD) {edge_response_d(p3->getba());}
								}
							}
							else if (p2nei>=2 && (el1r>1.8*BD || el2f>1.8*BD || el4f>1.8*BD || el3r>1.8*BD))
							{
								//cout<<"    -> edge breaks dismissed!"<<endl;
							}
						}
						else
						{
							cout<<"    -> burried case!"<<endl;
						}
					}
				}
			}
		}
		else
		{
			cout<<"  -> dead edge no action is needed!"<<endl;
		}
	}
}

void dbReader::point_response_b0(point3D *V)
{
	/******************************
	         \    /
			  \  /
			V  oo  V1
		C1    /  \     C2
		     /    \
	      
		    \    /
	C1	  V	 \  /  
		  o->-oo  V1  C2
		  o-<-o \
		     /   \
			/

			\    /
			 \  /  V1  C2
		   V  oo->-o
	C1		 / o-<-o
			/   \
			     \
	******************************/
	int Vmn = (int)V->getmp()->size();
	if (Vmn==1)
	{
		point3D *V1 = V->getmp()->at(0);
		cell *C1 = V->getab()->Cell();
		cell *C2 = V1->getab()->Cell();
		//cout<<"point breaks between cell "<<C1->id()<<" and cell "<<C2->id()<<" at point "<<V->id()<<endl;
		V->removemp(V1);
		V1->removemp(V);
		node_remove_0(V1);
		double V_abnx = V->getab()->getnormal(0);
		double V_abny = V->getab()->getnormal(1);
		double V_abnd = sqrt(V_abnx*V_abnx + V_abny*V_abny);
		double V_banx = V->getba()->getnormal(0);
		double V_bany = V->getba()->getnormal(1);
		double V_band = sqrt(V_banx*V_banx + V_bany*V_bany);
		double V_ab_normal_x = V_abnx/V_abnd;
		double V_ab_normal_y = V_abny/V_abnd;
		double V_ba_normal_x = V_banx/V_band;
		double V_ba_normal_y = V_bany/V_band;
		double V_conx = V_ab_normal_x + V_ba_normal_x;
		double V_cony = V_ab_normal_y + V_ba_normal_y;
		double V_cond = sqrt(V_conx*V_conx + V_cony*V_cony);
		double V_co_normal_x = V_conx/V_cond;
		double V_co_normal_y = V_cony/V_cond;
		double V1_abnx = V1->getab()->getnormal(0);
		double V1_abny = V1->getab()->getnormal(1);
		double V1_abnd = sqrt(V1_abnx*V1_abnx + V1_abny*V1_abny);
		double V1_banx = V1->getba()->getnormal(0);
		double V1_bany = V1->getba()->getnormal(1);
		double V1_band = sqrt(V1_banx*V1_banx + V1_bany*V1_bany);
		double V1_ab_normal_x = V1_abnx/V1_abnd;
		double V1_ab_normal_y = V1_abny/V1_abnd;
		double V1_ba_normal_x = V1_banx/V1_band;
		double V1_ba_normal_y = V1_bany/V1_band;
		double V1_conx = V1_ab_normal_x + V1_ba_normal_x;
		double V1_cony = V1_ab_normal_y + V1_ba_normal_y;
		double V1_cond = sqrt(V1_conx*V1_conx + V1_cony*V1_cony);
		double V1_co_normal_x = V1_conx/V1_cond;
		double V1_co_normal_y = V1_cony/V1_cond;
		double V_tension_x = V->getTension_ab(0) + V->getTension_ba(0);
		double V_tension_y = V->getTension_ab(1) + V->getTension_ba(1);
		double V1_tension_x = V1->getTension_ab(0) + V1->getTension_ba(0);
		double V1_tension_y = V1->getTension_ab(1) + V1->getTension_ba(1);
		double ld = 0.1;
		double V_dx = ld*V_tension_x + V->x();
		double V_dy = ld*V_tension_y + V->y();
		double V1_dx = ld*V1_tension_x + V1->x();
		double V1_dy = ld*V1_tension_y + V1->y();
		V->setx(V_dx);
		V->sety(V_dy);
		V1->setx(V1_dx);
		V1->sety(V1_dy);
		int c1pn = (int)C1->get_pairs()->size();
		for (int i=0;i<c1pn;i++)
		{
			interpair *Pl = C1->get_pairs()->at(i);
			if ((Pl->get_I1()==C1 && Pl->get_I2()==C2) ||
				(Pl->get_I2()==C1 && Pl->get_I1()==C2))
			{
				if (Pl->MP()->size()==0)
				{
					cout<<"    -> error: pair "<<Pl->id()<<": "<<C1->id()<<","<<C2->id()<<" was empty!"<<endl;
				}
				mergepair *Mp = Pl->MP()->at(0);
				double break_force = Pl->get_adhesion()/20;
				double Vbf_x = break_force*-V_co_normal_x;
				double Vbf_y = break_force*-V_co_normal_y;
				double V1bf_x = break_force*-V1_co_normal_x;
				double V1bf_y = break_force*-V1_co_normal_y;
				V->set_adhesion_break(Vbf_x, Vbf_y);
				V1->set_adhesion_break(V1bf_x, V1bf_y);
				Pl->remove_MP(Mp);
				break;
			}
		}
	}
	else
	{
		cell *C = V->getab()->Cell();
		vector<point3D*> Vmps;
		for (int i=0;i<Vmn;i++)
		{
			Vmps.push_back(V->getmp()->at(i));
			V->getmp()->at(i)->removemp(V);
		}
		V->clearmp();
		node_remove_0(V);
		double V_abnx = V->getab()->getnormal(0);
		double V_abny = V->getab()->getnormal(1);
		double V_abnd = sqrt(V_abnx*V_abnx + V_abny*V_abny);
		double V_banx = V->getba()->getnormal(0);
		double V_bany = V->getba()->getnormal(1);
		double V_band = sqrt(V_banx*V_banx + V_bany*V_bany);
		double V_ab_normal_x = V_abnx/V_abnd;
		double V_ab_normal_y = V_abny/V_abnd;
		double V_ba_normal_x = V_banx/V_band;
		double V_ba_normal_y = V_bany/V_band;
		double V_conx = V_ab_normal_x + V_ba_normal_x;
		double V_cony = V_ab_normal_y + V_ba_normal_y;
		double V_cond = sqrt(V_conx*V_conx + V_cony*V_cony);
		double V_co_normal_x = V_conx/V_cond;
		double V_co_normal_y = V_cony/V_cond;
		double V_tension_x = V->getTension_ab(0) + V->getTension_ba(0);
		double V_tension_y = V->getTension_ab(1) + V->getTension_ba(1);
		double ld = 0.1;
		double V_dx = ld*V_tension_x + V->x();
		double V_dy = ld*V_tension_y + V->y();
		V->setx(V_dx);
		V->sety(V_dy);
		double break_force = 0;
		for (int j=0;j<(int)Vmps.size();j++)
		{
			point3D *V1 = Vmps[j];
			cell *C1 = V1->getab()->Cell();
			int c1pn = (int)C1->get_pairs()->size();
			for (int i=0;i<c1pn;i++)
			{
				interpair *Pl = C1->get_pairs()->at(i);
				if (Pl->get_redundant()) continue;
				if ((Pl->get_I1()==C1 && Pl->get_I2()==C) ||
					(Pl->get_I2()==C1 && Pl->get_I1()==C))
				{
					break_force = Pl->get_adhesion()/20;
					if (Pl->MP()->size()==0)
					{
						cout<<"    -> error: pair "<<Pl->id()<<": "<<C1->id()<<","<<C->id()<<" was empty!"<<endl;
					}
					mergepair *Mp = Pl->MP()->at(0);
					double V1_abnx = V1->getab()->getnormal(0);
					double V1_abny = V1->getab()->getnormal(1);
					double V1_abnd = sqrt(V1_abnx*V1_abnx + V1_abny*V1_abny);
					double V1_banx = V1->getba()->getnormal(0);
					double V1_bany = V1->getba()->getnormal(1);
					double V1_band = sqrt(V1_banx*V1_banx + V1_bany*V1_bany);
					double V1_ab_normal_x = V1_abnx/V1_abnd;
					double V1_ab_normal_y = V1_abny/V1_abnd;
					double V1_ba_normal_x = V1_banx/V1_band;
					double V1_ba_normal_y = V1_bany/V1_band;
					double V1_conx = V1_ab_normal_x + V1_ba_normal_x;
					double V1_cony = V1_ab_normal_y + V1_ba_normal_y;
					double V1_cond = sqrt(V1_conx*V1_conx + V1_cony*V1_cony);
					double V1_co_normal_x = V1_conx/V1_cond;
					double V1_co_normal_y = V1_cony/V1_cond;
					double V1bf_x = break_force*-V1_co_normal_x;
					double V1bf_y = break_force*-V1_co_normal_y;
					V1->set_adhesion_break(V1bf_x, V1bf_y);
					Pl->remove_MP(Mp);
					break;
				}
			}
		}
		double Vbf_x = break_force*-V_co_normal_x;
		double Vbf_y = break_force*-V_co_normal_y;
		V->set_adhesion_break(Vbf_x, Vbf_y);
	}
}

void dbReader::point_response_b1(point3D *V)
{
	/*************************************

           o-<-o     o-<-o
          /     \   /     \
		 o    V1 ooo V2    o
		/       // \\       \
	   o     p2oo V oop3     o
		\ C2  //p1 p4\\   C3/
		 o->-oo       oo->-o
		      |       |
			  o   C1  o
			   \     /
			    o->-o

	*************************************/
	if (V->get_in_cell())
	{
		edge *V_ab = V->getab();
		edge *V_ba = V->getba();
		if (V_ab->get_attach())
		{
			edge_response_b0(V_ab);
		}
		if (V_ba->get_attach())
		{
			edge *V_ba_Nei = V_ba->get_Neighbor();
			edge_response_b0(V_ba_Nei);
		}
	}
}

void dbReader::point_response_b2(point3D *V)
{
	/***********************************
	             o-<-o
	          R /
	       -<-oo E1
		       \\ V E2 F
			 E3 \o-->-o->-o
			 V1  oo-<-o-<-o
				 || V2 E4
			 o->-oo->-o
	***********************************/
	if (V->get_in_cell() && V->getba()->get_attach() && V->getab()->get_attach() && 
		V->getba()->get_Pair()!=V->getab()->get_Pair())
	{
		point3D *R = V->getrp();
		point3D *F = V->getfp();
		edge *E1 = V->getba();
		edge *E2 = V->getab();
		edge *E3 = E1->get_Neighbor();
		edge *E4 = E2->get_Neighbor();
		point3D *V1 = E3->p1();
		point3D *V2 = E4->p2();
		//cout<<"point breaks at point "<<V->id()<<" that cell "<<E1->Cell()->id()<<" breaks away from cell "<<E3->Cell()->id()<<" and cell "<<E4->Cell()->id()<<endl;
		interpair *Pl1 = E1->get_Pair();
		interpair *Pl2 = E2->get_Pair();
		if (Pl1->MP()->at(0)->T1()==V)
		{
			Pl1->MP()->at(0)->setT1(R);
			Pl1->MP()->at(0)->setS2(V1->getfp());
			if (Pl1->MP()->at(0)->S1()==R)
			{
				Pl1->MP()->at(0)->set_F_t1(Pl1->MP()->at(0)->get_F_s1());
				Pl1->MP()->at(0)->set_F_s2(Pl1->MP()->at(0)->get_F_t2());
			}
			else
			{
				Pl1->MP()->at(0)->set_F_t1(R->get_stored_apart_force());
				Pl1->MP()->at(0)->set_F_s2(V1->getfp()->get_stored_apart_force());
			}
		}
		else if (Pl1->MP()->at(0)->T2()==V)
		{
			Pl1->MP()->at(0)->setT2(R);
			Pl1->MP()->at(0)->setS1(V1->getfp());
			if (Pl1->MP()->at(0)->S2()==R)
			{
				Pl1->MP()->at(0)->set_F_t2(Pl1->MP()->at(0)->get_F_s2());
				Pl1->MP()->at(0)->set_F_s1(Pl1->MP()->at(0)->get_F_t1());
			}
			else
			{
				Pl1->MP()->at(0)->set_F_t2(R->get_stored_apart_force());
				Pl1->MP()->at(0)->set_F_s1(V1->getfp()->get_stored_apart_force());
			}
		}
		E1->set_attach(0);E1->set_Neighbor(NULL);
		E3->set_attach(0);E3->set_Neighbor(NULL);
		if (Pl2->MP()->size()==0) {cout<<"    -> error: V-F pair was empty: "<<Pl2->get_I1()->id()<<","<<Pl2->get_I2()->id()<<endl;}
		if (Pl2->MP()->at(0)->S1()==V)
		{
			Pl2->MP()->at(0)->setS1(F);
			Pl2->MP()->at(0)->setT2(V2->getrp());
			if (Pl2->MP()->at(0)->T1()==F)
			{
				Pl2->MP()->at(0)->set_F_s1(Pl2->MP()->at(0)->get_F_t1());
				Pl2->MP()->at(0)->set_F_t2(Pl2->MP()->at(0)->get_F_s2());
			}
			else
			{
				Pl2->MP()->at(0)->set_F_s1(F->get_stored_apart_force());
				Pl2->MP()->at(0)->set_F_t2(V2->getrp()->get_stored_apart_force());
			}
		}
		else if (Pl2->MP()->at(0)->S2()==V)
		{
			Pl2->MP()->at(0)->setS2(F);
			Pl2->MP()->at(0)->setT1(V2->getrp());
			if (Pl2->MP()->at(0)->T2()==F)
			{
				Pl2->MP()->at(0)->set_F_s2(Pl2->MP()->at(0)->get_F_t2());
				Pl2->MP()->at(0)->set_F_t1(Pl2->MP()->at(0)->get_F_s1());
			}
			else
			{
				Pl2->MP()->at(0)->set_F_s2(F->get_stored_apart_force());
				Pl2->MP()->at(0)->set_F_t1(V2->getrp()->get_stored_apart_force());
			}
		}
		E2->set_attach(0);E2->set_Neighbor(NULL);
		E4->set_attach(0);E4->set_Neighbor(NULL);
		if (!V1->getba()->get_attach() && !V2->getab()->get_attach())
		{
			cell *C1 = V1->getba()->Cell();
			cell *C2 = V2->getab()->Cell();
			int pn = (int)C1->get_pairs()->size();
			for (int i=0;i<pn;i++)
			{
				if (C1->get_pairs()->at(i)->get_I1()==C1 &&
					C1->get_pairs()->at(i)->get_I2()==C2)
				{
					mergepair *Mp = C1->get_pairs()->at(i)->MP()->at(0);
					point3D *sC = Mp->S1();
					point3D *nC = Mp->S2();
					sC->removemp(nC);
					nC->removemp(sC);
					node_remove_0(sC);
					C1->get_pairs()->at(i)->remove_MP(Mp);
					break;
				}
				else if (C1->get_pairs()->at(i)->get_I2()==C1 &&
						 C1->get_pairs()->at(i)->get_I1()==C2)
				{
					mergepair *Mp = C1->get_pairs()->at(i)->MP()->at(0);
					point3D *sC = Mp->S2();
					point3D *nC = Mp->S1();
					sC->removemp(nC);
					nC->removemp(sC);
					node_remove_0(sC);
					C1->get_pairs()->at(i)->remove_MP(Mp);
					break;
				}
			}
		}
		//////// add the break force ///////
		double break_force = Pl1->get_adhesion()/20;
		double nbax = E1->getnormal(0);
		double nbay = E1->getnormal(1);
		double nabx = E2->getnormal(0);
		double naby = E2->getnormal(1);
		double nx = (nbax + nabx)*0.5;
		double ny = (naby + naby)*0.5;
		double nd = sqrt(nx*nx + ny*ny);
		double bf_x1 = -break_force*nx/nd + V->get_adhesion_break(0);
		double bf_y1 = -break_force*ny/nd + V->get_adhesion_break(1);
		double bf_x2 = break_force*nx/nd/2 + V1->get_adhesion_break(0);
		double bf_y2 = break_force*ny/nd/2 + V1->get_adhesion_break(1);
		double bf_x3 = break_force*nx/nd/2 + V2->get_adhesion_break(0);
		double bf_y3 = break_force*ny/nd/2 + V2->get_adhesion_break(1);
		double ld = 0.1;
		double V_dx = ld*-break_force*nx/nd + V->x();
		double V_dy = ld*-break_force*ny/nd + V->y();
		V->setx(V_dx);
		V->sety(V_dy);
		V->set_adhesion_break(bf_x1, bf_y1);
		V1->set_adhesion_break(bf_x2, bf_y2);
		V2->set_adhesion_break(bf_x3, bf_y3);
		V->removemp(V1);V->removemp(V2);
		V1->removemp(V);V2->removemp(V);
		node_remove_0(V);
		if (!R->getab()->get_attach() && !R->getba()->get_attach())
		{
			//cout<<"    -> additional ";
			point_response_b0(R);
		}
		if (!F->getab()->get_attach() && !F->getba()->get_attach())
		{
			//cout<<"    -> additional ";
			point_response_b0(F);
		}
	}
	else if (V->get_in_cell() && !V->getba()->get_attach() && V->getab()->get_attach())
	{
		//cout<<"    -> additional ";
		edge_response_b0(V->getab());
	}
	else if (V->get_in_cell() && V->getba()->get_attach() && !V->getab()->get_attach())
	{
		//cout<<"    -> additional ";
		edge_response_b0(V->getba()->get_Neighbor());
	}
}

void dbReader::cell_friction_response(cell *C, int time)
{
	/********************************
	 o       oo          o
	 |       ||          |
	 o->-o---oo--o->-o---o
	 o-<-o---o---o---o-<-o---o---o
								  \ t
								   o ==>
				C				   |     (dx,dy)
								   o ==>
								  / s
	 o---o->-o---o->-o->-o---o---o
	     o-<-o   o-<-o-<-o
		 |   |   |       |
		 o   o   o       o
	********************************/
	double x_angle = 110;
	C->center_refresh();
	edge *E = C->get_sides()->at(0);
	double xc = C->get_center(0);
	double yc = C->get_center(1);
	double dx = 0;
	double dy = 1;
	double angle = C->get_migrate_angle();
	double counter_angle = angle + 180;
	if (counter_angle>360) {counter_angle -= 360;}
	if (counter_angle==90 || counter_angle==270)
	{
		if (counter_angle==270) {dy = -1;}
		bool flagE = false;
		int sn = (int)C->get_sides()->size();
		for (int i=0;i<sn;i++)
		{
			double x1 = C->get_sides()->at(i)->p1()->x();
			double x2 = C->get_sides()->at(i)->p2()->x();
			if (x1<=xc && x2>xc)
			{
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = dy*ny;
				if (plot_dn>0)
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	else
	{
		dx = cos(counter_angle*PI/180);
		dy = sin(counter_angle*PI/180);
		double slope = tan(angle*PI/180);
		double intercept = yc - slope*xc;
		int sn = (int)C->get_sides()->size();
		bool flagE = false;
		for (int i=0;i<sn;i++)
		{
			point3D *p1 = C->get_sides()->at(i)->p1();
			point3D *p2 = C->get_sides()->at(i)->p2();
			double x1 = p1->x();
			double y1 = p1->y();
			double x2 = p2->x();
			double y2 = p2->y();
			double t1 = y1 - slope*x1 - intercept;
			double t2 = y2 - slope*x2 - intercept;
			if ((t1>=0 && t2<0) || (t1<0 && t2>=0))
			{
				double nx = C->get_sides()->at(i)->getnormal(0);
				double ny = C->get_sides()->at(i)->getnormal(1);
				double plot_dn = nx*dx + ny*dy;
				if (plot_dn>0) 
				{
					E = C->get_sides()->at(i);
					flagE = true;
					break;
				}
			}
			if (flagE) {break;}
		}
	}
	point3D *et = E->p1();
	point3D *es = E->p2();
	point3D *s = C->get_migrate_p(0);
	point3D *t = C->get_migrate_p(1);
	if (s==NULL || t==NULL) {cout<<"    -> error: s t are NULLs at migrating cell "<<C->id()<<endl;}
	bool front_b = true;
	vector<interpair*> front_single_pairs;
	int s_t_step = 0;
	if (s!=t)
	{
		point3D *ppp = s;
		while (ppp!=t)
		{
			if (ppp->getab()->get_attach()) {front_b = false;break;}
			if (ppp!=s && !ppp->getab()->get_attach() && !ppp->getba()->get_attach())
			{
				int pppmn = (int)ppp->getmp()->size();
				if (pppmn>0)
				{
					int pn = (int)C->get_pairs()->size();
					for (int i=0;i<pppmn;i++)
					{
						cell *Ct2 = ppp->getmp()->at(i)->getab()->Cell();
						for (int j=0;j<pn;j++)
						{
							if ((C->get_pairs()->at(j)->get_I1()==C &&
								 C->get_pairs()->at(j)->get_I2()==Ct2) ||
								(C->get_pairs()->at(j)->get_I2()==C &&
								 C->get_pairs()->at(j)->get_I1()==Ct2))
							{
								front_single_pairs.push_back(C->get_pairs()->at(j));
								break;
							}
						}
					}
				}
			}
			s_t_step++;
			ppp = ppp->getfp();
		}
	}
	////////////////
	point3D *ps = s;
	point3D *pt = t;
	double angle_s = 0;
	double angle_t = 0;
	vector<interpair*> Ips;
	vector<interpair*> Ipt;
	vector<int> Ips_type;
	vector<int> Ipt_type;
	vector<double> Ips_angle;
	vector<double> Ipt_angle;
	vector<int> Ips_final; // 0:select, 1:not select
	vector<int> Ipt_final; // 0:select, 1:not select
	vector<interpair*> Ips_single;
	vector<interpair*> Ipt_single;
	/*********************************
	 pair type:
	 ================================
	 type: 0
	     s       t
	 o---o->-o->-o---o      | 0   0 |
	     o-<-o-<-o          o->-o->-o
		 | 0   0 |      o-<-o-<-o-<-o-<-o
		                    t       s
	 type: 1
	     s       t
	 o---o->-o->-o---o     || 1   0 |
	 o--oo-<-o-<-o      o--oo->-o->-o
		|| 1   0 |      o---o-<-o-<-o---o
							t       s
	 type: 2
	     s       t
	 o---o->-o->-o---o      | 0   1 ||
	     o-<-o-<-oo--o      o->-o->-oo--o
		 | 0   1 ||     o---o-<-o-<-o---o
							t       s
	 type: 3
		 s       t
	 o---o->-o->-o---o     || 1   1 ||
	 o--oo-<-o-<-oo--o  o--oo->-o->-oo--o
		|| 1   1 ||     o---o-<-o-<-o---o
							t       s
	*********************************/
	if (!s->getba()->get_attach() || !s->getab()->get_attach())
	{
		s = s->getfp();
		while (angle_s<x_angle || s!=es)
		{
			if (s==es) break;
			s = s->getrp();
			edge *Es = s->getba();
			if (Es->get_attach())
			{
				interpair *PEs = Es->get_Pair();
				cell *C11 = PEs->get_I1();
				cell *C12 = PEs->get_I2();
				if (PEs->MP()->size()==0)
				{
					cout<<"      -> error at cell "<<C->id()<<" at point "<<s->id()<<":s_pair at edge "<<Es->id()<<":"<<C11->id()<<","<<C12->id()<<" was emtpy"<<endl;
				}
				point3D *PEs_s = NULL;
				point3D *PEs_t = NULL;
				if (C11==C)
				{
					PEs_s = PEs->MP()->at(0)->S1();
					PEs_t = PEs->MP()->at(0)->T1();
				}
				else
				{
					PEs_s = PEs->MP()->at(0)->S2();
					PEs_t = PEs->MP()->at(0)->T2();
				}
				s = PEs_s->getfp();
				double sabx = s->getab()->getnormal(0);
				double saby = s->getab()->getnormal(1);
				double sbax = s->getba()->getnormal(0);
				double sbay = s->getba()->getnormal(1);
				double sabd = sqrt(sabx*sabx + saby*saby);
				double sbad = sqrt(sbax*sbax + sbay*sbay);
				double snx = sabx/sabd + sbax/sbad;
				double sny = saby/sabd + sbay/sbad;
				double angle_sn = vector2angle(snx,sny);
				angle_s = angle - angle_sn;
				if (angle_s<0) {angle_s += 360;}
				if (angle_s>x_angle) {break;}
				int pair_type = -1;
				int s_mn = (int)PEs_s->getmp()->size();
				int t_mn = (int)PEs_t->getmp()->size();
				if      (t_mn==1 && s_mn==1) {pair_type = 0;}
				else if (t_mn==1 && s_mn>=2) {pair_type = 1;}
				else if (t_mn>=2 && s_mn==1) {pair_type = 2;}
				else if (t_mn>=2 && s_mn>=2) {pair_type = 3;}
				if (!E->get_attach())
				{
					Ips.push_back(PEs);
					Ips_angle.push_back(angle_s);
					Ips_type.push_back(pair_type);
					Ips_final.push_back(0);
				}
				else
				{
					if (PEs!=E->get_Pair())
					{
						Ips.push_back(PEs);
						Ips_angle.push_back(angle_s);
						Ips_type.push_back(pair_type);
						Ips_final.push_back(0);
					}
					else
					{
						angle_s = 180;
						s = es;
						break;
					}
				}
			}
			else
			{
				if ((int)s->getmp()->size()==1 &&
				        !s->getab()->get_attach() &&
				        !s->getba()->get_attach())
				{
					cell *Cns = s->getmp()->at(0)->getab()->Cell();
					int cpn = (int)C->get_pairs()->size();
					for (int j=0;j<cpn;j++)
					{
						if ((C->get_pairs()->at(j)->get_I1()==C &&
							 C->get_pairs()->at(j)->get_I2()==Cns) ||
							(C->get_pairs()->at(j)->get_I2()==C &&
							 C->get_pairs()->at(j)->get_I1()==Cns))
						{
							Ips_single.push_back(C->get_pairs()->at(j));
							break;
						}
					}
				}
				double sabx = s->getab()->getnormal(0);
				double saby = s->getab()->getnormal(1);
				double sbax = s->getba()->getnormal(0);
				double sbay = s->getba()->getnormal(1);
				double sabd = sqrt(sabx*sabx + saby*saby);
				double sbad = sqrt(sbax*sbax + sbay*sbay);
				double snx = sabx/sabd + sbax/sbad;
				double sny = saby/sabd + sbay/sbad;
				double angle_sn = 0;
				if (s==ps) {angle_sn = vector2angle(sbax,sbay);}
				else       {angle_sn = vector2angle(snx,sny);}
				angle_s = angle - angle_sn;
				if (angle_s<0) {angle_s += 360;}
				if (angle_s>x_angle) {break;}
			}
		}
	}
	if (!t->getab()->get_attach() || !t->getba()->get_attach())
	{
		t = t->getrp();
		while (angle_t<96 || t!=et)
		{
			if (t==et) break;
			t = t->getfp();
			edge *Et = t->getab();
			if (Et->get_attach())
			{
				interpair *PEt = Et->get_Pair();
				cell *C11 = PEt->get_I1();
				cell *C12 = PEt->get_I2();
				if (PEt->MP()->size()==0)
				{
					cout<<"      -> error at cell "<<C->id()<<" at point "<<t->id()<<":t_pair at edge "<<Et->id()<<":"<<C11->id()<<","<<C12->id()<<" was emtpy"<<endl;
				}
				point3D *PEt_s = NULL;
				point3D *PEt_t = NULL;
				if (C11==C)
				{
					PEt_s = PEt->MP()->at(0)->S1();
					PEt_t = PEt->MP()->at(0)->T1();
				}
				else
				{
					PEt_s = PEt->MP()->at(0)->S2();
					PEt_t = PEt->MP()->at(0)->T2();
				}
				t = PEt_t->getrp();
				double tabx = t->getab()->getnormal(0);
				double taby = t->getab()->getnormal(1);
				double tbax = t->getba()->getnormal(0);
				double tbay = t->getba()->getnormal(1);
				double tabd = sqrt(tabx*tabx + taby*taby);
				double tbad = sqrt(tbax*tbax + tbay*tbay);
				double tnx = tabx/tabd + tbax/tbad;
				double tny = taby/tabd + tbay/tbad;
				double angle_tn = vector2angle(tnx,tny);
				angle_t = angle_tn - angle;
				if (angle_t<0) {angle_t += 360;}
				if (angle_t>x_angle) {break;}
				int pair_type = -1;
				int s_mn = (int)PEt_s->getmp()->size();
				int t_mn = (int)PEt_t->getmp()->size();
				if      (s_mn==1 && t_mn==1) {pair_type = 0;}
				else if (s_mn==1 && t_mn>=2) {pair_type = 1;}
				else if (s_mn>=2 && t_mn==1) {pair_type = 2;}
				else if (s_mn>=2 && t_mn>=2) {pair_type = 3;}
				if (!E->get_attach())
				{
					Ipt.push_back(PEt);
					Ipt_angle.push_back(angle_t);
					Ipt_type.push_back(pair_type);
					Ipt_final.push_back(0);
				}
				else
				{
					if (E->get_Pair()!=PEt)
					{
						Ipt.push_back(PEt);
						Ipt_angle.push_back(angle_t);
						Ipt_type.push_back(pair_type);
						Ipt_final.push_back(0);
					}
					else
					{
						angle_t = 180;
						t = et;
						break;
					}
				}
			}
			else
			{
				if ((int)t->getmp()->size()==1 &&
					    !t->getab()->get_attach() &&
						!t->getba()->get_attach())
				{
					cell *Cnt = t->getmp()->at(0)->getab()->Cell();
					int cpn = (int)C->get_pairs()->size();
					for (int j=0;j<cpn;j++)
					{
						if ((C->get_pairs()->at(j)->get_I1()==C &&
							 C->get_pairs()->at(j)->get_I2()==Cnt) ||
							(C->get_pairs()->at(j)->get_I2()==C &&
							 C->get_pairs()->at(j)->get_I1()==Cnt))
						{
							Ipt_single.push_back(C->get_pairs()->at(j));
							break;
						}
					}
				}
				double tabx = t->getab()->getnormal(0);
				double taby = t->getab()->getnormal(1);
				double tbax = t->getba()->getnormal(0);
				double tbay = t->getba()->getnormal(1);
				double tabd = sqrt(tabx*tabx + taby*taby);
				double tbad = sqrt(tbax*tbax + tbay*tbay);
				double tnx = tabx/tabd + tbax/tbad;
				double tny = taby/tabd + tbay/tbad;
				double angle_tn = 0;
				if (t==pt) {angle_tn = vector2angle(tabx,taby);}
				else       {angle_tn = vector2angle(tnx,tny);}
				angle_t = angle_tn - angle;
				if (angle_t<0) {angle_t += 360;}
				if (angle_t>x_angle) {break;}
			}
		}
	}
	int Ipsn = (int)Ips.size();
	int Iptn = (int)Ipt.size();
	dx = cos(angle*PI/180);
	dy = sin(angle*PI/180);
	if      (Ipsn==0 && Iptn==0)
	{
		point3D *s_end_p = NULL;
		point3D *t_end_p = NULL;
		if (t->getrp()==s->getfp())
		{
			s_end_p = t;
			t_end_p = s;
		}
		else if (t->getrp()==s)
		{
			s_end_p = t;
			t_end_p = s;
		}
		else
		{
			s_end_p = t->getrp();
			t_end_p = s->getfp();
		}
		bool run_f = true;
		while (s_end_p!=t_end_p)
		{
			if (s_end_p->getab()->get_attach()) {run_f = false; break;}
			s_end_p = s_end_p->getfp();
		}
		if (run_f)
		{
			int Ipsn_s = (int)Ips_single.size();
			int Iptn_s = (int)Ipt_single.size();
			int Ipfs_s = (int)front_single_pairs.size();
			if (Ipfs_s>0 && front_b)
			{
				for (int j=0;j<Ipfs_s;j++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (front_single_pairs[j]->get_I1()==C)
					{
						mergepair *Mp = front_single_pairs[j]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[j]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = front_single_pairs[j]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[j]->remove_MP(Mp);
					}
				}
			}
			if ((Ipsn_s>0 || Iptn_s>0) && front_b)
			{
				for (int j=0;j<Ipsn_s;j++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ips_single[j]->get_I1()==C)
					{
						mergepair *Mp = Ips_single[j]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[j]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ips_single[j]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[j]->remove_MP(Mp);
					}
				}
				for (int j=0;j<Iptn_s;j++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ipt_single[j]->get_I1()==C)
					{
						mergepair *Mp = Ipt_single[j]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[j]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ipt_single[j]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[j]->remove_MP(Mp);
					}
				}
				//////// position refresh /////////
				point3D *pmp = ps->getfp();
				double dmpx = BD*dx;
				double dmpy = BD*dy;
				while (pmp!=ps)
				{
					double px = pmp->x();
					double py = pmp->y();
					pmp->setx(px + dmpx);
					pmp->sety(py + dmpy);
					pmp = pmp->getfp();
				}
				double pxt = ps->x();
				double pyt = ps->y();
				ps->setx(pxt + dmpx);
				ps->sety(pyt + dmpy);
			}
			else if ((Ipsn_s==0 || Iptn_s==0) && front_b)
			{
				point3D *pmp = ps->getfp();
				double dmpx = BD*dx;
				double dmpy = BD*dy;
				while (pmp!=ps)
				{
					double px = pmp->x();
					double py = pmp->y();
					pmp->setx(px + dmpx);
					pmp->sety(py + dmpy);
					pmp = pmp->getfp();
				}
				double pxt = ps->x();
				double pyt = ps->y();
				ps->setx(pxt + dmpx);
				ps->sety(pyt + dmpy);
			}
		}
	}
	else if (Ipsn==0 && Iptn>0)
	{
		point3D *t_end_p = NULL;
		if (Ipt[Iptn-1]->get_I1()==C) {t_end_p = Ipt[Iptn-1]->MP()->at(0)->T1();}
		else                          {t_end_p = Ipt[Iptn-1]->MP()->at(0)->T2();}
		bool run_f = true;
		if (!t_end_p->getab()->get_attach())
		{
			point3D *pmp = t_end_p->getfp();
			while (pmp!=ps)
			{
				if (pmp->getba()->get_attach()) {run_f = false;break;}
				pmp = pmp->getfp();
			}
		}
		else {run_f = false;}
		if (run_f)
		{
			double Ipt_Ff_max = 0;
			double Ipt_ab_max = 0;
			for (int i=0;i<Iptn;i++)
			{
				point3D *s1 = Ipt[i]->MP()->at(0)->S1();
				point3D *t1 = Ipt[i]->MP()->at(0)->T1();
				point3D *s2 = Ipt[i]->MP()->at(0)->S2();
				point3D *t2 = Ipt[i]->MP()->at(0)->T2();
				double Ff_s = Ipt[i]->MP()->at(0)->get_Fr_s1();
				double Ff_t = Ipt[i]->MP()->at(0)->get_Fr_t1();
				int Ff_d = Ipt[i]->MP()->at(0)->get_Fr_dir(); // -1:null, 1:s->t, 0:s<-t
				double Ipt_ab = Ipt[i]->get_adhesion();
				double e_force_x_s1 = s1->get_elastic_force(0);
				double e_force_y_s1 = s1->get_elastic_force(1);
				double e_force_x_t1 = t1->get_elastic_force(0);
				double e_force_y_t1 = t1->get_elastic_force(1);
				double e_force_x_s2 = s2->get_elastic_force(0);
				double e_force_y_s2 = s2->get_elastic_force(1);
				double e_force_x_t2 = t2->get_elastic_force(0);
				double e_force_y_t2 = t2->get_elastic_force(1);
				double s1ldx = 0;double s1ldy = 0;double t1ldx = 0;double t1ldy = 0;
				double s2ldx = 0;double s2ldy = 0;double t2ldx = 0;double t2ldy = 0;
				if (Ipt[i]->get_I1()==C)
				{
					double s1lx = s1->getfp()->x() - s1->x();
					double s1ly = s1->getfp()->y() - s1->y();
					double s1ld = sqrt(s1lx*s1lx + s1ly*s1ly);
					s1ldx = s1lx/s1ld;
					s1ldy = s1ly/s1ld;
					t2ldx = s1ldx;
					t2ldy = s1ldy;
					double t1lx = t1->getfp()->x() - t1->x();
					double t1ly = t1->getfp()->y() - t1->y();
					double t1ld = sqrt(t1lx*t1lx + t1ly*t1ly);
					t1ldx = t1lx/t1ld;
					t1ldy = t1ly/t1ld;
					s2ldx = t1ldx;
					s2ldy = t1ldy;
				}
				else
				{
					double s2lx = s2->getfp()->x() - s2->x();
					double s2ly = s2->getfp()->y() - s2->y();
					double s2ld = sqrt(s2lx*s2lx + s2ly*s2ly);
					s2ldx = s2lx/s2ld;
					s2ldy = s2ly/s2ld;
					t1ldx = s2ldx;
					t1ldy = s2ldy;
					double t2lx = t2->getfp()->x() - t2->x();
					double t2ly = t2->getfp()->y() - t2->y();
					double t2ld = sqrt(t2lx*t2lx + t2ly*t2ly);
					t2ldx = t2lx/t2ld;
					t2ldy = t2ly/t2ld;
					s1ldx = t2ldx;
					s1ldy = t2ldy;
				}
				double Ff_s1_inc = abs(e_force_x_s1*s1ldx + e_force_y_s1*s1ldy);
				double Ff_t2_inc = abs(e_force_x_t2*t2ldx + e_force_y_t2*t2ldy);
				double Ff_s2_inc = abs(e_force_x_s2*s2ldx + e_force_y_s2*s2ldy);
				double Ff_t1_inc = abs(e_force_x_t1*t1ldx + e_force_y_t1*t1ldy);
				double Ff_s1_t2_inc = 0;
				if (Ff_s1_inc>Ff_t2_inc) {Ff_s1_t2_inc = Ff_s1_inc;}
				else                     {Ff_s1_t2_inc = Ff_t2_inc;}
				double Ff_s2_t1_inc = 0;
				if (Ff_s2_inc>Ff_t1_inc) {Ff_s2_t1_inc = Ff_s2_inc;}
				else                     {Ff_s2_t1_inc = Ff_t1_inc;}
				if (Ipt_ab_max<Ipt_ab) {Ipt_ab_max = Ipt_ab;}
				if (Ipt[i]->get_I1()==C)
				{
					if      (Ff_d==-1) // null
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(1);
					}
					else if (Ff_d==1)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					else if (Ff_d==0)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(1);
					}
					double C_Fr_s1 = Ipt[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ipt[i]->MP()->at(0)->get_Fr_t1();
					if (Ipt_Ff_max<C_Fr_s1) {Ipt_Ff_max = C_Fr_s1;}
					if (Ipt_Ff_max<C_Fr_t1) {Ipt_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = s1->getfp();
					while (ptC!=t1)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getfp()->x() - ptC->x();
						double ply = ptC->getfp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==1)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Ipt_Ff_max<Fr_tC) {Ipt_Ff_max = Fr_tC;}
						if (Ipt_Ff_max<Fr_mC) {Ipt_Ff_max = Fr_mC;}
						ptC = ptC->getfp();
					}
				}
				else
				{
					if      (Ff_d==-1)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==1)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==0)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					double C_Fr_s1 = Ipt[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ipt[i]->MP()->at(0)->get_Fr_t1();
					if (Ipt_Ff_max<C_Fr_s1) {Ipt_Ff_max = C_Fr_s1;}
					if (Ipt_Ff_max<C_Fr_t1) {Ipt_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = s2->getfp();
					while (ptC!=t2)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getfp()->x() - ptC->x();
						double ply = ptC->getfp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==0)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Ipt_Ff_max<Fr_tC) {Ipt_Ff_max = Fr_tC;}
						if (Ipt_Ff_max<Fr_mC) {Ipt_Ff_max = Fr_mC;}
						ptC = ptC->getfp();
					}
				}
			}
			/////// boundary edges sliding ///////
			if (Ipt_Ff_max>Ipt_ab_max && front_b)
			{
				point3D *t_s1 = NULL;
				point3D *t_t1 = NULL;
				vector<point3D*> Ipt_Pt;
				vector<point3D*> Ipt_NePt;
				vector<interpair*> Ipt_Ip;
				vector<edge*> Ipt_NeEg;
				vector<node*> Ipt_Nd;
				if (Ipt[0]->get_I1()==C) {t_s1 = Ipt[0]->MP()->at(0)->S1();}
				else                     {t_s1 = Ipt[0]->MP()->at(0)->S2();}
				if (Ipt[Iptn-1]->get_I1()==C) {t_t1 = Ipt[Iptn-1]->MP()->at(0)->T1();}
				else                          {t_t1 = Ipt[Iptn-1]->MP()->at(0)->T2();}
				double m_l = 10000;
				double t_t1_ab_nx = t_t1->getab()->getnormal(0);
				double t_t1_ab_ny = t_t1->getab()->getnormal(1);
				double t_t1_angle = vector2angle(t_t1_ab_nx,t_t1_ab_ny) - angle;
				if (t_t1_angle<0) {t_t1_angle += 360;}
				//////////// heat up ////////////
				point3D *t_tm = t_s1;
				int tspn = 0;
				while (t_tm!=t_t1->getfp())
				{
					Ipt_Pt.push_back(t_tm);
					node* tNd = t_tm->get_node_id();
					Ipt_Nd.push_back(tNd);
					tNd->removevertex(t_tm);
					if ((int)t_tm->getmp()->size()==0)
					{
						Ipt_NePt.push_back(NULL);
					}
					else
					{
						Ipt_NePt.push_back(t_tm->getmp()->at(0));
						for (int i=0;i<(int)t_tm->getmp()->size();i++)
						{
							t_tm->getmp()->at(i)->removemp(t_tm);
						}
						t_tm->clearmp(); // clear the neighbors
					}
					if (t_tm->getab()->get_attach())
					{
						double neegl = t_tm->getab()->getlength();if (neegl<m_l) {m_l = neegl;}
						Ipt_NeEg.push_back(t_tm->getab()->get_Neighbor());
						Ipt_Ip.push_back(t_tm->getab()->get_Pair());
					}
					else
					{
						Ipt_NeEg.push_back(NULL);
						Ipt_Ip.push_back(NULL);
					}
					t_tm = t_tm->getfp();
					tspn++;
				}
				//////////// single point pair removal /////////////
				int Iptn_s = (int)Ipt_single.size();
				for (int i=0;i<Iptn_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ipt_single[i]->get_I1()==C)
					{
						mergepair *Mp = Ipt_single[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ipt_single[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[i]->remove_MP(Mp);
					}
				}
				int Ipsn_s = (int)Ips_single.size();
				for (int i=0;i<Ipsn_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ips_single[i]->get_I1()==C)
					{
						mergepair *Mp = Ips_single[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ips_single[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[i]->remove_MP(Mp);
					}
				}
				int Ipfs_s = (int)front_single_pairs.size();
				for (int i=0;i<Ipfs_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (front_single_pairs[i]->get_I1()==C)
					{
						mergepair *Mp = front_single_pairs[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = front_single_pairs[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[i]->remove_MP(Mp);
					}
				}
				//////////// position refresh ////////////
				for (int i=tspn-1;i>=0;i--)
				{
					if (i==0)
					{
						int nd_size = (int)nodeList.size();
						node *A = new node(nd_size);
						nodeList.push_back(A);
						nodeList[nd_size]->pushvertex(Ipt_Pt[i]);
						nodeList[nd_size]->set_in_global(1);
						Ipt_Pt[i]->set_node_id(nodeList[nd_size]);
						Ipt_Pt[i]->getab()->set_attach(0);
					}
					else
					{
						if (tspn>0)
						{
							/////////////// node //////////////
							Ipt_Nd[i-1]->pushvertex(Ipt_Pt[i]);
							Ipt_Pt[i]->set_node_id(Ipt_Nd[i-1]);
							double rpx = 0;
							double rpy = 0;
							if (Ipt_NePt[i-1]!=NULL)
							{
								rpx = Ipt_NePt[i-1]->x();
								rpy = Ipt_NePt[i-1]->y();
							}
							else
							{
								rpx = Ipt_Pt[i-1]->x();
								rpy = Ipt_Pt[i-1]->y();
							}
							Ipt_Pt[i]->setx(rpx);
							Ipt_Pt[i]->sety(rpy);
							if (Ipt_NePt[i-1]!=NULL)
							{
								vector<point3D*> rpmp;
								rpmp.push_back(Ipt_NePt[i-1]);
								Ipt_Pt[i]->pushmp(Ipt_NePt[i-1]);
								for (int j=0;j<(int)Ipt_NePt[i-1]->getmp()->size();j++)
								{
									rpmp.push_back(Ipt_NePt[i-1]->getmp()->at(j));
									Ipt_Pt[i]->pushmp(Ipt_NePt[i-1]->getmp()->at(j));
								}
								for (int j=0;j<(int)rpmp.size();j++)
								{
									rpmp[j]->pushmp(Ipt_Pt[i]);
									edge *rpmp_ab = rpmp[j]->getab();
									edge *rpmp_ba = rpmp[j]->getba();
									cell *Crp = rpmp_ab->Cell();
									cell *C1r = NULL;
									cell *C2r = NULL;
									if (rpmp_ab->get_attach())
									{
										if (rpmp_ab->get_Pair()->get_I1()==Crp) {C1r = rpmp_ab->get_Pair()->get_I2();}
										else {C1r = rpmp_ab->get_Pair()->get_I1();}
									}
									if (rpmp_ba->get_attach())
									{
										if (rpmp_ba->get_Pair()->get_I1()==Crp) {C2r = rpmp_ba->get_Pair()->get_I2();}
										else {C2r = rpmp_ba->get_Pair()->get_I1();}
									}
									if (C1r!=C && C2r!=C)
									{
										int crpn = (int)Crp->get_pairs()->size();
										for (int k=0;k<crpn;k++)
										{
											if ((Crp->get_pairs()->at(k)->get_I1()==Crp &&
												Crp->get_pairs()->at(k)->get_I2()==C) ||
												(Crp->get_pairs()->at(k)->get_I2()==Crp &&
												Crp->get_pairs()->at(k)->get_I1()==C))
											{
												if ((int)Crp->get_pairs()->at(k)->MP()->size()>0)
												{
													if      (Crp->get_pairs()->at(k)->MP()->at(0)->S1()==rpmp[j])
													{
														Crp->get_pairs()->at(k)->MP()->at(0)->setS2(Ipt_Pt[i]);
														Crp->get_pairs()->at(k)->MP()->at(0)->setT2(Ipt_Pt[i]);
													}
													else if (Crp->get_pairs()->at(k)->MP()->at(0)->S2()==rpmp[j])
													{
														Crp->get_pairs()->at(k)->MP()->at(0)->setS1(Ipt_Pt[i]);
														Crp->get_pairs()->at(k)->MP()->at(0)->setT1(Ipt_Pt[i]);
													}
												}
												break;
											}
										}
									}
								}
							}
							if (i<tspn-1)
							{
								//////////// edge ////////////
								if (Ipt_NeEg[i-1]==NULL)
								{
									Ipt_Pt[i]->getab()->set_attach(0);
								}
								else
								{
									Ipt_Pt[i]->getab()->set_attach(1);
									Ipt_NeEg[i-1]->set_attach(1);
									Ipt_Pt[i]->getab()->set_Neighbor(Ipt_NeEg[i-1]);
									Ipt_NeEg[i-1]->set_Neighbor(Ipt_Pt[i]->getab());
								}
								//////////// pair ////////////
								Ipt_Pt[i]->getab()->set_Pair(Ipt_Ip[i-1]);
							}
							else
							{
								//////////// edge ////////////
								Ipt_Pt[i]->getab()->set_attach(0);
								//////////// pair ////////////
								Ipt_Pt[i]->getab()->set_Pair(NULL);
							}
						}
					}
				}
				/////// points other than points from t_s->t_t ///////
				double dpmx = dx*m_l;
				double dpmy = dy*m_l;
				if      (t_t1_angle>96)
				{
					vector<point3D*> tmpt;
					tmpt.push_back(Ipt_NePt[tspn-1]);
					for (int j=0;j<(int)Ipt_NePt[tspn-1]->getmp()->size();j++)
					{tmpt.push_back(Ipt_NePt[tspn-1]->getmp()->at(j));}
					for (int j=0;j<(int)tmpt.size();j++)
					{
						point3D *mppt = tmpt[j];
						edge *mpab = mppt->getab();
						edge *mpba = mppt->getba();
						cell *Cmp = mpab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (mpab->get_attach())
						{
							if (mpab->get_Pair()->get_I1()==Cmp) {C1m = mpab->get_Pair()->get_I2();}
							else {C1m = mpab->get_Pair()->get_I1();}
						}
						if (mpba->get_attach())
						{
							if (mpba->get_Pair()->get_I1()==Cmp) {C2m = mpba->get_Pair()->get_I2();}
							else {C2m = mpba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Cmp->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Cmp->get_pairs()->at(k)->get_I1()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I2()==C) ||
									(Cmp->get_pairs()->at(k)->get_I2()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I1()==C))
								{
									mergepair *Mp = Cmp->get_pairs()->at(k)->MP()->at(0);
									Cmp->get_pairs()->at(k)->remove_MP(Mp);
									break;
								}
							}
						}
					}
					Ipt_NeEg[tspn-2]->set_attach(0);
					point3D *ptp = t_t1->getfp();
					while (ptp!=t_s1->getfp())
					{
						double px = ptp->x();
						double py = ptp->y();
						if ((int)ptp->getmp()->size()==0)
						{
							ptp->setx(px + dpmx*1.0);
							ptp->sety(py + dpmy*1.0);
						}
						else
						{
							ptp->setx(px + dpmx*0.5);
							ptp->sety(py + dpmy*0.5);
							for (int j=0;j<(int)ptp->getmp()->size();j++)
							{
								ptp->getmp()->at(j)->setx(px + dpmx*0.5);
								ptp->getmp()->at(j)->sety(py + dpmy*0.5);
							}
						}
						ptp = ptp->getfp();
					}
				}
				else if (t_t1_angle<=96)
				{
					double tx = Ipt_NePt[tspn-1]->x();
					double ty = Ipt_NePt[tspn-1]->y();
					t_t1->getfp()->setx(tx);
					t_t1->getfp()->sety(ty);
					t_t1->getab()->set_attach(1);
					t_t1->getab()->set_Neighbor(Ipt_NeEg[tspn-2]);
					Ipt_NeEg[tspn-2]->set_Neighbor(t_t1->getab());
					t_t1->getab()->set_Pair(Ipt_NeEg[tspn-2]->get_Pair());
					node_merge(Ipt_NePt[tspn-1],t_t1->getfp());
					vector<point3D*> t1mp;
					t1mp.push_back(Ipt_NePt[tspn-1]);
					t_t1->getfp()->pushmp(Ipt_NePt[tspn-1]);
					for (int j=0;j<(int)Ipt_NePt[tspn-1]->getmp()->size();j++)
					{
						t1mp.push_back(Ipt_NePt[tspn-1]->getmp()->at(j));
						t_t1->getfp()->pushmp(Ipt_NePt[tspn-1]->getmp()->at(j));
					}
					for (int j=0;j<(int)t1mp.size();j++)
					{
						t1mp[j]->pushmp(t_t1->getfp());
						edge *tmab = t1mp[j]->getab();
						edge *tmba = t1mp[j]->getba();
						cell *Ctm = tmab->Cell();
						cell *C1t = NULL;
						cell *C2t = NULL;
						if (tmab->get_attach())
						{
							if (tmab->get_Pair()->get_I1()==Ctm) {C1t = tmab->get_Pair()->get_I2();}
							else {C1t = tmab->get_Pair()->get_I1();}
						}
						if (tmba->get_attach())
						{
							if (tmba->get_Pair()->get_I1()==Ctm) {C2t = tmba->get_Pair()->get_I2();}
							else {C2t = tmba->get_Pair()->get_I1();}
						}
						if (C1t!=C && C2t!=C)
						{
							int cmpn = (int)Ctm->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Ctm->get_pairs()->at(k)->get_I1()==Ctm &&
									 Ctm->get_pairs()->at(k)->get_I2()==C) ||
									(Ctm->get_pairs()->at(k)->get_I2()==Ctm &&
									 Ctm->get_pairs()->at(k)->get_I1()==C))
								{
									if      (Ctm->get_pairs()->at(k)->MP()->at(0)->S1()==t1mp[j])
									{
										Ctm->get_pairs()->at(k)->MP()->at(0)->setS2(t_t1->getfp());
										Ctm->get_pairs()->at(k)->MP()->at(0)->setT2(t_t1->getfp());
									}
									else if (Ctm->get_pairs()->at(k)->MP()->at(0)->S2()==t1mp[j])
									{
										Ctm->get_pairs()->at(k)->MP()->at(0)->setS1(t_t1->getfp());
										Ctm->get_pairs()->at(k)->MP()->at(0)->setT1(t_t1->getfp());
									}
									break;
								}
							}
						}
					}
					point3D *ptp = t_t1->getfp()->getfp();
					while (ptp!=t_s1->getfp())
					{
						double px = ptp->x();
						double py = ptp->y();
						if ((int)ptp->getmp()->size()==0)
						{
							ptp->setx(px + dpmx*1.0);
							ptp->sety(py + dpmy*1.0);
						}
						else
						{
							ptp->setx(px + dpmx*0.5);
							ptp->sety(py + dpmy*0.5);
							for (int j=0;j<(int)ptp->getmp()->size();j++)
							{
								ptp->getmp()->at(j)->setx(px + dpmx*0.5);
								ptp->getmp()->at(j)->sety(py + dpmy*0.5);
							}
						}
						ptp = ptp->getfp();
					}
				}
				////////// interpair end-points update //////////
				if (Iptn>0)
				{
					for (int i=0;i<Iptn;i++)
					{
						if (Ipt[i]->get_I1()==C)
						{
							if (i==Iptn-1)
							{
								point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
								point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
								Ipt[i]->MP()->at(0)->setS1(ps1->getfp());
								if (t_t1->getab()->get_attach())
								{
									Ipt[i]->MP()->at(0)->setT1(t_t1->getfp());
								}
								else
								{
									Ipt[i]->MP()->at(0)->setS2(ps2->getfp());
								}
								ps2 = Ipt[i]->MP()->at(0)->S2();
								point3D *pt2 = Ipt[i]->MP()->at(0)->T2();
								if (ps2!=pt2)
								{
									if (pt2->getrp()==ps2)
									{
										double Fr_m = Ipt[i]->MP()->at(0)->get_Fr_t1();
										Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_m);
									}
									else
									{
										double Fr_m = pt2->getrp()->get_stored_friction_force();
										Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_m);
										point3D *mpt = pt2->getrp();
										while (mpt!=ps2)
										{
											point3D *mpmt = mpt->getmp()->at(0);
											double Fr_m1 = mpt->get_stored_friction_force();
											double Fr_m2 = mpmt->get_stored_friction_force();
											double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
											mpt->set_stored_friction_force(Fr_m0);
											mpmt->set_stored_friction_force(Fr_m0);
											mpt = mpt->getrp();
										}
									}
								}
							}
							else
							{
								point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
								point3D *pt1 = Ipt[i]->MP()->at(0)->T1();
								Ipt[i]->MP()->at(0)->setS1(ps1->getfp());
								Ipt[i]->MP()->at(0)->setT1(pt1->getfp());
								if (ps1!=pt1)
								{
									if (ps1->getfp()==pt1)
									{
										double Fr_t1 = Ipt[i]->MP()->at(0)->get_Fr_t1();
										Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_t1);
									}
									else
									{
										point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
										point3D *pt2 = Ipt[i]->MP()->at(0)->T2();
										double Fr_mt2 = pt2->getrp()->get_stored_friction_force();
										double Fr_ms1 = pt2->getrp()->getmp()->at(0)->get_stored_friction_force();
										double Fr_m = Fr_mt2;if (Fr_mt2<Fr_ms1) {Fr_m = Fr_ms1;}
										Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_m);
										point3D *mpt = pt2->getrp();
										while (mpt!=ps2)
										{
											point3D *mpmt = mpt->getmp()->at(0);
											double Fr_m1 = mpt->get_stored_friction_force();
											double Fr_m2 = mpmt->get_stored_friction_force();
											double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
											mpt->set_stored_friction_force(Fr_m0);
											mpmt->set_stored_friction_force(Fr_m0);
											mpt = mpt->getrp();
										}
									}
								}
							}
						}
						else
						{
							if (i==Iptn-1)
							{
								point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
								point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
								Ipt[i]->MP()->at(0)->setS2(ps2->getfp());
								if (t_t1->getab()->get_attach())
								{
									Ipt[i]->MP()->at(0)->setT2(t_t1->getfp());
								}
								else
								{
									Ipt[i]->MP()->at(0)->setS1(ps1->getfp());
								}
								ps1 = Ipt[i]->MP()->at(0)->S1();
								point3D *pt1 = Ipt[i]->MP()->at(0)->T1();
								if (ps1!=pt1)
								{
									if (pt1->getrp()==ps1)
									{
										double Fr_m = Ipt[i]->MP()->at(0)->get_Fr_s1();
										Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
									}
									else
									{
										double Fr_m = pt1->getrp()->get_stored_friction_force();
										Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
										point3D *mpt = pt1->getrp();
										while (mpt!=ps1)
										{
											point3D *mpmt = mpt->getmp()->at(0);
											double Fr_m1 = mpt->get_stored_friction_force();
											double Fr_m2 = mpmt->get_stored_friction_force();
											double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
											mpt->set_stored_friction_force(Fr_m0);
											mpmt->set_stored_friction_force(Fr_m0);
											mpt = mpt->getrp();
										}
									}
								}
							}
							else
							{
								point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
								point3D *pt2 = Ipt[i]->MP()->at(0)->T2();
								Ipt[i]->MP()->at(0)->setS2(ps2->getfp());
								Ipt[i]->MP()->at(0)->setT2(pt2->getfp());
								point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
								point3D *pt1 = Ipt[i]->MP()->at(0)->T1();
								if (ps1!=pt1)
								{
									if (ps1->getfp()==pt1)
									{
										double Fr_m = Ipt[i]->MP()->at(0)->get_Fr_s1();
										Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
									}
									else
									{
										double Fr_m1 = pt1->getrp()->get_stored_friction_force();
										double Fr_m2 = pt1->getrp()->getmp()->at(0)->get_stored_friction_force();
										double Fr_m = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m = Fr_m2;}
										Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
										point3D *mpt = pt1->getrp();
										while (mpt!=ps1)
										{
											point3D *mpmt = mpt->getmp()->at(0);
											double Fr_m1 = mpt->get_stored_friction_force();
											double Fr_m2 = mpmt->get_stored_friction_force();
											double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
											mpt->set_stored_friction_force(Fr_m0);
											mpmt->set_stored_friction_force(Fr_m0);
											mpt = mpt->getrp();
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else if (Ipsn>0 && Iptn==0)
	{
		point3D *s_end_p = NULL;
		if (Ips[Ipsn-1]->get_I1()==C) {s_end_p = Ips[Ipsn-1]->MP()->at(0)->S1();}
		else                          {s_end_p = Ips[Ipsn-1]->MP()->at(0)->S2();}
		bool run_f = true;
		if (!s_end_p->getba()->get_attach())
		{
			point3D *pmp = s_end_p->getrp();
			while (pmp!=pt)
			{
				if (pmp->getab()->get_attach()) {run_f = false;break;}
				pmp = pmp->getrp();
			}
		}
		else {run_f = false;}
		if (run_f)
		{
			double Ips_Ff_max = 0;
			double Ips_ab_max = 0;
			for (int i=0;i<Ipsn;i++)
			{
				point3D *s1 = Ips[i]->MP()->at(0)->S1();
				point3D *t1 = Ips[i]->MP()->at(0)->T1();
				point3D *s2 = Ips[i]->MP()->at(0)->S2();
				point3D *t2 = Ips[i]->MP()->at(0)->T2();
				double Ff_s = Ips[i]->MP()->at(0)->get_Fr_s1();
				double Ff_t = Ips[i]->MP()->at(0)->get_Fr_t1();
				int Ff_d = Ips[i]->MP()->at(0)->get_Fr_dir();
				double Ips_ab = Ips[i]->get_adhesion();
				double e_force_x_s1 = s1->get_elastic_force(0);
				double e_force_y_s1 = s1->get_elastic_force(1);
				double e_force_x_t1 = t1->get_elastic_force(0);
				double e_force_y_t1 = t1->get_elastic_force(1);
				double e_force_x_s2 = s2->get_elastic_force(0);
				double e_force_y_s2 = s2->get_elastic_force(1);
				double e_force_x_t2 = t2->get_elastic_force(0);
				double e_force_y_t2 = t2->get_elastic_force(1);
				double s1ldx = 0;double s1ldy = 0;double t1ldx = 0;double t1ldy = 0;
				double s2ldx = 0;double s2ldy = 0;double t2ldx = 0;double t2ldy = 0;
				if (Ips[i]->get_I1()==C)
				{
					double s1lx = s1->x() - s1->getrp()->x();
					double s1ly = s1->y() - s1->getrp()->y();
					double s1ld = sqrt(s1lx*s1lx + s1ly*s1ly);
					s1ldx = s1lx/s1ld;
					s1ldy = s1ly/s1ld;
					t2ldx = s1ldx;
					t2ldy = s1ldy;
					double t1lx = t1->x() - t1->getrp()->x();
					double t1ly = t1->y() - t1->getrp()->y();
					double t1ld = sqrt(t1lx*t1lx + t1ly*t1ly);
					t1ldx = t1lx/t1ld;
					t1ldy = t1ly/t1ld;
					s2ldx = t1ldx;
					s2ldy = t1ldy;
				}
				else
				{
					double s2lx = s2->x() - s2->getrp()->x();
					double s2ly = s2->y() - s2->getrp()->y();
					double s2ld = sqrt(s2lx*s2lx + s2ly*s2ly);
					s2ldx = s2lx/s2ld;
					s2ldy = s2ly/s2ld;
					t1ldx = s2ldx;
					t1ldy = s2ldy;
					double t2lx = t2->x() - t2->getrp()->x();
					double t2ly = t2->y() - t2->getrp()->y();
					double t2ld = sqrt(t2lx*t2lx + t2ly*t2ly);
					t2ldx = t2lx/t2ld;
					t2ldy = t2ly/t2ld;
					s1ldx = t2ldx;
					s1ldy = t2ldy;
				}
				double Ff_s1_inc = abs(e_force_x_s1*s1ldx + e_force_y_s1*s1ldy);
				double Ff_t2_inc = abs(e_force_x_t2*t2ldx + e_force_y_t2*t2ldy);
				double Ff_s2_inc = abs(e_force_x_s2*s2ldx + e_force_y_s2*s2ldy);
				double Ff_t1_inc = abs(e_force_x_t1*t1ldx + e_force_y_t1*t1ldy);
				double Ff_s1_t2_inc = 0;
				if (Ff_s1_inc>Ff_t2_inc) {Ff_s1_t2_inc = Ff_s1_inc;}
				else                     {Ff_s1_t2_inc = Ff_t2_inc;}
				double Ff_s2_t1_inc = 0;
				if (Ff_s2_inc>Ff_t1_inc) {Ff_s2_t1_inc = Ff_s2_inc;}
				else                     {Ff_s2_t1_inc = Ff_t1_inc;}
				if (Ips_ab_max<Ips_ab) {Ips_ab_max = Ips_ab;}
				if (Ips[i]->get_I1()==C)
				{
					if      (Ff_d==-1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==0)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					double C_Fr_s1 = Ips[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ips[i]->MP()->at(0)->get_Fr_t1();
					if (C_Fr_s1>Ips_Ff_max) {Ips_Ff_max = C_Fr_s1;}
					if (C_Fr_t1>Ips_Ff_max) {Ips_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = t1->getrp();
					while (ptC!=s1)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getrp()->x() - ptC->x();
						double ply = ptC->getrp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==0)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Fr_tC>Ips_Ff_max) {Ips_Ff_max = Fr_tC;}
						if (Fr_mC>Ips_Ff_max) {Ips_Ff_max = Fr_mC;}
						ptC = ptC->getrp();
					}
				}
				else
				{
					if      (Ff_d==-1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(1);
					}
					else if (Ff_d==1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					else if (Ff_d==0)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(1);
					}
					double C_Fr_s1 = Ips[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ips[i]->MP()->at(0)->get_Fr_t1();
					if (C_Fr_s1>Ips_Ff_max) {Ips_Ff_max = C_Fr_s1;}
					if (C_Fr_t1>Ips_Ff_max) {Ips_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = t2->getrp();
					while (ptC!=s2)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getrp()->x() - ptC->x();
						double ply = ptC->getrp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==1)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Fr_tC>Ips_Ff_max) {Ips_Ff_max = Fr_tC;}
						if (Fr_mC>Ips_Ff_max) {Ips_Ff_max = Fr_mC;}
						ptC = ptC->getrp();
					}
				}
			}
			/////// boundary edges sliding ///////
			if (Ips_Ff_max>Ips_ab_max && front_b)
			{
				point3D *s_t1 = NULL;
				point3D *s_s1 = NULL;
				vector<point3D*> Ips_Pt;
				vector<point3D*> Ips_NePt;
				vector<interpair*> Ips_Ip;
				vector<edge*> Ips_NeEg;
				vector<node*> Ips_Nd;
				if (Ips[0]->get_I1()==C) {s_t1 = Ips[0]->MP()->at(0)->T1();}
				else                     {s_t1 = Ips[0]->MP()->at(0)->T2();}
				if (Ips[Ipsn-1]->get_I1()==C) {s_s1 = Ips[Ipsn-1]->MP()->at(0)->S1();}
				else                          {s_s1 = Ips[Ipsn-1]->MP()->at(0)->S2();}
				double m_l = 10000;
				double s_s1_ba_nx = s_s1->getba()->getnormal(0);
				double s_s1_ba_ny = s_s1->getba()->getnormal(1);
				double s_s1_angle = angle - vector2angle(s_s1_ba_nx,s_s1_ba_ny);
				if (s_s1_angle<0) {s_s1_angle += 360;}
				////////////// heat up //////////////
				point3D *s_tm = s_t1;
				int sspn = 0;
				while (s_tm!=s_s1->getrp())
				{
					Ips_Pt.push_back(s_tm);
					node* sNd = s_tm->get_node_id();
					Ips_Nd.push_back(sNd);
					sNd->removevertex(s_tm);
					if ((int)s_tm->getmp()->size()==0)
					{
						Ips_NePt.push_back(NULL);
					}
					else
					{
						Ips_NePt.push_back(s_tm->getmp()->at(0));
						for (int i=0;i<(int)s_tm->getmp()->size();i++)
						{
							s_tm->getmp()->at(i)->removemp(s_tm);
						}
						s_tm->clearmp(); // clear the neighbors
					}
					if (s_tm->getba()->get_attach())
					{
						double neegl = s_tm->getba()->getlength();if (neegl<m_l) {m_l = neegl;}
						Ips_NeEg.push_back(s_tm->getba()->get_Neighbor());
						Ips_Ip.push_back(s_tm->getba()->get_Pair());
					}
					else
					{
						Ips_NeEg.push_back(NULL);
						Ips_Ip.push_back(NULL);
					}
					s_tm = s_tm->getrp();
					sspn++;
				}
				////////////// single point pair removal //////////////
				int Iptn_s = (int)Ipt_single.size();
				for (int i=0;i<Iptn_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ipt_single[i]->get_I1()==C)
					{
						mergepair *Mp = Ipt_single[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ipt_single[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[i]->remove_MP(Mp);
					}
				}
				int Ipsn_s = (int)Ips_single.size();
				for (int i=0;i<Ipsn_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ips_single[i]->get_I1()==C)
					{
						mergepair *Mp = Ips_single[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ips_single[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[i]->remove_MP(Mp);
					}
				}
				int Ipfs_s = (int)front_single_pairs.size();
				for (int i=0;i<Ipfs_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (front_single_pairs[i]->get_I1()==C)
					{
						mergepair *Mp = front_single_pairs[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = front_single_pairs[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[i]->remove_MP(Mp);
					}
				}
				////////////// position refresh /////////////
				for (int i=sspn-1;i>=0;i--)
				{
					if (i==0)
					{
						int nd_size = (int)nodeList.size();
						node *A = new node(nd_size);
						nodeList.push_back(A);
						nodeList[nd_size]->pushvertex(Ips_Pt[i]);
						nodeList[nd_size]->set_in_global(1);
						Ips_Pt[i]->set_node_id(nodeList[nd_size]);
						Ips_Pt[i]->getba()->set_attach(0);
					}
					else
					{
						/////////// node ////////////
						Ips_Nd[i-1]->pushvertex(Ips_Pt[i]);
						Ips_Pt[i]->set_node_id(Ips_Nd[i-1]);
						double fpx = 0;
						double fpy = 0;
						if (Ips_NePt[i-1]!=NULL)
						{
							fpx = Ips_NePt[i-1]->x();
							fpy = Ips_NePt[i-1]->y();
						}
						else
						{
							fpx = Ips_Pt[i-1]->x();
							fpy = Ips_Pt[i-1]->y();
						}
						Ips_Pt[i]->setx(fpx);
						Ips_Pt[i]->sety(fpy);
						if (Ips_NePt[i-1]!=NULL)
						{
							vector<point3D*> fpmp;
							fpmp.push_back(Ips_NePt[i-1]);
							Ips_Pt[i]->pushmp(Ips_NePt[i-1]);
							for (int j=0;j<(int)Ips_NePt[i-1]->getmp()->size();j++)
							{
								fpmp.push_back(Ips_NePt[i-1]->getmp()->at(j));
								Ips_Pt[i]->pushmp(Ips_NePt[i-1]->getmp()->at(j));
							}
							for (int j=0;j<(int)fpmp.size();j++)
							{
								fpmp[j]->pushmp(Ips_Pt[i]);
								edge *fpmp_ab = fpmp[j]->getab();
								edge *fpmp_ba = fpmp[j]->getba();
								cell *Cfp = fpmp_ab->Cell();
								cell *C1f = NULL;
								cell *C2f = NULL;
								if (fpmp_ab->get_attach())
								{
									if (fpmp_ab->get_Pair()->get_I1()==Cfp) {C1f = fpmp_ab->get_Pair()->get_I2();}
									else {C1f = fpmp_ab->get_Pair()->get_I1();}
								}
								if (fpmp_ba->get_attach())
								{
									if (fpmp_ba->get_Pair()->get_I1()==Cfp) {C2f = fpmp_ba->get_Pair()->get_I2();}
									else {C2f = fpmp_ba->get_Pair()->get_I1();}
								}
								if (C1f!=C && C2f!=C)
								{
									int cfpn = (int)Cfp->get_pairs()->size();
									for (int k=0;k<cfpn;k++)
									{
										if ((Cfp->get_pairs()->at(k)->get_I1()==Cfp &&
											 Cfp->get_pairs()->at(k)->get_I2()==C) ||
											(Cfp->get_pairs()->at(k)->get_I2()==Cfp &&
											 Cfp->get_pairs()->at(k)->get_I1()==C))
										{
											if ((int)Cfp->get_pairs()->at(k)->MP()->size()>0)
											{
												if      (Cfp->get_pairs()->at(k)->MP()->at(0)->S1()==fpmp[j])
												{
													Cfp->get_pairs()->at(k)->MP()->at(0)->setS2(Ips_Pt[i]);
													Cfp->get_pairs()->at(k)->MP()->at(0)->setT2(Ips_Pt[i]);
												}
												else if (Cfp->get_pairs()->at(k)->MP()->at(0)->S2()==fpmp[j])
												{
													Cfp->get_pairs()->at(k)->MP()->at(0)->setS1(Ips_Pt[i]);
													Cfp->get_pairs()->at(k)->MP()->at(0)->setT1(Ips_Pt[i]);
												}
											}
											break;
										}
									}
								}
							}
							fpmp.clear();
						}
						if (i<sspn-1)
						{
							///////// edge ///////////
							if (Ips_NeEg[i-1]==NULL)
							{
								Ips_Pt[i]->getba()->set_attach(0);
							}
							else
							{
								Ips_Pt[i]->getba()->set_attach(1);
								Ips_NeEg[i-1]->set_attach(1);
								Ips_Pt[i]->getba()->set_Neighbor(Ips_NeEg[i-1]);
								Ips_NeEg[i-1]->set_Neighbor(Ips_Pt[i]->getba());
							}
							////////// pair //////////
							Ips_Pt[i]->getba()->set_Pair(Ips_Ip[i-1]);
						}
						else
						{
							///////// edge ///////////
							Ips_Pt[i]->getba()->set_attach(0);
							////////// pair //////////
							Ips_Pt[i]->getba()->set_Pair(NULL);
						}
					}
				}
				/////// points other than points from s_s1->s_t1 ///////
				double dpmx = dx*m_l;
				double dpmy = dy*m_l;
				if      (s_s1_angle>96)
				{
					vector<point3D*> smpt;
					smpt.push_back(Ips_NePt[sspn-1]);
					for (int j=0;j<(int)Ips_NePt[sspn-1]->getmp()->size();j++)
					{smpt.push_back(Ips_NePt[sspn-1]->getmp()->at(j));}
					for (int j=0;j<(int)smpt.size();j++)
					{
						point3D *mppt = smpt[j];
						edge *mpab = mppt->getab();
						edge *mpba = mppt->getba();
						cell *Cmp = mpab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (mpab->get_attach())
						{
							if (mpab->get_Pair()->get_I1()==Cmp) {C1m = mpab->get_Pair()->get_I2();}
							else {C1m = mpab->get_Pair()->get_I1();}
						}
						if (mpba->get_attach())
						{
							if (mpba->get_Pair()->get_I1()==Cmp) {C2m = mpba->get_Pair()->get_I2();}
							else {C2m = mpba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Cmp->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Cmp->get_pairs()->at(k)->get_I1()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I2()==C) ||
									(Cmp->get_pairs()->at(k)->get_I2()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I1()==C))
								{
									mergepair *Mp = Cmp->get_pairs()->at(k)->MP()->at(0);
									Cmp->get_pairs()->at(k)->remove_MP(Mp);
									break;
								}
							}
						}
					}
					Ips_NeEg[sspn-2]->set_attach(0);
					point3D *ptp = s_s1->getrp();
					while (ptp!=s_t1->getrp())
					{
						ptp = ptp->getrp();
						double px = ptp->x();
						double py = ptp->y();
						if ((int)ptp->getmp()->size()==0)
						{
							ptp->setx(px + dpmx*1.0);
							ptp->sety(py + dpmy*1.0);
						}
						else
						{
							ptp->setx(px + dpmx*0.5);
							ptp->sety(py + dpmy*0.5);
							for (int j=0;j<(int)ptp->getmp()->size();j++)
							{
								ptp->getmp()->at(j)->setx(px + dpmx*0.5);
								ptp->getmp()->at(j)->sety(py + dpmy*0.5);
							}
						}
					}
				}
				else if (s_s1_angle<=96)
				{
					double sx = Ips_NePt[sspn-1]->x();
					double sy = Ips_NePt[sspn-1]->y();
					s_s1->getrp()->setx(sx);
					s_s1->getrp()->sety(sy);
					s_s1->getba()->set_attach(1);
					s_s1->getba()->set_Neighbor(Ips_NeEg[sspn-2]);
					Ips_NeEg[sspn-2]->set_Neighbor(s_s1->getba());
					s_s1->getba()->set_Pair(Ips_NeEg[sspn-2]->get_Pair());
					node_merge(Ips_NePt[sspn-1],s_s1->getrp());
					vector<point3D*> s1mp;
					s1mp.push_back(Ips_NePt[sspn-1]);
					s_s1->getrp()->pushmp(Ips_NePt[sspn-1]);
					for (int j=0;j<(int)Ips_NePt[sspn-1]->getmp()->size();j++)
					{
						s1mp.push_back(Ips_NePt[sspn-1]->getmp()->at(j));
						s_s1->getrp()->pushmp(Ips_NePt[sspn-1]->getmp()->at(j));
					}
					for (int j=0;j<(int)s1mp.size();j++)
					{
						s1mp[j]->pushmp(s_s1->getrp());
						edge *smab = s1mp[j]->getab();
						edge *smba = s1mp[j]->getba();
						cell *Csm = smab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (smab->get_attach())
						{
							if (smab->get_Pair()->get_I1()==Csm) {C1m = smab->get_Pair()->get_I2();}
							else {C1m = smab->get_Pair()->get_I1();}
						}
						if (smba->get_attach())
						{
							if (smba->get_Pair()->get_I1()==Csm) {C2m = smba->get_Pair()->get_I2();}
							else {C2m = smba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Csm->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Csm->get_pairs()->at(k)->get_I1()==Csm &&
									 Csm->get_pairs()->at(k)->get_I2()==C) ||
									(Csm->get_pairs()->at(k)->get_I2()==Csm &&
									 Csm->get_pairs()->at(k)->get_I1()==C))
								{
									if      (Csm->get_pairs()->at(k)->MP()->at(0)->S1()==s1mp[j])
									{
										Csm->get_pairs()->at(k)->MP()->at(0)->setS2(s_s1->getrp());
										Csm->get_pairs()->at(k)->MP()->at(0)->setT2(s_s1->getrp());
									}
									else if (Csm->get_pairs()->at(k)->MP()->at(0)->S2()==s1mp[j])
									{
										Csm->get_pairs()->at(k)->MP()->at(0)->setS1(s_s1->getrp());
										Csm->get_pairs()->at(k)->MP()->at(0)->setT1(s_s1->getrp());
									}
									break;
								}
							}
						}
					}
					point3D *ptp = s_s1->getrp()->getrp();
					while (ptp!=s_t1->getrp())
					{
						double px = ptp->x();
						double py = ptp->y();
						if ((int)ptp->getmp()->size()==0)
						{
							ptp->setx(px + dpmx*1.0);
							ptp->sety(py + dpmy*1.0);
						}
						else
						{
							ptp->setx(px + dpmx*0.5);
							ptp->sety(py + dpmy*0.5);
							for (int j=0;j<(int)ptp->getmp()->size();j++)
							{
								ptp->getmp()->at(j)->setx(px + dpmx*0.5);
								ptp->getmp()->at(j)->sety(py + dpmy*0.5);
							}
						}
						ptp = ptp->getrp();
					}
				}
				////////// interpair end-points update //////////
				for (int i=0;i<Ipsn;i++)
				{
					if (Ips[i]->get_I1()==C)
					{
						if (i==Ipsn-1)
						{
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							Ips[i]->MP()->at(0)->setT1(pt1->getrp());
							if (s_s1->getba()->get_attach())
							{
								Ips[i]->MP()->at(0)->setS1(s_s1->getrp());
							}
							else
							{
								Ips[i]->MP()->at(0)->setT2(pt2->getrp());
							}
							pt2 = Ips[i]->MP()->at(0)->T2();
							point3D *ps2 = Ips[i]->MP()->at(0)->S2();
							if (ps2!=pt2)
							{
								if (ps2->getfp()==pt2)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_s1();
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps2->getfp()->get_stored_friction_force();
									double Fr_m2 = ps2->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m0);
									point3D *mpt = ps2->getfp();
									while (mpt!=pt2)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
						else
						{
							point3D *ps1 = Ips[i]->MP()->at(0)->S1();
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							Ips[i]->MP()->at(0)->setS1(ps1->getrp());
							Ips[i]->MP()->at(0)->setT1(pt1->getrp());
							point3D *ps2 = Ips[i]->MP()->at(0)->S2();
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							if (ps2!=pt2)
							{
								if (ps2->getfp()==pt2)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_s1();
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps2->getfp()->get_stored_friction_force();
									double Fr_m2 = ps2->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m0);
									point3D *mpt = ps2->getfp();
									while (mpt!=pt2)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
					}
					else
					{
						if (i==Ipsn-1)
						{
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							Ips[i]->MP()->at(0)->setT2(pt2->getrp());
							if (s_s1->getba()->get_attach())
							{
								Ips[i]->MP()->at(0)->setS2(s_s1->getrp());
							}
							else
							{
								Ips[i]->MP()->at(0)->setT1(pt1->getrp());
							}
							pt1 = Ips[i]->MP()->at(0)->T1();
							point3D *ps1 = Ips[i]->MP()->at(0)->S1();
							if (ps1!=pt1)
							{
								if (ps1->getfp()==pt1)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_t1();
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps1->getfp()->get_stored_friction_force();
									double Fr_m2 = ps1->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m0);
									point3D *mpt = ps1->getfp();
									while (mpt!=pt1)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
						else
						{
							point3D *ps2 = Ips[i]->MP()->at(0)->S2();
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							Ips[i]->MP()->at(0)->setS2(ps2->getrp());
							Ips[i]->MP()->at(0)->setT2(pt2->getrp());
							point3D *ps1 = Ips[i]->MP()->at(0)->S1();
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							if (ps1!=pt1)
							{
								if (ps1->getfp()==pt1)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_t1();
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps1->getfp()->get_stored_friction_force();
									double Fr_m2 = ps1->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m0);
									point3D *mpt = ps1->getfp();
									while (mpt!=pt1)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else if (Ipsn>0 && Iptn>0)
	{
		////////// check validity //////////
		bool run_f = true;
		point3D *t_end_p = NULL;
		point3D *s_end_p = NULL;
		if (Ipt[Iptn-1]->get_I1()==C) {t_end_p = Ipt[Iptn-1]->MP()->at(0)->T1();}
		else                          {t_end_p = Ipt[Iptn-1]->MP()->at(0)->T2();}
		if (Ips[Ipsn-1]->get_I1()==C) {s_end_p = Ips[Ipsn-1]->MP()->at(0)->S1();}
		else                          {s_end_p = Ips[Ipsn-1]->MP()->at(0)->S2();}
		if (!t_end_p->getab()->get_attach() && !s_end_p->getba()->get_attach())
		{
			if (t_end_p->getfp()!=s_end_p)
			{
				point3D *pmp = t_end_p->getfp();
				while (pmp!=s_end_p)
				{
					int mpn = (int)pmp->getmp()->size();
					if (mpn>0) {run_f = false;break;}
					pmp = pmp->getfp();
				}
			}
		}
		else {run_f = false;}
		if (run_f)
		{
			double Ipt_Ff_max = 0;
			double Ips_Ff_max = 0;
			double Ipt_ab_max = 0;
			double Ips_ab_max = 0;
			for (int i=0;i<Iptn;i++)
			{
				point3D *s1 = Ipt[i]->MP()->at(0)->S1();
				point3D *t1 = Ipt[i]->MP()->at(0)->T1();
				point3D *s2 = Ipt[i]->MP()->at(0)->S2();
				point3D *t2 = Ipt[i]->MP()->at(0)->T2();
				double Ff_s = Ipt[i]->MP()->at(0)->get_Fr_s1();
				double Ff_t = Ipt[i]->MP()->at(0)->get_Fr_t1();
				int Ff_d = Ipt[i]->MP()->at(0)->get_Fr_dir(); // -1:null, 1:s->t, 0:s<-t
				double Ipt_ab = Ipt[i]->get_adhesion();
				/**********************************************
				|               |  Ipt
				o---o---o---o---o
				o-<-o-<-o-<-o-<-o---o t
				t1              s1   \
									  o ==>
						C			  |
									  o ==>
									 /
				o---o---o---o---o---o s

				Friction direction should be opposite to the 
				direction of migration: 1 s1->t1
				**********************************************/
				double e_force_x_s1 = s1->get_elastic_force(0);
				double e_force_y_s1 = s1->get_elastic_force(1);
				double e_force_x_t1 = t1->get_elastic_force(0);
				double e_force_y_t1 = t1->get_elastic_force(1);
				double e_force_x_s2 = s2->get_elastic_force(0);
				double e_force_y_s2 = s2->get_elastic_force(1);
				double e_force_x_t2 = t2->get_elastic_force(0);
				double e_force_y_t2 = t2->get_elastic_force(1);
				double s1ldx = 0;double s1ldy = 0;double t1ldx = 0;double t1ldy = 0;
				double s2ldx = 0;double s2ldy = 0;double t2ldx = 0;double t2ldy = 0;
				if (Ipt[i]->get_I1()==C)
				{
					double s1lx = s1->getfp()->x() - s1->x();
					double s1ly = s1->getfp()->y() - s1->y();
					double s1ld = sqrt(s1lx*s1lx + s1ly*s1ly);
					s1ldx = s1lx/s1ld;
					s1ldy = s1ly/s1ld;
					t2ldx = s1ldx;
					t2ldy = s1ldy;
					double t1lx = t1->x() - t1->getrp()->x();
					double t1ly = t1->y() - t1->getrp()->y();
					double t1ld = sqrt(t1lx*t1lx + t1ly*t1ly);
					t1ldx = t1lx/t1ld;
					t1ldy = t1ly/t1ld;
					s2ldx = t1ldx;
					s2ldy = t1ldy;
				}
				else
				{
					double s2lx = s2->getfp()->x() - s2->x();
					double s2ly = s2->getfp()->y() - s2->y();
					double s2ld = sqrt(s2lx*s2lx + s2ly*s2ly);
					s2ldx = s2lx/s2ld;
					s2ldy = s2ly/s2ld;
					t1ldx = s2ldx;
					t1ldy = s2ldy;
					double t2lx = t2->x() - t2->getrp()->x();
					double t2ly = t2->y() - t2->getrp()->y();
					double t2ld = sqrt(t2lx*t2lx + t2ly*t2ly);
					t2ldx = t2lx/t2ld;
					t2ldy = t2ly/t2ld;
					s1ldx = t2ldx;
					s1ldy = t2ldy;
				}
				double Ff_s1_inc = abs(e_force_x_s1*s1ldx + e_force_y_s1*s1ldy);
				double Ff_t2_inc = abs(e_force_x_t2*t2ldx + e_force_y_t2*t2ldy);
				double Ff_s2_inc = abs(e_force_x_s2*s2ldx + e_force_y_s2*s2ldy);
				double Ff_t1_inc = abs(e_force_x_t1*t1ldx + e_force_y_t1*t1ldy);
				double Ff_s1_t2_inc = 0;
				if (Ff_s1_inc>Ff_t2_inc) {Ff_s1_t2_inc = Ff_s1_inc;}
				else                     {Ff_s1_t2_inc = Ff_t2_inc;}
				double Ff_s2_t1_inc = 0;
				if (Ff_s2_inc>Ff_t1_inc) {Ff_s2_t1_inc = Ff_s2_inc;}
				else                     {Ff_s2_t1_inc = Ff_t1_inc;}
				if (Ipt_ab_max<Ipt_ab) {Ipt_ab_max = Ipt_ab;}
				if (Ipt[i]->get_I1()==C)
				{
					if      (Ff_d==-1) // null
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(1);
					}
					else if (Ff_d==1)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					else if (Ff_d==0)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(1);
					}
					double C_Fr_s1 = Ipt[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ipt[i]->MP()->at(0)->get_Fr_t1();
					if (Ipt_Ff_max<C_Fr_s1) {Ipt_Ff_max = C_Fr_s1;}
					if (Ipt_Ff_max<C_Fr_t1) {Ipt_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = s1->getfp();
					while (ptC!=t1)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getfp()->x() - ptC->x();
						double ply = ptC->getfp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==1)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Ipt_Ff_max<Fr_tC) {Ipt_Ff_max = Fr_tC;}
						if (Ipt_Ff_max<Fr_mC) {Ipt_Ff_max = Fr_mC;}
						ptC = ptC->getfp();
					}
				}
				else
				{
					if      (Ff_d==-1)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==1)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ipt[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==0)
					{
						Ipt[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ipt[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					double C_Fr_s1 = Ipt[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ipt[i]->MP()->at(0)->get_Fr_t1();
					if (Ipt_Ff_max<C_Fr_s1) {Ipt_Ff_max = C_Fr_s1;}
					if (Ipt_Ff_max<C_Fr_t1) {Ipt_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = s2->getfp();
					while (ptC!=t2)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getfp()->x() - ptC->x();
						double ply = ptC->getfp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==0)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Ipt_Ff_max<Fr_tC) {Ipt_Ff_max = Fr_tC;}
						if (Ipt_Ff_max<Fr_mC) {Ipt_Ff_max = Fr_mC;}
						ptC = ptC->getfp();
					}
				}
			}
			for (int i=0;i<Ipsn;i++)
			{
				point3D *s1 = Ips[i]->MP()->at(0)->S1();
				point3D *t1 = Ips[i]->MP()->at(0)->T1();
				point3D *s2 = Ips[i]->MP()->at(0)->S2();
				point3D *t2 = Ips[i]->MP()->at(0)->T2();
				double Ff_s = Ips[i]->MP()->at(0)->get_Fr_s1();
				double Ff_t = Ips[i]->MP()->at(0)->get_Fr_t1();
				int Ff_d = Ips[i]->MP()->at(0)->get_Fr_dir();
				double Ips_ab = Ips[i]->get_adhesion();
				double e_force_x_s1 = s1->get_elastic_force(0);
				double e_force_y_s1 = s1->get_elastic_force(1);
				double e_force_x_t1 = t1->get_elastic_force(0);
				double e_force_y_t1 = t1->get_elastic_force(1);
				double e_force_x_s2 = s2->get_elastic_force(0);
				double e_force_y_s2 = s2->get_elastic_force(1);
				double e_force_x_t2 = t2->get_elastic_force(0);
				double e_force_y_t2 = t2->get_elastic_force(1);
				double s1ldx = 0;double s1ldy = 0;double t1ldx = 0;double t1ldy = 0;
				double s2ldx = 0;double s2ldy = 0;double t2ldx = 0;double t2ldy = 0;
				if (Ips[i]->get_I1()==C)
				{
					double s1lx = s1->x() - s1->getfp()->x();
					double s1ly = s1->y() - s1->getfp()->y();
					double s1ld = sqrt(s1lx*s1lx + s1ly*s1ly);
					s1ldx = s1lx/s1ld;
					s1ldy = s1ly/s1ld;
					t2ldx = s1ldx;
					t2ldy = s1ldy;
					double t1lx = t1->getrp()->x() - t1->x();
					double t1ly = t1->getrp()->y() - t1->y();
					double t1ld = sqrt(t1lx*t1lx + t1ly*t1ly);
					t1ldx = t1lx/t1ld;
					t1ldy = t1ly/t1ld;
					s2ldx = t1ldx;
					s2ldy = t1ldy;
				}
				else
				{
					double s2lx = s2->x() - s2->getfp()->x();
					double s2ly = s2->y() - s2->getfp()->y();
					double s2ld = sqrt(s2lx*s2lx + s2ly*s2ly);
					s2ldx = s2lx/s2ld;
					s2ldy = s2ly/s2ld;
					t1ldx = s2ldx;
					t1ldy = s2ldy;
					double t2lx = t2->getrp()->x() - t2->x();
					double t2ly = t2->getrp()->y() - t2->y();
					double t2ld = sqrt(t2lx*t2lx + t2ly*t2ly);
					t2ldx = t2lx/t2ld;
					t2ldy = t2ly/t2ld;
					s1ldx = t2ldx;
					s1ldy = t2ldy;
				}
				double Ff_s1_inc = abs(e_force_x_s1*s1ldx + e_force_y_s1*s1ldy);
				double Ff_t2_inc = abs(e_force_x_t2*t2ldx + e_force_y_t2*t2ldy);
				double Ff_s2_inc = abs(e_force_x_s2*s2ldx + e_force_y_s2*s2ldy);
				double Ff_t1_inc = abs(e_force_x_t1*t1ldx + e_force_y_t1*t1ldy);
				double Ff_s1_t2_inc = 0;
				if (Ff_s1_inc>Ff_t2_inc) {Ff_s1_t2_inc = Ff_s1_inc;}
				else                     {Ff_s1_t2_inc = Ff_t2_inc;}
				double Ff_s2_t1_inc = 0;
				if (Ff_s2_inc>Ff_t1_inc) {Ff_s2_t1_inc = Ff_s2_inc;}
				else                     {Ff_s2_t1_inc = Ff_t1_inc;}
				if (Ips_ab_max<Ips_ab) {Ips_ab_max = Ips_ab;}
				if (Ips[i]->get_I1()==C)
				{
					if      (Ff_d==-1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(0);
					}
					else if (Ff_d==0)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					double C_Fr_s1 = Ips[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ips[i]->MP()->at(0)->get_Fr_t1();
					if (C_Fr_s1>Ips_Ff_max) {Ips_Ff_max = C_Fr_s1;}
					if (C_Fr_t1>Ips_Ff_max) {Ips_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = t1->getrp();
					while (ptC!=s1)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getrp()->x() - ptC->x();
						double ply = ptC->getrp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==0)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Fr_tC>Ips_Ff_max) {Ips_Ff_max = Fr_tC;}
						if (Fr_mC>Ips_Ff_max) {Ips_Ff_max = Fr_mC;}
						ptC = ptC->getrp();
					}
				}
				else
				{
					if      (Ff_d==-1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(1);
					}
					else if (Ff_d==1)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc + Ff_s);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc + Ff_t);
					}
					else if (Ff_d==0)
					{
						Ips[i]->MP()->at(0)->set_Fr_s1(Ff_s1_t2_inc);
						Ips[i]->MP()->at(0)->set_Fr_t1(Ff_s2_t1_inc);
						Ips[i]->MP()->at(0)->set_Fr_dir(1);
					}
					double C_Fr_s1 = Ips[i]->MP()->at(0)->get_Fr_s1();
					double C_Fr_t1 = Ips[i]->MP()->at(0)->get_Fr_t1();
					if (C_Fr_s1>Ips_Ff_max) {Ips_Ff_max = C_Fr_s1;}
					if (C_Fr_t1>Ips_Ff_max) {Ips_Ff_max = C_Fr_t1;}
					/////// for middle points ///////
					point3D *ptC = t2->getrp();
					while (ptC!=s2)
					{
						point3D *pmC = ptC->getmp()->at(0);
						double plx = ptC->getrp()->x() - ptC->x();
						double ply = ptC->getrp()->y() - ptC->y();
						double pld = sqrt(plx*plx + ply*ply);
						double pldx = plx/pld;
						double pldy = ply/pld;
						double e_pt_x = ptC->get_elastic_force(0);
						double e_pt_y = ptC->get_elastic_force(1);
						double e_pm_x = pmC->get_elastic_force(0);
						double e_pm_y = pmC->get_elastic_force(1);
						double Ff_pt_inc = abs(e_pt_x*pldx + e_pt_y*pldy);
						double Ff_pm_inc = abs(e_pm_x*pldx + e_pm_y*pldy);
						double Ff_p_inc = 0;
						if (Ff_pt_inc>Ff_pm_inc) {Ff_p_inc = Ff_pt_inc;}
						else                     {Ff_p_inc = Ff_pm_inc;}
						if (Ff_d==1)
						{
							double Ff_pt = ptC->get_stored_friction_force();
							double Ff_pm = pmC->get_stored_friction_force();
							ptC->set_stored_friction_force(Ff_pt + Ff_p_inc);
							pmC->set_stored_friction_force(Ff_pm + Ff_p_inc);
						}
						else
						{
							ptC->set_stored_friction_force(Ff_p_inc);
							pmC->set_stored_friction_force(Ff_p_inc);
						}
						double Fr_tC = ptC->get_stored_friction_force();
						double Fr_mC = pmC->get_stored_friction_force();
						if (Fr_tC>Ips_Ff_max) {Ips_Ff_max = Fr_tC;}
						if (Fr_mC>Ips_Ff_max) {Ips_Ff_max = Fr_mC;}
						ptC = ptC->getrp();
					}
				}
			}
			/////// boundary edges sliding ///////
			if ((Ipt_Ff_max>Ipt_ab_max || Ips_Ff_max>Ips_ab_max) && front_b)
			{
				point3D *t_s1 = NULL;
				point3D *t_t1 = NULL;
				point3D *s_t1 = NULL;
				point3D *s_s1 = NULL;
				vector<point3D*> Ipt_Pt;
				vector<point3D*> Ipt_NePt;
				vector<interpair*> Ipt_Ip;
				vector<edge*> Ipt_NeEg;
				vector<node*> Ipt_Nd;
				vector<point3D*> Ips_Pt;
				vector<point3D*> Ips_NePt;
				vector<interpair*> Ips_Ip;
				vector<edge*> Ips_NeEg;
				vector<node*> Ips_Nd;
				if (Ipt[0]->get_I1()==C) {t_s1 = Ipt[0]->MP()->at(0)->S1();}
				else                     {t_s1 = Ipt[0]->MP()->at(0)->S2();}
				if (Ipt[Iptn-1]->get_I1()==C) {t_t1 = Ipt[Iptn-1]->MP()->at(0)->T1();}
				else                          {t_t1 = Ipt[Iptn-1]->MP()->at(0)->T2();}
				if (Ips[0]->get_I1()==C) {s_t1 = Ips[0]->MP()->at(0)->T1();}
				else                     {s_t1 = Ips[0]->MP()->at(0)->T2();}
				if (Ips[Ipsn-1]->get_I1()==C) {s_s1 = Ips[Ipsn-1]->MP()->at(0)->S1();}
				else                          {s_s1 = Ips[Ipsn-1]->MP()->at(0)->S2();}
				double m_l = 10000;
				double t_t1_ab_nx = t_t1->getab()->getnormal(0);
				double t_t1_ab_ny = t_t1->getab()->getnormal(1);
				double t_t1_angle = vector2angle(t_t1_ab_nx,t_t1_ab_ny) - angle;
				double s_s1_ba_nx = s_s1->getba()->getnormal(0);
				double s_s1_ba_ny = s_s1->getba()->getnormal(1);
				double s_s1_angle = angle - vector2angle(s_s1_ba_nx,s_s1_ba_ny);
				if (t_t1_angle<0) {t_t1_angle += 360;}
				if (s_s1_angle<0) {s_s1_angle += 360;}
				point3D *t_tm = t_s1;
				int tspn = 0;
				while (t_tm!=t_t1->getfp())
				{
					Ipt_Pt.push_back(t_tm);
					node* tNd = t_tm->get_node_id();
					Ipt_Nd.push_back(tNd);
					tNd->removevertex(t_tm);
					if ((int)t_tm->getmp()->size()==0)
					{
						Ipt_NePt.push_back(NULL);
					}
					else
					{
						Ipt_NePt.push_back(t_tm->getmp()->at(0));
						for (int i=0;i<(int)t_tm->getmp()->size();i++)
						{
							t_tm->getmp()->at(i)->removemp(t_tm);
						}
						t_tm->clearmp(); // clear the neighbors
					}
					if (t_tm->getab()->get_attach())
					{
						double neegl = t_tm->getab()->getlength();if (neegl<m_l) {m_l = neegl;}
						Ipt_NeEg.push_back(t_tm->getab()->get_Neighbor());
						Ipt_Ip.push_back(t_tm->getab()->get_Pair());
					}
					else
					{
						Ipt_NeEg.push_back(NULL);
						Ipt_Ip.push_back(NULL);
					}
					t_tm = t_tm->getfp();
					tspn++;
				}
				point3D *s_tm = s_t1;
				int sspn = 0;
				while (s_tm!=s_s1->getrp())
				{
					Ips_Pt.push_back(s_tm);
					node* sNd = s_tm->get_node_id();
					Ips_Nd.push_back(sNd);
					sNd->removevertex(s_tm);
					if ((int)s_tm->getmp()->size()==0)
					{
						Ips_NePt.push_back(NULL);
					}
					else
					{
						Ips_NePt.push_back(s_tm->getmp()->at(0));
						for (int i=0;i<(int)s_tm->getmp()->size();i++)
						{
							s_tm->getmp()->at(i)->removemp(s_tm);
						}
						s_tm->clearmp(); // clear the neighbors
					}
					if (s_tm->getba()->get_attach())
					{
						double neegl = s_tm->getba()->getlength();if (neegl<m_l) {m_l = neegl;}
						Ips_NeEg.push_back(s_tm->getba()->get_Neighbor());
						Ips_Ip.push_back(s_tm->getba()->get_Pair());
					}
					else
					{
						Ips_NeEg.push_back(NULL);
						Ips_Ip.push_back(NULL);
					}
					s_tm = s_tm->getrp();
					sspn++;
				}
				/////////////// single point pair removal ////////////////
				int Ipsn_s = (int)Ips_single.size();
				for (int i=0;i<Ipsn_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ips_single[i]->get_I1()==C)
					{
						mergepair *Mp = Ips_single[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ips_single[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ips_single[i]->remove_MP(Mp);
					}
				}
				int Iptn_s = (int)Ipt_single.size();
				for (int i=0;i<Iptn_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (Ipt_single[i]->get_I1()==C)
					{
						mergepair *Mp = Ipt_single[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = Ipt_single[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						Ipt_single[i]->remove_MP(Mp);
					}
				}
				int Ipfs_s = (int)front_single_pairs.size();
				for (int i=0;i<Ipfs_s;i++)
				{
					point3D *sC = NULL;
					point3D *nC = NULL;
					if (front_single_pairs[i]->get_I1()==C)
					{
						mergepair *Mp = front_single_pairs[i]->MP()->at(0);
						sC = Mp->S1();
						nC = Mp->S2();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[i]->remove_MP(Mp);
					}
					else
					{
						mergepair *Mp = front_single_pairs[i]->MP()->at(0);
						sC = Mp->S2();
						nC = Mp->S1();
						sC->removemp(nC);
						nC->removemp(sC);
						node_remove_0(sC);
						front_single_pairs[i]->remove_MP(Mp);
					}
				}
				for (int i=tspn-1;i>=0;i--)
				{
					if (i==0)
					{
						int nd_size = (int)nodeList.size();
						node *A = new node(nd_size);
						nodeList.push_back(A);
						nodeList[nd_size]->pushvertex(Ipt_Pt[i]);
						nodeList[nd_size]->set_in_global(1);
						Ipt_Pt[i]->set_node_id(nodeList[nd_size]);
						Ipt_Pt[i]->getab()->set_attach(0);
					}
					else
					{
						/////////////// node //////////////
						Ipt_Nd[i-1]->pushvertex(Ipt_Pt[i]);
						Ipt_Pt[i]->set_node_id(Ipt_Nd[i-1]);
						double rpx = 0;
						double rpy = 0;
						if (Ipt_NePt[i-1]!=NULL)
						{
							rpx = Ipt_NePt[i-1]->x();
							rpy = Ipt_NePt[i-1]->y();
						}
						else
						{
							rpx = Ipt_Pt[i-1]->x();
							rpy = Ipt_Pt[i-1]->y();
						}
						Ipt_Pt[i]->setx(rpx);
						Ipt_Pt[i]->sety(rpy);
						if (Ipt_NePt[i-1]!=NULL)
						{
							vector<point3D*> rpmp;
							rpmp.push_back(Ipt_NePt[i-1]);
							Ipt_Pt[i]->pushmp(Ipt_NePt[i-1]);
							for (int j=0;j<(int)Ipt_NePt[i-1]->getmp()->size();j++)
							{
								rpmp.push_back(Ipt_NePt[i-1]->getmp()->at(j));
								Ipt_Pt[i]->pushmp(Ipt_NePt[i-1]->getmp()->at(j));
							}
							for (int j=0;j<(int)rpmp.size();j++)
							{
								rpmp[j]->pushmp(Ipt_Pt[i]);
								edge *rpmp_ab = rpmp[j]->getab();
								edge *rpmp_ba = rpmp[j]->getba();
								cell *Crp = rpmp_ab->Cell();
								cell *C1r = NULL;
								cell *C2r = NULL;
								if (rpmp_ab->get_attach())
								{
									if (rpmp_ab->get_Pair()->get_I1()==Crp) {C1r = rpmp_ab->get_Pair()->get_I2();}
									else {C1r = rpmp_ab->get_Pair()->get_I1();}
								}
								if (rpmp_ba->get_attach())
								{
									if (rpmp_ba->get_Pair()->get_I1()==Crp) {C2r = rpmp_ba->get_Pair()->get_I2();}
									else {C2r = rpmp_ba->get_Pair()->get_I1();}
								}
								if (C1r!=C && C2r!=C)
								{
									int crpn = (int)Crp->get_pairs()->size();
									for (int k=0;k<crpn;k++)
									{
										if ((Crp->get_pairs()->at(k)->get_I1()==Crp &&
											 Crp->get_pairs()->at(k)->get_I2()==C) ||
											(Crp->get_pairs()->at(k)->get_I2()==Crp &&
											 Crp->get_pairs()->at(k)->get_I1()==C))
										{
											if ((int)Crp->get_pairs()->at(k)->MP()->size()>0)
											{
												if      (Crp->get_pairs()->at(k)->MP()->at(0)->S1()==rpmp[j])
												{
													Crp->get_pairs()->at(k)->MP()->at(0)->setS2(Ipt_Pt[i]);
													Crp->get_pairs()->at(k)->MP()->at(0)->setT2(Ipt_Pt[i]);
												}
												else if (Crp->get_pairs()->at(k)->MP()->at(0)->S2()==rpmp[j])
												{
													Crp->get_pairs()->at(k)->MP()->at(0)->setS1(Ipt_Pt[i]);
													Crp->get_pairs()->at(k)->MP()->at(0)->setT1(Ipt_Pt[i]);
												}
											}
											break;
										}
									}
								}
							}
						}
						if (i<tspn-1)
						{
							//////////// edge ////////////
							if (Ipt_NeEg[i-1]==NULL)
							{
								Ipt_Pt[i]->getab()->set_attach(0);
							}
							else
							{
								Ipt_Pt[i]->getab()->set_attach(1);
								Ipt_NeEg[i-1]->set_attach(1);
								Ipt_Pt[i]->getab()->set_Neighbor(Ipt_NeEg[i-1]);
								Ipt_NeEg[i-1]->set_Neighbor(Ipt_Pt[i]->getab());
							}
							//////////// pair ////////////
							Ipt_Pt[i]->getab()->set_Pair(Ipt_Ip[i-1]);
						}
						else
						{
							//////////// edge ////////////
							Ipt_Pt[i]->getab()->set_attach(0);
							//////////// pair ////////////
							Ipt_Pt[i]->getab()->set_Pair(NULL);
						}
					}
				}
				for (int i=sspn-1;i>=0;i--)
				{
					if (i==0)
					{
						int nd_size = (int)nodeList.size();
						node *A = new node(nd_size);
						nodeList.push_back(A);
						nodeList[nd_size]->pushvertex(Ips_Pt[i]);
						nodeList[nd_size]->set_in_global(1);
						Ips_Pt[i]->set_node_id(nodeList[nd_size]);
						Ips_Pt[i]->getba()->set_attach(0);
					}
					else
					{
						/////////// node ////////////
						Ips_Nd[i-1]->pushvertex(Ips_Pt[i]);
						Ips_Pt[i]->set_node_id(Ips_Nd[i-1]);
						double fpx = 0;
						double fpy = 0;
						if (Ips_NePt[i-1]!=NULL)
						{
							fpx = Ips_NePt[i-1]->x();
							fpy = Ips_NePt[i-1]->y();
						}
						else
						{
							fpx = Ips_Pt[i-1]->x();
							fpy = Ips_Pt[i-1]->y();
						}
						Ips_Pt[i]->setx(fpx);
						Ips_Pt[i]->sety(fpy);
						if (Ips_NePt[i-1]!=NULL)
						{
							vector<point3D*> fpmp;
							fpmp.push_back(Ips_NePt[i-1]);
							Ips_Pt[i]->pushmp(Ips_NePt[i-1]);
							for (int j=0;j<(int)Ips_NePt[i-1]->getmp()->size();j++)
							{
								fpmp.push_back(Ips_NePt[i-1]->getmp()->at(j));
								Ips_Pt[i]->pushmp(Ips_NePt[i-1]->getmp()->at(j));
							}
							for (int j=0;j<(int)fpmp.size();j++)
							{
								fpmp[j]->pushmp(Ips_Pt[i]);
								edge *fpmp_ab = fpmp[j]->getab();
								edge *fpmp_ba = fpmp[j]->getba();
								cell *Cfp = fpmp_ab->Cell();
								cell *C1f = NULL;
								cell *C2f = NULL;
								if (fpmp_ab->get_attach())
								{
									if (fpmp_ab->get_Pair()->get_I1()==Cfp) {C1f = fpmp_ab->get_Pair()->get_I2();}
									else {C1f = fpmp_ab->get_Pair()->get_I1();}
								}
								if (fpmp_ba->get_attach())
								{
									if (fpmp_ba->get_Pair()->get_I1()==Cfp) {C2f = fpmp_ba->get_Pair()->get_I2();}
									else {C2f = fpmp_ba->get_Pair()->get_I1();}
								}
								if (C1f!=C && C2f!=C)
								{
									int cfpn = (int)Cfp->get_pairs()->size();
									for (int k=0;k<cfpn;k++)
									{
										if ((Cfp->get_pairs()->at(k)->get_I1()==Cfp &&
											 Cfp->get_pairs()->at(k)->get_I2()==C) ||
											(Cfp->get_pairs()->at(k)->get_I2()==Cfp &&
											 Cfp->get_pairs()->at(k)->get_I1()==C))
										{
											if ((int)Cfp->get_pairs()->at(k)->MP()->size()>0)
											{
												if      (Cfp->get_pairs()->at(k)->MP()->at(0)->S1()==fpmp[j])
												{
													Cfp->get_pairs()->at(k)->MP()->at(0)->setS2(Ips_Pt[i]);
													Cfp->get_pairs()->at(k)->MP()->at(0)->setT2(Ips_Pt[i]);
												}
												else if (Cfp->get_pairs()->at(k)->MP()->at(0)->S2()==fpmp[j])
												{
													Cfp->get_pairs()->at(k)->MP()->at(0)->setS1(Ips_Pt[i]);
													Cfp->get_pairs()->at(k)->MP()->at(0)->setT1(Ips_Pt[i]);
												}
											}
											break;
										}
									}
								}
							}
						}
						if (i<sspn-1)
						{
							///////// edge ///////////
							if (Ips_NeEg[i-1]==NULL)
							{
								Ips_Pt[i]->getba()->set_attach(0);
							}
							else
							{
								Ips_Pt[i]->getba()->set_attach(1);
								Ips_NeEg[i-1]->set_attach(1);
								Ips_Pt[i]->getba()->set_Neighbor(Ips_NeEg[i-1]);
								Ips_NeEg[i-1]->set_Neighbor(Ips_Pt[i]->getba());
							}
							////////// pair //////////
							Ips_Pt[i]->getba()->set_Pair(Ips_Ip[i-1]);
						}
						else
						{
							///////// edge ///////////
							Ips_Pt[i]->getba()->set_attach(0);
							////////// pair //////////
							Ips_Pt[i]->getba()->set_Pair(NULL);
						}
					}
				}
				/////// points behind t_t and s_s ///////
				if      (t_t1_angle>96 && s_s1_angle>96)
				{
					vector<point3D*> tmpt;
					tmpt.push_back(Ipt_NePt[tspn-1]);
					for (int j=0;j<(int)Ipt_NePt[tspn-1]->getmp()->size();j++)
					{tmpt.push_back(Ipt_NePt[tspn-1]->getmp()->at(j));}
					for (int j=0;j<(int)tmpt.size();j++)
					{
						point3D *mppt = tmpt[j];
						edge *mpab = mppt->getab();
						edge *mpba = mppt->getba();
						cell *Cmp = mpab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (mpab->get_attach())
						{
							if (mpab->get_Pair()->get_I1()==Cmp) {C1m = mpab->get_Pair()->get_I2();}
							else {C1m = mpab->get_Pair()->get_I1();}
						}
						if (mpba->get_attach())
						{
							if (mpba->get_Pair()->get_I1()==Cmp) {C2m = mpba->get_Pair()->get_I2();}
							else {C2m = mpba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Cmp->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Cmp->get_pairs()->at(k)->get_I1()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I2()==C) ||
									(Cmp->get_pairs()->at(k)->get_I2()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I1()==C))
								{
									mergepair *Mp = Cmp->get_pairs()->at(k)->MP()->at(0);
									Cmp->get_pairs()->at(k)->remove_MP(Mp);
									break;
								}
							}
						}
					}
					vector<point3D*> smpt;
					smpt.push_back(Ips_NePt[sspn-1]);
					for (int j=0;j<(int)Ips_NePt[sspn-1]->getmp()->size();j++)
					{smpt.push_back(Ips_NePt[sspn-1]->getmp()->at(j));}
					for (int j=0;j<(int)smpt.size();j++)
					{
						point3D *mppt = smpt[j];
						edge *mpab = mppt->getab();
						edge *mpba = mppt->getba();
						cell *Cmp = mpab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (mpab->get_attach())
						{
							if (mpab->get_Pair()->get_I1()==Cmp) {C1m = mpab->get_Pair()->get_I2();}
							else {C1m = mpab->get_Pair()->get_I1();}
						}
						if (mpba->get_attach())
						{
							if (mpba->get_Pair()->get_I1()==Cmp) {C2m = mpba->get_Pair()->get_I2();}
							else {C2m = mpba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Cmp->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Cmp->get_pairs()->at(k)->get_I1()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I2()==C) ||
									(Cmp->get_pairs()->at(k)->get_I2()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I1()==C))
								{
									mergepair *Mp = Cmp->get_pairs()->at(k)->MP()->at(0);
									Cmp->get_pairs()->at(k)->remove_MP(Mp);
									break;
								}
							}
						}
					}
					double cx = (t_t1->x() + s_s1->x())*0.5;
					double cy = (t_t1->y() + s_s1->y())*0.5;
					Ipt_NeEg[tspn-2]->set_attach(0);
					Ips_NeEg[sspn-2]->set_attach(0);
					if (t_t1->getfp()!=s_s1)
					{
						point3D *ptp = t_t1->getfp();
						while (ptp!=s_s1)
						{
							double px = ptp->x();
							double py = ptp->y();
							double dpx = cx - px;
							double dpy = cy - py;
							double dpd = sqrt(dpx*dpx + dpy*dpy);
							double dpmx = dpx/dpd*m_l*0.5;
							double dpmy = dpy/dpd*m_l*0.5;
							ptp->setx(px + dpmx);
							ptp->sety(py + dpmy);
							ptp = ptp->getfp();
						}
					}
				}
				else if (t_t1_angle<=96 && s_s1_angle<=96) // t_t1->fp and s_s1->rp 
				{
					double cx = (Ipt_NePt[tspn-1]->x() + Ips_NePt[sspn-1]->x())*0.5;
					double cy = (Ipt_NePt[tspn-1]->y() + Ips_NePt[sspn-1]->y())*0.5;
					double tx = Ipt_NePt[tspn-1]->x();
					double ty = Ipt_NePt[tspn-1]->y();
					double sx = Ips_NePt[sspn-1]->x();
					double sy = Ips_NePt[sspn-1]->y();
					t_t1->getfp()->setx(tx);
					t_t1->getfp()->sety(ty);
					s_s1->getrp()->setx(sx);
					s_s1->getrp()->sety(sy);
					t_t1->getab()->set_attach(1);
					t_t1->getab()->set_Neighbor(Ipt_NeEg[tspn-2]);
					Ipt_NeEg[tspn-2]->set_Neighbor(t_t1->getab());
					t_t1->getab()->set_Pair(Ipt_NeEg[tspn-2]->get_Pair());
					s_s1->getba()->set_attach(1);
					s_s1->getba()->set_Neighbor(Ips_NeEg[sspn-2]);
					Ips_NeEg[sspn-2]->set_Neighbor(s_s1->getba());
					s_s1->getba()->set_Pair(Ips_NeEg[sspn-2]->get_Pair());
					node_merge(Ipt_NePt[tspn-1],t_t1->getfp());
					node_merge(Ips_NePt[sspn-1],s_s1->getrp());
					vector<point3D*> t1mp;
					t1mp.push_back(Ipt_NePt[tspn-1]);
					t_t1->getfp()->pushmp(Ipt_NePt[tspn-1]);
					for (int j=0;j<(int)Ipt_NePt[tspn-1]->getmp()->size();j++)
					{
						t1mp.push_back(Ipt_NePt[tspn-1]->getmp()->at(j));
						t_t1->getfp()->pushmp(Ipt_NePt[tspn-1]->getmp()->at(j));
					}
					for (int j=0;j<(int)t1mp.size();j++)
					{
						t1mp[j]->pushmp(t_t1->getfp());
						edge *tmab = t1mp[j]->getab();
						edge *tmba = t1mp[j]->getba();
						cell *Ctm = tmab->Cell();
						cell *C1t = NULL;
						cell *C2t = NULL;
						if (tmab->get_attach())
						{
							if (tmab->get_Pair()->get_I1()==Ctm) {C1t = tmab->get_Pair()->get_I2();}
							else {C1t = tmab->get_Pair()->get_I1();}
						}
						if (tmba->get_attach())
						{
							if (tmba->get_Pair()->get_I1()==Ctm) {C2t = tmba->get_Pair()->get_I2();}
							else {C2t = tmba->get_Pair()->get_I1();}
						}
						if (C1t!=C && C2t!=C)
						{
							int cmpn = (int)Ctm->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Ctm->get_pairs()->at(k)->get_I1()==Ctm &&
									 Ctm->get_pairs()->at(k)->get_I2()==C) ||
									(Ctm->get_pairs()->at(k)->get_I2()==Ctm &&
									 Ctm->get_pairs()->at(k)->get_I1()==C))
								{
									if      (Ctm->get_pairs()->at(k)->MP()->at(0)->S1()==t1mp[j])
									{
										Ctm->get_pairs()->at(k)->MP()->at(0)->setS2(t_t1->getfp());
										Ctm->get_pairs()->at(k)->MP()->at(0)->setT2(t_t1->getfp());
									}
									else if (Ctm->get_pairs()->at(k)->MP()->at(0)->S2()==t1mp[j])
									{
										Ctm->get_pairs()->at(k)->MP()->at(0)->setS1(t_t1->getfp());
										Ctm->get_pairs()->at(k)->MP()->at(0)->setT1(t_t1->getfp());
									}
									break;
								}
							}
						}
					}
					vector<point3D*> s1mp;
					s1mp.push_back(Ips_NePt[sspn-1]);
					s_s1->getrp()->pushmp(Ips_NePt[sspn-1]);
					for (int j=0;j<(int)Ips_NePt[sspn-1]->getmp()->size();j++)
					{
						s1mp.push_back(Ips_NePt[sspn-1]->getmp()->at(j));
						s_s1->getrp()->pushmp(Ips_NePt[sspn-1]->getmp()->at(j));
					}
					for (int j=0;j<(int)s1mp.size();j++)
					{
						s1mp[j]->pushmp(s_s1->getrp());
						edge *smab = s1mp[j]->getab();
						edge *smba = s1mp[j]->getba();
						cell *Csm = smab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (smab->get_attach())
						{
							if (smab->get_Pair()->get_I1()==Csm) {C1m = smab->get_Pair()->get_I2();}
							else {C1m = smab->get_Pair()->get_I1();}
						}
						if (smba->get_attach())
						{
							if (smba->get_Pair()->get_I1()==Csm) {C2m = smba->get_Pair()->get_I2();}
							else {C2m = smba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Csm->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Csm->get_pairs()->at(k)->get_I1()==Csm &&
									 Csm->get_pairs()->at(k)->get_I2()==C) ||
									(Csm->get_pairs()->at(k)->get_I2()==Csm &&
									 Csm->get_pairs()->at(k)->get_I1()==C))
								{
									if      (Csm->get_pairs()->at(k)->MP()->at(0)->S1()==s1mp[j])
									{
										Csm->get_pairs()->at(k)->MP()->at(0)->setS2(s_s1->getrp());
										Csm->get_pairs()->at(k)->MP()->at(0)->setT2(s_s1->getrp());
									}
									else if (Csm->get_pairs()->at(k)->MP()->at(0)->S2()==s1mp[j])
									{
										Csm->get_pairs()->at(k)->MP()->at(0)->setS1(s_s1->getrp());
										Csm->get_pairs()->at(k)->MP()->at(0)->setT1(s_s1->getrp());
									}
									break;
								}
							}
						}
					}
					if (t_t1->getfp()->getfp()!=s_s1->getrp())
					{
						point3D *ptp = t_t1->getfp()->getfp();
						while (ptp!=s_s1->getrp())
						{
							double px = ptp->x();
							double py = ptp->y();
							double dpx = cx - px;
							double dpy = cy - py;
							double dpd = sqrt(dpx*dpx + dpy*dpy);
							double dpmx = dpx/dpd*m_l*0.5;
							double dpmy = dpy/dpd*m_l*0.5;
							ptp->setx(px + dpmx);
							ptp->sety(py + dpmy);
							ptp = ptp->getfp();
						}
					}
				}
				else if (t_t1_angle>96 && s_s1_angle<=96) // s_s1->rp
				{
					vector<point3D*> tmpt;
					tmpt.push_back(Ipt_NePt[tspn-1]);
					for (int j=0;j<(int)Ipt_NePt[tspn-1]->getmp()->size();j++)
					{tmpt.push_back(Ipt_NePt[tspn-1]->getmp()->at(j));}
					for (int j=0;j<(int)tmpt.size();j++)
					{
						point3D *mppt = tmpt[j];
						edge *mpab = mppt->getab();
						edge *mpba = mppt->getba();
						cell *Cmp = mpab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (mpab->get_attach())
						{
							if (mpab->get_Pair()->get_I1()==Cmp) {C1m = mpab->get_Pair()->get_I2();}
							else {C1m = mpab->get_Pair()->get_I1();}
						}
						if (mpba->get_attach())
						{
							if (mpba->get_Pair()->get_I1()==Cmp) {C2m = mpba->get_Pair()->get_I2();}
							else {C2m = mpba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Cmp->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Cmp->get_pairs()->at(k)->get_I1()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I2()==C) ||
									(Cmp->get_pairs()->at(k)->get_I2()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I1()==C))
								{
									mergepair *Mp = Cmp->get_pairs()->at(k)->MP()->at(0);
									Cmp->get_pairs()->at(k)->remove_MP(Mp);
									break;
								}
							}
						}
					}
					double cx = (t_t1->x() + Ips_NePt[sspn-1]->x())*0.5;
					double cy = (t_t1->y() + Ips_NePt[sspn-1]->y())*0.5;
					double sx = Ips_NePt[sspn-1]->x();
					double sy = Ips_NePt[sspn-1]->y();
					Ipt_NeEg[tspn-2]->set_attach(0);
					s_s1->getrp()->setx(sx);
					s_s1->getrp()->sety(sy);
					s_s1->getba()->set_attach(1);
					s_s1->getba()->set_Neighbor(Ips_NeEg[sspn-2]);
					Ips_NeEg[sspn-2]->set_Neighbor(s_s1->getba());
					s_s1->getba()->set_Pair(Ips_NeEg[sspn-2]->get_Pair());
					node_merge(Ips_NePt[sspn-1],s_s1->getrp());
					vector<point3D*> s1mp;
					s1mp.push_back(Ips_NePt[sspn-1]);
					s_s1->getrp()->pushmp(Ips_NePt[sspn-1]);
					for (int j=0;j<(int)Ips_NePt[sspn-1]->getmp()->size();j++)
					{
						s1mp.push_back(Ips_NePt[sspn-1]->getmp()->at(j));
						s_s1->getrp()->pushmp(Ips_NePt[sspn-1]->getmp()->at(j));
					}
					for (int j=0;j<(int)s1mp.size();j++)
					{
						s1mp[j]->pushmp(s_s1->getrp());
						edge *smab = s1mp[j]->getab();
						edge *smba = s1mp[j]->getba();
						cell *Csm = smab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (smab->get_attach())
						{
							if (smab->get_Pair()->get_I1()==Csm) {C1m = smab->get_Pair()->get_I2();}
							else {C1m = smab->get_Pair()->get_I1();}
						}
						if (smba->get_attach())
						{
							if (smba->get_Pair()->get_I1()==Csm) {C2m = smba->get_Pair()->get_I2();}
							else {C2m = smba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Csm->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Csm->get_pairs()->at(k)->get_I1()==Csm &&
									 Csm->get_pairs()->at(k)->get_I2()==C) ||
									(Csm->get_pairs()->at(k)->get_I2()==Csm &&
									 Csm->get_pairs()->at(k)->get_I1()==C))
								{
									if      (Csm->get_pairs()->at(k)->MP()->at(0)->S1()==s1mp[j])
									{
										Csm->get_pairs()->at(k)->MP()->at(0)->setS2(s_s1->getrp());
										Csm->get_pairs()->at(k)->MP()->at(0)->setT2(s_s1->getrp());
									}
									else if (Csm->get_pairs()->at(k)->MP()->at(0)->S2()==s1mp[j])
									{
										Csm->get_pairs()->at(k)->MP()->at(0)->setS1(s_s1->getrp());
										Csm->get_pairs()->at(k)->MP()->at(0)->setT1(s_s1->getrp());
									}
									break;
								}
							}
						}
					}
					if (t_t1->getfp()!=s_s1->getrp())
					{
						point3D *ptp = t_t1->getfp();
						while (ptp!=s_s1->getrp())
						{
							double px = ptp->x();
							double py = ptp->y();
							double dpx = cx - px;
							double dpy = cy - py;
							double dpd = sqrt(dpx*dpx + dpy*dpy);
							double dpmx = dpx/dpd*m_l*0.5;
							double dpmy = dpy/dpd*m_l*0.5;
							ptp->setx(px + dpmx);
							ptp->sety(py + dpmy);
							ptp = ptp->getfp();
						}
					}
				}
				else if (t_t1_angle<=96 && s_s1_angle>96) // t_t1->fp
				{
					vector<point3D*> smpt;
					smpt.push_back(Ips_NePt[sspn-1]);
					for (int j=0;j<(int)Ips_NePt[sspn-1]->getmp()->size();j++)
					{smpt.push_back(Ips_NePt[sspn-1]->getmp()->at(j));}
					for (int j=0;j<(int)smpt.size();j++)
					{
						point3D *mppt = smpt[j];
						edge *mpab = mppt->getab();
						edge *mpba = mppt->getba();
						cell *Cmp = mpab->Cell();
						cell *C1m = NULL;
						cell *C2m = NULL;
						if (mpab->get_attach())
						{
							if (mpab->get_Pair()->get_I1()==Cmp) {C1m = mpab->get_Pair()->get_I2();}
							else {C1m = mpab->get_Pair()->get_I1();}
						}
						if (mpba->get_attach())
						{
							if (mpba->get_Pair()->get_I1()==Cmp) {C2m = mpba->get_Pair()->get_I2();}
							else {C2m = mpba->get_Pair()->get_I1();}
						}
						if (C1m!=C && C2m!=C)
						{
							int cmpn = (int)Cmp->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Cmp->get_pairs()->at(k)->get_I1()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I2()==C) ||
									(Cmp->get_pairs()->at(k)->get_I2()==Cmp &&
									 Cmp->get_pairs()->at(k)->get_I1()==C))
								{
									mergepair *Mp = Cmp->get_pairs()->at(k)->MP()->at(0);
									Cmp->get_pairs()->at(k)->remove_MP(Mp);
									break;
								}
							}
						}
					}
					double cx = (Ipt_NePt[tspn-1]->x() + s_s1->x())*0.5;
					double cy = (Ipt_NePt[tspn-1]->y() + s_s1->y())*0.5;
					Ips_NeEg[sspn-2]->set_attach(0);
					double tx = Ipt_NePt[tspn-1]->x();
					double ty = Ipt_NePt[tspn-1]->y();
					t_t1->getfp()->setx(tx);
					t_t1->getfp()->sety(ty);
					t_t1->getab()->set_attach(1);
					t_t1->getab()->set_Neighbor(Ipt_NeEg[tspn-2]);
					Ipt_NeEg[tspn-2]->set_Neighbor(t_t1->getab());
					t_t1->getab()->set_Pair(Ipt_NeEg[tspn-2]->get_Pair());
					node_merge(Ipt_NePt[tspn-1],t_t1->getfp());
					vector<point3D*> t1mp;
					t1mp.push_back(Ipt_NePt[tspn-1]);
					t_t1->getfp()->pushmp(Ipt_NePt[tspn-1]);
					for (int j=0;j<(int)Ipt_NePt[tspn-1]->getmp()->size();j++)
					{
						t1mp.push_back(Ipt_NePt[tspn-1]->getmp()->at(j));
						t_t1->getfp()->pushmp(Ipt_NePt[tspn-1]->getmp()->at(j));
					}
					for (int j=0;j<(int)t1mp.size();j++)
					{
						t1mp[j]->pushmp(t_t1->getfp());
						edge *tmab = t1mp[j]->getab();
						edge *tmba = t1mp[j]->getba();
						cell *Ctm = tmab->Cell();
						cell *C1t = NULL;
						cell *C2t = NULL;
						if (tmab->get_attach())
						{
							if (tmab->get_Pair()->get_I1()==Ctm) {C1t = tmab->get_Pair()->get_I2();}
							else {C1t = tmab->get_Pair()->get_I1();}
						}
						if (tmba->get_attach())
						{
							if (tmba->get_Pair()->get_I1()==Ctm) {C2t = tmba->get_Pair()->get_I2();}
							else {C2t = tmba->get_Pair()->get_I1();}
						}
						if (C1t!=C && C2t!=C)
						{
							int cmpn = (int)Ctm->get_pairs()->size();
							for (int k=0;k<cmpn;k++)
							{
								if ((Ctm->get_pairs()->at(k)->get_I1()==Ctm &&
									 Ctm->get_pairs()->at(k)->get_I2()==C) ||
									(Ctm->get_pairs()->at(k)->get_I2()==Ctm &&
									 Ctm->get_pairs()->at(k)->get_I1()==C))
								{
									if      (Ctm->get_pairs()->at(k)->MP()->at(0)->S1()==t1mp[j])
									{
										Ctm->get_pairs()->at(k)->MP()->at(0)->setS2(t_t1->getfp());
										Ctm->get_pairs()->at(k)->MP()->at(0)->setT2(t_t1->getfp());
									}
									else if (Ctm->get_pairs()->at(k)->MP()->at(0)->S2()==t1mp[j])
									{
										Ctm->get_pairs()->at(k)->MP()->at(0)->setS1(t_t1->getfp());
										Ctm->get_pairs()->at(k)->MP()->at(0)->setT1(t_t1->getfp());
									}
									break;
								}
							}
						}
					}
					if (t_t1->getfp()->getfp()!=s_s1)
					{
						point3D *ptp = t_t1->getfp()->getfp();
						while (ptp!=s_s1)
						{
							double px = ptp->x();
							double py = ptp->y();
							double dpx = cx - px;
							double dpy = cy - py;
							double dpd = sqrt(dpx*dpx + dpy*dpy);
							double dpmx = dpx/dpd*m_l*0.5;
							double dpmy = dpy/dpd*m_l*0.5;
							ptp->setx(px + dpmx);
							ptp->sety(py + dpmy);
							ptp = ptp->getfp();
						}
					}
				}
				////////// points before t_s1 and s_t1 ///////////
				point3D *ptp = s_t1;
				double dpmx = dx*m_l;
				double dpmy = dy*m_l;
				while (ptp!=t_s1->getfp())
				{
					double px = ptp->x();
					double py = ptp->y();
					if ((int)ptp->getmp()->size()==0)
					{
						ptp->setx(px + dpmx*1.0);
						ptp->sety(py + dpmy*1.0);
					}
					else
					{
						ptp->setx(px + dpmx*0.5);
						ptp->sety(py + dpmy*0.5);
						for (int j=0;j<(int)ptp->getmp()->size();j++)
						{
							ptp->getmp()->at(j)->setx(px + dpmx*0.5);
							ptp->getmp()->at(j)->sety(py + dpmy*0.5);
						}
					}
					ptp = ptp->getfp();
				}
				////////// interpair end-points update //////////
				for (int i=0;i<Iptn;i++)
				{
					if (Ipt[i]->get_I1()==C)
					{
						if (i==Iptn-1)
						{
							point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
							point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
							Ipt[i]->MP()->at(0)->setS1(ps1->getfp());
							if (t_t1->getab()->get_attach())
							{
								Ipt[i]->MP()->at(0)->setT1(t_t1->getfp());
							}
							else
							{
								Ipt[i]->MP()->at(0)->setS2(ps2->getfp());
							}
							ps2 = Ipt[i]->MP()->at(0)->S2();
							point3D *pt2 = Ipt[i]->MP()->at(0)->T2();
							if (ps2!=pt2)
							{
								if (pt2->getrp()==ps2)
								{
									double Fr_m = Ipt[i]->MP()->at(0)->get_Fr_t1();
									Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_m);
								}
								else
								{
									double Fr_m = pt2->getrp()->get_stored_friction_force();
									Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_m);
									point3D *mpt = pt2->getrp();
									while (mpt!=ps2)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getrp();
									}
								}
							}
						}
						else
						{
							point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
							point3D *pt1 = Ipt[i]->MP()->at(0)->T1();
							Ipt[i]->MP()->at(0)->setS1(ps1->getfp());
							Ipt[i]->MP()->at(0)->setT1(pt1->getfp());
							if (ps1!=pt1)
							{
								if (ps1->getfp()==pt1)
								{
									double Fr_t1 = Ipt[i]->MP()->at(0)->get_Fr_t1();
									Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_t1);
								}
								else
								{
									point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
									point3D *pt2 = Ipt[i]->MP()->at(0)->T2();
									double Fr_mt2 = pt2->getrp()->get_stored_friction_force();
									double Fr_ms1 = pt2->getrp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m = Fr_mt2;if (Fr_mt2<Fr_ms1) {Fr_m = Fr_ms1;}
									Ipt[i]->MP()->at(0)->set_Fr_s1(Fr_m);
									point3D *mpt = pt2->getrp();
									while (mpt!=ps2)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getrp();
									}
								}
							}
						}
					}
					else
					{
						if (i==Iptn-1)
						{
							point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
							point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
							Ipt[i]->MP()->at(0)->setS2(ps2->getfp());
							if (t_t1->getab()->get_attach())
							{
								Ipt[i]->MP()->at(0)->setT2(t_t1->getfp());
							}
							else
							{
								Ipt[i]->MP()->at(0)->setS1(ps1->getfp());
							}
							ps1 = Ipt[i]->MP()->at(0)->S1();
							point3D *pt1 = Ipt[i]->MP()->at(0)->T1();
							if (ps1!=pt1)
							{
								if (pt1->getrp()==ps1)
								{
									double Fr_m = Ipt[i]->MP()->at(0)->get_Fr_s1();
									Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
								}
								else
								{
									double Fr_m = pt1->getrp()->get_stored_friction_force();
									Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
									point3D *mpt = pt1->getrp();
									while (mpt!=ps1)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getrp();
									}
								}
							}
						}
						else
						{
							point3D *ps2 = Ipt[i]->MP()->at(0)->S2();
							point3D *pt2 = Ipt[i]->MP()->at(0)->T2();
							Ipt[i]->MP()->at(0)->setS2(ps2->getfp());
							Ipt[i]->MP()->at(0)->setT2(pt2->getfp());
							point3D *ps1 = Ipt[i]->MP()->at(0)->S1();
							point3D *pt1 = Ipt[i]->MP()->at(0)->T1();
							if (ps1!=pt1)
							{
								if (ps1->getfp()==pt1)
								{
									double Fr_m = Ipt[i]->MP()->at(0)->get_Fr_s1();
									Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
								}
								else
								{
									double Fr_m1 = pt1->getrp()->get_stored_friction_force();
									double Fr_m2 = pt1->getrp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m = Fr_m2;}
									Ipt[i]->MP()->at(0)->set_Fr_t1(Fr_m);
									point3D *mpt = pt1->getrp();
									while (mpt!=ps1)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getrp();
									}
								}
							}
						}
					}
				}
				for (int i=0;i<Ipsn;i++)
				{
					if (Ips[i]->get_I1()==C)
					{
						if (i==Ipsn-1)
						{
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							Ips[i]->MP()->at(0)->setT1(pt1->getrp());
							if (s_s1->getba()->get_attach())
							{
								Ips[i]->MP()->at(0)->setS1(s_s1->getrp());
							}
							else
							{
								Ips[i]->MP()->at(0)->setT2(pt2->getrp());
							}
							pt2 = Ips[i]->MP()->at(0)->T2();
							point3D *ps2 = Ips[i]->MP()->at(0)->S2();
							if (ps2!=pt2)
							{
								if (ps2->getfp()==pt2)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_s1();
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps2->getfp()->get_stored_friction_force();
									double Fr_m2 = ps2->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m0);
									point3D *mpt = ps2->getfp();
									while (mpt!=pt2)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
						else
						{
							point3D *ps1 = Ips[i]->MP()->at(0)->S1();
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							Ips[i]->MP()->at(0)->setS1(ps1->getrp());
							Ips[i]->MP()->at(0)->setT1(pt1->getrp());
							point3D *ps2 = Ips[i]->MP()->at(0)->S2();
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							if (ps2!=pt2)
							{
								if (ps2->getfp()==pt2)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_s1();
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps2->getfp()->get_stored_friction_force();
									double Fr_m2 = ps2->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_t1(Fr_m0);
									point3D *mpt = ps2->getfp();
									while (mpt!=pt2)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
					}
					else
					{
						if (i==Ipsn-1)
						{
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							Ips[i]->MP()->at(0)->setT2(pt2->getrp());
							if (s_s1->getba()->get_attach())
							{
								Ips[i]->MP()->at(0)->setS2(s_s1->getrp());
							}
							else
							{
								Ips[i]->MP()->at(0)->setT1(pt1->getrp());
							}
							pt1 = Ips[i]->MP()->at(0)->T1();
							point3D *ps1 = Ips[i]->MP()->at(0)->S1();
							if (ps1!=pt1)
							{
								if (ps1->getfp()==pt1)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_t1();
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps1->getfp()->get_stored_friction_force();
									double Fr_m2 = ps1->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m0);
									point3D *mpt = ps1->getfp();
									while (mpt!=pt1)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
						else
						{
							point3D *ps2 = Ips[i]->MP()->at(0)->S2();
							point3D *pt2 = Ips[i]->MP()->at(0)->T2();
							Ips[i]->MP()->at(0)->setS2(ps2->getrp());
							Ips[i]->MP()->at(0)->setT2(pt2->getrp());
							point3D *ps1 = Ips[i]->MP()->at(0)->S1();
							point3D *pt1 = Ips[i]->MP()->at(0)->T1();
							if (ps1!=pt1)
							{
								if (ps1->getfp()==pt1)
								{
									double Fr_m = Ips[i]->MP()->at(0)->get_Fr_t1();
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m);
								}
								else
								{
									double Fr_m1 = ps1->getfp()->get_stored_friction_force();
									double Fr_m2 = ps1->getfp()->getmp()->at(0)->get_stored_friction_force();
									double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
									Ips[i]->MP()->at(0)->set_Fr_s1(Fr_m0);
									point3D *mpt = ps1->getfp();
									while (mpt!=pt1)
									{
										point3D *mpmt = mpt->getmp()->at(0);
										double Fr_m1 = mpt->get_stored_friction_force();
										double Fr_m2 = mpmt->get_stored_friction_force();
										double Fr_m0 = Fr_m1;if (Fr_m1<Fr_m2) {Fr_m0 = Fr_m2;}
										mpt->set_stored_friction_force(Fr_m0);
										mpmt->set_stored_friction_force(Fr_m0);
										mpt = mpt->getfp();
									}
								}
							}
						}
					}
				}
			}
		}
	}
	/*if (time==28 && C->id()==957)
	{
		cout<<"Pair 35851: "<<collisionpairList[35851]->MP()->size()<<endl;
		if (collisionpairList[35851]->MP()->size()>0)
		{
			cout<<"S1: "<<collisionpairList[35851]->MP()->at(0)->S1()->id()<<endl;
			cout<<"T1: "<<collisionpairList[35851]->MP()->at(0)->T1()->id()<<endl;
			cout<<"S2: "<<collisionpairList[35851]->MP()->at(0)->S2()->id()<<endl;
			cout<<"T2: "<<collisionpairList[35851]->MP()->at(0)->T2()->id()<<endl;
			point3D *s1tm = collisionpairList[35851]->MP()->at(0)->S1();
			point3D *t1tm = collisionpairList[35851]->MP()->at(0)->T1();
			int stn = 0;
			while (s1tm!=t1tm)
			{
				stn++;
				cout<<"p id:"<<s1tm->id()<<","<<s1tm->getmp()->at(0)->id()<<endl;
				cout<<"p x: "<<s1tm->x()<<","<<s1tm->getmp()->at(0)->x()<<endl;
				cout<<"p mn: "<<s1tm->getmp()->size()<<","<<s1tm->getmp()->at(0)->getmp()->size()<<endl;
				s1tm = s1tm->getfp();
			}
			cout<<"p id:"<<s1tm->id()<<","<<s1tm->getmp()->at(0)->id()<<endl;
			cout<<"p x: "<<s1tm->x()<<","<<s1tm->getmp()->at(0)->x()<<endl;
			cout<<"p mn: "<<s1tm->getmp()->size()<<","<<s1tm->getmp()->at(0)->getmp()->size()<<endl;
			cout<<"p num: "<<stn<<endl;
		}
	}*/
}

void dbReader::cell_dynamic(int number)
{
}

int dbReader::cell_dynamic_FEM(int number)
{
	int w_id1 = -1;
	int cn = (int)cellList.size();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		cellList[i]->center_refresh();
		/////// update the area info ////////
		for (int j=9;j>=1;j--)
		{
			cellList[i]->set_area_p(j,cellList[i]->get_area_p(j-1));
		}
		double area0 = cellList[i]->get_area();
		cellList[i]->set_area_p(0,area0);
		cellList[i]->set_area();
		////// kill cells which are too small //////
		area0 = cellList[i]->get_area();
		if (area0<AR*0.1) {cellList[i]->set_set_dead(1);}
	}
	/**************************************
			Cell death response
	**************************************/
	time_t t_b_dead, t_e_dead;
	t_b_dead = clock();
	cell_death_correction();
	t_e_dead = clock();
	cout<<"cell death done: "<<(t_e_dead-t_b_dead)*0.001<<" seconds!\n";

	double default_area = 6*GAR/6;
	cn = (int)cellList.size();
	/*************************************
	 The way to store matrix: COO format
	 (Coordinate list)
	 IA: row index
	 JA: column index
	 VA:  value
	 ----------------
	 | A11 A12 |
	 | A21 A22 |
	*************************************/
	vector<double> F; // the force-x-y-coordinate on node i
	vector<double> F_p; // the force-x-y-coordinate on node i: without area refinement force
	vector<le_sm*> sms;

	//int F_mode = 0; // 0: second newton mode; 1: elastic mode;
	// build up triangle stiffness matrix first //
	// elastic mode: estimate the force magnitude //
	/*if (F_mode==1)
	{
		for (int i=0;i<cn;i++)
		{
			int smi_id = 0;
			int sn = (int)cellList[i]->get_sides()->size();
			for (int j=0;j<sn;j++)
			{
				edge *E = cellList[i]->get_sides()->at(j);
				E->set_virtual_stiffness_matrix_value_only();
				E->p1()->set_SM_ID_single(smi_id);
				E->p1()->clear_node_sm_single();
				E->p1()->setA_ij_single(0,0,0,0);
				smi_id++;
			}
			int ipn = (int)cellList[i]->get_inner_p()->size();
			for (int j=0;j<ipn;j++)
			{
				cellList[i]->get_inner_p()->at(j)->set_SM_ID_single(smi_id);
				cellList[i]->get_inner_p()->at(j)->clear_node_sm_single();
				cellList[i]->get_inner_p()->at(j)->setA_ij_single(0,0,0,0);
				smi_id++;
			}
			int tn = (int)cellList[i]->get_inner_t()->size();
			vector<le_sm*> sms_single;
			double *SOR_II_V_single = new double[smi_id*2];
			int *SOR_N_single = new int[smi_id*2];
			int *SOR_Rn_single = new int[smi_id*2];

			for (int j=0;j<tn;j++)
			{
				triangle *T1 = cellList[i]->get_inner_t()->at(j);
				T1->set_Stiffness_matrix();
				point3D* TA = T1->getA();
				point3D* TB = T1->getB();
				point3D* TC = T1->getC();
				double v00 = T1->get_Stiffness_matrix(0,0);
				double v01 = T1->get_Stiffness_matrix(0,1);
				double v02 = T1->get_Stiffness_matrix(0,2);
				double v03 = T1->get_Stiffness_matrix(0,3);
				double v04 = T1->get_Stiffness_matrix(0,4);
				double v05 = T1->get_Stiffness_matrix(0,5);
				double v11 = T1->get_Stiffness_matrix(1,1);
				double v12 = T1->get_Stiffness_matrix(1,2);
				double v13 = T1->get_Stiffness_matrix(1,3);
				double v14 = T1->get_Stiffness_matrix(1,4);
				double v15 = T1->get_Stiffness_matrix(1,5);
				double v22 = T1->get_Stiffness_matrix(2,2);
				double v23 = T1->get_Stiffness_matrix(2,3);
				double v24 = T1->get_Stiffness_matrix(2,4);
				double v25 = T1->get_Stiffness_matrix(2,5);
				double v33 = T1->get_Stiffness_matrix(3,3);
				double v34 = T1->get_Stiffness_matrix(3,4);
				double v35 = T1->get_Stiffness_matrix(3,5);
				double v44 = T1->get_Stiffness_matrix(4,4);
				double v45 = T1->get_Stiffness_matrix(4,5);
				double v55 = T1->get_Stiffness_matrix(5,5);
				int smid_1 = TA->get_SM_ID_single();
				int smid_2 = TB->get_SM_ID_single();
				int smid_3 = TC->get_SM_ID_single();
				bool n12 = false;
				bool n13 = false;
				bool n23 = false;
				int N1smn = (int)TA->get_node_sm_single()->size();
				for (int k=0;k<N1smn;k++)
				{
					le_sm *SM1 = TA->get_node_sm_single()->at(k);
					point3D *SM11 = SM1->getN1();
					point3D *SM12 = SM1->getN2();
					if (SM11==TA && SM12==TB)
					{
						n12 = true;
						double f11 = SM1->get_stiffness_matrix_12(0) + v02;
						double f12 = SM1->get_stiffness_matrix_12(1) + v03;
						double f21 = SM1->get_stiffness_matrix_12(2) + v12;
						double f22 = SM1->get_stiffness_matrix_12(3) + v13;
						SM1->set_sm_12(f11,f12,f21,f22);
						SM1->set_sm_21(f11,f21,f12,f22);
					}
					else if (SM11==TB && SM12==TA)
					{
						n12 = true;
						double f11 = SM1->get_stiffness_matrix_21(0) + v02;
						double f12 = SM1->get_stiffness_matrix_21(1) + v03;
						double f21 = SM1->get_stiffness_matrix_21(2) + v12;
						double f22 = SM1->get_stiffness_matrix_21(3) + v13;
						SM1->set_sm_21(f11,f12,f21,f22);
						SM1->set_sm_12(f11,f21,f12,f22);
					}
					if (SM11==TA && SM12==TC)
					{
						n13 = true;
						double f11 = SM1->get_stiffness_matrix_12(0) + v04;
						double f12 = SM1->get_stiffness_matrix_12(1) + v05;
						double f21 = SM1->get_stiffness_matrix_12(2) + v14;
						double f22 = SM1->get_stiffness_matrix_12(3) + v15;
						SM1->set_sm_12(f11,f12,f21,f22);
						SM1->set_sm_21(f11,f21,f12,f22);
					}
					else if (SM11==TC && SM12==TA)
					{
						n13 = true;
						double f11 = SM1->get_stiffness_matrix_21(0) + v04;
						double f12 = SM1->get_stiffness_matrix_21(1) + v05;
						double f21 = SM1->get_stiffness_matrix_21(2) + v14;
						double f22 = SM1->get_stiffness_matrix_21(3) + v15;
						SM1->set_sm_21(f11,f12,f21,f22);
						SM1->set_sm_12(f11,f21,f12,f22);
					}
				}
				int N2smn = (int)TB->get_node_sm_single()->size();
				for (int k=0;k<N2smn;k++)
				{
					le_sm *SM2 = TB->get_node_sm_single()->at(k);
					point3D *SM21 = SM2->getN1();
					point3D *SM22 = SM2->getN2();
					if (SM21==TB && SM22==TC)
					{
						n23 = true;
						double f11 = SM2->get_stiffness_matrix_12(0) + v24;
						double f12 = SM2->get_stiffness_matrix_12(1) + v25;
						double f21 = SM2->get_stiffness_matrix_12(2) + v34;
						double f22 = SM2->get_stiffness_matrix_12(3) + v35;
						SM2->set_sm_12(f11,f12,f21,f22);
						SM2->set_sm_21(f11,f21,f12,f22);
						break;
					}
					else if (SM21==TC && SM22==TB)
					{
						n23 = true;
						double f11 = SM2->get_stiffness_matrix_21(0) + v24;
						double f12 = SM2->get_stiffness_matrix_21(1) + v25;
						double f21 = SM2->get_stiffness_matrix_21(2) + v34;
						double f22 = SM2->get_stiffness_matrix_21(3) + v35;
						SM2->set_sm_21(f11,f12,f21,f22);
						SM2->set_sm_12(f11,f21,f12,f22);
						break;
					}
				}
				if (!n12)
				{
					int l_n = (int)TA->get_node_sm_single()->size();
					le_sm *L = new le_sm(l_n);
					sms_single.push_back(L);
					L->setN1(TA);
					L->setN2(TB);
					L->set_sm_12(v02,v03,v12,v13);
					L->set_sm_21(v02,v12,v03,v13);
					TA->push_node_sm_single(L);
					TB->push_node_sm_single(L);
				}
				if (!n13)
				{
					int l_n = (int)TA->get_node_sm_single()->size();
					le_sm *L = new le_sm(l_n);
					sms_single.push_back(L);
					L->setN1(TA);
					L->setN2(TC);
					L->set_sm_12(v04,v05,v14,v15);
					L->set_sm_21(v04,v14,v05,v15);
					TA->push_node_sm_single(L);
					TC->push_node_sm_single(L);
				}
				if (!n23)
				{
					int l_n = (int)TB->get_node_sm_single()->size();
					le_sm *L = new le_sm(l_n);
					sms_single.push_back(L);
					L->setN1(TB);
					L->setN2(TC);
					L->set_sm_12(v24,v25,v34,v35);
					L->set_sm_21(v24,v34,v25,v35);
					TB->push_node_sm_single(L);
					TC->push_node_sm_single(L);
				}
				double v1_11 = TA->getA_ij_single(0) + v00;
				double v1_12 = TA->getA_ij_single(1) + v01;
				double v1_21 = TA->getA_ij_single(2) + v01;
				double v1_22 = TA->getA_ij_single(3) + v11;
				double v2_11 = TB->getA_ij_single(0) + v22;
				double v2_12 = TB->getA_ij_single(1) + v23;
				double v2_21 = TB->getA_ij_single(2) + v23;
				double v2_22 = TB->getA_ij_single(3) + v33;
				double v3_11 = TC->getA_ij_single(0) + v44;
				double v3_12 = TC->getA_ij_single(1) + v45;
				double v3_21 = TC->getA_ij_single(2) + v45;
				double v3_22 = TC->getA_ij_single(3) + v55;
				TA->setA_ij_single(v1_11,v1_12,v1_21,v1_22);
				TB->setA_ij_single(v2_11,v2_12,v2_21,v2_22);
				TC->setA_ij_single(v3_11,v3_12,v3_21,v3_22);
				SOR_II_V_single[smid_1*2] = v1_11;
				SOR_II_V_single[smid_1*2+1] = v1_22;
				SOR_II_V_single[smid_2*2] = v2_11;
				SOR_II_V_single[smid_2*2+1] = v2_22;
				SOR_II_V_single[smid_3*2] = v3_11;
				SOR_II_V_single[smid_3*2+1] = v3_22;
			}
			for (int j=0;j<sn;j++)
			{
				edge *E = cellList[i]->get_sides()->at(j);
				point3D *p1 = E->p1();
				point3D *p2 = E->p2();
				double K02 = E->get_virtual_triangle_stiffness_matrix(0,2);
				double K03 = E->get_virtual_triangle_stiffness_matrix(0,3);
				double K12 = E->get_virtual_triangle_stiffness_matrix(1,2);
				double K13 = E->get_virtual_triangle_stiffness_matrix(1,3);
				double K00 = E->get_virtual_triangle_stiffness_matrix(0,0);
				double K01 = E->get_virtual_triangle_stiffness_matrix(0,1);
				double K11 = E->get_virtual_triangle_stiffness_matrix(1,1);
				double K22 = E->get_virtual_triangle_stiffness_matrix(2,2);
				double K23 = E->get_virtual_triangle_stiffness_matrix(2,3);
				double K33 = E->get_virtual_triangle_stiffness_matrix(3,3);
				for (int k=0;k<(int)p1->get_node_sm_single()->size();k++)
				{
					if (p1->get_node_sm_single()->at(k)->getN1()==p1 &&
						p1->get_node_sm_single()->at(k)->getN2()==p2)
					{
						double a11 = K02 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_12(0);
						double a12 = K03 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_12(1);
						double a21 = K12 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_12(2);
						double a22 = K13 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_12(3);
						p1->get_node_sm_single()->at(k)->set_sm_12(a11,a12,a21,a22);
						p1->get_node_sm_single()->at(k)->set_sm_21(a11,a21,a12,a22);
						break;
					}
					else if (p1->get_node_sm_single()->at(k)->getN1()==p2 &&
						p1->get_node_sm_single()->at(k)->getN2()==p1)
					{
						double a11 = K02 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_21(0);
						double a12 = K03 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_21(1);
						double a21 = K12 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_21(2);
						double a22 = K13 + p1->get_node_sm_single()->at(k)->get_stiffness_matrix_21(3);
						p1->get_node_sm_single()->at(k)->set_sm_21(a11,a12,a21,a22);
						p1->get_node_sm_single()->at(k)->set_sm_12(a11,a21,a12,a22);
						break;
					}
				}
				double a11_1 = K00 + p1->getA_ij_single(0);
				double a12_1 = K01 + p1->getA_ij_single(1);
				double a21_1 = K01 + p1->getA_ij_single(2);
				double a22_1 = K11 + p1->getA_ij_single(3);
				double a11_2 = K22 + p2->getA_ij_single(0);
				double a12_2 = K23 + p2->getA_ij_single(1);
				double a21_2 = K23 + p2->getA_ij_single(2);
				double a22_2 = K33 + p2->getA_ij_single(3);
				p1->setA_ij_single(a11_1,a12_1,a21_1,a22_1);
				p2->setA_ij_single(a11_2,a12_2,a21_2,a22_2);
				int smid_1 = p1->get_SM_ID_single();
				int smid_2 = p2->get_SM_ID_single();
				double v1_11 = p1->getA_ij_single(0);
				double v1_22 = p1->getA_ij_single(3);
				double v2_11 = p2->getA_ij_single(0);
				double v2_22 = p2->getA_ij_single(3);
				SOR_II_V_single[smid_1*2] = v1_11;
				SOR_II_V_single[smid_1*2+1] = v1_22;
				SOR_II_V_single[smid_2*2] = v2_11;
				SOR_II_V_single[smid_2*2+1] = v2_22;
			}
			double *U_single = new double[smi_id*2];
			double *Un_single = new double[smi_id*2];
			double *F_single = new double[smi_id*2];
			double inc_areaT = cellList[i]->get_GR();
			double P_scale = 0;
			for (int j=0;j<sn;j++)
			{
				point3D *pt = cellList[i]->get_sides()->at(j)->p1();
				double inc_areaTM = pt->get_migration_rate();
				double inc_areaS = pt->get_single_growth_rate();
				double inc_area = 0;
				if (inc_areaT!=0)
				{
					inc_area += inc_areaT/sn;
				}
				else
				{
					double cell_area = cellList[i]->get_area();
					double cell_initial_area = cellList[i]->get_initial_area();
					double area_diff = cell_initial_area - cell_area;
					if      (area_diff>default_area)  {area_diff = default_area;}
					else if (area_diff<-default_area) {area_diff = -default_area;}
					area_diff = area_diff/sn;
					inc_area += area_diff;
				}
				if (inc_areaS!=0)
				{
					inc_area += inc_areaS/sn;
				}
				double ab_l = pt->getab()->getlength();
				double ba_l = pt->getba()->getlength();
				P_scale = 2*inc_area/(ab_l + ba_l);
				double ab_normal_x = pt->getab()->getnormal(0);
				double ab_normal_y = pt->getab()->getnormal(1);
				double ba_normal_x = pt->getba()->getnormal(0);
				double ba_normal_y = pt->getba()->getnormal(1);
				double ab_d = sqrt(ab_normal_x*ab_normal_x + ab_normal_y*ab_normal_y);
				double ba_d = sqrt(ba_normal_x*ba_normal_x + ba_normal_y*ba_normal_y);
				double pt_normal_x = ab_normal_x/ab_d + ba_normal_x/ba_d;
				double pt_normal_y = ab_normal_y/ba_d + ba_normal_y/ba_d;
				double pt_normal_d = sqrt(pt_normal_x*pt_normal_x + pt_normal_y*pt_normal_y);

				double mag_dp_x = P_scale*pt_normal_x/pt_normal_d;
				double mag_dp_y = P_scale*pt_normal_y/pt_normal_d;

				/// migration rate added ///
				double P_migrate_scale = 0;
				if (inc_areaTM!=0)
				{
					double inc_aream = inc_areaTM/sn;
					P_migrate_scale = 2*inc_aream/(ab_l + ba_l);
					double angle_migrate = cellList[i]->get_migrate_angle();
					double x_scale = cos(angle_migrate*PI/180);
					double y_scale = sin(angle_migrate*PI/180);
					mag_dp_x += x_scale*P_migrate_scale;
					mag_dp_y += y_scale*P_migrate_scale;
				}
				F_single[j*2] = mag_dp_x;
				F_single[j*2+1] = mag_dp_y;
				pt->set_pre_dp_mag(P_scale);
				pt->setP_in(mag_dp_x,mag_dp_y);
			}
			for (int j=0;j<ipn;j++)
			{
				F_single[(j+sn)*2] = 0;
				F_single[(j+sn)*2+1] = 0;
			}
			for (int j=0;j<sn;j++)
			{
				int sm_ind = cellList[i]->get_sides()->at(j)->p1()->get_SM_ID_single();
				int nnzr = cellList[i]->get_sides()->at(j)->p1()->get_node_sm_single()->size();
				SOR_Rn_single[sm_ind*2] = nnzr*2 + 1;
				SOR_Rn_single[sm_ind*2+1] = nnzr*2 + 1;
				U_single[sm_ind*2] = 0.1;
				U_single[sm_ind*2+1] = 0.1;
				Un_single[sm_ind*2] = 0.1;
				Un_single[sm_ind*2] = 0.1;
			}
			for (int j=0;j<ipn;j++)
			{
				int sm_ind = cellList[i]->get_inner_p()->at(j)->get_SM_ID_single();
				int nnzr = cellList[i]->get_inner_p()->at(j)->get_node_sm_single()->size();
				SOR_Rn_single[sm_ind*2] = nnzr*2 + 1;
				SOR_Rn_single[sm_ind*2+1] = nnzr*2 + 1;
				U_single[sm_ind*2] = 0.1;
				U_single[sm_ind*2+1] = 0.1;
				Un_single[sm_ind*2] = 0.1;
				Un_single[sm_ind*2] = 0.1;
			}
			int nnz = SOR_Rn_single[0];
			SOR_N_single[0] = nnz;
			int smi_idt = smi_id*2;
			for (int j=1;j<smi_idt;j++)
			{
				nnz = SOR_Rn_single[j] + SOR_N_single[j-1];
				SOR_N_single[j] = nnz;
			}
			delete[] SOR_Rn_single;
			int *JA = new int[nnz];
			double *VA = new double[nnz];
			/// j = 0 ///
			int sm_ind0 = cellList[i]->get_sides()->at(0)->p1()->get_SM_ID_single();
			int nnzr0 = cellList[i]->get_sides()->at(0)->p1()->get_node_sm_single()->size();
			for (int j=0;j<nnzr0;j++)
			{
				le_sm *SMt0 = cellList[i]->get_sides()->at(0)->p1()->get_node_sm_single()->at(j);
				int sm_id10 = SMt0->getN1()->get_SM_ID_single();
				int sm_id20 = SMt0->getN2()->get_SM_ID_single();
				if (sm_id10==sm_ind0)
				{
					int i10 = sm_id20*2;
					int i20 = sm_id20*2+1;
					JA[j*2] = i10;
					JA[j*2 + 1] = i20;
					JA[nnzr0*2 + 1 + j*2] = i10;
					JA[nnzr0*2 + 1 + j*2 + 1] = i20;
					VA[j*2] = SMt0->get_stiffness_matrix_12(0);
					VA[j*2 + 1] = SMt0->get_stiffness_matrix_12(1);
					VA[nnzr0*2 + 1 + j*2] = SMt0->get_stiffness_matrix_12(2);
					VA[nnzr0*2 + 1 + j*2 + 1] = SMt0->get_stiffness_matrix_12(3);
				}
				else
				{
					int i10 = sm_id10*2;
					int i20 = sm_id10*2+1;
					JA[j*2] = i10;
					JA[j*2 + 1] = i20;
					JA[nnzr0*2 + 1 + j*2] = i10;
					JA[nnzr0*2 + 1 + j*2 + 1] = i20;
					VA[j*2] = SMt0->get_stiffness_matrix_21(0);
					VA[j*2 + 1] = SMt0->get_stiffness_matrix_21(1);
					VA[nnzr0*2 + 1 + j*2] = SMt0->get_stiffness_matrix_21(2);
					VA[nnzr0*2 + 1 + j*2 + 1] = SMt0->get_stiffness_matrix_21(3);
				}
			}
			JA[2*nnzr0] = sm_ind0*2 + 1;
			JA[4*nnzr0 + 1] = sm_ind0*2;
			VA[2*nnzr0] = cellList[i]->get_sides()->at(0)->p1()->getA_ij_single(1);
			VA[4*nnzr0 + 1] = cellList[i]->get_sides()->at(0)->p1()->getA_ij_single(2);
			/// 0<j<smi_id ///
			for (int j=1;j<smi_id;j++)
			{
				int sm_ind = 0;
				int nnzr = 0;
				point3D *pt = NULL;
				if (j<sn) 
				{
					pt = cellList[i]->get_sides()->at(j)->p1();
					sm_ind = pt->get_SM_ID_single();
					nnzr = pt->get_node_sm_single()->size();
				}
				else
				{
					pt = cellList[i]->get_inner_p()->at(j-sn);
					sm_ind = pt->get_SM_ID_single();
					nnzr = pt->get_node_sm_single()->size();
				}
				int ind = SOR_N_single[sm_ind*2 - 1];
				for (int k=0;k<nnzr;k++)
				{
					le_sm *SMl = pt->get_node_sm_single()->at(k);
					int sm_id1 = SMl->getN1()->get_SM_ID_single();
					int sm_id2 = SMl->getN2()->get_SM_ID_single();
					if (sm_id1==sm_ind)
					{
						int i1 = sm_id2*2;
						int i2 = sm_id2*2+1;
						JA[ind + k*2] = i1;
						JA[ind + k*2 + 1] = i2;
						JA[ind + nnzr*2 + 1 + k*2] = i1;
						JA[ind + nnzr*2 + 1 + k*2 + 1] = i2;
						VA[ind + k*2] = SMl->get_stiffness_matrix_12(0);
						VA[ind + k*2 + 1] = SMl->get_stiffness_matrix_12(1);
						VA[ind + nnzr*2 + 1 + k*2] = SMl->get_stiffness_matrix_12(2);
						VA[ind + nnzr*2 + 1 + k*2 + 1] = SMl->get_stiffness_matrix_12(3);
					}
					else
					{
						int i1 = sm_id1*2;
						int i2 = sm_id1*2+1;
						JA[ind + k*2] = i1;
						JA[ind + k*2 + 1] = i2;
						JA[ind + nnzr*2 + 1 + k*2] = i1;
						JA[ind + nnzr*2 + 1 + k*2 + 1] = i2;
						VA[ind + k*2] = SMl->get_stiffness_matrix_21(0);
						VA[ind + k*2 + 1] = SMl->get_stiffness_matrix_21(1);
						VA[ind + nnzr*2 + 1 + k*2] = SMl->get_stiffness_matrix_21(2);
						VA[ind + nnzr*2 + 1 + k*2 + 1] = SMl->get_stiffness_matrix_21(3);
					}
				}
				JA[ind + nnzr*2] = sm_ind*2 + 1;
				JA[ind + nnzr*4 + 1] = sm_ind*2;
				VA[ind + nnzr*2] = pt->getA_ij_single(1);
				VA[ind + nnzr*4 + 1] = pt->getA_ij_single(2);
			}
			////// from force to displacement ///////
			double converge = 5.0e-3;
			double w = 1.5;
			double diff = 100;
			int flagtime = 0;
			int runningstep = 500;
			diff = SOR_solver(smi_idt, converge, w, runningstep, flagtime, SOR_N_single, JA, VA, F_single, Un_single, U_single, SOR_II_V_single);
			//cout<<"steps: "<<flagtime<<"; diff: "<<diff<<";"<<endl;
			//////////////
			double control_scale_g = 0; // growth control
			int control_scale_gn = 0;
			double control_scale_m = 0; // migration control
			int control_scale_mn = 0;
			double scale_est = (cellList[i]->get_lame(0)*2 + cellList[i]->get_lame(1))/BD;
			for (int j=0;j<sn;j++)
			{
				point3D *pt = cellList[i]->get_sides()->at(j)->p1();
				if (pt->get_migration_rate()!=0) // migration points
				{
					int ida = pt->get_SM_ID_single();
					double xn = Un_single[ida*2];
					double yn = Un_single[ida*2+1];
					double c_dp = sqrt(xn*xn + yn*yn);
					double r_dp = pt->get_pre_dp_mag();
					double dp_scale = 1;if (c_dp!=0) {dp_scale = r_dp/c_dp;}
					control_scale_m += dp_scale;
					control_scale_mn++;
				}
				else // static points
				{
					int ida = pt->get_SM_ID_single();
					double xn = Un_single[ida*2];
					double yn = Un_single[ida*2+1];
					double c_dp = sqrt(xn*xn + yn*yn);
					double r_dp = pt->get_pre_dp_mag();
					double dp_scale = 1;if (c_dp!=0) {dp_scale = r_dp/c_dp;}
					control_scale_g += dp_scale;
					control_scale_gn++;
					if (i==2)
					{
						cout<<"F,U_"<<pt->id()<<": "<<F_single[ida*2]<<","<<F_single[ida*2+1]<<"; "<<xn<<","<<yn<<"; dp: "<<r_dp<<"; dp_scale: "<<dp_scale<<endl;
					}
				}
			}
			if (control_scale_gn==0) {control_scale_g = 1;}
			else                     {control_scale_g = control_scale_g/control_scale_gn;}
			if (control_scale_mn==0) {control_scale_m = 1;}
			else                     {control_scale_m = control_scale_m/control_scale_mn;}
			control_scale_g = scale_est;
			control_scale_m = scale_est;
			if (i==2) cout<<"controal_scale_"<<i<<": "<<control_scale_g<<","<<control_scale_m<<endl;
			for (int j=0;j<sn;j++)
			{
				point3D *pt = cellList[i]->get_sides()->at(j)->p1();
				double inPx = pt->getP_in(0);
				double inPy = pt->getP_in(1);
				if (pt->get_migration_rate()!=0)
				{
					inPx = inPx*control_scale_m;
					inPy = inPy*control_scale_m;
					pt->setP_in(inPx,inPy);
				}
				else
				{
					inPx = inPx*control_scale_g;
					inPy = inPy*control_scale_g;
					pt->setP_in(inPx,inPy);
				}
			}
			/////////////////////////
			delete[] SOR_II_V_single;
			delete[] SOR_N_single;
			delete[] JA;
			delete[] VA;
			delete[] U_single;
			delete[] Un_single;
			delete[] F_single;
			sms_single.clear();
		}
	}*/

	int sm_id = 0; // matrix id
	/********************************
	 To set up id and pressure force 
	 exerted on the boundary nodes
	********************************/
	time_t t_b_F, t_e_F;
	t_b_F = clock();
	int nn = (int)nodeList.size();
	vector<point3D*> nodepoints;
	for (int i=0;i<nn;i++)
	{
		if (nodeList[i]->get_in_global())
		{
			//////////// clear the old nodes ////////////
			if (nodeList[i]->getvertex()->size()==0) {cout<<"  -> error: node "<<i<<" was empty!"<<endl;}
			nodeList[i]->getvertex()->at(0)->clear_node_sm();
			double fx = 0, fy = 0;
			double f_px = 0, f_py = 0; // without area refinement force
			double a11 = 0, a12 = 0, a21 = 0, a22 = 0;
			nodeList[i]->set_SM_ID(sm_id);
			nodeList[i]->getvertex()->at(0)->set_SM_ID(sm_id);
			nodepoints.push_back(nodeList[i]->getvertex()->at(0));
			int vn = (int)nodeList[i]->getvertex()->size();
			for (int j=0;j<vn;j++)
			{
				point3D *pt = nodeList[i]->getvertex()->at(j);
				double ab_l = pt->getab()->getlength();
				double ba_l = pt->getba()->getlength();
				cell *Cpt = pt->getab()->Cell();
				//if (Cpt->Soften()) continue;
				int csn = (int)Cpt->get_sides()->size();
				double P_scale = 0;
				double P_scale_p = 0; // without area refinement force
				double P_migrate_scale = 0;
				double P_migrate_scale_x = 0;
				double P_migrate_scale_y = 0;
				double P_migrate_response_x = 0;
				double P_migrate_response_y = 0;
				double sFx = pt->get_single_Force(0);
				double sFy = pt->get_single_Force(1);
				///////// incremental pressure due to migration ///////////
				double inc_areaTM = pt->get_migration_rate();
				if (inc_areaTM!=0)
				{
					double inc_area = inc_areaTM/csn;
					//P_migrate_scale = 4*inc_area*Cpt->get_mass()/TM/TM;
					P_migrate_scale = inc_area/BD*(Cpt->get_lame(0)*2 + Cpt->get_lame(1))/BD;
					double angle_migrate = Cpt->get_migrate_angle();
					double x_scale = cos(angle_migrate*PI/180);
					double y_scale = sin(angle_migrate*PI/180);
					P_migrate_scale_x = x_scale*P_migrate_scale;
					P_migrate_scale_y = y_scale*P_migrate_scale; // second newton mode
				}
				if (pt->get_migrate_response_force(0)!=0 ||
					pt->get_migrate_response_force(1)!=0)
				{
					P_migrate_response_x = pt->get_migrate_response_force(0);
					P_migrate_response_y = pt->get_migrate_response_force(1);
					pt->set_migrate_response_force(0,0);
					pt->set_pop_out(1);
				}
				///////// incremental pressure due to incremental area /////////
				double inc_areaT = Cpt->get_GR();
				double inc_areaS = pt->get_single_growth_rate();
				if (inc_areaT!=0 || inc_areaS!=0)
				{
					double inc_area = (inc_areaT + inc_areaS)/csn;
					//if (F_mode==0)
					//{
						//P_scale = 4*inc_area*Cpt->get_mass()/TM/TM; // the second-newton-law mode
						P_scale = inc_area/BD*(Cpt->get_lame(0)*2 + Cpt->get_lame(1))/BD; // virtual length: 1 micrometer
					//}
					//else if (F_mode==1)
					//{
					//	P_scale = 2*inc_area/(ab_l + ba_l); // the elastic mode 
					//}
						if (!Cpt->get_migrate_relax()) {P_scale_p = P_scale;} // cells only grow but not migrate
				}
				else
				{
					if (!Cpt->Soften())
					{
						double cell_area = Cpt->get_area();
						double cell_initial_area = Cpt->get_initial_area();
						double area_diff = cell_initial_area - cell_area;
						if      (area_diff>default_area)  {area_diff = default_area;}
						else if (area_diff<-default_area) {area_diff = -default_area;}
						area_diff = area_diff/csn;
						//if (F_mode==0)
						//{
						//	P_scale = 4*area_diff*Cpt->get_mass()/TM/TM; // the second-newton-law mode
							P_scale = area_diff/BD*(Cpt->get_lame(0)*2 + Cpt->get_lame(1))/BD; // virtual length: 1 micrometer
						//}
						//else if (F_mode==1)
						//{
						//	P_scale = 2*area_diff/(ab_l + ba_l); // the elastic mode
						//}
					}
				}
				/////// important to determine internal force direction ////////
				if      (inc_areaT==0) {pt->set_pop_out(0);}
				else if (inc_areaT<0)  {pt->set_pop_out(-1);}
				else if (inc_areaT>0)  {pt->set_pop_out(1);}
				if      (inc_areaTM>0) {pt->set_pop_out(1);}
				double ab_normal_x = pt->getab()->getnormal(0);
				double ab_normal_y = pt->getab()->getnormal(1);
				double ba_normal_x = pt->getba()->getnormal(0);
				double ba_normal_y = pt->getba()->getnormal(1);
				double ab_d = sqrt(ab_normal_x*ab_normal_x + ab_normal_y*ab_normal_y);
				double ba_d = sqrt(ba_normal_x*ba_normal_x + ba_normal_y*ba_normal_y);
				double pt_normal_x = ab_normal_x/ab_d + ba_normal_x/ba_d;
				double pt_normal_y = ab_normal_y/ba_d + ba_normal_y/ba_d;
				double pt_normal_d = sqrt(pt_normal_x*pt_normal_x + pt_normal_y*pt_normal_y);
				double inc_P_ab_x = 0;
				double inc_P_ab_y = 0;
				double inc_P_ba_x = 0;
				double inc_P_ba_y = 0;
				//if (F_mode==0)
				//{
					inc_P_ab_x = ab_normal_x/ab_d*(P_scale) + P_migrate_scale_x*0.5;
					inc_P_ab_y = ab_normal_y/ab_d*(P_scale) + P_migrate_scale_y*0.5;
					inc_P_ba_x = ba_normal_x/ba_d*(P_scale) + P_migrate_scale_x*0.5;
					inc_P_ba_y = ba_normal_y/ba_d*(P_scale) + P_migrate_scale_y*0.5; // the second-newton-law mode
				//}
				///
				double inc_P_p_ab_x = ab_normal_x/ab_d*(P_scale_p) + P_migrate_scale_x*0.5;
				double inc_P_p_ab_y = ab_normal_y/ab_d*(P_scale_p) + P_migrate_scale_y*0.5;
				double inc_P_p_ba_x = ba_normal_x/ba_d*(P_scale_p) + P_migrate_scale_x*0.5;
				double inc_P_p_ba_y = ba_normal_y/ba_d*(P_scale_p) + P_migrate_scale_y*0.5;
				///
				double inc_P_x = 0;
				double inc_P_y = 0;
				//if (F_mode==1)
				//{
				//	inc_P_x = pt->getP_in(0);
				//	inc_P_y = pt->getP_in(1); // the elastic mode, get the force directly from the database
				//}
				//////////////// the basic pressure and tension ///////////////
				double T_ab_x = (pt->getfp()->x() - pt->x())*Cpt->get_tcoef();
				double T_ab_y = (pt->getfp()->y() - pt->y())*Cpt->get_tcoef();
				double T_ba_x = (pt->getrp()->x() - pt->x())*Cpt->get_tcoef();
				double T_ba_y = (pt->getrp()->y() - pt->y())*Cpt->get_tcoef();
				double P_ab_x = 0.5*ab_normal_x/ab_d*ab_l*Cpt->get_pcoef();
				double P_ab_y = 0.5*ab_normal_y/ab_d*ab_l*Cpt->get_pcoef();
				double P_ba_x = 0.5*ba_normal_x/ab_d*ba_l*Cpt->get_pcoef();
				double P_ba_y = 0.5*ba_normal_y/ab_d*ba_l*Cpt->get_pcoef();
				double BF_x = pt->get_adhesion_break(0);
				double BF_y = pt->get_adhesion_break(1);
				double Force_x = T_ab_x + T_ba_x + P_ab_x + P_ba_x + BF_x;
				double Force_y = T_ab_y + T_ba_y + P_ab_y + P_ba_y + BF_y;
				Force_x += P_migrate_response_x;
				Force_y += P_migrate_response_y;
				Force_x += sFx;
				Force_y += sFy;
				pt->setTension_ab(T_ab_x, T_ab_y);
				pt->setTension_ba(T_ba_x, T_ba_y);
				pt->setPressure_ab(P_ab_x, P_ab_y);
				pt->setPressure_ba(P_ba_x, P_ba_y);
				pt->setForce(Force_x, Force_y);
				//if (F_mode==0)
				//{
					pt->setP_in_ab(inc_P_ab_x, inc_P_ab_y);
					pt->setP_in_ba(inc_P_ba_x, inc_P_ba_y);
				//}
				//else if (F_mode==1)
				//{
				//	pt->setP_in_ab(inc_P_x/2, inc_P_y/2);
				//	pt->setP_in_ba(inc_P_x/2, inc_P_y/2);
				//}
				pt->set_adhesion_break(0,0); // the break adhesion releases instantly after the edge breaks
				fx += pt->getForce(0);
				fy += pt->getForce(1);
				fx += pt->getP_in_ab(0);
				fx += pt->getP_in_ba(0);
				fy += pt->getP_in_ab(1);
				fy += pt->getP_in_ba(1);
				///
				f_px += inc_P_p_ab_x + inc_P_p_ba_x + P_migrate_response_x + sFx;
				f_py += inc_P_p_ab_y + inc_P_p_ba_y + P_migrate_response_y + sFy;
				/////////// clear the old diagnal value ///////////
				pt->setA_ij(0,0,0,0);
				//////// clear the old single growth rate /////////
				pt->set_single_growth_rate(0);
			}
			F.push_back(fx);
			F.push_back(fy);
			F_p.push_back(f_px);
			F_p.push_back(f_py);
			sm_id++;
		}
	}
	/***********************************
	 To set up id for the interior nodes 
	 and calculate the stiffness matrices 
	 for the interior triangles
	***********************************/
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (cellList[i]->Soften()) continue;
		//////// clear the old mechanical signal info ///////
		cellList[i]->set_sense_mech(false);
		cellList[i]->set_sense_mech_vec(0, 0);
		cellList[i]->set_sense_migr(false);
		cellList[i]->set_sense_migr_vec(0, 0);
		int pn = (int)cellList[i]->get_inner_p()->size();
		if (pn>0)
		{
			for (int j=0;j<pn;j++)
			{
				point3D *pin = cellList[i]->get_inner_p()->at(j);
				pin->setid(sm_id); // node id assigned to the point only for linear equation
				pin->set_SM_ID(sm_id); // node id assigned to the point only for linear equation
				nodepoints.push_back(pin);
				////////// clear the old diagonal value  and nodes ///////////
				pin->setA_ij(0,0,0,0);
				pin->clear_node_sm();
				F.push_back(0); 
				F.push_back(0); // assuming interior forces on interior points equal to 0
				F_p.push_back(0);
				F_p.push_back(0);
				sm_id++;
			}
		}
	}
	t_e_F = clock();
	cout<<"Force vector build up: "<<(t_e_F-t_b_F)*0.001<<" seconds!"<<endl;

	/***********************************
	 Set up global Stiffness matrix
	 -----------------------------------
	 volume: sm_id X sm_id

	 * 0 * 0 0 * 0 0 | 0 | 
	 0 * 0 * 0 0 0 0 | 0 | 
	 0 0 * 0 0 * 0 0 | 0 | 
	 0 0 0 * 0 * 0 0 | 0 | 
	 0 0 0 0 * 0 * * | 0 | 
	 0 0 0 0 0 * * 0 | 0 | 
	 0 0 0 0 0 0 * 0 | 0 | 
	 0 0 0 0 0 0 0 * | 0 | 
	 ----------------|---|
	 0 0 0 0 0 0 0 0 | * | A
	 ----------------|---|
					   A   
	 ----------------------------------
	 One virtual point A is added to
	 ensure the matrix is invertible:
				A
			   / \
		    ==1=<=0==
	 0->-1 is one single edge without 
	 neighbor edge attached to it. The 
	 virtual triangle 01A is assigned 
	 with soft material property
	***********************************/
	double *SOR_II_V = new double[sm_id*2]; // int i, diagonal elements
	int *SOR_N = new int[sm_id*2];          // number of nonzero elements
	int *SOR_Rn = new int[sm_id*2];         // number of nonzero elements of each row
	
	time_t t_begin2, t_end2;
	t_begin2 = clock();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (!cellList[i]->Soften())
		{
			int tn = (int)cellList[i]->get_inner_t()->size();
			for (int j=0;j<tn;j++)
			{
				triangle *T1 = cellList[i]->get_inner_t()->at(j);
				//if (F_mode==0) {T1->set_Stiffness_matrix();}
				T1->set_Stiffness_matrix();
				point3D* TA = T1->getA();
				point3D* TB = T1->getB();
				point3D* TC = T1->getC();
				double v00 = T1->get_Stiffness_matrix(0,0);
				double v01 = T1->get_Stiffness_matrix(0,1);
				double v02 = T1->get_Stiffness_matrix(0,2);
				double v03 = T1->get_Stiffness_matrix(0,3);
				double v04 = T1->get_Stiffness_matrix(0,4);
				double v05 = T1->get_Stiffness_matrix(0,5);
				double v11 = T1->get_Stiffness_matrix(1,1);
				double v12 = T1->get_Stiffness_matrix(1,2);
				double v13 = T1->get_Stiffness_matrix(1,3);
				double v14 = T1->get_Stiffness_matrix(1,4);
				double v15 = T1->get_Stiffness_matrix(1,5);
				double v22 = T1->get_Stiffness_matrix(2,2);
				double v23 = T1->get_Stiffness_matrix(2,3);
				double v24 = T1->get_Stiffness_matrix(2,4);
				double v25 = T1->get_Stiffness_matrix(2,5);
				double v33 = T1->get_Stiffness_matrix(3,3);
				double v34 = T1->get_Stiffness_matrix(3,4);
				double v35 = T1->get_Stiffness_matrix(3,5);
				double v44 = T1->get_Stiffness_matrix(4,4);
				double v45 = T1->get_Stiffness_matrix(4,5);
				double v55 = T1->get_Stiffness_matrix(5,5);
				/*********************************************************************
				To construct the nodes interaction
					N1      N2      N3
					---------------------
				N1  00 01 | 02 03 | 04 05
					10 11 | 12 13 | 14 15
					---------------------
				N2  20 21 | 22 23 | 24 25
					30 31 | 32 33 | 34 35
					---------------------
				N3  40 41 | 42 43 | 44 45
					50 51 | 52 53 | 54 55
					---------------------
				*********************************************************************/
				point3D *N1, *N2, *N3;
				int smid_1 = 0;
				int smid_2 = 0;
				int smid_3 = 0;
				if (TA->get_cell_boundary()) 
				{
					if (TA->get_node_id()==NULL)
					{
						cout<<"  -> error: TA "<<TA->id()<<" from cell "<<i<<": node is NULL!"<<endl;
					}
					N1 = TA->get_node_id()->getvertex()->at(0);
					smid_1 = TA->get_node_id()->SM_ID();
				}
				else 
				{
					N1 = TA;
					smid_1 = TA->get_SM_ID();
				}
				if (TB->get_cell_boundary()) 
				{
					if (TB->get_node_id()==NULL)
					{
						cout<<"  -> error: TB "<<TB->id()<<" from cell "<<i<<": node is NULL!"<<endl;
					}
					N2 = TB->get_node_id()->getvertex()->at(0);
					smid_2 = TB->get_node_id()->SM_ID();
				}
				else 
				{
					N2 = TB;
					smid_2 = TB->get_SM_ID();
				}
				if (TC->get_cell_boundary()) 
				{
					if (TC->get_node_id()==NULL)
					{
						cout<<"  -> error: TC "<<TC->id()<<" from cell "<<i<<": node is NULL!"<<endl;
					}
					N3 = TC->get_node_id()->getvertex()->at(0);
					smid_3 = TC->get_node_id()->SM_ID();
				}
				else 
				{
					N3 = TC;
					smid_3 = TC->get_SM_ID();
				}
				/*if (number==21 && (i==3271 || i==5178 || i==3357))
				{
					if (smid_1==smid_2 || smid_1==smid_3 || smid_2==smid_3)
					{
						cout<<" cell_"<<i<<":triangle sm_id:"<<j<<":"<<smid_1<<","<<smid_2<<","<<smid_3<<endl;
						cout<<" cell_"<<i<<":triangle N id:"<<j<<":"<<N1->id()<<","<<N2->id()<<","<<N3->id()<<endl;
						cout<<" cell_"<<i<<":triangle N node id:"<<j<<":";
						if (N1->get_cell_boundary()) {cout<<N1->get_node_id()->id()<<",";} else {cout<<"N1"<<",";}
						if (N2->get_cell_boundary()) {cout<<N2->get_node_id()->id()<<",";} else {cout<<"N2"<<",";}
						if (N3->get_cell_boundary()) {cout<<N3->get_node_id()->id()<<endl;} else {cout<<"N3"<<endl;}
						cout<<" cell_"<<i<<":triangle p id:"<<j<<":"<<TA->id()<<","<<TB->id()<<","<<TC->id()<<endl;
						cout<<" cell_"<<i<<":triangle p x:"<<j<<":"<<TA->x()<<","<<TB->x()<<","<<TC->x()<<endl;
						cout<<" cell_"<<i<<":triangle v:"<<j<<v00<<","<<v01<<","<<v02<<","<<v03<<","<<v04<<","<<v05<<",";
						cout<<v11<<","<<v12<<","<<v13<<","<<v14<<","<<v15<<",";
						cout<<v22<<","<<v23<<","<<v24<<","<<v25<<",";
						cout<<v33<<","<<v34<<","<<v35<<",";
						cout<<v44<<","<<v45<<",";
						cout<<v55<<endl;
					}
				}*/
				bool n12 = false;
				bool n13 = false;
				bool n23 = false;
				int N1smn = (int)N1->get_node_sm()->size();
				for (int k=0;k<N1smn;k++)
				{
					le_sm *SM1 = N1->get_node_sm()->at(k);
					point3D *SM11 = SM1->getN1();
					point3D *SM12 = SM1->getN2();
					if (SM11==N1 && SM12==N2)
					{
						//		N1 N2
						//	N1  00 00
						//	N2  00 00
						//	12: N1-N2
						//	21: N2-N1
						//diagonal element value 
						//independent from index
						n12 = true;
						double f11 = SM1->get_stiffness_matrix_12(0) + v02;//T1->get_Stiffness_matrix(0,2);
						double f12 = SM1->get_stiffness_matrix_12(1) + v03;//T1->get_Stiffness_matrix(0,3);
						double f21 = SM1->get_stiffness_matrix_12(2) + v12;//T1->get_Stiffness_matrix(1,2);
						double f22 = SM1->get_stiffness_matrix_12(3) + v13;//T1->get_Stiffness_matrix(1,3);
						SM1->set_sm_12(f11,f12,f21,f22);
						SM1->set_sm_21(f11,f21,f12,f22);
					}
					else if (SM11==N2 && SM12==N1)
					{
						//	   N2 N1
						//	N2 00 00
						//	N1 00 00
						//	12: N2-N1
						//	21: N1-N2
						//diagonal element value 
						//independent from index
						n12 = true;
						double f11 = SM1->get_stiffness_matrix_21(0) + v02;//T1->get_Stiffness_matrix(0,2);
						double f12 = SM1->get_stiffness_matrix_21(1) + v03;//T1->get_Stiffness_matrix(0,3);
						double f21 = SM1->get_stiffness_matrix_21(2) + v12;//T1->get_Stiffness_matrix(1,2);
						double f22 = SM1->get_stiffness_matrix_21(3) + v13;//T1->get_Stiffness_matrix(1,3);
						SM1->set_sm_21(f11,f12,f21,f22); // f11: N1-N1, f12: N1-N2
						SM1->set_sm_12(f11,f21,f12,f22); // f22: N2-N2, f21: N2-N1
					}
					if (SM11==N1 && SM12==N3)
					{
						n13 = true;
						double f11 = SM1->get_stiffness_matrix_12(0) + v04;//T1->get_Stiffness_matrix(0,4);
						double f12 = SM1->get_stiffness_matrix_12(1) + v05;//T1->get_Stiffness_matrix(0,5);
						double f21 = SM1->get_stiffness_matrix_12(2) + v14;//T1->get_Stiffness_matrix(1,4);
						double f22 = SM1->get_stiffness_matrix_12(3) + v15;//T1->get_Stiffness_matrix(1,5);
						SM1->set_sm_12(f11,f12,f21,f22); // f11: N1-N1, f12: N1-N3
						SM1->set_sm_21(f11,f21,f12,f22); // f22: N3-N3, f21: N3-N1
					}
					else if (SM11==N3 && SM12==N1)
					{
						n13 = true;
						double f11 = SM1->get_stiffness_matrix_21(0) + v04;//T1->get_Stiffness_matrix(0,4);
						double f12 = SM1->get_stiffness_matrix_21(1) + v05;//T1->get_Stiffness_matrix(0,5);
						double f21 = SM1->get_stiffness_matrix_21(2) + v14;//T1->get_Stiffness_matrix(1,4);
						double f22 = SM1->get_stiffness_matrix_21(3) + v15;//T1->get_Stiffness_matrix(1,5);
						SM1->set_sm_21(f11,f12,f21,f22); // f11: N1-N1, f12: N1-N3
						SM1->set_sm_12(f11,f21,f12,f22); // f22: N3-N3, f21: N3-N1
					}
				}
				int N2smn = (int)N2->get_node_sm()->size();
				for (int k=0;k<N2smn;k++)
				{
					le_sm *SM2 = N2->get_node_sm()->at(k);
					point3D *SM21 = SM2->getN1();
					point3D *SM22 = SM2->getN2();
					if (SM21==N2 && SM22==N3)
					{
						n23 = true;
						double f11 = SM2->get_stiffness_matrix_12(0) + v24;//T1->get_Stiffness_matrix(2,4);
						double f12 = SM2->get_stiffness_matrix_12(1) + v25;//T1->get_Stiffness_matrix(2,5);
						double f21 = SM2->get_stiffness_matrix_12(2) + v34;//T1->get_Stiffness_matrix(3,4);
						double f22 = SM2->get_stiffness_matrix_12(3) + v35;//T1->get_Stiffness_matrix(3,5);
						SM2->set_sm_12(f11,f12,f21,f22); // f11: N2-N2, f12: N2-N3
						SM2->set_sm_21(f11,f21,f12,f22); // f22: N3-N3, f21: N3-N2
					}
					else if (SM21==N3 && SM22==N2)
					{
						n23 = true;
						double f11 = SM2->get_stiffness_matrix_21(0) + v24;//T1->get_Stiffness_matrix(2,4);
						double f12 = SM2->get_stiffness_matrix_21(1) + v25;//T1->get_Stiffness_matrix(2,5);
						double f21 = SM2->get_stiffness_matrix_21(2) + v34;//T1->get_Stiffness_matrix(3,4);
						double f22 = SM2->get_stiffness_matrix_21(3) + v35;//T1->get_Stiffness_matrix(3,5);
						SM2->set_sm_21(f11,f12,f21,f22); // f11: N2-N2, f12: N2-N3
						SM2->set_sm_12(f11,f21,f12,f22); // f22: N3-N3, f21: N3-N2
					}
				}
				if (!n12)
				{
					int l_n = (int)N1->get_node_sm()->size();
					le_sm *L = new le_sm(l_n);
					sms.push_back(L);
					L->setN1(N1);
					L->setN2(N2);
					L->set_sm_12(v02,v03,v12,v13);
					L->set_sm_21(v02,v12,v03,v13);
					N1->push_node_sm(L);
					N2->push_node_sm(L);
					N1->set_node_sm_n();
					N2->set_node_sm_n();
				}
				if (!n13)
				{
					int l_n = (int)N1->get_node_sm()->size();
					le_sm *L = new le_sm(l_n);
					sms.push_back(L);
					L->setN1(N1);
					L->setN2(N3);
					L->set_sm_12(v04,v05,v14,v15);
					L->set_sm_21(v04,v14,v05,v15);
					N1->push_node_sm(L);
					N3->push_node_sm(L);
					N1->set_node_sm_n();
					N3->set_node_sm_n();
				}
				if (!n23)
				{
					int l_n = (int)N2->get_node_sm()->size();
					le_sm *L = new le_sm(l_n);
					sms.push_back(L);
					L->setN1(N2);
					L->setN2(N3);
					L->set_sm_12(v24,v25,v34,v35);
					L->set_sm_21(v24,v34,v25,v35);
					N2->push_node_sm(L);
					N3->push_node_sm(L);
					N2->set_node_sm_n();
					N3->set_node_sm_n();
				}
				double v1_11 = N1->getA_ij(0) + v00;//T1->get_Stiffness_matrix(0,0);
				double v1_12 = N1->getA_ij(1) + v01;//T1->get_Stiffness_matrix(0,1);
				double v1_21 = N1->getA_ij(2) + v01;//T1->get_Stiffness_matrix(1,0);
				double v1_22 = N1->getA_ij(3) + v11;//T1->get_Stiffness_matrix(1,1);
				double v2_11 = N2->getA_ij(0) + v22;//T1->get_Stiffness_matrix(2,2);
				double v2_12 = N2->getA_ij(1) + v23;//T1->get_Stiffness_matrix(2,3);
				double v2_21 = N2->getA_ij(2) + v23;//T1->get_Stiffness_matrix(3,2);
				double v2_22 = N2->getA_ij(3) + v33;//T1->get_Stiffness_matrix(3,3);
				double v3_11 = N3->getA_ij(0) + v44;//T1->get_Stiffness_matrix(4,4);
				double v3_12 = N3->getA_ij(1) + v45;//T1->get_Stiffness_matrix(4,5);
				double v3_21 = N3->getA_ij(2) + v45;//T1->get_Stiffness_matrix(5,4);
				double v3_22 = N3->getA_ij(3) + v55;//T1->get_Stiffness_matrix(5,5);
				N1->setA_ij(v1_11,v1_12,v1_21,v1_22);
				N2->setA_ij(v2_11,v2_12,v2_21,v2_22);
				N3->setA_ij(v3_11,v3_12,v3_21,v3_22);
				SOR_II_V[smid_1*2] = v1_11;
				SOR_II_V[smid_1*2+1] = v1_22;
				SOR_II_V[smid_2*2] = v2_11;
				SOR_II_V[smid_2*2+1] = v2_22;
				SOR_II_V[smid_3*2] = v3_11;
				SOR_II_V[smid_3*2+1] = v3_22;
			}
		}
		else
		{
			int csn = (int)cellList[i]->get_sides()->size();
			for (int j=0;j<csn;j++)
			{
				edge *E = cellList[i]->get_sides()->at(j);
				E->set_soften_stiffness_matrix();
				point3D *TA = E->p1();
				point3D *TB = E->p2();
				double v00 = E->get_soften_triangle_stiffness_matrix(0,0);
				double v01 = E->get_soften_triangle_stiffness_matrix(0,1);
				double v02 = E->get_soften_triangle_stiffness_matrix(0,2);
				double v03 = E->get_soften_triangle_stiffness_matrix(0,3);
				double v11 = E->get_soften_triangle_stiffness_matrix(1,1);
				double v12 = E->get_soften_triangle_stiffness_matrix(1,2);
				double v13 = E->get_soften_triangle_stiffness_matrix(1,3);
				double v22 = E->get_soften_triangle_stiffness_matrix(2,2);
				double v23 = E->get_soften_triangle_stiffness_matrix(2,3);
				double v33 = E->get_soften_triangle_stiffness_matrix(3,3);
				/*********************************
				To construct the nodes interaction
					N1      N2
					-------------
				N1  00 01 | 02 03
					10 11 | 12 13
					-------------
				N2  20 21 | 22 23
					30 31 | 32 33
					-------------
				*********************************/
				point3D *N1, *N2;
				int smid_1 = 0;
				int smid_2 = 0;
				N1 = TA->get_node_id()->getvertex()->at(0);
				smid_1 = TA->get_node_id()->SM_ID();
				N2 = TB->get_node_id()->getvertex()->at(0);
				smid_2 = TB->get_node_id()->SM_ID();
				bool n12 = false;
				int N1smn = (int)N1->get_node_sm()->size();
				for (int k=0;k<N1smn;k++)
				{
					le_sm *SM1 = N1->get_node_sm()->at(k);
					point3D *SM11 = SM1->getN1();
					point3D *SM12 = SM1->getN2();
					if (SM11==N1 && SM12==N2)
					{
						n12 = true;
						double f11 = SM1->get_stiffness_matrix_12(0) + v02;
						double f12 = SM1->get_stiffness_matrix_12(1) + v03;
						double f21 = SM1->get_stiffness_matrix_12(2) + v12;
						double f22 = SM1->get_stiffness_matrix_12(3) + v13;
						SM1->set_sm_12(f11,f12,f21,f22);
						SM1->set_sm_21(f11,f21,f12,f22);
					}
					else if (SM11==N2 && SM12==N1)
					{
						n12 = true;
						double f11 = SM1->get_stiffness_matrix_21(0) + v02;
						double f12 = SM1->get_stiffness_matrix_21(1) + v03;
						double f21 = SM1->get_stiffness_matrix_21(2) + v12;
						double f22 = SM1->get_stiffness_matrix_21(3) + v13;
						SM1->set_sm_21(f11,f12,f21,f22);
						SM1->set_sm_12(f11,f21,f12,f22);
					}
				}
				if (!n12)
				{
					int l_n = (int)N1->get_node_sm()->size();
					le_sm *L = new le_sm(l_n);
					sms.push_back(L);
					L->setN1(N1);
					L->setN2(N2);
					L->set_sm_12(v02,v03,v12,v13);
					L->set_sm_21(v02,v12,v03,v13);
					N1->push_node_sm(L);
					N2->push_node_sm(L);
					N1->set_node_sm_n();
					N2->set_node_sm_n();
				}
				double v1_11 = N1->getA_ij(0) + v00;
				double v1_12 = N1->getA_ij(1) + v01;
				double v1_21 = N1->getA_ij(2) + v01;
				double v1_22 = N1->getA_ij(3) + v11;
				double v2_11 = N2->getA_ij(0) + v22;
				double v2_12 = N2->getA_ij(1) + v23;
				double v2_21 = N2->getA_ij(2) + v23;
				double v2_22 = N2->getA_ij(3) + v33;
				N1->setA_ij(v1_11,v1_12,v1_21,v1_22);
				N2->setA_ij(v2_11,v2_12,v2_21,v2_22);
				SOR_II_V[smid_1*2] = v1_11;
				SOR_II_V[smid_1*2+1] = v1_22;
				SOR_II_V[smid_2*2] = v2_11;
				SOR_II_V[smid_2*2+1] = v2_22;
			}
		}
		if (!cellList[i]->Soften())
		{
			int sn = (int)cellList[i]->get_sides()->size();
			for (int j=0;j<sn;j++)
			{
				edge *E = cellList[i]->get_sides()->at(j);
				point3D *p1 = E->p1()->get_node_id()->getvertex()->at(0);
				// update x coordinates of last 3 steps
				double uo0 = p1->x_old(0);
				double uo1 = p1->x_old(1);
				double vo0 = p1->y_old(0);
				double vo1 = p1->y_old(1);
				p1->setx_old(uo1,2);
				p1->setx_old(uo0,1);
				p1->setx_old(p1->x(),0);
				p1->sety_old(vo1,2);
				p1->sety_old(vo0,1);
				p1->sety_old(p1->y(),0);
				//
				if (E->get_attach()) continue;
				E->set_virtual_stiffness_matrix();
				point3D *p2 = E->p2()->get_node_id()->getvertex()->at(0);
				double v1_11 = p1->getA_ij(0);
				double v1_12 = p1->getA_ij(1);
				double v1_21 = p1->getA_ij(2);
				double v1_22 = p1->getA_ij(3);
				double v2_11 = p2->getA_ij(0);
				double v2_12 = p2->getA_ij(1);
				double v2_21 = p2->getA_ij(2);
				double v2_22 = p2->getA_ij(3);
				double f11 = cellList[i]->get_sides()->at(j)->get_virtual_stiffness_matrix(0);
				double f12 = cellList[i]->get_sides()->at(j)->get_virtual_stiffness_matrix(1);
				double f21 = cellList[i]->get_sides()->at(j)->get_virtual_stiffness_matrix(2);
				double f22 = cellList[i]->get_sides()->at(j)->get_virtual_stiffness_matrix(3);
				int smid_1 = p1->get_node_id()->SM_ID();
				int smid_2 = p2->get_node_id()->SM_ID();
				SOR_II_V[smid_1*2] = v1_11;
				SOR_II_V[smid_1*2+1] = v1_22;
				SOR_II_V[smid_2*2] = v2_11;
				SOR_II_V[smid_2*2+1] = v2_22;
			}
		}
	}
	
	double *Xn = new double[sm_id*2];
	double *X = new double[sm_id*2];
	double *Fn = new double[sm_id*2];
	/// without area refinement force
	double *Xn_p = new double [sm_id*2];
	double *X_p = new double [sm_id*2];
	double *Fn_p = new double [sm_id*2];
	for (int i=0;i<sm_id;i++)
	{
		int sm_ind = nodepoints[i]->get_SM_ID();
		int nnzr = nodepoints[i]->get_node_sm()->size();
		SOR_Rn[sm_ind*2] = nnzr*2 + 1;
		SOR_Rn[sm_ind*2+1] = nnzr*2 + 1;
		Xn[sm_ind*2] = 0.1;
		Xn[sm_ind*2+1] = 0.1;
		X[sm_ind*2] = 0.1;
		X[sm_ind*2+1] = 0.1;
		Fn[sm_ind*2] = F[sm_ind*2];
		Fn[sm_ind*2+1] = F[sm_ind*2+1];
		///
		Xn_p[sm_ind*2] = 0.1;
		Xn_p[sm_ind*2+1] = 0.1;
		X_p[sm_ind*2] = 0.1;
		X_p[sm_ind*2+1] = 0.1;
		Fn_p[sm_ind*2] = F_p[sm_ind*2];
		Fn_p[sm_ind*2+1] = F_p[sm_ind*2+1];
	}
	int nnz = SOR_Rn[0];
	SOR_N[0] = nnz;
	int sm_idt = sm_id*2;
	for (int i=1;i<sm_idt;i++)
	{
		nnz = SOR_Rn[i] + SOR_N[i-1];
		SOR_N[i] = nnz;
	}
	delete SOR_Rn;
	int *JA = new int[nnz];
	double *VA = new double[nnz];
	///////// i==0 //////////
	int sm_ind0 = nodepoints[0]->get_SM_ID();
	int nnzr0 = nodepoints[0]->get_node_sm()->size();
	for (int j=0;j<nnzr0;j++)
	{
		le_sm *SMt0 = nodepoints[0]->get_node_sm()->at(j);
		int sm_id10 = SMt0->getN1()->get_SM_ID();
		int sm_id20 = SMt0->getN2()->get_SM_ID();
		if (sm_id10==sm_ind0)
		{
			int i10 = sm_id20*2;
			int i20 = sm_id20*2+1;
			JA[j*2] = i10;
			JA[j*2 + 1] = i20;
			JA[nnzr0*2 + 1 + j*2] = i10;
			JA[nnzr0*2 + 1 + j*2 + 1] = i20;
			VA[j*2] = SMt0->get_stiffness_matrix_12(0);
			VA[j*2 + 1] = SMt0->get_stiffness_matrix_12(1);
			VA[nnzr0*2 + 1 + j*2] = SMt0->get_stiffness_matrix_12(2);
			VA[nnzr0*2 + 1 + j*2 + 1] = SMt0->get_stiffness_matrix_12(3);
		}
		else
		{
			int i10 = sm_id10*2;
			int i20 = sm_id10*2+1;
			JA[j*2] = i10;
			JA[j*2 + 1] = i20;
			JA[nnzr0*2 + 1 + j*2] = i10;
			JA[nnzr0*2 + 1 + j*2 + 1] = i20;
			VA[j*2] = SMt0->get_stiffness_matrix_21(0);
			VA[j*2 + 1] = SMt0->get_stiffness_matrix_21(1);
			VA[nnzr0*2 + 1 + j*2] = SMt0->get_stiffness_matrix_21(2);
			VA[nnzr0*2 + 1 + j*2 + 1] = SMt0->get_stiffness_matrix_21(3);
		}
	}
	JA[2*nnzr0] = sm_ind0*2 + 1;
	JA[4*nnzr0 + 1] = sm_ind0*2;
	VA[2*nnzr0] = nodepoints[0]->getA_ij(1);
	VA[4*nnzr0 + 1] = nodepoints[0]->getA_ij(2);
	/////// 0<i<sm_id ///////
	for (int i=1;i<sm_id;i++)
	{
		int sm_ind = nodepoints[i]->get_SM_ID();
		int nnzr = nodepoints[i]->get_node_sm()->size();
		int ind = SOR_N[sm_ind*2 - 1];
		for (int j=0;j<nnzr;j++)
		{
			le_sm *SMl = nodepoints[i]->get_node_sm()->at(j);
			int sm_id1 = SMl->getN1()->get_SM_ID();
			int sm_id2 = SMl->getN2()->get_SM_ID();
			if (sm_id1==sm_ind)
			{
				int i1 = sm_id2*2;
				int i2 = sm_id2*2+1;
				JA[ind + j*2] = i1;
				JA[ind + j*2 + 1] = i2;
				JA[ind + nnzr*2 + 1 + j*2] = i1;
				JA[ind + nnzr*2 + 1 + j*2 + 1] = i2;
				VA[ind + j*2] = SMl->get_stiffness_matrix_12(0);
				VA[ind + j*2 + 1] = SMl->get_stiffness_matrix_12(1);
				VA[ind + nnzr*2 + 1 + j*2] = SMl->get_stiffness_matrix_12(2);;
				VA[ind + nnzr*2 + 1 + j*2 + 1] = SMl->get_stiffness_matrix_12(3);
				/*if (number==21)
				{
					if (VA[ind + j*2]>1000000 || VA[ind + j*2]<-1000000 ||
						VA[ind + j*2 + 1]>1000000 || VA[ind + j*2 + 1]<-1000000 || 
						VA[ind + nnzr*2 + 1 + j*2]>1000000 || VA[ind + nnzr*2 + 1 + j*2]<-1000000 ||
						VA[ind + nnzr*2 + 1 + j*2 + 1]>1000000 || VA[ind + nnzr*2 + 1 + j*2 + 1]<-1000000)
					{
						if (SMl->getN1()->get_cell_boundary())
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN1()->getab()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						else
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_SM_ID()<<" of cell "<<SMl->getN1()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						if (SMl->getN2()->get_cell_boundary())
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN2()->getab()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<": ";
							if (SMl->getN2()->getmp()->size()>0) {cout<<SMl->getN2()->getmp()->at(0)->getab()->Cell()->id()<<endl;}
							else {cout<<endl;}
						}
						else
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_SM_ID()<<" of cell "<<SMl->getN2()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<endl;
						}
						cout<<" triangle NULL N1-N2: 00: "<<SMl->get_stiffness_matrix_12(0)<<endl;
						cout<<" triangle NULL N1-N2: 01: "<<SMl->get_stiffness_matrix_12(1)<<endl;
						cout<<" triangle NULL N1-N2: 10: "<<SMl->get_stiffness_matrix_12(2)<<endl;
						cout<<" triangle NULL N1-N2: 11: "<<SMl->get_stiffness_matrix_12(3)<<endl;
					}
				}*/
			}
			else
			{
				int i1 = sm_id1*2;
				int i2 = sm_id1*2+1;
				JA[ind + j*2] = i1;
				JA[ind + j*2 + 1] = i2;
				JA[ind + nnzr*2 + 1 + j*2] = i1;
				JA[ind + nnzr*2 + 1 + j*2 + 1] = i2;
				VA[ind + j*2] = SMl->get_stiffness_matrix_21(0);
				VA[ind + j*2 + 1] = SMl->get_stiffness_matrix_21(1);
				VA[ind + nnzr*2 + 1 + j*2] = SMl->get_stiffness_matrix_21(2);;
				VA[ind + nnzr*2 + 1 + j*2 + 1] = SMl->get_stiffness_matrix_21(3);
				/*if (number==21)
				{
					if (VA[ind + j*2]>1000000 || VA[ind + j*2]<-1000000 ||
						VA[ind + j*2 + 1]>1000000 || VA[ind + j*2 + 1]<-1000000 || 
						VA[ind + nnzr*2 + 1 + j*2]>1000000 || VA[ind + nnzr*2 + 1 + j*2]<-1000000 ||
						VA[ind + nnzr*2 + 1 + j*2 + 1]>1000000 || VA[ind + nnzr*2 + 1 + j*2 + 1]<-1000000)
					{
						if (SMl->getN1()->get_cell_boundary())
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN1()->getab()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						else
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_SM_ID()<<" of cell "<<SMl->getN1()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						if (SMl->getN2()->get_cell_boundary())
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN2()->getab()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<": ";
							if (SMl->getN2()->getmp()->size()>0) {cout<<SMl->getN2()->getmp()->at(0)->getab()->Cell()->id()<<endl;}
							else {cout<<endl;}
						}
						else
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_SM_ID()<<" of cell "<<SMl->getN2()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<endl;							
						}
						cout<<" triangle NULL N1-N2: 00: "<<SMl->get_stiffness_matrix_12(0)<<endl;
						cout<<" triangle NULL N1-N2: 01: "<<SMl->get_stiffness_matrix_12(1)<<endl;
						cout<<" triangle NULL N1-N2: 10: "<<SMl->get_stiffness_matrix_12(2)<<endl;
						cout<<" triangle NULL N1-N2: 11: "<<SMl->get_stiffness_matrix_12(3)<<endl;
					}
				}*/
			}
		}
		JA[ind + nnzr*2] = sm_ind*2 + 1;
		JA[ind + nnzr*4 + 1] = sm_ind*2;
		VA[ind + nnzr*2] = nodepoints[i]->getA_ij(1);
		VA[ind + nnzr*4 + 1] = nodepoints[i]->getA_ij(2);
	}
	t_end2 = clock();
	cout<<"Stiffness matrix build up: "<<(int)cellList.size()<<": "<<(t_end2-t_begin2)*0.001<<" seconds!"<<endl;
	////////// SOR linear solver //////////
	time_t t_begin, t_end, t_total=0;
	t_begin = clock();

	double converge = 5.0e-3;
	double w = 1.5;
	double diff = 100;
	int flagtime = 0; // iteration time
	int runningstep = 500;
	diff = SOR_solver(sm_idt, converge, w, runningstep, flagtime, SOR_N, JA, VA, Fn, Xn, X, SOR_II_V);
	t_end = clock();
	cout<<"Running SOR: Steps: "<<flagtime<<", Time: "<<(t_end - t_begin)*0.001<<" seconds! Maximum difference: "<<diff<<endl;
	/// without area refinement force
	diff = 100;
	flagtime = 0;
	diff = SOR_solver(sm_idt, converge, w, runningstep, flagtime, SOR_N, JA, VA, Fn_p, Xn_p, X_p, SOR_II_V);
	/**********************************************
	 update node position:
	 1. first update cell boundary vertices
	 2. second update the interior vertices
	**********************************************/
	for (int i=0;i<sm_id;i++)
	{
		if (nodepoints[i]->get_cell_boundary())
		{
			node *Nod = nodepoints[i]->get_node_id();
			int nvn = (int)Nod->getvertex()->size();
			int ida = nodepoints[i]->get_SM_ID();
			double xn = nodepoints[i]->x() + Xn[ida*2];
			double yn = nodepoints[i]->y() + Xn[ida*2+1];
			for (int j=0;j<nvn;j++)
			{
				point3D *pnj = Nod->getvertex()->at(j);
				//if (pnj->getab()->Cell()->id()==0)
				//{
				//	cout<<"Fn, Xn "<<pnj->id()<<": "<<Fn[ida*2]<<","<<Fn[ida*2+1]<<";"<<Xn[ida*2]<<","<<Xn[ida*2+1]<<endl;
				//}
				pnj->setx(xn);
				pnj->sety(yn);
			}
		}
		else
		{
			int ncellt = nodepoints[i]->Cell()->get_cell_type();
			int ida = nodepoints[i]->get_SM_ID();
			double xn = nodepoints[i]->x() + Xn[ida*2];
			double yn = nodepoints[i]->y() + Xn[ida*2+1];
			nodepoints[i]->setx(xn);
			nodepoints[i]->sety(yn);
		}
	}
	/********************************************
	 Displacement Relaxation
	 -------------------------------------------
	 Slight refinement for the system. To ensure 
	 each cell still preserves approximatedly the 
	 same area. Accuracy depends on number of 
	 running steps. The default number of steps 
	 is set to 3.
	********************************************/
	time_t t_b_rel, t_e_rel;
	t_b_rel = clock();
	int flag_relax = 0;
	default_area = 6*GAR/12;
	double *F_relax = new double[sm_id*2];
	double *X_relax = new double[sm_id*2];
	double *X_relax_n = new double[sm_id*2];
	double retrov[2];retrov[0] = 0;retrov[1] = 0;
	int retrovn = 0;
	while (flag_relax<3)
	{
		flag_relax++;
		for (int i=0;i<cn;i++)
		{
			if (cellList[i]->Dead()) continue;
			cellList[i]->set_area();
		}
		for (int i=0;i<sm_id;i++)
		{
			X_relax[i*2] = 0.1;
			X_relax_n[i*2] = 0.1;
			X_relax[i*2+1] = 0.1;
			X_relax_n[i*2+1] = 0.1;
			if (nodepoints[i]->get_cell_boundary())
			{
				double fx = 0;
				double fy = 0;
				node *Nod = nodepoints[i]->get_node_id();
				int nvn = (int)Nod->getvertex()->size();
				for (int j=0;j<nvn;j++)
				{
					point3D *pt = Nod->getvertex()->at(j);
					cell *Cpt = pt->getab()->Cell();
					if (Cpt->Soften()) continue;
					if (Cpt->Migrate_response()) continue;
					int csn = (int)Cpt->get_sides()->size();
					double P_scale = 0;
					double inc_areaT = Cpt->get_GR();
					if (inc_areaT==0 || pt->get_migration_rate()==0)
					{
						double cell_area = Cpt->get_area();
						double cell_initial_area = Cpt->get_initial_area();
						double area_diff = cell_initial_area - cell_area;
						if      (area_diff>default_area)  {area_diff = default_area;}
						else if (area_diff<-default_area) {area_diff = -default_area;}
						area_diff = area_diff/csn;
						//P_scale = 4*area_diff*Cpt->get_mass()/TM/TM;
						P_scale = area_diff/BD*(Cpt->get_lame(0)*2 + Cpt->get_lame(1))/BD;
						double ab_normal_x = pt->getab()->getnormal(0);
						double ab_normal_y = pt->getab()->getnormal(1);
						double ba_normal_x = pt->getba()->getnormal(0);
						double ba_normal_y = pt->getba()->getnormal(1);
						double ab_d = sqrt(ab_normal_x*ab_normal_x + ab_normal_y*ab_normal_y);
						double ba_d = sqrt(ba_normal_x*ba_normal_x + ba_normal_y*ba_normal_y);
						double inc_P_ab_x = ab_normal_x/ab_d*P_scale;
						double inc_P_ab_y = ab_normal_y/ab_d*P_scale;
						double inc_P_ba_x = ba_normal_x/ba_d*P_scale;
						double inc_P_ba_y = ba_normal_y/ba_d*P_scale;
						fx += inc_P_ab_x;
						fx += inc_P_ba_x;
						fy += inc_P_ab_y;
						fy += inc_P_ba_y;
					}
					if (flag_relax==2)
					{
						pt->set_migration_rate(0); // reset to NULL
					}
				}
				F_relax[i*2] = fx;
				F_relax[i*2+1] = fy;
			}
			else
			{
				F_relax[i*2] = 0;
				F_relax[i*2+1] = 0;
			}
		}
		diff = 100;
		flagtime = 0;
		converge = 5.0e-3;
		w = 1.0;
		flagtime = 0;
		runningstep = 500;
		diff = SOR_solver(sm_idt, converge, w, runningstep, flagtime, SOR_N, JA, VA, F_relax, X_relax_n, X_relax, SOR_II_V);
		cout<<"Relaxation step "<<flag_relax<<": "<<flagtime<<" steps; maxium difference: "<<diff<<endl;
		for (int i=0;i<sm_id;i++)
		{
			if (nodepoints[i]->get_cell_boundary())
			{
				node *Nod = nodepoints[i]->get_node_id();
				int nvn = (int)Nod->getvertex()->size();
				int ida = nodepoints[i]->get_SM_ID();
				double xn = nodepoints[i]->x() + X_relax_n[ida*2];
				double yn = nodepoints[i]->y() + X_relax_n[ida*2+1];
				for (int j=0;j<nvn;j++)
				{
					point3D *pnj = Nod->getvertex()->at(j);
					pnj->setx(xn);
					pnj->sety(yn);
				}
				if (flag_relax==2)
				{
					bool testfix = false;
					for (int j=0;j<nvn;j++)
					{
						int ncellt = Nod->getvertex()->at(j)->getab()->Cell()->get_cell_type();
						int v_fix = fixed_cells[ncellt];
						if (v_fix==1) {testfix = true;break;}
					}
					if (testfix)
					{
						retrov[0] += xn - nodepoints[i]->x_old(3);
						retrov[1] += yn - nodepoints[i]->y_old(3);
						retrovn++;
					}
				}
			}
			else
			{
				int ncellt = nodepoints[i]->Cell()->get_cell_type();
				int v_fix = fixed_cells[ncellt];
				int ida = nodepoints[i]->get_SM_ID();
				double xn = nodepoints[i]->x() + X_relax_n[ida*2];
				double yn = nodepoints[i]->y() + X_relax_n[ida*2+1];
				nodepoints[i]->setx(xn);
				nodepoints[i]->sety(yn);
			}
		}
	}
	if (retrovn>0)
	{
		retrov[0] = retrov[0]/retrovn;
		retrov[1] = retrov[1]/retrovn;
		cout<<"  -> average displacement vector of static cells: "<<retrov[0]<<","<<retrov[1]<<endl;
	}
	///// update the tension and pressure /////
	for (int i=0;i<sm_id;i++)
	{
		if (nodepoints[i]->get_cell_boundary())
		{
			node *Nod = nodepoints[i]->get_node_id();
			int nvn = (int)Nod->getvertex()->size();
			for (int j=0;j<nvn;j++)
			{
				point3D *pt = Nod->getvertex()->at(j);
				cell *Cpt = pt->getab()->Cell();
				double ab_normal_x = pt->getab()->getnormal(0);
				double ab_normal_y = pt->getab()->getnormal(1);
				double ba_normal_x = pt->getba()->getnormal(0);
				double ba_normal_y = pt->getba()->getnormal(1);
				double ab_d = sqrt(ab_normal_x*ab_normal_x + ab_normal_y*ab_normal_y);
				double ba_d = sqrt(ba_normal_x*ba_normal_x + ba_normal_y*ba_normal_y);
				double ab_l = pt->getab()->getlength();
				double ba_l = pt->getba()->getlength();
				double T_ab_x = (pt->getfp()->x() - pt->x())*Cpt->get_tcoef();
				double T_ab_y = (pt->getfp()->y() - pt->y())*Cpt->get_tcoef();
				double T_ba_x = (pt->getrp()->x() - pt->x())*Cpt->get_tcoef();
				double T_ba_y = (pt->getrp()->y() - pt->y())*Cpt->get_tcoef();
				double P_ab_x = 0.5*ab_normal_x/ab_d*ab_l*Cpt->get_pcoef();
				double P_ab_y = 0.5*ab_normal_y/ab_d*ab_l*Cpt->get_pcoef();
				double P_ba_x = 0.5*ba_normal_x/ab_d*ba_l*Cpt->get_pcoef();
				double P_ba_y = 0.5*ba_normal_y/ab_d*ba_l*Cpt->get_pcoef();
				pt->setTension_ab(T_ab_x, T_ab_y);
				pt->setTension_ba(T_ba_x, T_ba_y);
				pt->setPressure_ab(P_ab_x, P_ab_y);
				pt->setPressure_ba(P_ba_x, P_ba_y);
				// add correction vector
				double ptx = pt->x();
				double pty = pt->y();
				pt->setx(ptx - retrov[0]);
				pt->sety(pty - retrov[1]);
			}
		}
	}
	t_e_rel = clock();
	cout<<"displacement relaxation done: "<<(t_e_rel-t_b_rel)*0.001<<" seconds!\n";
	/********************************************
	////////// internal force recovery //////////
	********************************************/
	/////////////////////////////
	int sense_migration_flag = 0; // 0: regardless of adhesion angle between cells; 1: regard to adhesion angle between cells;
	/////////////////////////////
	time_t t_b_re, t_e_re;
	t_b_re = clock();
	vector<edge*> edgebreak0;
	vector<edge*> edgebreak1;
	vector<point3D*> pointbreak0;
	vector<point3D*> pointbreak1;
	vector<point3D*> pointbreak2;
	int copln = (int)collisionpairList.size();
	for (int i=0;i<copln;i++)
	{
		if (collisionpairList[i]->get_redundant()) continue;
		if (collisionpairList[i]->MP()->size()==0) continue;
		bool keeprun = false;
		cell *C11 = collisionpairList[i]->get_I1();
		cell *C12 = collisionpairList[i]->get_I2();
		if (C11->Migrate() ||
			C12->Migrate() ||
			C11->Soften() ||
			C12->Soften()) {keeprun = true;}
		if (!keeprun) continue;
		//cout<<"collisionpair "<<i<<":"<<collisionpairList[i]->get_I1()->id()<<" "<<collisionpairList[i]->get_I2()->id()<<": "<<endl;
		mergepair *Mp = collisionpairList[i]->MP()->at(0);
		double F_adhesion = collisionpairList[i]->get_adhesion();
		point3D *s11 = Mp->S1();
		point3D *t11 = Mp->T1();
		point3D *s21 = Mp->S2();
		point3D *t21 = Mp->T2();
		double Fs11x = 0;double Fs11y = 0;
		double Ft11x = 0;double Ft11y = 0;
		double Fs21x = 0;double Fs21y = 0;
		double Ft21x = 0;double Ft21y = 0;
		double Fs11x_m = 0;double Fs11y_m = 0;
		double Ft11x_m = 0;double Ft11y_m = 0;
		double Fs21x_m = 0;double Fs21y_m = 0;
		double Ft21x_m = 0;double Ft21y_m = 0;
		/*****************************************
		case 1: Only 2 points join together. Do 
		        force recovery within the loop

			o-<-oo-<-o
			  t1||s2
				oo
			  s1||t2
			o->-oo->-o

		-----------------------------------------
		case 2: More than 2 points join together, 
				but not fill up the around angle.
				Do force recovery later after the 
				loop

		    o-<-o-<-ooo-<-o-<-o
				t1 // \\ s2  
		          oo   oo
		      s1 //     \\ t2
			o->-oo       oo->-o

		-----------------------------------------
		case 3: More than 2 points join together,
		        and fill up the around angle.
				Do force recovery later after the 
				loop

              o--<<--o
              |      |
		 o-<-oo      oo-<-o
			  \\    //
			   oo  oo
			    \\//
			  t1 oo s2
				 ||
				 oo
			  s1 || t2
			 o->-oo->-o
		*****************************************/
		if      (s11==t11)
		{
			//cout<<"s11==t11: collisionpair:"<<collisionpairList[i]->get_I1()->id()<<" "<<collisionpairList[i]->get_I2()->id()<<": "<<endl;
			edge *s11_ba = s11->getba();
			if (!s11_ba->get_attach())
			{
				double vs20 = s11_ba->get_virtual_triangle_stiffness_matrix(2,0);
				double vs21 = s11_ba->get_virtual_triangle_stiffness_matrix(2,1);
				double vs22 = s11_ba->get_virtual_triangle_stiffness_matrix(2,2);
				double vs23 = s11_ba->get_virtual_triangle_stiffness_matrix(2,3);
				double vs30 = s11_ba->get_virtual_triangle_stiffness_matrix(3,0);
				double vs31 = s11_ba->get_virtual_triangle_stiffness_matrix(3,1);
				double vs32 = s11_ba->get_virtual_triangle_stiffness_matrix(3,2);
				double vs33 = s11_ba->get_virtual_triangle_stiffness_matrix(3,3);
				int p1s_sm_id = s11->getrp()->get_node_id()->SM_ID();
				int p2s_sm_id = s11->get_node_id()->SM_ID();
				double us0 = Xn[p1s_sm_id*2];
				double us1 = Xn[p1s_sm_id*2+1];
				double us2 = Xn[p2s_sm_id*2];
				double us3 = Xn[p2s_sm_id*2+1];
				double us0_m = Xn_p[p1s_sm_id*2];
				double us1_m = Xn_p[p1s_sm_id*2+1];
				double us2_m = Xn_p[p2s_sm_id*2];
				double us3_m = Xn_p[p2s_sm_id*2+1];
				Fs11x += vs20*us0 + vs21*us1 + vs22*us2 + vs23*us3;
				Fs11y += vs30*us0 + vs31*us1 + vs32*us2 + vs33*us3;
				Fs11x_m += vs20*us0_m + vs21*us1_m + vs22*us2_m + vs23*us3_m;
				Fs11y_m += vs30*us0_m + vs31*us1_m + vs32*us2_m + vs33*us3_m;
			}
			if (!C11->Soften())
			{
				int at_ns = (int)s11->get_inner_at()->size();
				for (int j=0;j<at_ns;j++)
				{
					triangle *Tj = s11->get_inner_at()->at(j);
					point3D *TAj = Tj->getA();
					point3D *TBj = Tj->getB();
					point3D *TCj = Tj->getC();
					int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
					int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
					int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
					double u0 = Xn[TAj_sm_id*2];
					double u1 = Xn[TAj_sm_id*2+1];
					double u2 = Xn[TBj_sm_id*2];
					double u3 = Xn[TBj_sm_id*2+1];
					double u4 = Xn[TCj_sm_id*2];
					double u5 = Xn[TCj_sm_id*2+1];
					double u0_m = Xn_p[TAj_sm_id*2];
					double u1_m = Xn_p[TAj_sm_id*2+1];
					double u2_m = Xn_p[TBj_sm_id*2];
					double u3_m = Xn_p[TBj_sm_id*2+1];
					double u4_m = Xn_p[TCj_sm_id*2];
					double u5_m = Xn_p[TCj_sm_id*2+1];
					if      (TAj==s11)
					{
						double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
						double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
						double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
						double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
						double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
						double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
						Fs11x += FAx;
						Fs11y += FAy;
						Fs11x_m += FAx_m;
						Fs11y_m += FAy_m;
					}
					else if (TBj==s11)
					{
						double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
						double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
						double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
						double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
						double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
						double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
						Fs11x += FBx;
						Fs11y += FBy;
						Fs11x_m += FBx_m;
						Fs11y_m += FBy_m;
					}
					else if (TCj==s11)
					{
						double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
						double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
						double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
						double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
						double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
						double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
						Fs11x += FCx;
						Fs11y += FCy;
						Fs11x_m += FCx_m;
						Fs11y_m += FCy_m;
					}
				}
			}
			else
			{
				int N1_sm_id = s11->get_node_id()->SM_ID();
				int N2_sm_id = s11->getrp()->get_node_id()->SM_ID();
				int N3_sm_id = s11->getfp()->get_node_id()->SM_ID();
				double u0 = Xn[N1_sm_id*2];
				double u1 = Xn[N1_sm_id*2+1];
				double u2 = Xn[N2_sm_id*2];
				double u3 = Xn[N2_sm_id*2+1];
				double u4 = Xn[N3_sm_id*2];
				double u5 = Xn[N3_sm_id*2+1];
				double u0_m = Xn_p[N1_sm_id*2];
				double u1_m = Xn_p[N1_sm_id*2+1];
				double u2_m = Xn_p[N2_sm_id*2];
				double u3_m = Xn_p[N2_sm_id*2+1];
				double u4_m = Xn_p[N3_sm_id*2];
				double u5_m = Xn_p[N3_sm_id*2+1];
				double v00 = s11->getab()->get_soften_triangle_stiffness_matrix(0,0);
				double v01 = s11->getab()->get_soften_triangle_stiffness_matrix(0,1);
				double v02 = s11->getab()->get_soften_triangle_stiffness_matrix(0,2);
				double v03 = s11->getab()->get_soften_triangle_stiffness_matrix(0,3);
				double v10 = s11->getab()->get_soften_triangle_stiffness_matrix(1,0);
				double v11 = s11->getab()->get_soften_triangle_stiffness_matrix(1,1);
				double v12 = s11->getab()->get_soften_triangle_stiffness_matrix(1,2);
				double v13 = s11->getab()->get_soften_triangle_stiffness_matrix(1,3);
				double v20 = s11->getba()->get_soften_triangle_stiffness_matrix(2,0);
				double v21 = s11->getba()->get_soften_triangle_stiffness_matrix(2,1);
				double v22 = s11->getba()->get_soften_triangle_stiffness_matrix(2,2);
				double v23 = s11->getba()->get_soften_triangle_stiffness_matrix(2,3);
				double v30 = s11->getba()->get_soften_triangle_stiffness_matrix(3,0);
				double v31 = s11->getba()->get_soften_triangle_stiffness_matrix(3,1);
				double v32 = s11->getba()->get_soften_triangle_stiffness_matrix(3,2);
				double v33 = s11->getba()->get_soften_triangle_stiffness_matrix(3,3);
				double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
				double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
				double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
				double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
				double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
				double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
				double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
				double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
				Fs11x += (Frx + Ffx);
				Fs11y += (Fry + Ffy);
				Fs11x_m += (Frx_m + Ffx_m);
				Fs11y_m += (Fry_m + Ffy_m);
			}
			edge *t21_ab = t21->getab();
			if (!t21_ab->get_attach())
			{
				double vt00 = t21_ab->get_virtual_triangle_stiffness_matrix(0,0);
				double vt01 = t21_ab->get_virtual_triangle_stiffness_matrix(0,1);
				double vt02 = t21_ab->get_virtual_triangle_stiffness_matrix(0,2);
				double vt03 = t21_ab->get_virtual_triangle_stiffness_matrix(0,3);
				double vt10 = t21_ab->get_virtual_triangle_stiffness_matrix(1,0);
				double vt11 = t21_ab->get_virtual_triangle_stiffness_matrix(1,1);
				double vt12 = t21_ab->get_virtual_triangle_stiffness_matrix(1,2);
				double vt13 = t21_ab->get_virtual_triangle_stiffness_matrix(1,3);
				int p1t_sm_id = t21->get_node_id()->SM_ID();
				int p2t_sm_id = t21->getfp()->get_node_id()->SM_ID();
				double ut0 = Xn[p1t_sm_id*2];
				double ut1 = Xn[p1t_sm_id*2+1];
				double ut2 = Xn[p2t_sm_id*2];
				double ut3 = Xn[p2t_sm_id*2+1];
				double ut0_m = Xn_p[p1t_sm_id*2];
				double ut1_m = Xn_p[p1t_sm_id*2+1];
				double ut2_m = Xn_p[p2t_sm_id*2];
				double ut3_m = Xn_p[p2t_sm_id*2+1];
				Ft21x += vt00*ut0 + vt01*ut1 + vt02*ut2 + vt03*ut3;
				Ft21y += vt10*ut0 + vt11*ut1 + vt12*ut2 + vt13*ut3;
				Ft21x_m += vt00*ut0_m + vt01*ut1_m + vt02*ut2_m + vt03*ut3_m;
				Ft21y_m += vt10*ut0_m + vt11*ut1_m + vt12*ut2_m + vt13*ut3_m;
			}
			if (!C12->Soften())
			{
				int at_nt = (int)t21->get_inner_at()->size();
				for (int j=0;j<at_nt;j++)
				{
					triangle *Tj = t21->get_inner_at()->at(j);
					point3D *TAj = Tj->getA();
					point3D *TBj = Tj->getB();
					point3D *TCj = Tj->getC();
					//cout<<"Triangle check: "<<t21->id()<<" endpoints: "<<TAj->id()<<" "<<TBj->id()<<" "<<TCj->id()<<endl;
					int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
					int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
					int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
					double u0 = Xn[TAj_sm_id*2];
					double u1 = Xn[TAj_sm_id*2+1];
					double u2 = Xn[TBj_sm_id*2];
					double u3 = Xn[TBj_sm_id*2+1];
					double u4 = Xn[TCj_sm_id*2];
					double u5 = Xn[TCj_sm_id*2+1];
					double u0_m = Xn_p[TAj_sm_id*2];
					double u1_m = Xn_p[TAj_sm_id*2+1];
					double u2_m = Xn_p[TBj_sm_id*2];
					double u3_m = Xn_p[TBj_sm_id*2+1];
					double u4_m = Xn_p[TCj_sm_id*2];
					double u5_m = Xn_p[TCj_sm_id*2+1];
					if      (TAj==t21)
					{
						double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
						double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
						double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
						double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
						double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
						double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
						Ft21x += FAx;
						Ft21y += FAy;
						Ft21x_m += FAx_m;
						Ft21y_m += FAy_m;
					}
					else if (TBj==t21)
					{
						double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
						double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
						double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
						double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
						double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
						double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
						Ft21x += FBx;
						Ft21y += FBy;
						Ft21x_m += FBx_m;
						Ft21y_m += FBy_m;
					}
					else if (TCj==t21)
					{
						double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
						double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
						double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
						double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
						double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
						double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
						Ft21x += FCx;
						Ft21y += FCy;
						Ft21x_m += FCx_m;
						Ft21y_m += FCy_m;
					}
				}
			}
			else
			{
				int N1_sm_id = t21->get_node_id()->SM_ID();
				int N2_sm_id = t21->getrp()->get_node_id()->SM_ID();
				int N3_sm_id = t21->getfp()->get_node_id()->SM_ID();
				double u0 = Xn[N1_sm_id*2];
				double u1 = Xn[N1_sm_id*2+1]; // p
				double u2 = Xn[N2_sm_id*2];
				double u3 = Xn[N2_sm_id*2+1]; // rp
				double u4 = Xn[N3_sm_id*2];
				double u5 = Xn[N3_sm_id*2+1]; // fp
				double u0_m = Xn_p[N1_sm_id*2];
				double u1_m = Xn_p[N1_sm_id*2+1]; // p
				double u2_m = Xn_p[N2_sm_id*2];
				double u3_m = Xn_p[N2_sm_id*2+1]; // rp
				double u4_m = Xn_p[N3_sm_id*2];
				double u5_m = Xn_p[N3_sm_id*2+1]; // fp
				double v00 = t21->getab()->get_soften_triangle_stiffness_matrix(0,0);
				double v01 = t21->getab()->get_soften_triangle_stiffness_matrix(0,1);
				double v02 = t21->getab()->get_soften_triangle_stiffness_matrix(0,2);
				double v03 = t21->getab()->get_soften_triangle_stiffness_matrix(0,3);
				double v10 = t21->getab()->get_soften_triangle_stiffness_matrix(1,0);
				double v11 = t21->getab()->get_soften_triangle_stiffness_matrix(1,1);
				double v12 = t21->getab()->get_soften_triangle_stiffness_matrix(1,2);
				double v13 = t21->getab()->get_soften_triangle_stiffness_matrix(1,3);
				double v20 = t21->getba()->get_soften_triangle_stiffness_matrix(2,0);
				double v21 = t21->getba()->get_soften_triangle_stiffness_matrix(2,1);
				double v22 = t21->getba()->get_soften_triangle_stiffness_matrix(2,2);
				double v23 = t21->getba()->get_soften_triangle_stiffness_matrix(2,3);
				double v30 = t21->getba()->get_soften_triangle_stiffness_matrix(3,0);
				double v31 = t21->getba()->get_soften_triangle_stiffness_matrix(3,1);
				double v32 = t21->getba()->get_soften_triangle_stiffness_matrix(3,2);
				double v33 = t21->getba()->get_soften_triangle_stiffness_matrix(3,3);
				double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
				double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
				double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
				double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
				double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
				double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
				double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
				double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
				Ft21x += (Frx + Ffx);
				Ft21y += (Fry + Ffy);
				Ft21x_m += (Frx_m + Ffx_m);
				Ft21y_m += (Fry_m + Ffy_m);
			}
			edge *s21_ba = s21->getba();
			if (!s21_ba->get_attach())
			{
				double vs20 = s21_ba->get_virtual_triangle_stiffness_matrix(2,0);
				double vs21 = s21_ba->get_virtual_triangle_stiffness_matrix(2,1);
				double vs22 = s21_ba->get_virtual_triangle_stiffness_matrix(2,2);
				double vs23 = s21_ba->get_virtual_triangle_stiffness_matrix(2,3);
				double vs30 = s21_ba->get_virtual_triangle_stiffness_matrix(3,0);
				double vs31 = s21_ba->get_virtual_triangle_stiffness_matrix(3,1);
				double vs32 = s21_ba->get_virtual_triangle_stiffness_matrix(3,2);
				double vs33 = s21_ba->get_virtual_triangle_stiffness_matrix(3,3);
				int p1s_sm_id = s21->getrp()->get_node_id()->SM_ID();
				int p2s_sm_id = s21->get_node_id()->SM_ID();
				double us0 = Xn[p1s_sm_id*2];
				double us1 = Xn[p1s_sm_id*2+1];
				double us2 = Xn[p2s_sm_id*2];
				double us3 = Xn[p2s_sm_id*2+1];
				double us0_m = Xn_p[p1s_sm_id*2];
				double us1_m = Xn_p[p1s_sm_id*2+1];
				double us2_m = Xn_p[p2s_sm_id*2];
				double us3_m = Xn_p[p2s_sm_id*2+1];
				Fs21x += vs20*us0 + vs21*us1 + vs22*us2 + vs23*us3;
				Fs21y += vs30*us0 + vs31*us1 + vs32*us2 + vs33*us3;
				Fs21x_m += vs20*us0_m + vs21*us1_m + vs22*us2_m + vs23*us3_m;
				Fs21y_m += vs30*us0_m + vs31*us1_m + vs32*us2_m + vs33*us3_m;
			}
			edge *t11_ab = t11->getab();
			if (!t11_ab->get_attach())
			{
				double vt00 = t11_ab->get_virtual_triangle_stiffness_matrix(0,0);
				double vt01 = t11_ab->get_virtual_triangle_stiffness_matrix(0,1);
				double vt02 = t11_ab->get_virtual_triangle_stiffness_matrix(0,2);
				double vt03 = t11_ab->get_virtual_triangle_stiffness_matrix(0,3);
				double vt10 = t11_ab->get_virtual_triangle_stiffness_matrix(1,0);
				double vt11 = t11_ab->get_virtual_triangle_stiffness_matrix(1,1);
				double vt12 = t11_ab->get_virtual_triangle_stiffness_matrix(1,2);
				double vt13 = t11_ab->get_virtual_triangle_stiffness_matrix(1,3);
				int p1t_sm_id = t11->get_node_id()->SM_ID();
				int p2t_sm_id = t11->getfp()->get_node_id()->SM_ID();
				double ut0 = Xn[p1t_sm_id*2];
				double ut1 = Xn[p1t_sm_id*2+1];
				double ut2 = Xn[p2t_sm_id*2];
				double ut3 = Xn[p2t_sm_id*2+1];
				double ut0_m = Xn_p[p1t_sm_id*2];
				double ut1_m = Xn_p[p1t_sm_id*2+1];
				double ut2_m = Xn_p[p2t_sm_id*2];
				double ut3_m = Xn_p[p2t_sm_id*2+1];
				Ft11x += vt00*ut0 + vt01*ut1 + vt02*ut2 + vt03*ut3;
				Ft11y += vt10*ut0 + vt11*ut1 + vt12*ut2 + vt13*ut3;
				Ft11x_m += vt00*ut0_m + vt01*ut1_m + vt02*ut2_m + vt03*ut3_m;
				Ft11y_m += vt10*ut0_m + vt11*ut1_m + vt12*ut2_m + vt13*ut3_m;
			}
			//cout<<s11->id()<<" "<<s21->id()<<": "<<Fs11x+Ft11x<<", "<<Fs11y+Ft11y<<"; "<<Fs21x+Ft21x<<", "<<Fs21y+Ft21y<<endl;
			s11->set_elastic_force(Fs11x+Ft11x, Fs11y+Ft11y);
			s21->set_elastic_force(Fs21x+Ft21x, Fs21y+Ft21y);
			int s11mn = (int)s11->getmp()->size();
			if (s11mn==1)
			{
				double ns1_abx = s11->getab()->getnormal(0);
				double ns1_aby = s11->getab()->getnormal(1);
				double ns1_abd = sqrt(ns1_abx*ns1_abx + ns1_aby*ns1_aby);
				double ns1_ab_normal_x = ns1_abx/ns1_abd;
				double ns1_ab_normal_y = ns1_aby/ns1_abd;
				double ns1_bax = s11->getba()->getnormal(0);
				double ns1_bay = s11->getba()->getnormal(1);
				double ns1_bad = sqrt(ns1_bax*ns1_bax + ns1_bay*ns1_bay);
				double ns1_ba_normal_x = ns1_bax/ns1_bad;
				double ns1_ba_normal_y = ns1_bay/ns1_bad;
				double ns1_cox = ns1_ab_normal_x + ns1_ba_normal_x;
				double ns1_coy = ns1_ab_normal_y + ns1_ba_normal_y;
				double ns1_cod = sqrt(ns1_cox*ns1_cox + ns1_coy*ns1_coy);
				double ns1_co_normal_x = ns1_cox/ns1_cod;
				double ns1_co_normal_y = ns1_coy/ns1_cod;
				double ns2_abx = s21->getab()->getnormal(0);
				double ns2_aby = s21->getab()->getnormal(1);
				double ns2_abd = sqrt(ns2_abx*ns2_abx + ns2_aby*ns2_aby);
				double ns2_ab_normal_x = ns2_abx/ns2_abd;
				double ns2_ab_normal_y = ns2_aby/ns2_abd;
				double ns2_bax = s21->getba()->getnormal(0);
				double ns2_bay = s21->getba()->getnormal(1);
				double ns2_bad = sqrt(ns2_bax*ns2_bax + ns2_bay*ns2_bay);
				double ns2_ba_normal_x = ns2_bax/ns2_bad;
				double ns2_ba_normal_y = ns2_bay/ns2_bad;
				double ns2_cox = ns2_ab_normal_x + ns2_ba_normal_x;
				double ns2_coy = ns2_ab_normal_y + ns2_ba_normal_y;
				double ns2_cod = sqrt(ns2_cox*ns2_cox + ns2_coy*ns2_coy);
				double ns2_co_normal_x = ns2_cox/ns2_cod;
				double ns2_co_normal_y = ns2_coy/ns2_cod;
				bool s11_vert = false;
				if      ((s11->get_pop_out()>0 && t21->get_pop_out()>0)  ||
						 (s11->get_pop_out()>0 && t21->get_pop_out()==0) ||
						 (s11->get_pop_out()==0 && t21->get_pop_out()>0))
				{
					Mp->set_F_s1(0);Mp->set_F_t2(0);
					Mp->set_F_s2(0);Mp->set_F_t1(0);
				}
				else if (s11->get_pop_out()>0 && t21->get_pop_out()<0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double plot_s1 = abs((Fs11x + Ft11x)*ns1_co_normal_x + (Fs11y + Ft11y)*ns1_co_normal_y);
					double plot_t2 = abs((Ft21x + Fs21x)*ns2_co_normal_x + (Ft21y + Fs21y)*ns2_co_normal_y);
					double sF_s1_new = sF_s1 + plot_s1;
					double sF_t2_new = sF_t2 + plot_t2;
					double diff_t2_s1 = sF_t2_new - sF_s1_new;
					if (diff_t2_s1>0) {sF_t2_new = diff_t2_s1;sF_s1_new = diff_t2_s1;}
					else {sF_t2_new = 0;sF_s1_new = 0;}
					Mp->set_F_s1(sF_s1_new);Mp->set_F_t1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);Mp->set_F_s2(sF_t2_new);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_vert = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					if (diff_t2_s1>0)
					{
						double Ft21_d = sqrt((Ft21x + Fs21x)*(Ft21x + Fs21x) + (Ft21y + Fs21y)*(Ft21y + Fs21y));
						C1_sense_mech_vec_x += diff_t2_s1*(Ft21x + Fs21x)/Ft21_d;
						C1_sense_mech_vec_y += diff_t2_s1*(Ft21y + Fs21y)/Ft21_d;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()<0 && t21->get_pop_out()>0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double plot_s1 = abs((Fs11x + Ft11x)*ns1_co_normal_x + (Fs11y + Ft11y)*ns1_co_normal_y);
					double plot_t2 = abs((Ft21x + Fs21x)*ns2_co_normal_x + (Ft21y + Fs21y)*ns2_co_normal_y);
					double sF_s1_new = sF_s1 + plot_s1;
					double sF_t2_new = sF_t2 + plot_t2;
					double diff_s1_t2 = sF_s1_new - sF_t2_new;
					if (diff_s1_t2>0) {sF_t2_new = diff_s1_t2;sF_s1_new = diff_s1_t2;}
					else {sF_t2_new = 0;sF_s1_new = 0;}
					Mp->set_F_s1(sF_s1_new);Mp->set_F_t1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);Mp->set_F_s2(sF_t2_new);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_vert = true;}
					/////
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (diff_s1_t2>0)
					{
						double Fs11_d = sqrt((Fs11x + Ft11x)*(Fs11x + Ft11x) + (Fs11y + Ft11y)*(Fs11y + Ft11y));
						C2_sense_mech_vec_x += diff_s1_t2*(Fs11x + Ft11x)/Fs11_d;
						C2_sense_mech_vec_y += diff_s1_t2*(Fs11y + Ft11y)/Fs11_d;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C12->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()==0 && t21->get_pop_out()==0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double sF_s1_new = sF_s1;// + plot_s1;
					double sF_t2_new = sF_t2;// + plot_t2;
					int Dir_s1 = -1;
					int Dir_t2 = -1; // direction: -1,opposite; 1, same;
					double plot_s1 = (Fs11x + Ft11x)*ns1_co_normal_x + (Fs11y + Ft11y)*ns1_co_normal_y;
					double plot_t2 = (Ft21x + Fs21x)*ns2_co_normal_x + (Ft21y + Fs21y)*ns2_co_normal_y;
					if (plot_s1>0) {Dir_s1 = 1;} // same as normal
					if (plot_t2>0) {Dir_t2 = 1;} // same as normal
					if (Dir_s1>0) {sF_s1_new += plot_s1;} else {sF_s1_new = 0;}
					if (Dir_t2>0) {sF_t2_new += plot_t2;} else {sF_t2_new = 0;}
					Mp->set_F_s1(sF_s1_new);Mp->set_F_t1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);Mp->set_F_s2(sF_t2_new);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_vert = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (Dir_s1>0)
					{
						C1_sense_mech_vec_x += Fs11x + Ft11x;
						C1_sense_mech_vec_y += Fs11y + Ft11y;
						C2_sense_mech_vec_x += Ft21x + Fs21x;
						C2_sense_mech_vec_y += Ft21y + Fs21y;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C11->set_sense_mech(true);
						C12->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()<0 && t21->get_pop_out()==0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double sF_s1_new = sF_s1;
					double sF_t2_new = sF_t2;
					int Dir_t2 = -1; // direction: -1,opposite; 1, same;
					double plot_s1 = abs((Fs11x + Ft11x)*ns1_co_normal_x + (Fs11y + Ft11y)*ns1_co_normal_y);
					double plot_t2 = (Ft21x + Fs21x)*ns2_co_normal_x + (Ft21y + Fs21y)*ns2_co_normal_y;
					if (plot_t2>0) {Dir_t2 = 1;} // same as normal
					sF_s1_new += plot_s1;
					if (Dir_t2>0) {sF_t2_new += plot_t2;} else {sF_t2_new = 0;}
					Mp->set_F_s1(sF_s1_new);Mp->set_F_t1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);Mp->set_F_s2(sF_t2_new);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_vert = true;}
					/////
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (Dir_t2>0)
					{
						C2_sense_mech_vec_x += Ft21x + Fs21x;
						C2_sense_mech_vec_y += Ft21y + Fs21y;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C12->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()==0 && t21->get_pop_out()<0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double sF_s1_new = sF_s1;// + plot_s1;
					double sF_t2_new = sF_t2;// + plot_t2;
					int Dir_s1 = -1;
					double plot_s1 = (Fs11x + Ft11x)*ns1_co_normal_x + (Fs11y + Ft11y)*ns1_co_normal_y;
					double plot_t2 = abs((Ft21x + Fs21x)*ns2_co_normal_x + (Ft21y + Fs21y)*ns2_co_normal_y);
					if (plot_s1>0) {Dir_s1 = 1;} // same as normal
					if (Dir_s1>0) {sF_s1_new += plot_s1;} else {sF_s1_new = 0;}
					sF_t2_new += plot_t2;
					Mp->set_F_s1(sF_s1_new);Mp->set_F_t1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);Mp->set_F_s2(sF_t2_new);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_vert = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					if (Dir_s1>0)
					{
						C1_sense_mech_vec_x += Fs11x + Ft11x;
						C1_sense_mech_vec_y += Fs11y + Ft11y;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()<0 && t21->get_pop_out()<0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double plot_s1 = abs((Fs11x + Ft11x)*ns1_co_normal_x + (Fs11y + Ft11y)*ns1_co_normal_y);
					double plot_t2 = abs((Ft21x + Fs21x)*ns2_co_normal_x + (Ft21y + Fs21y)*ns2_co_normal_y);
					double sF_s1_new = sF_s1 + plot_s1;
					double sF_t2_new = sF_t2 + plot_t2;
					Mp->set_F_s1(sF_s1_new);Mp->set_F_t1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);Mp->set_F_s2(sF_t2_new);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_vert = true;}
				}
				if (s11_vert)
				{
					pointbreak0.push_back(s11);
				}
				///
				if (!s11->get_migrate_mark() && !t21->get_migrate_mark())
				{
					if      (sense_migration_flag==0)
					{
						double C1_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						if (C11->Migrate())
						{
							double C11cx = C11->get_center(0);
							double C11cy = C11->get_center(1);
							double C12cx = C12->get_center(0);
							double C12cy = C12->get_center(1);
							double m21x = C11cx - C12cx;
							double m21y = C11cy - C12cy;
							double m21d = sqrt(m21x*m21x + m21y*m21y);
							C2_sense_migr_vec_x += m21x/m21d;
							C2_sense_migr_vec_y += m21y/m21d;
							C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
							C12->set_sense_migr(true);
						}
						if (C12->Migrate())
						{
							double C11cx = C11->get_center(0);
							double C11cy = C11->get_center(1);
							double C12cx = C12->get_center(0);
							double C12cy = C12->get_center(1);
							double m12x = C12cx - C11cx;
							double m12y = C12cy - C11cy;
							double m12d = sqrt(m12x*m12x + m12y*m12y);
							C1_sense_migr_vec_x += m12x/m12d;
							C1_sense_migr_vec_y += m12y/m12d;
							C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
							C11->set_sense_migr(true);
						}
					}
					else if (sense_migration_flag==1)
					{
						int Dir_s1 = -1;
						int Dir_t2 = -1; // direction: -1, opposite; 1: same;
						double plot_s1 = (Fs11x_m + Ft11x_m)*ns1_co_normal_x + (Fs11y_m + Ft11y_m)*ns1_co_normal_y;
						double plot_t2 = (Ft21x_m + Fs21x_m)*ns2_co_normal_x + (Ft21y_m + Fs21y_m)*ns2_co_normal_y;
						if (plot_s1>0) {Dir_s1 = 1;} // same as normal
						if (plot_t2>0) {Dir_t2 = 1;} // same as normal
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						if (Dir_s1>0)
						{
							/// strategy 2:
							if (C11->Migrate())
							{
								double C11_ma = C11->get_migrate_angle();
								double C11_ma_x = cos(C11_ma*PI/180);
								double C11_ma_y = sin(C11_ma*PI/180);
								double plot_C11_ma_p = -(Fs11x_m + Ft11x_m)*C11_ma_x -(Fs11y_m + Ft11y_m)*C11_ma_y;
								double plot_C11_ma_p_d = sqrt((Fs11x_m + Ft11x_m)*(Fs11x_m + Ft11x_m) + (Fs11y_m + Ft11y_m)*(Fs11y_m + Ft11y_m));
								if (plot_C11_ma_p>0.71*plot_C11_ma_p_d) // 45 degree
								{
									C2_sense_migr_vec_x -= Fs11x_m + Ft11x_m;
									C2_sense_migr_vec_y -= Fs11y_m + Ft11y_m;
									C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
									C12->set_sense_migr(true);
								}
							}
						}
						if (Dir_t2>0)
						{
							/// strategy 2:
							if (C12->Migrate())
							{
								double C12_ma = C12->get_migrate_angle();
								double C12_ma_x = cos(C12_ma*PI/180);
								double C12_ma_y = sin(C12_ma*PI/180);
								double plot_C12_ma_m = -(Ft21x_m + Fs21x_m)*C12_ma_x -(Ft21y_m + Fs21y_m)*C12_ma_y;
								double plot_C12_ma_m_d = sqrt((Ft21x_m + Fs21x_m)*(Ft21x_m + Fs21x_m) + (Ft21y_m + Fs21y_m)*(Ft21y_m + Fs21y_m));
								if (plot_C12_ma_m>0.71*plot_C12_ma_m_d)
								{
									C1_sense_migr_vec_x -= Ft21x_m + Fs21x_m;
									C1_sense_migr_vec_y -= Ft21y_m + Fs21y_m;						
									C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);						
									C11->set_sense_migr(true);
								}
							}
						}
					}
				}
				else if (s11->get_migrate_mark() && !t21->get_migrate_mark())
				{
					if (sense_migration_flag==1)
					{
						double Dir_s1 = -1;
						double plot_s1 = (Fs11x_m + Ft11x_m)*ns1_co_normal_x + (Fs11y_m + Ft11y_m)*ns1_co_normal_y;
						double F_m_d = sqrt((Fs11x_m + Ft11x_m)*(Fs11x_m + Ft11x_m) + (Fs11y_m + Ft11y_m)*(Fs11y_m + Ft11y_m));
						Dir_s1 = plot_s1/F_m_d;
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						double C11_ma = C11->get_migrate_angle();
						double C11_ma_x = cos(C11_ma*PI/180);
						double C11_ma_y = sin(C11_ma*PI/180);
						double C11_ma_2 = (Ft21x_m + Fs21x_m)*C11_ma_x + (Ft21y_m + Fs21y_m)*C11_ma_y;
						double C11_ma_2_d = sqrt((Ft21x_m + Fs21x_m)*(Ft21x_m + Fs21x_m) + (Ft21y_m + Fs21y_m)*(Ft21y_m + Fs21y_m));
						if (C11_ma_2>0.71*C11_ma_2_d && Dir_s1<0.25) // 15 degree
						{
							C2_sense_migr_vec_x += Ft21x_m + Fs21x_m;
							C2_sense_migr_vec_y += Ft21y_m + Fs21y_m;
							C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
							C12->set_sense_migr(true);
						}
					}
				}
				else if (!s11->get_migrate_mark() && t21->get_migrate_mark())
				{
					if (sense_migration_flag==1)
					{
						double Dir_t2 = -1;
						double plot_t2 = (Ft21x_m + Fs21x_m)*ns2_co_normal_x + (Ft21y_m + Fs21y_m)*ns2_co_normal_y;
						double F_m_d = sqrt((Ft21x_m + Fs21x_m)*(Ft21x_m + Fs21x_m) + (Ft21y_m + Fs21y_m)*(Ft21y_m + Fs21y_m));
						Dir_t2 = plot_t2/F_m_d;
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C12_ma = C12->get_migrate_angle();
						double C12_ma_x = cos(C12_ma*PI/180);
						double C12_ma_y = sin(C12_ma*PI/180);
						double C12_ma_1 = (Fs11x_m + Ft11x_m)*C12_ma_x + (Fs11y_m + Ft11y_m)*C12_ma_y;
						double C12_ma_1_d = sqrt((Fs11x_m + Ft11x_m)*(Fs11x_m + Ft11x_m) + (Fs11y_m + Ft11y_m)*(Fs11y_m + Ft11y_m));
						if (C12_ma_1>0.71*C12_ma_1_d && Dir_t2<0.25) // 15 degree
						{
							C1_sense_migr_vec_x += Fs11x_m + Ft11x_m;
							C1_sense_migr_vec_y += Fs11y_m + Ft11y_m;
							C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
							C11->set_sense_migr(true);
						}
					}
				}
			}	
		}
		else if (s11!=t11)
		{
			//cout<<"s11!=t11: collisionpair:"<<collisionpairList[i]->get_I1()->id()<<" "<<collisionpairList[i]->get_I2()->id()<<": "<<endl;
			if      (!s11->getba()->get_attach() || !t21->getab()->get_attach())
			{
				if (!s11->getba()->get_attach())
				{
					edge *s11_ba = s11->getba();
					double vs20 = s11_ba->get_virtual_triangle_stiffness_matrix(2,0);
					double vs21 = s11_ba->get_virtual_triangle_stiffness_matrix(2,1);
					double vs22 = s11_ba->get_virtual_triangle_stiffness_matrix(2,2);
					double vs23 = s11_ba->get_virtual_triangle_stiffness_matrix(2,3);
					double vs30 = s11_ba->get_virtual_triangle_stiffness_matrix(3,0);
					double vs31 = s11_ba->get_virtual_triangle_stiffness_matrix(3,1);
					double vs32 = s11_ba->get_virtual_triangle_stiffness_matrix(3,2);
					double vs33 = s11_ba->get_virtual_triangle_stiffness_matrix(3,3);
					int p1s_sm_id = s11->getrp()->get_node_id()->SM_ID();
					int p2s_sm_id = s11->get_node_id()->SM_ID();
					double us0 = Xn[p1s_sm_id*2];
					double us1 = Xn[p1s_sm_id*2+1];
					double us2 = Xn[p2s_sm_id*2];
					double us3 = Xn[p2s_sm_id*2+1];
					double us0_m = Xn_p[p1s_sm_id*2];
					double us1_m = Xn_p[p1s_sm_id*2+1];
					double us2_m = Xn_p[p2s_sm_id*2];
					double us3_m = Xn_p[p2s_sm_id*2+1];
					Fs11x += vs20*us0 + vs21*us1 + vs22*us2 + vs23*us3;
					Fs11y += vs30*us0 + vs31*us1 + vs32*us2 + vs33*us3;
					Fs11x_m += vs20*us0_m + vs21*us1_m + vs22*us2_m + vs23*us3_m;
					Fs11y_m += vs30*us0_m + vs31*us1_m + vs32*us2_m + vs33*us3_m;
				}
				if (!C11->Soften())
				{
					int at_ns = (int)s11->get_inner_at()->size();
					for (int j=0;j<at_ns;j++)
					{
						triangle *Tj = s11->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==s11)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Fs11x += FAx;
							Fs11y += FAy;
							Fs11x_m += FAx_m;
							Fs11y_m += FAy_m;
						}
						else if (TBj==s11)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Fs11x += FBx;
							Fs11y += FBy;
							Fs11x_m += FBx_m;
							Fs11y_m += FBy_m;
						}
						else if (TCj==s11)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Fs11x += FCx;
							Fs11y += FCy;
							Fs11x_m += FCx_m;
							Fs11y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = s11->get_node_id()->SM_ID();
					int N2_sm_id = s11->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = s11->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1];
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1];
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1];
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1];
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1];
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1];
					double v00 = s11->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = s11->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = s11->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = s11->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = s11->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = s11->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = s11->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = s11->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = s11->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = s11->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = s11->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = s11->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = s11->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = s11->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = s11->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = s11->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Fs11x += (Frx + Ffx);
					Fs11y += (Fry + Ffy);
					Fs11x_m += (Frx_m + Ffx_m);
					Fs11y_m += (Fry_m + Ffy_m);
				}
				if (!t21->getab()->get_attach())
				{
					edge *t21_ab = t21->getab();
					double vt00 = t21_ab->get_virtual_triangle_stiffness_matrix(0,0);
					double vt01 = t21_ab->get_virtual_triangle_stiffness_matrix(0,1);
					double vt02 = t21_ab->get_virtual_triangle_stiffness_matrix(0,2);
					double vt03 = t21_ab->get_virtual_triangle_stiffness_matrix(0,3);
					double vt10 = t21_ab->get_virtual_triangle_stiffness_matrix(1,0);
					double vt11 = t21_ab->get_virtual_triangle_stiffness_matrix(1,1);
					double vt12 = t21_ab->get_virtual_triangle_stiffness_matrix(1,2);
					double vt13 = t21_ab->get_virtual_triangle_stiffness_matrix(1,3);
					int p1t_sm_id = t21->get_node_id()->SM_ID();
					int p2t_sm_id = t21->getfp()->get_node_id()->SM_ID();
					double ut0 = Xn[p1t_sm_id*2];
					double ut1 = Xn[p1t_sm_id*2+1];
					double ut2 = Xn[p2t_sm_id*2];
					double ut3 = Xn[p2t_sm_id*2+1];
					double ut0_m = Xn_p[p1t_sm_id*2];
					double ut1_m = Xn_p[p1t_sm_id*2+1];
					double ut2_m = Xn_p[p2t_sm_id*2];
					double ut3_m = Xn_p[p2t_sm_id*2+1];
					Ft21x += vt00*ut0 + vt01*ut1 + vt02*ut2 + vt03*ut3;
					Ft21y += vt10*ut0 + vt11*ut1 + vt12*ut2 + vt13*ut3;
					Ft21x_m += vt00*ut0_m + vt01*ut1_m + vt02*ut2_m + vt03*ut3_m;
					Ft21y_m += vt10*ut0_m + vt11*ut1_m + vt12*ut2_m + vt13*ut3_m;
				}
				if (!C12->Soften())
				{
					int at_nt = (int)t21->get_inner_at()->size();
					for (int j=0;j<at_nt;j++)
					{
						triangle *Tj = t21->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==t21)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Ft21x += FAx;
							Ft21y += FAy;
							Ft21x_m += FAx_m;
							Ft21y_m += FAy_m;
						}
						else if (TBj==t21)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Ft21x += FBx;
							Ft21y += FBy;
							Ft21x_m += FBx_m;
							Ft21y_m += FBy_m;
						}
						else if (TCj==t21)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Ft21x += FCx;
							Ft21y += FCy;
							Ft21x_m += FCx_m;
							Ft21y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = t21->get_node_id()->SM_ID();
					int N2_sm_id = t21->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = t21->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1]; // p
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1]; // rp
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1]; // fp
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1]; // p
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1]; // rp
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1]; // fp
					double v00 = t21->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = t21->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = t21->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = t21->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = t21->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = t21->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = t21->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = t21->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = t21->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = t21->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = t21->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = t21->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = t21->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = t21->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = t21->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = t21->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Ft21x += (Frx + Ffx);
					Ft21y += (Fry + Ffy);
					Ft21x_m += (Frx_m + Ffx_m);
					Ft21y_m += (Fry_m + Ffy_m);
				}
			}
			else if (s11->getba()->get_attach() && t21->getab()->get_attach())
			{
				if (!C11->Soften())
				{
					int at_ns = (int)s11->get_inner_at()->size();
					for (int j=0;j<at_ns;j++)
					{
						triangle *Tj = s11->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==s11)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Fs11x += FAx;
							Fs11y += FAy;
							Fs11x_m += FAx_m;
							Fs11y_m += FAy_m;
						}
						else if (TBj==s11)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Fs11x += FBx;
							Fs11y += FBy;
							Fs11x_m += FBx_m;
							Fs11y_m += FBy_m;
						}
						else if (TCj==s11)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Fs11x += FCx;
							Fs11y += FCy;
							Fs11x_m += FCx_m;
							Fs11y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = s11->get_node_id()->SM_ID();
					int N2_sm_id = s11->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = s11->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1];
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1];
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1];
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1];
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1];
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1];
					double v00 = s11->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = s11->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = s11->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = s11->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = s11->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = s11->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = s11->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = s11->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = s11->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = s11->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = s11->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = s11->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = s11->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = s11->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = s11->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = s11->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Fs11x += (Frx + Ffx);
					Fs11y += (Fry + Ffy);
					Fs11x_m += (Frx_m + Ffx_m);
					Fs11y_m += (Fry_m + Ffy_m);
				}
				if (!C12->Soften())
				{
					int at_nt = (int)t21->get_inner_at()->size();
					for (int j=0;j<at_nt;j++)
					{
						triangle *Tj = t21->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==t21)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Ft21x += FAx;
							Ft21y += FAy;
							Ft21x_m += FAx_m;
							Ft21y_m += FAy_m;
						}
						else if (TBj==t21)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Ft21x += FBx;
							Ft21y += FBy;
							Ft21x_m += FBx_m;
							Ft21y_m += FBy_m;
						}
						else if (TCj==t21)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Ft21x += FCx;
							Ft21y += FCy;
							Ft21x_m += FCx_m;
							Ft21y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = t21->get_node_id()->SM_ID();
					int N2_sm_id = t21->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = t21->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1]; // p
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1]; // rp
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1]; // fp
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1]; // p
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1]; // rp
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1]; // fp
					double v00 = t21->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = t21->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = t21->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = t21->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = t21->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = t21->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = t21->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = t21->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = t21->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = t21->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = t21->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = t21->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = t21->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = t21->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = t21->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = t21->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Ft21x += (Frx + Ffx);
					Ft21y += (Fry + Ffy);
					Ft21x_m += (Frx_m + Ffx_m);
					Ft21y_m += (Fry_m + Ffy_m);
				}
			}
			if      (!s21->getba()->get_attach() || !t11->getab()->get_attach())
			{
				if (!s21->getba()->get_attach())
				{
					edge *s21_ba = s21->getba();
					double vs20 = s21_ba->get_virtual_triangle_stiffness_matrix(2,0);
					double vs21 = s21_ba->get_virtual_triangle_stiffness_matrix(2,1);
					double vs22 = s21_ba->get_virtual_triangle_stiffness_matrix(2,2);
					double vs23 = s21_ba->get_virtual_triangle_stiffness_matrix(2,3);
					double vs30 = s21_ba->get_virtual_triangle_stiffness_matrix(3,0);
					double vs31 = s21_ba->get_virtual_triangle_stiffness_matrix(3,1);
					double vs32 = s21_ba->get_virtual_triangle_stiffness_matrix(3,2);
					double vs33 = s21_ba->get_virtual_triangle_stiffness_matrix(3,3);
					int p1s_sm_id = s21->getrp()->get_node_id()->SM_ID();
					int p2s_sm_id = s21->get_node_id()->SM_ID();
					double us0 = Xn[p1s_sm_id*2];
					double us1 = Xn[p1s_sm_id*2+1];
					double us2 = Xn[p2s_sm_id*2];
					double us3 = Xn[p2s_sm_id*2+1];
					double us0_m = Xn_p[p1s_sm_id*2];
					double us1_m = Xn_p[p1s_sm_id*2+1];
					double us2_m = Xn_p[p2s_sm_id*2];
					double us3_m = Xn_p[p2s_sm_id*2+1];
					Fs21x += vs20*us0 + vs21*us1 + vs22*us2 + vs23*us3;
					Fs21y += vs30*us0 + vs31*us1 + vs32*us2 + vs33*us3;
					Fs21x_m += vs20*us0_m + vs21*us1_m + vs22*us2_m + vs23*us3_m;
					Fs21y_m += vs30*us0_m + vs31*us1_m + vs32*us2_m + vs33*us3_m;
				}
				if (!C12->Soften())
				{
					int at_ns = (int)s21->get_inner_at()->size();
					for (int j=0;j<at_ns;j++)
					{
						triangle *Tj = s21->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==s21)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Fs21x += FAx;
							Fs21y += FAy;
							Fs21x_m += FAx_m;
							Fs21y_m += FAy_m;
						}
						else if (TBj==s21)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Fs21x += FBx;
							Fs21y += FBy;
							Fs21x_m += FBx_m;
							Fs21y_m += FBy_m;
						}
						else if (TCj==s21)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Fs21x += FCx;
							Fs21y += FCy;
							Fs21x_m += FCx_m;
							Fs21y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = s21->get_node_id()->SM_ID();
					int N2_sm_id = s21->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = s21->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1]; // p
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1]; // rp
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1]; // fp
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1]; // p
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1]; // rp
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1]; // fp
					double v00 = s21->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = s21->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = s21->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = s21->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = s21->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = s21->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = s21->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = s21->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = s21->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = s21->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = s21->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = s21->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = s21->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = s21->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = s21->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = s21->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Fs21x += (Frx + Ffx);
					Fs21y += (Fry + Ffy);
					Fs21x_m += (Frx_m + Ffx_m);
					Fs21y_m += (Fry_m + Ffy_m);
				}
				if (!t11->getab()->get_attach())
				{
					edge *t11_ab = t11->getab();
					double vt00 = t11_ab->get_virtual_triangle_stiffness_matrix(0,0);
					double vt01 = t11_ab->get_virtual_triangle_stiffness_matrix(0,1);
					double vt02 = t11_ab->get_virtual_triangle_stiffness_matrix(0,2);
					double vt03 = t11_ab->get_virtual_triangle_stiffness_matrix(0,3);
					double vt10 = t11_ab->get_virtual_triangle_stiffness_matrix(1,0);
					double vt11 = t11_ab->get_virtual_triangle_stiffness_matrix(1,1);
					double vt12 = t11_ab->get_virtual_triangle_stiffness_matrix(1,2);
					double vt13 = t11_ab->get_virtual_triangle_stiffness_matrix(1,3);
					int p1t_sm_id = t11->get_node_id()->SM_ID();
					int p2t_sm_id = t11->getfp()->get_node_id()->SM_ID();
					double ut0 = Xn[p1t_sm_id*2];
					double ut1 = Xn[p1t_sm_id*2+1];
					double ut2 = Xn[p2t_sm_id*2];
					double ut3 = Xn[p2t_sm_id*2+1];
					double ut0_m = Xn_p[p1t_sm_id*2];
					double ut1_m = Xn_p[p1t_sm_id*2+1];
					double ut2_m = Xn_p[p2t_sm_id*2];
					double ut3_m = Xn_p[p2t_sm_id*2+1];
					Ft11x += vt00*ut0 + vt01*ut1 + vt02*ut2 + vt03*ut3;
					Ft11y += vt10*ut0 + vt11*ut1 + vt12*ut2 + vt13*ut3;
					Ft11x_m += vt00*ut0_m + vt01*ut1_m + vt02*ut2_m + vt03*ut3_m;
					Ft11y_m += vt10*ut0_m + vt11*ut1_m + vt12*ut2_m + vt13*ut3_m;
				}
				if (!C11->Soften())
				{
					int at_nt = (int)t11->get_inner_at()->size();
					for (int j=0;j<at_nt;j++)
					{
						triangle *Tj = t11->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==t11)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Ft11x += FAx;
							Ft11y += FAy;
							Ft11x_m += FAx_m;
							Ft11y_m += FAy_m;
						}
						else if (TBj==t11)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Ft11x += FBx;
							Ft11y += FBy;
							Ft11x_m += FBx_m;
							Ft11y_m += FBy_m;
						}
						else if (TCj==t11)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Ft11x += FCx;
							Ft11y += FCy;
							Ft11x_m += FCx_m;
							Ft11y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = t11->get_node_id()->SM_ID();
					int N2_sm_id = t11->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = t11->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1];
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1];
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1];
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1];
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1];
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1];
					double v00 = t11->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = t11->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = t11->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = t11->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = t11->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = t11->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = t11->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = t11->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = t11->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = t11->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = t11->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = t11->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = t11->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = t11->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = t11->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = t11->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Ft11x += (Frx + Ffx);
					Ft11y += (Fry + Ffy);
					Ft11x_m += (Frx_m + Ffx_m);
					Ft11y_m += (Fry_m + Ffy_m);
				}
			}
			else if (s21->getba()->get_attach() && t11->getab()->get_attach())
			{
				if (!C12->Soften())
				{
					int at_ns = (int)s21->get_inner_at()->size();
					for (int j=0;j<at_ns;j++)
					{
						triangle *Tj = s21->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==s21)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Fs21x += FAx;
							Fs21y += FAy;
							Fs21x_m += FAx_m;
							Fs21y_m += FAy_m;
						}
						else if (TBj==s21)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Fs21x += FBx;
							Fs21y += FBy;
							Fs21x_m += FBx_m;
							Fs21y_m += FBy_m;
						}
						else if (TCj==s21)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Fs21x += FCx;
							Fs21y += FCy;
							Fs21x_m += FCx_m;
							Fs21y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = s21->get_node_id()->SM_ID();
					int N2_sm_id = s21->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = s21->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1]; // p
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1]; // rp
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1]; // fp
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1]; // p
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1]; // rp
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1]; // fp
					double v00 = s21->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = s21->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = s21->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = s21->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = s21->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = s21->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = s21->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = s21->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = s21->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = s21->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = s21->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = s21->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = s21->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = s21->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = s21->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = s21->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Fs21x += (Frx + Ffx);
					Fs21y += (Fry + Ffy);
					Fs21x_m += (Frx_m + Ffx_m);
					Fs21y_m += (Fry_m + Ffy_m);
				}
				if (!C11->Soften())
				{
					int at_nt = (int)t11->get_inner_at()->size();
					for (int j=0;j<at_nt;j++)
					{
						triangle *Tj = t11->get_inner_at()->at(j);
						point3D *TAj = Tj->getA();
						point3D *TBj = Tj->getB();
						point3D *TCj = Tj->getC();
						int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
						int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
						int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
						double u0 = Xn[TAj_sm_id*2];
						double u1 = Xn[TAj_sm_id*2+1];
						double u2 = Xn[TBj_sm_id*2];
						double u3 = Xn[TBj_sm_id*2+1];
						double u4 = Xn[TCj_sm_id*2];
						double u5 = Xn[TCj_sm_id*2+1];
						double u0_m = Xn_p[TAj_sm_id*2];
						double u1_m = Xn_p[TAj_sm_id*2+1];
						double u2_m = Xn_p[TBj_sm_id*2];
						double u3_m = Xn_p[TBj_sm_id*2+1];
						double u4_m = Xn_p[TCj_sm_id*2];
						double u5_m = Xn_p[TCj_sm_id*2+1];
						if      (TAj==t11)
						{
							double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
							double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
							double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
							double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
							double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
							double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
							Ft11x += FAx;
							Ft11y += FAy;
							Ft11x_m += FAx_m;
							Ft11y_m += FAy_m;
						}
						else if (TBj==t11)
						{
							double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
							double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
							double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
							double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
							double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
							double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
							Ft11x += FBx;
							Ft11y += FBy;
							Ft11x_m += FBx_m;
							Ft11y_m += FBy_m;
						}
						else if (TCj==t11)
						{
							double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
							double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
							double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
							double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
							double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
							double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
							Ft11x += FCx;
							Ft11y += FCy;
							Ft11x_m += FCx_m;
							Ft11y_m += FCy_m;
						}
					}
				}
				else
				{
					int N1_sm_id = t11->get_node_id()->SM_ID();
					int N2_sm_id = t11->getrp()->get_node_id()->SM_ID();
					int N3_sm_id = t11->getfp()->get_node_id()->SM_ID();
					double u0 = Xn[N1_sm_id*2];
					double u1 = Xn[N1_sm_id*2+1];
					double u2 = Xn[N2_sm_id*2];
					double u3 = Xn[N2_sm_id*2+1];
					double u4 = Xn[N3_sm_id*2];
					double u5 = Xn[N3_sm_id*2+1];
					double u0_m = Xn_p[N1_sm_id*2];
					double u1_m = Xn_p[N1_sm_id*2+1];
					double u2_m = Xn_p[N2_sm_id*2];
					double u3_m = Xn_p[N2_sm_id*2+1];
					double u4_m = Xn_p[N3_sm_id*2];
					double u5_m = Xn_p[N3_sm_id*2+1];
					double v00 = t11->getab()->get_soften_triangle_stiffness_matrix(0,0);
					double v01 = t11->getab()->get_soften_triangle_stiffness_matrix(0,1);
					double v02 = t11->getab()->get_soften_triangle_stiffness_matrix(0,2);
					double v03 = t11->getab()->get_soften_triangle_stiffness_matrix(0,3);
					double v10 = t11->getab()->get_soften_triangle_stiffness_matrix(1,0);
					double v11 = t11->getab()->get_soften_triangle_stiffness_matrix(1,1);
					double v12 = t11->getab()->get_soften_triangle_stiffness_matrix(1,2);
					double v13 = t11->getab()->get_soften_triangle_stiffness_matrix(1,3);
					double v20 = t11->getba()->get_soften_triangle_stiffness_matrix(2,0);
					double v21 = t11->getba()->get_soften_triangle_stiffness_matrix(2,1);
					double v22 = t11->getba()->get_soften_triangle_stiffness_matrix(2,2);
					double v23 = t11->getba()->get_soften_triangle_stiffness_matrix(2,3);
					double v30 = t11->getba()->get_soften_triangle_stiffness_matrix(3,0);
					double v31 = t11->getba()->get_soften_triangle_stiffness_matrix(3,1);
					double v32 = t11->getba()->get_soften_triangle_stiffness_matrix(3,2);
					double v33 = t11->getba()->get_soften_triangle_stiffness_matrix(3,3);
					double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
					double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
					double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
					double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
					double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
					double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
					double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
					double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
					Ft11x += (Frx + Ffx);
					Ft11y += (Fry + Ffy);
					Ft11x_m += (Frx_m + Ffx_m);
					Ft11y_m += (Fry_m + Ffy_m);
				}
			}
			/******************************************************
			 case 1: all static
			 -----------------------------------------------------
					 status: static
					 direction of recovered force: opposite to normal
					 ==> attach together
					 ---------------------------------------------
					 status: static
					 direction of recovered force: same as normal
					 ==> break apart
					 ---------------------------------------------
				   o---o---o
				  /  cell 0 \  break
				 o      -->  oo---o---o
				 | A        // <--     \
				 | |   --> //           o
		 attach  o====o===oo <-- cell 2 |  ===>  dragging force
				 | |   --> \\           o  
				 | V        \\ <--     /
				 o      -->  oo---o---o
				  \  cell 1 /  break
				   o---o---o
			 =====================================================
			 case 2: all static
			 -----------------------------------------------------
					 status: static
					 direction of recovered force: opposite to normal
					 ==> attach together
					 ---------------------------------------------
					 status: static
					 direction of recovered force: same as normal
					 ==> break apart
					 ---------------------------------------------
				   o---o---o
				  /  cell 0 \  attach
				 o      <--  oo---o---o
				 | |     |  // -->     \
				 | V     V //           o
		  break  o====o===oo --> cell 2 |  <===  pushing force
				 | A     A \\           o  
				 | |     |  \\ -->     /
				 o      <--  oo---o---o
				  \  cell 1 /  attach 
				   o---o---o
			 =====================================================
			 case 3: all grow
			 -----------------------------------------------------
			         status: grow
					 direction of recovered force: same as normal
					 ==> attach together
					 ---------------------------------------------
				  o---o---o
				  /  cell 0 \  attach
				 o      -->  oo---o---o
				 | |     |  // <--     \
				 | V     V //           o
		 attach  o====o===oo <-- cell 2 |  
				 | A     A \\           o  
				 | |     |  \\ <--     /
				 o      -->  oo---o---o
				  \  cell 1 /  attach 
				   o---o---o
			 =====================================================
			 case 4: all shrink
			 -----------------------------------------------------
			         status: shrink
					 direction of recovered force: opposite to normal
					 ==> attach together
					 ---------------------------------------------
				   o---o---o
				  /  cell 0 \  attach
				 o      <--  oo---o---o
				 | A     A  // -->     \
				 | |     | //           o
		 attach  o====o===oo --> cell 2 |
				 | |     | \\           o  
				 | V     V  \\ -->     /
				 o      <--  oo---o---o
				  \  cell 1 /  attach
				   o---o---o
			=====================================================
			 case 5: cell 2: grow; cell 0,1: static 
			 -----------------------------------------------------
					 status: static
					 direction of recovered force: opposite to normal
					 ==> attach together
					 ---------------------------------------------
					 status: static
					 direction of recovered force: same as normal
					 ==> break apart
					 ---------------------------------------------
					 status: grow
					 direction of recovered force: same as normal
					 ==> attach together
					 ---------------------------------------------
				   o---o---o
				  /  cell 0 \  attach
				 o      <--  oo---o---o
				 | |     |  // <--     \
				 | V     V //           o
		  break  o====o===oo <-- cell 2 |
				 | A     A \\           o
				 | |     |  \\ <--     /
				 o      <--  oo---o---o
				  \  cell 1 /  attach
				   o---o---o
			=====================================================
			 case 6: cell 2: shrink; cell 0,1,3: static 
			 -----------------------------------------------------
			         status: static
					 direction of recovered force: opposite to normal
					 ==> attach together
					 ---------------------------------------------
					 status: static
					 direction of recovered force: same as normal
					 ==> break apart
					 ---------------------------------------------
					 status: shrink
					 direction of recovered force: opposite to normal
					 ==> break apart
					 ---------------------------------------------
				   o---o---o
				  /  cell 3 \
				 o           o
				 |           |
				 o  A     |  o
				  \ |     V /
		attach     o===o===o  break
				  / |     A \
				 o  V     |  o
				 |           |  break
				 o      -->  oo---o---o
				 | A cell 0 // -->     \
				 | |   --> //           o
		 attach  o====o===oo --> cell 2 |  
				 | |   --> \\           o  
				 | V        \\ -->     /
				 o      -->  oo---o---o
				  \  cell 1 /  break 
				   o---o---o
			******************************************************/
			s11->set_elastic_force(Fs11x, Fs11y);
			t11->set_elastic_force(Ft11x, Ft11y);
			s21->set_elastic_force(Fs21x, Fs21y);
			t21->set_elastic_force(Ft21x, Ft21y);
			bool s11_ab = false;
			bool s21_ab = false;
			int s1mn = (int)s11->getmp()->size();
			int t1mn = (int)t11->getmp()->size();
			if (s1mn==1)
			{
				double ns1x = s11->getab()->getnormal(0);
				double ns1y = s11->getab()->getnormal(1);
				double ns1d = sqrt(ns1x*ns1x + ns1y*ns1y);
				double normal_s1_x = ns1x/ns1d;
				double normal_s1_y = ns1y/ns1d;
				if      ((s11->get_pop_out()>0 && t21->get_pop_out()>0)  ||
						 (s11->get_pop_out()>0 && t21->get_pop_out()==0) ||
						 (s11->get_pop_out()==0 && t21->get_pop_out()>0)) 
				{
					Mp->set_F_s1(0);
					Mp->set_F_t2(0); // stored-apart-force refreshed to NULL
					s11->set_stored_force_vis(0,0);
					t21->set_stored_force_vis(0,0);
				}
				else if (s11->get_pop_out()==0 && t21->get_pop_out()==0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double sF_s1_new = sF_s1;// + plot_s1;
					double sF_t2_new = sF_t2;// + plot_t2;
					int Dir_s1 = -1;
					int Dir_t2 = -1;
					double plot_s1 = Fs11x*normal_s1_x + Fs11y*normal_s1_y;
					double plot_t2 = -Ft21x*normal_s1_x -Ft21y*normal_s1_y;
					if (plot_s1>0) {Dir_s1 = 1;}
					if (plot_t2>0) {Dir_t2 = 1;}
					if (Dir_s1>0) {sF_s1_new += plot_s1;} else {sF_s1_new = 0;}
					if (Dir_t2>0) {sF_t2_new += plot_t2;} else {sF_t2_new = 0;}
					Mp->set_F_s1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);
					double s_f_s1_x = sF_s1_new*normal_s1_x;
					double s_f_s1_y = sF_s1_new*normal_s1_y;
					double s_f_t2_x = -sF_t2_new*normal_s1_x;
					double s_f_t2_y = -sF_t2_new*normal_s1_y;
					s11->set_stored_force_vis(s_f_s1_x,s_f_s1_y);
					t21->set_stored_force_vis(s_f_t2_x,s_f_t2_y);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_ab = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (Dir_s1>0)
					{
						C1_sense_mech_vec_x += Fs11x;
						C1_sense_mech_vec_y += Fs11y;
						C2_sense_mech_vec_x += Ft21x;
						C2_sense_mech_vec_y += Ft21y;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C11->set_sense_mech(true);
						C12->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()<0 && t21->get_pop_out()>0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double plot_s1 = abs(Fs11x*normal_s1_x + Fs11y*normal_s1_y);
					double plot_t2 = abs(-Ft21x*normal_s1_x -Ft21y*normal_s1_y);
					double sF_s1_new = sF_s1 + plot_s1;
					double sF_t2_new = sF_t2 + plot_t2;
					double diff_s1_t2 = sF_s1_new - sF_t2_new;
					if (diff_s1_t2>0) {sF_t2_new = diff_s1_t2;sF_s1_new = diff_s1_t2;}
					else {sF_t2_new = 0;sF_s1_new = 0;}
					Mp->set_F_s1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);
					double s_f_s1_x = -sF_s1_new*normal_s1_x;
					double s_f_s1_y = -sF_s1_new*normal_s1_y;
					s11->set_stored_force_vis(s_f_s1_x,s_f_s1_y);
					t21->set_stored_force_vis(0,0);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_ab = true;}
					/////
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (diff_s1_t2>0)
					{
						double Fs11_d = sqrt(Fs11x*Fs11x + Fs11y*Fs11y);
						C2_sense_mech_vec_x += diff_s1_t2*Fs11x/Fs11_d;
						C2_sense_mech_vec_y += diff_s1_t2*Fs11y/Fs11_d;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C12->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()>0 && t21->get_pop_out()<0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double plot_s1 = abs(Fs11x*normal_s1_x + Fs11y*normal_s1_y);
					double plot_t2 = abs(-Ft21x*normal_s1_x -Ft21y*normal_s1_y);
					double sF_s1_new = sF_s1 + plot_s1;
					double sF_t2_new = sF_t2 + plot_t2;
					double diff_t2_s1 = sF_t2_new - sF_s1_new;
					if (diff_t2_s1>0) {sF_t2_new = diff_t2_s1;sF_s1_new = diff_t2_s1;}
					else {sF_t2_new = 0;sF_s1_new = 0;}
					Mp->set_F_s1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);
					double s_f_t2_x = sF_t2_new*normal_s1_x;
					double s_f_t2_y = sF_t2_new*normal_s1_y;
					s11->set_stored_force_vis(0,0);
					t21->set_stored_force_vis(s_f_t2_x,s_f_t2_y);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_ab = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					if (diff_t2_s1>0)
					{
						double Ft21_d = sqrt(Ft21x*Ft21x + Ft21y*Ft21y);
						C1_sense_mech_vec_x += diff_t2_s1*Ft21x/Ft21_d;
						C1_sense_mech_vec_y += diff_t2_s1*Ft21y/Ft21_d;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()<0 && t21->get_pop_out()==0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double sF_s1_new = sF_s1;// + plot_s1;
					double sF_t2_new = sF_t2;// + plot_t2;
					int Dir_t2 = -1;
					double plot_s1 = abs(Fs11x*normal_s1_x + Fs11y*normal_s1_y);
					double plot_t2 = -Ft21x*normal_s1_x -Ft21y*normal_s1_y;
					if (plot_t2>0) {Dir_t2 = 1;}
					sF_s1_new += plot_s1;
					if (Dir_t2>0) {sF_t2_new += plot_t2;} else {sF_t2_new = 0;}
					Mp->set_F_s1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);
					double s_f_s1_x = -sF_s1_new*normal_s1_x;
					double s_f_s1_y = -sF_s1_new*normal_s1_y;
					double s_f_t2_x = -sF_t2_new*normal_s1_x;
					double s_f_t2_y = -sF_t2_new*normal_s1_y;
					s11->set_stored_force_vis(s_f_s1_x,s_f_s1_y);
					t21->set_stored_force_vis(s_f_t2_x,s_f_t2_y);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_ab = true;}
					/////
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (Dir_t2>0)
					{
						C2_sense_mech_vec_x += Ft21x;
						C2_sense_mech_vec_y += Ft21y;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C12->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()==0 && t21->get_pop_out()<0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double sF_s1_new = sF_s1;// + plot_s1;
					double sF_t2_new = sF_t2;// + plot_t2;
					int Dir_s1 = -1;
					double plot_s1 = Fs11x*normal_s1_x + Fs11y*normal_s1_y;
					double plot_t2 = abs(-Ft21x*normal_s1_x -Ft21y*normal_s1_y);
					if (plot_s1>0) {Dir_s1 = 1;}
					if (Dir_s1>0) {sF_s1_new += plot_s1;} else {sF_s1_new = 0;}
					sF_t2_new += plot_t2;
					Mp->set_F_s1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);
					double s_f_s1_x = sF_s1_new*normal_s1_x;
					double s_f_s1_y = sF_s1_new*normal_s1_y;
					double s_f_t2_x = sF_t2_new*normal_s1_x;
					double s_f_t2_y = sF_t2_new*normal_s1_y;
					s11->set_stored_force_vis(s_f_s1_x,s_f_s1_y);
					t21->set_stored_force_vis(s_f_t2_x,s_f_t2_y);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_ab = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					if (Dir_s1>0)
					{
						C1_sense_mech_vec_x += Fs11x;
						C1_sense_mech_vec_y += Fs11y;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
				}
				else if (s11->get_pop_out()<0 && t21->get_pop_out()<0)
				{
					double sF_s1 = Mp->get_F_s1();
					double sF_t2 = Mp->get_F_t2();
					double plot_s1 = abs(Fs11x*normal_s1_x + Fs11y*normal_s1_y);
					double plot_t2 = abs(-Ft21x*normal_s1_x -Ft21y*normal_s1_y);
					double sF_s1_new = sF_s1 + plot_s1;
					double sF_t2_new = sF_t2 + plot_t2;
					Mp->set_F_s1(sF_s1_new);
					Mp->set_F_t2(sF_t2_new);
					double s_f_s1_x = -sF_s1_new*normal_s1_x;
					double s_f_s1_y = -sF_s1_new*normal_s1_y;
					double s_f_t2_x = sF_t2_new*normal_s1_x;
					double s_f_t2_y = sF_t2_new*normal_s1_y;
					s11->set_stored_force_vis(s_f_s1_x,s_f_s1_y);
					t21->set_stored_force_vis(s_f_t2_x,s_f_t2_y);
					if (sF_s1_new>F_adhesion || sF_t2_new>F_adhesion) {s11_ab = true;}
				}
				///
				if (!s11->get_migrate_mark() && !t21->get_migrate_mark())
				{
					if      (sense_migration_flag==0)
					{
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						if (C11->Migrate())
						{
							double C11cx = C11->get_center(0);
							double C11cy = C11->get_center(1);
							double C12cx = C12->get_center(0);
							double C12cy = C12->get_center(1);
							double m21x = C11cx - C12cx;
							double m21y = C11cy - C12cy;
							double m21d = sqrt(m21x*m21x + m21y*m21y);
							C2_sense_migr_vec_x += m21x/m21d;
							C2_sense_migr_vec_y += m21y/m21d;
							C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
							C12->set_sense_migr(true);
						}
						if (C12->Migrate())
						{
							double C11cx = C11->get_center(0);
							double C11cy = C11->get_center(1);
							double C12cx = C12->get_center(0);
							double C12cy = C12->get_center(1);
							double m12x = C12cx - C11cx;
							double m12y = C12cy - C11cy;
							double m12d = sqrt(m12x*m12x + m12y*m12y);
							C1_sense_migr_vec_x += m12x/m12d;
							C1_sense_migr_vec_y += m12y/m12d;
							C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
							C11->set_sense_migr(true);
						}
					}
					else if (sense_migration_flag==1)
					{
						int Dir_s1 = -1;
						int Dir_t2 = -1;
						double plot_s1 = Fs11x_m*normal_s1_x + Fs11y_m*normal_s1_y;
						double plot_t2 = -Ft21x_m*normal_s1_x -Ft21y_m*normal_s1_y;
						if (plot_s1>0) {Dir_s1 = 1;}
						if (plot_t2>0) {Dir_t2 = 1;}
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						if (Dir_s1>0)
						{
							/// strategy 2:
							if (C11->Migrate())
							{
								double C11_ma = C11->get_migrate_angle();
								double C11_ma_x = cos(C11_ma*PI/180);
								double C11_ma_y = sin(C11_ma*PI/180);
								double plot_C11_ma_1 = -Fs11x_m*C11_ma_x -Fs11y_m*C11_ma_y;
								double plot_C11_ma_1_d = sqrt(Fs11x_m*Fs11x_m + Fs11y_m*Fs11y_m);
								if (plot_C11_ma_1>0.71*plot_C11_ma_1_d)
								{
									C2_sense_migr_vec_x -= Fs11x_m;
									C2_sense_migr_vec_y -= Fs11y_m;
									C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
									C12->set_sense_migr(true);
									if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<"; Dir_s1:"<<Dir_s1<<"; "<<Fs11x_m<<","<<Fs11y_m<<"; "<<Ft21x_m<<","<<Ft21y_m<<endl;}
								}
							}
						}
						if (Dir_t2>0)
						{
							/// strategy 2:
							if (C12->Migrate())
							{
								double C12_ma = C12->get_migrate_angle();
								double C12_ma_x = cos(C12_ma*PI/180);
								double C12_ma_y = sin(C12_ma*PI/180);
								double plot_C12_ma_2 = -Ft21x_m*C12_ma_x -Ft21y_m*C12_ma_y;
								double plot_C12_ma_2_d = sqrt(Ft21x_m*Ft21x_m + Ft21y_m*Ft21y_m);
								if (plot_C12_ma_2>0.71*plot_C12_ma_2_d) // 45 degrees
								{
									C1_sense_migr_vec_x -= Ft21x_m;
									C1_sense_migr_vec_y -= Ft21y_m;					
									C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);					
									C11->set_sense_migr(true);
									if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<"; Dir_t2:"<<Dir_t2<<"; "<<Fs11x_m<<","<<Fs11y_m<<"; "<<Ft21x_m<<","<<Ft21y_m<<endl;}
								}
							}
						}
					}
				}
				else if (s11->get_migrate_mark() && !t21->get_migrate_mark())
				{
					if (sense_migration_flag==1)
					{
						double Dir_s1 = -1;
						double plot_s1 = Fs11x_m*normal_s1_x + Fs11y_m*normal_s1_y;
						double plot_t2 = -Ft21x_m*normal_s1_x -Ft21y_m*normal_s1_y;
						double F_m_d = sqrt(Fs11x_m*Fs11x_m + Fs11y_m*Fs11y_m);
						Dir_s1 = plot_s1/F_m_d;
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						double C11_ma = C11->get_migrate_angle();
						double C11_ma_x = cos(C11_ma*PI/180);
						double C11_ma_y = sin(C11_ma*PI/180);
						double C11_ma_2 = Ft21x_m*C11_ma_x + Ft21y_m*C11_ma_y;
						double C11_ma_2_d = sqrt(Ft21x_m*Ft21x_m + Ft21y_m*Ft21y_m);
						if (C11_ma_2>0.71*C11_ma_2_d && Dir_s1<0.25 && plot_t2>0) // 15 degree
						{
							C2_sense_migr_vec_x += Ft21x_m;
							C2_sense_migr_vec_y += Ft21y_m;
							C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
							C12->set_sense_migr(true);
							if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<": m_angle:"<<C11->get_migrate_angle()<<"; Dir_s1:"<<Dir_s1<<"; "<<Fs11x_m<<","<<Fs11y_m<<"; "<<Ft21x_m<<","<<Ft21y_m<<endl;}
						}
					}
				}
				else if (!s11->get_migrate_mark() && t21->get_migrate_mark())
				{
					if (sense_migration_flag==1)
					{
						double Dir_t2 = -1;
						double plot_s1 = Fs11x_m*normal_s1_x + Fs11y_m*normal_s1_y;
						double plot_t2 = -Ft21x_m*normal_s1_x -Ft21y_m*normal_s1_y;
						double F_m_d = sqrt(Ft21x_m*Ft21x_m + Ft21y_m*Ft21y_m);
						Dir_t2 = plot_t2/F_m_d;
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C12_ma = C12->get_migrate_angle();
						double C12_ma_x = cos(C12_ma*PI/180);
						double C12_ma_y = sin(C12_ma*PI/180);
						double C12_ma_1 = Fs11x_m*C12_ma_x + Fs11y_m*C12_ma_y;
						double C12_ma_1_d = sqrt(Fs11x_m*Fs11x_m + Fs11y_m*Fs11y_m);
						if (C12_ma_1>0.71*C12_ma_1_d && Dir_t2<0.25 && plot_s1>0) // 15 degree
						{
							C1_sense_migr_vec_x += Fs11x_m;
							C1_sense_migr_vec_y += Fs11y_m;
							C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
							C11->set_sense_migr(true);
							if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<": m_angle:"<<C12->get_migrate_angle()<<"; Dir_t2:"<<Dir_t2<<"; "<<Fs11x_m<<","<<Fs11y_m<<"; "<<Ft21x_m<<","<<Ft21y_m<<endl;}
						}
					}
				}
			}
			if (t1mn==1)
			{
				double nt1x = t11->getba()->getnormal(0);
				double nt1y = t11->getba()->getnormal(1);
				double nt1d = sqrt(nt1x*nt1x + nt1y*nt1y);
				double normal_t1_x = nt1x/nt1d;
				double normal_t1_y = nt1y/nt1d;
				if      ((s21->get_pop_out()>0 && t11->get_pop_out()>0)  ||
						 (s21->get_pop_out()>0 && t11->get_pop_out()==0) ||
						 (s21->get_pop_out()==0 && t11->get_pop_out()>0)) 
				{
					Mp->set_F_s2(0);
					Mp->set_F_t1(0);
					s21->set_stored_force_vis(0,0);
					t11->set_stored_force_vis(0,0);
				}
				else if (s21->get_pop_out()==0 && t11->get_pop_out()==0)
				{
					double sF_t1 = Mp->get_F_t1();
					double sF_s2 = Mp->get_F_s2();
					double sF_t1_new = sF_t1;// + plot_t1;
					double sF_s2_new = sF_s2;// + plot_s2;
					int Dir_t1 = -1;
					int Dir_s2 = -1;
					double plot_t1 = Ft11x*normal_t1_x + Ft11y*normal_t1_y;
					double plot_s2 = -Fs21x*normal_t1_x -Fs21y*normal_t1_y;
					if (plot_t1>0) {Dir_t1 = 1;}
					if (plot_s2>0) {Dir_s2 = 1;}
					if (Dir_t1>0) {sF_t1_new += plot_t1;} else {sF_t1_new = 0;}
					if (Dir_s2>0) {sF_s2_new += plot_s2;} else {sF_s2_new = 0;}
					Mp->set_F_s2(sF_s2_new);
					Mp->set_F_t1(sF_t1_new);
					double s_f_s2_x = -sF_s2_new*normal_t1_x;
					double s_f_s2_y = -sF_s2_new*normal_t1_y;
					double s_f_t1_x = sF_t1_new*normal_t1_x;
					double s_f_t1_y = sF_t1_new*normal_t1_y;
					s21->set_stored_force_vis(s_f_s2_x,s_f_s2_y);
					t11->set_stored_force_vis(s_f_t1_x,s_f_t1_y);
					if (sF_s2_new>F_adhesion || sF_t1_new>F_adhesion) {s21_ab = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (Dir_t1>0)
					{
						C1_sense_mech_vec_x += Ft11x;
						C1_sense_mech_vec_y += Ft11y;
						C2_sense_mech_vec_x += Fs21x;
						C2_sense_mech_vec_y += Fs21y;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C11->set_sense_mech(true);
						C12->set_sense_mech(true);
					}
				}
				else if (s21->get_pop_out()<0 && t11->get_pop_out()>0)
				{
					double sF_t1 = Mp->get_F_t1();
					double sF_s2 = Mp->get_F_s2();
					double plot_t1 = abs(Ft11x*normal_t1_x + Ft11y*normal_t1_y);
					double plot_s2 = abs(-Fs21x*normal_t1_x -Fs21y*normal_t1_y);
					double sF_t1_new = sF_t1 + plot_t1;
					double sF_s2_new = sF_s2 + plot_s2;
					double diff_s2_t1 = sF_s2_new - sF_t1_new;
					if (diff_s2_t1>0) {sF_s2_new = diff_s2_t1;sF_t1_new = diff_s2_t1;}
					else {sF_s2_new = 0;sF_t1_new = 0;}
					Mp->set_F_s2(sF_s2_new);
					Mp->set_F_t1(sF_t1_new);
					double s_f_s2_x = sF_s2_new*normal_t1_x;
					double s_f_s2_y = sF_s2_new*normal_t1_y;
					s21->set_stored_force_vis(s_f_s2_x,s_f_s2_y);
					t11->set_stored_force_vis(0,0);
					if (sF_s2_new>F_adhesion || sF_t1_new>F_adhesion) {s21_ab = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					if (diff_s2_t1>0)
					{
						double Fs21_d = sqrt(Fs21x*Fs21x + Fs21y*Fs21y);
						C1_sense_mech_vec_x += diff_s2_t1*Fs21x/Fs21_d;
						C1_sense_mech_vec_y += diff_s2_t1*Fs21y/Fs21_d;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
				}
				else if (s21->get_pop_out()>0 && t11->get_pop_out()<0)
				{
					double sF_t1 = Mp->get_F_t1();
					double sF_s2 = Mp->get_F_s2();
					double plot_t1 = abs(Ft11x*normal_t1_x + Ft11y*normal_t1_y);
					double plot_s2 = abs(-Fs21x*normal_t1_x -Fs21y*normal_t1_y);
					double sF_t1_new = sF_t1 + plot_t1;
					double sF_s2_new = sF_s2 + plot_s2;
					double diff_t1_s2 = sF_t1_new - sF_s2_new;
					if (diff_t1_s2>0) {sF_t1_new = diff_t1_s2;sF_s2_new = diff_t1_s2;}
					else {sF_t1_new = 0;sF_s2_new = 0;}
					Mp->set_F_s2(sF_s2_new);
					Mp->set_F_t1(sF_t1_new);
					double s_f_t1_x = -sF_t1_new*normal_t1_x;
					double s_f_t1_y = -sF_t1_new*normal_t1_y;
					s21->set_stored_force_vis(0,0);
					t11->set_stored_force_vis(s_f_t1_x,s_f_t1_y);
					if (sF_s2_new>F_adhesion || sF_t1_new>F_adhesion) {s21_ab = true;}
					/////
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (diff_t1_s2>0)
					{
						double Ft11_d = sqrt(Ft11x*Ft11x + Ft11y*Ft11y);
						C2_sense_mech_vec_x += diff_t1_s2*Ft11x/Ft11_d;
						C2_sense_mech_vec_y += diff_t1_s2*Ft11y/Ft11_d;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C12->set_sense_mech(true);
					}
				}
				else if (s21->get_pop_out()<0 && t11->get_pop_out()==0)
				{
					double sF_t1 = Mp->get_F_t1();
					double sF_s2 = Mp->get_F_s2();
					double sF_t1_new = sF_t1;// + plot_t1;
					double sF_s2_new = sF_s2;// + plot_s2;
					int Dir_s2 = -1;
					double plot_t1 = abs(Ft11x*normal_t1_x + Ft11y*normal_t1_y);
					double plot_s2 = -Fs21x*normal_t1_x -Fs21y*normal_t1_y;
					if (plot_s2>0) {Dir_s2 = 1;}
					sF_t1_new += plot_t1;
					if (Dir_s2>0) {sF_s2_new += plot_s2;} else {sF_s2_new = 0;}
					Mp->set_F_s2(sF_s2_new);
					Mp->set_F_t1(sF_t1_new);
					double s_f_s2_x = sF_s2_new*normal_t1_x;
					double s_f_s2_y = sF_s2_new*normal_t1_y;
					double s_f_t1_x = sF_t1_new*normal_t1_x;
					double s_f_t1_y = sF_t1_new*normal_t1_y;
					s21->set_stored_force_vis(s_f_s2_x,s_f_s2_y);
					t11->set_stored_force_vis(s_f_t1_x,s_f_t1_y);
					if (sF_s2_new>F_adhesion || sF_t1_new>F_adhesion) {s21_ab = true;}
					/////
					double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
					double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
					if (Dir_s2>0)
					{
						C1_sense_mech_vec_x += Ft11x;
						C1_sense_mech_vec_y += Ft11y;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
				}
				else if (s21->get_pop_out()==0 && t11->get_pop_out()<0)
				{
					double sF_t1 = Mp->get_F_t1();
					double sF_s2 = Mp->get_F_s2();
					double sF_t1_new = sF_t1;// + plot_t1;
					double sF_s2_new = sF_s2;// + plot_s2;
					int Dir_t1 = -1;
					double plot_t1 = Ft11x*normal_t1_x + Ft11y*normal_t1_y;
					double plot_s2 = abs(-Fs21x*normal_t1_x -Fs21y*normal_t1_y);
					if (plot_t1>0) {Dir_t1 = 1;}
					if (Dir_t1>0) {sF_t1_new += plot_t1;} else {sF_t1_new = 0;}
					sF_s2_new += plot_s2;
					Mp->set_F_s2(sF_s2_new);
					Mp->set_F_t1(sF_t1_new);
					double s_f_s2_x = -sF_s2_new*normal_t1_x;
					double s_f_s2_y = -sF_s2_new*normal_t1_y;
					double s_f_t1_x = -sF_t1_new*normal_t1_x;
					double s_f_t1_y = -sF_t1_new*normal_t1_y;
					s21->set_stored_force_vis(s_f_s2_x,s_f_s2_y);
					t11->set_stored_force_vis(s_f_t1_x,s_f_t1_y);
					if (sF_s2_new>F_adhesion || sF_t1_new>F_adhesion) {s21_ab = true;}
					/////
					double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
					double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
					if (Dir_t1>0)
					{
						C2_sense_mech_vec_x += Fs21x;
						C2_sense_mech_vec_y += Fs21y;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
						C12->set_sense_mech(true);
					}
				}
				else if (s21->get_pop_out()<0 && t11->get_pop_out()<0)
				{
					double sF_t1 = Mp->get_F_t1();
					double sF_s2 = Mp->get_F_s2();
					double plot_t1 = abs(Ft11x*normal_t1_x + Ft11y*normal_t1_y);
					double plot_s2 = abs(-Fs21x*normal_t1_x -Fs21y*normal_t1_y);
					double sF_t1_new = sF_t1 + plot_t1;
					double sF_s2_new = sF_s2 + plot_s2;
					Mp->set_F_s2(sF_s2_new);
					Mp->set_F_t1(sF_t1_new);
					double s_f_s2_x = sF_s2_new*normal_t1_x;
					double s_f_s2_y = sF_s2_new*normal_t1_y;
					double s_f_t1_x = -sF_t1_new*normal_t1_x;
					double s_f_t1_y = -sF_t1_new*normal_t1_y;
					s21->set_stored_force_vis(s_f_s2_x,s_f_s2_y);
					t11->set_stored_force_vis(s_f_t1_x,s_f_t1_y);
					if (sF_s2_new>F_adhesion || sF_t1_new>F_adhesion) {s21_ab = true;}
				}
				///
				if (!s21->get_migrate_mark() && !t11->get_migrate_mark())
				{
					if      (sense_migration_flag==0)
					{
						double C1_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						if (C11->Migrate())
						{
							double C11cx = C11->get_center(0);
							double C11cy = C11->get_center(1);
							double C12cx = C12->get_center(0);
							double C12cy = C12->get_center(1);
							double m21x = C11cx - C12cx;
							double m21y = C11cy - C12cy;
							double m21d = sqrt(m21x*m21x + m21y*m21y);
							C2_sense_migr_vec_x += m21x/m21d;
							C2_sense_migr_vec_y += m21y/m21d;
							C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
							C12->set_sense_migr(true);
						}
						if (C12->Migrate())
						{
							double C11cx = C11->get_center(0);
							double C11cy = C11->get_center(1);
							double C12cx = C12->get_center(0);
							double C12cy = C12->get_center(1);
							double m12x = C12cx - C11cx;
							double m12y = C12cy - C11cy;
							double m12d = sqrt(m12x*m12x + m12y*m12y);
							C1_sense_migr_vec_x += m12x/m12d;
							C1_sense_migr_vec_y += m12y/m12d;
							C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
							C11->set_sense_migr(true);
						}
					}
					else if (sense_migration_flag==1)
					{
						int Dir_s2 = -1;
						int Dir_t1 = -1;
						double plot_t1 = Ft11x_m*normal_t1_x + Ft11y_m*normal_t1_y;
						double plot_s2 = -Fs21x_m*normal_t1_x -Fs21y_m*normal_t1_y;
						if (plot_t1>0) {Dir_t1 = 1;}
						if (plot_s2>0) {Dir_s2 = 1;}
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						if (Dir_s2>0)
						{
							/// strategy 2:
							if (C11->Migrate())
							{
								double C11_ma = C11->get_migrate_angle();
								double C11_ma_x = cos(C11_ma*PI/180);
								double C11_ma_y = sin(C11_ma*PI/180);
								double C11_ma_1 = -Ft11x_m*C11_ma_x -Ft11y_m*C11_ma_y;
								double C11_ma_1_d = sqrt(Ft11x_m*Ft11x_m + Ft11y_m*Ft11y_m);
								if (C11_ma_1>0.71*C11_ma_1_d) // 45 degrees
								{
									C2_sense_migr_vec_x -= Ft11x_m;
									C2_sense_migr_vec_y -= Ft11y_m;
									C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
									C12->set_sense_migr(true);
									if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<"; Dir_s2:"<<Dir_s2<<"; "<<Ft11x_m<<","<<Ft11y_m<<"; "<<Fs21x_m<<","<<Fs21y_m<<endl;}
								}
							}
						}
						if (Dir_t1>0)
						{
							/// strategy 2:
							if (C12->Migrate())
							{
								double C12_ma = C12->get_migrate_angle();
								double C12_ma_x = cos(C12_ma*PI/180);
								double C12_ma_y = sin(C12_ma*PI/180);
								double C12_ma_2 = -Fs21x_m*C12_ma_x -Fs21y_m*C12_ma_y;
								double C12_ma_2_d = sqrt(Fs21x_m*Fs21x_m + Fs21y_m*Fs21y_m);
								if (C12_ma_2>0.71*C12_ma_2_d) // 45 degrees
								{
									C1_sense_migr_vec_x -= Fs21x_m;
									C1_sense_migr_vec_y -= Fs21y_m;						
									C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);						
									C11->set_sense_migr(true);
									if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<"; Dir_t1:"<<Dir_t1<<"; "<<Ft11x_m<<","<<Ft11y_m<<"; "<<Fs21x_m<<","<<Fs21y_m<<endl;}
								}
							}
						}
					}
				}
				else if (s21->get_migrate_mark() && !t11->get_migrate_mark())
				{
					if (sense_migration_flag==1)
					{
						double Dir_s2 = -1;
						double plot_t1 = Ft11x_m*normal_t1_x + Ft11y_m*normal_t1_y;
						double plot_s2 = -Fs21x_m*normal_t1_x -Fs21y_m*normal_t1_y;
						double F_m_d = sqrt(Fs21x_m*Fs21x_m + Fs21y_m*Fs21y_m);
						Dir_s2 = plot_s2/F_m_d;
						double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
						double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
						double C12_ma = C12->get_migrate_angle();
						double C12_ma_x = cos(C12_ma*PI/180);
						double C12_ma_y = sin(C12_ma*PI/180);
						double C12_ma_1 = Ft11x_m*C12_ma_x + Ft11y_m*C12_ma_y;
						double C12_ma_1_d = sqrt(Ft11x_m*Ft11x_m + Ft11y_m*Ft11y_m);
						if (C12_ma_1>0.71*C12_ma_1_d && Dir_s2<0.25 && plot_t1>0)
						{
							C1_sense_migr_vec_x += Ft11x_m;
							C1_sense_migr_vec_y += Ft11y_m;
							C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
							C11->set_sense_migr(true);
							if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<"; m_angle: "<<C12->get_migrate_angle()<<"; Dir_s2: "<<Dir_s2<<"; "<<Ft11x_m<<","<<Ft11y_m<<"; "<<Fs21x_m<<","<<Fs21y_m<<endl;}
						}
					}
				}
				else if (!s21->get_migrate_mark() && t11->get_migrate_mark())
				{
					if (sense_migration_flag==1)
					{
						double Dir_t1 = -1;
						double plot_t1 = Ft11x_m*normal_t1_x + Ft11y_m*normal_t1_y;
						double plot_s2 = -Fs21x_m*normal_t1_x -Fs21y_m*normal_t1_y;
						double F_m_d = sqrt(Ft11x_m*Ft11x_m + Ft11y_m*Ft11y_m);
						Dir_t1 = plot_t1/F_m_d;
						double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
						double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
						double C11_ma = C11->get_migrate_angle();
						double C11_ma_x = cos(C11_ma*PI/180);
						double C11_ma_y = sin(C11_ma*PI/180);
						double C11_ma_2 = Fs21x_m*C11_ma_x + Fs21y_m*C11_ma_y;
						double C11_ma_2_d = sqrt(Fs21x_m*Fs21x_m + Fs21y_m*Fs21y_m);
						if (C11_ma_2>0.71*C11_ma_2_d && Dir_t1<0.25 && plot_s2>0) // 15 degree
						{
							C2_sense_migr_vec_x += Fs21x_m;
							C2_sense_migr_vec_y += Fs21y_m;
							C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
							C12->set_sense_migr(true);
							if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> cells "<<C11->id()<<","<<C12->id()<<"; m_angle: "<<C11->get_migrate_angle()<<" ;Dir_t1: "<<Dir_t1<<"; "<<Ft11x_m<<","<<Ft11y_m<<"; "<<Fs21x_m<<","<<Fs21y_m<<endl;}
						}
					}
				}
			}
			if (s11->getfp()==t11)
			{
				if      (s11_ab || s21_ab && s1mn==1 && t1mn==1) {edgebreak0.push_back(s11->getab());}
				else if (s11_ab && !s21_ab && s1mn==1)           {edgebreak0.push_back(s11->getab());}
				else if (!s11_ab && s21_ab && t1mn==1)           {edgebreak0.push_back(s21->getab());}
			}
			else if (s11->getfp()!=t11)
			{
				if (s1mn==1 && s11_ab) {edgebreak0.push_back(s11->getab());}
				if (t1mn==1 && s21_ab) {edgebreak0.push_back(s21->getab());}
				/////// internal forces recovered for points between endpoints should also be calculated //////
				point3D *ptemp = s11->getfp();
				point3D *ptemm = ptemp->getmp()->at(0);
				while (ptemp!=t11)
				{
					int at_np = (int)ptemp->get_inner_at()->size();
					int at_nm = (int)ptemm->get_inner_at()->size();
					double Fptempx = 0;
					double Fptempy = 0;
					double Fptemmx = 0;
					double Fptemmy = 0;
					double Fptempx_m = 0;
					double Fptempy_m = 0;
					double Fptemmx_m = 0;
					double Fptemmy_m = 0;
					if (!C11->Soften())
					{
						for (int j=0;j<at_np;j++)
						{
							triangle *Tj = ptemp->get_inner_at()->at(j);
							point3D *TAj = Tj->getA();
							point3D *TBj = Tj->getB();
							point3D *TCj = Tj->getC();
							int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
							int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
							int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
							double u0 = Xn[TAj_sm_id*2];
							double u1 = Xn[TAj_sm_id*2+1];
							double u2 = Xn[TBj_sm_id*2];
							double u3 = Xn[TBj_sm_id*2+1];
							double u4 = Xn[TCj_sm_id*2];
							double u5 = Xn[TCj_sm_id*2+1];
							double u0_m = Xn_p[TAj_sm_id*2];
							double u1_m = Xn_p[TAj_sm_id*2+1];
							double u2_m = Xn_p[TBj_sm_id*2];
							double u3_m = Xn_p[TBj_sm_id*2+1];
							double u4_m = Xn_p[TCj_sm_id*2];
							double u5_m = Xn_p[TCj_sm_id*2+1];
							if      (TAj==ptemp)
							{
								double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
								double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
								double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
								double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
								double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
								double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
								Fptempx += FAx;
								Fptempy += FAy;
								Fptempx_m += FAx_m;
								Fptempy_m += FAy_m;
							}
							else if (TBj==ptemp)
							{
								double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
								double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
								double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
								double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
								double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
								double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
								Fptempx += FBx;
								Fptempy += FBy;
								Fptempx_m += FBx_m;
								Fptempy_m += FBy_m;
							}
							else if (TCj==ptemp)
							{
								double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
								double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
								double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
								double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
								double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
								double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
								Fptempx += FCx;
								Fptempy += FCy;
								Fptempx_m += FCx_m;
								Fptempy_m += FCy_m;
							}
						}
					}
					else
					{
						int N1_sm_id = ptemp->get_node_id()->SM_ID();
						int N2_sm_id = ptemp->getrp()->get_node_id()->SM_ID();
						int N3_sm_id = ptemp->getfp()->get_node_id()->SM_ID();
						double u0 = Xn[N1_sm_id*2];
						double u1 = Xn[N1_sm_id*2+1];
						double u2 = Xn[N2_sm_id*2];
						double u3 = Xn[N2_sm_id*2+1];
						double u4 = Xn[N3_sm_id*2];
						double u5 = Xn[N3_sm_id*2+1];
						double u0_m = Xn_p[N1_sm_id*2];
						double u1_m = Xn_p[N1_sm_id*2+1];
						double u2_m = Xn_p[N2_sm_id*2];
						double u3_m = Xn_p[N2_sm_id*2+1];
						double u4_m = Xn_p[N3_sm_id*2];
						double u5_m = Xn_p[N3_sm_id*2+1];
						double v00 = ptemp->getab()->get_soften_triangle_stiffness_matrix(0,0);
						double v01 = ptemp->getab()->get_soften_triangle_stiffness_matrix(0,1);
						double v02 = ptemp->getab()->get_soften_triangle_stiffness_matrix(0,2);
						double v03 = ptemp->getab()->get_soften_triangle_stiffness_matrix(0,3);
						double v10 = ptemp->getab()->get_soften_triangle_stiffness_matrix(1,0);
						double v11 = ptemp->getab()->get_soften_triangle_stiffness_matrix(1,1);
						double v12 = ptemp->getab()->get_soften_triangle_stiffness_matrix(1,2);
						double v13 = ptemp->getab()->get_soften_triangle_stiffness_matrix(1,3);
						double v20 = ptemp->getba()->get_soften_triangle_stiffness_matrix(2,0);
						double v21 = ptemp->getba()->get_soften_triangle_stiffness_matrix(2,1);
						double v22 = ptemp->getba()->get_soften_triangle_stiffness_matrix(2,2);
						double v23 = ptemp->getba()->get_soften_triangle_stiffness_matrix(2,3);
						double v30 = ptemp->getba()->get_soften_triangle_stiffness_matrix(3,0);
						double v31 = ptemp->getba()->get_soften_triangle_stiffness_matrix(3,1);
						double v32 = ptemp->getba()->get_soften_triangle_stiffness_matrix(3,2);
						double v33 = ptemp->getba()->get_soften_triangle_stiffness_matrix(3,3);
						double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
						double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
						double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
						double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
						double Frx_m = v20*u2_m + v21*u3_m + v22*u0_m + v23*u1_m;
						double Fry_m = v30*u2_m + v31*u3_m + v32*u0_m + v33*u1_m;
						double Ffx_m = v00*u0_m + v01*u1_m + v02*u4_m + v03*u5_m;
						double Ffy_m = v10*u0_m + v11*u1_m + v12*u4_m + v13*u5_m;
						Fptempx += (Frx + Ffx);
						Fptempy += (Fry + Ffy);
						Fptempx_m += (Frx_m + Ffx_m);
						Fptempy_m += (Fry_m + Ffy_m);
					}
					if (!C12->Soften())
					{
						for (int j=0;j<at_nm;j++)
						{
							triangle *Tj = ptemm->get_inner_at()->at(j);
							point3D *TAj = Tj->getA();
							point3D *TBj = Tj->getB();
							point3D *TCj = Tj->getC();
							int TAj_sm_id = 0;if (TAj->get_cell_boundary()) {TAj_sm_id = TAj->get_node_id()->SM_ID();} else {TAj_sm_id = TAj->get_SM_ID();}
							int TBj_sm_id = 0;if (TBj->get_cell_boundary()) {TBj_sm_id = TBj->get_node_id()->SM_ID();} else {TBj_sm_id = TBj->get_SM_ID();}
							int TCj_sm_id = 0;if (TCj->get_cell_boundary()) {TCj_sm_id = TCj->get_node_id()->SM_ID();} else {TCj_sm_id = TCj->get_SM_ID();}
							double u0 = Xn[TAj_sm_id*2];
							double u1 = Xn[TAj_sm_id*2+1];
							double u2 = Xn[TBj_sm_id*2];
							double u3 = Xn[TBj_sm_id*2+1];
							double u4 = Xn[TCj_sm_id*2];
							double u5 = Xn[TCj_sm_id*2+1];
							double u0_m = Xn_p[TAj_sm_id*2];
							double u1_m = Xn_p[TAj_sm_id*2+1];
							double u2_m = Xn_p[TBj_sm_id*2];
							double u3_m = Xn_p[TBj_sm_id*2+1];
							double u4_m = Xn_p[TCj_sm_id*2];
							double u5_m = Xn_p[TCj_sm_id*2+1];
							if      (TAj==ptemm)
							{
								double v00 = Tj->get_Stiffness_matrix(0,0);double v01 = Tj->get_Stiffness_matrix(0,1);double v02 = Tj->get_Stiffness_matrix(0,2);double v03 = Tj->get_Stiffness_matrix(0,3);double v04 = Tj->get_Stiffness_matrix(0,4);double v05 = Tj->get_Stiffness_matrix(0,5);
								double v10 = Tj->get_Stiffness_matrix(1,0);double v11 = Tj->get_Stiffness_matrix(1,1);double v12 = Tj->get_Stiffness_matrix(1,2);double v13 = Tj->get_Stiffness_matrix(1,3);double v14 = Tj->get_Stiffness_matrix(1,4);double v15 = Tj->get_Stiffness_matrix(1,5);
								double FAx = v00*u0 + v01*u1 + v02*u2 + v03*u3 + v04*u4 + v05*u5;
								double FAy = v10*u0 + v11*u1 + v12*u2 + v13*u3 + v14*u4 + v15*u5;
								double FAx_m = v00*u0_m + v01*u1_m + v02*u2_m + v03*u3_m + v04*u4_m + v05*u5_m;
								double FAy_m = v10*u0_m + v11*u1_m + v12*u2_m + v13*u3_m + v14*u4_m + v15*u5_m;
								Fptemmx += FAx;
								Fptemmy += FAy;
								Fptemmx_m += FAx_m;
								Fptemmy_m += FAy_m;
							}
							else if (TBj==ptemm)
							{
								double v20 = Tj->get_Stiffness_matrix(2,0);double v21 = Tj->get_Stiffness_matrix(2,1);double v22 = Tj->get_Stiffness_matrix(2,2);double v23 = Tj->get_Stiffness_matrix(2,3);double v24 = Tj->get_Stiffness_matrix(2,4);double v25 = Tj->get_Stiffness_matrix(2,5);
								double v30 = Tj->get_Stiffness_matrix(3,0);double v31 = Tj->get_Stiffness_matrix(3,1);double v32 = Tj->get_Stiffness_matrix(3,2);double v33 = Tj->get_Stiffness_matrix(3,3);double v34 = Tj->get_Stiffness_matrix(3,4);double v35 = Tj->get_Stiffness_matrix(3,5);
								double FBx = v20*u0 + v21*u1 + v22*u2 + v23*u3 + v24*u4 + v25*u5;
								double FBy = v30*u0 + v31*u1 + v32*u2 + v33*u3 + v34*u4 + v35*u5;
								double FBx_m = v20*u0_m + v21*u1_m + v22*u2_m + v23*u3_m + v24*u4_m + v25*u5_m;
								double FBy_m = v30*u0_m + v31*u1_m + v32*u2_m + v33*u3_m + v34*u4_m + v35*u5_m;
								Fptemmx += FBx;
								Fptemmy += FBy;
								Fptemmx_m += FBx_m;
								Fptemmy_m += FBy_m;
							}
							else if (TCj==ptemp)
							{
								double v40 = Tj->get_Stiffness_matrix(4,0);double v41 = Tj->get_Stiffness_matrix(4,1);double v42 = Tj->get_Stiffness_matrix(4,2);double v43 = Tj->get_Stiffness_matrix(4,3);double v44 = Tj->get_Stiffness_matrix(4,4);double v45 = Tj->get_Stiffness_matrix(4,5);
								double v50 = Tj->get_Stiffness_matrix(5,0);double v51 = Tj->get_Stiffness_matrix(5,1);double v52 = Tj->get_Stiffness_matrix(5,2);double v53 = Tj->get_Stiffness_matrix(5,3);double v54 = Tj->get_Stiffness_matrix(5,4);double v55 = Tj->get_Stiffness_matrix(5,5);
								double FCx = v40*u0 + v41*u1 + v42*u2 + v43*u3 + v44*u4 + v45*u5;
								double FCy = v50*u0 + v51*u1 + v52*u2 + v53*u3 + v54*u4 + v55*u5;
								double FCx_m = v40*u0_m + v41*u1_m + v42*u2_m + v43*u3_m + v44*u4_m + v45*u5_m;
								double FCy_m = v50*u0_m + v51*u1_m + v52*u2_m + v53*u3_m + v54*u4_m + v55*u5_m;
								Fptemmx += FCx;
								Fptemmy += FCy;
								Fptemmx_m += FCx_m;
								Fptemmy_m += FCy_m;
							}
						}
					}
					else
					{
						int N1_sm_id = ptemm->get_node_id()->SM_ID();
						int N2_sm_id = ptemm->getrp()->get_node_id()->SM_ID();
						int N3_sm_id = ptemm->getfp()->get_node_id()->SM_ID();
						double u0 = Xn[N1_sm_id*2];
						double u1 = Xn[N1_sm_id*2+1];
						double u2 = Xn[N2_sm_id*2];
						double u3 = Xn[N2_sm_id*2+1];
						double u4 = Xn[N3_sm_id*2];
						double u5 = Xn[N3_sm_id*2+1];
						double u0_m = Xn_p[N1_sm_id*2];
						double u1_m = Xn_p[N1_sm_id*2+1];
						double u2_m = Xn_p[N2_sm_id*2];
						double u3_m = Xn_p[N2_sm_id*2+1];
						double u4_m = Xn_p[N3_sm_id*2];
						double u5_m = Xn_p[N3_sm_id*2+1];
						double v00 = ptemm->getab()->get_soften_triangle_stiffness_matrix(0,0);
						double v01 = ptemm->getab()->get_soften_triangle_stiffness_matrix(0,1);
						double v02 = ptemm->getab()->get_soften_triangle_stiffness_matrix(0,2);
						double v03 = ptemm->getab()->get_soften_triangle_stiffness_matrix(0,3);
						double v10 = ptemm->getab()->get_soften_triangle_stiffness_matrix(1,0);
						double v11 = ptemm->getab()->get_soften_triangle_stiffness_matrix(1,1);
						double v12 = ptemm->getab()->get_soften_triangle_stiffness_matrix(1,2);
						double v13 = ptemm->getab()->get_soften_triangle_stiffness_matrix(1,3);
						double v20 = ptemm->getba()->get_soften_triangle_stiffness_matrix(2,0);
						double v21 = ptemm->getba()->get_soften_triangle_stiffness_matrix(2,1);
						double v22 = ptemm->getba()->get_soften_triangle_stiffness_matrix(2,2);
						double v23 = ptemm->getba()->get_soften_triangle_stiffness_matrix(2,3);
						double v30 = ptemm->getba()->get_soften_triangle_stiffness_matrix(3,0);
						double v31 = ptemm->getba()->get_soften_triangle_stiffness_matrix(3,1);
						double v32 = ptemm->getba()->get_soften_triangle_stiffness_matrix(3,2);
						double v33 = ptemm->getba()->get_soften_triangle_stiffness_matrix(3,3);
						double Frx = v20*u2 + v21*u3 + v22*u0 + v23*u1;
						double Fry = v30*u2 + v31*u3 + v32*u0 + v33*u1;
						double Ffx = v00*u0 + v01*u1 + v02*u4 + v03*u5;
						double Ffy = v10*u0 + v11*u1 + v12*u4 + v13*u5;
						double Frx_m = v20*u2 + v21*u3 + v22*u0 + v23*u1;
						double Fry_m = v30*u2 + v31*u3 + v32*u0 + v33*u1;
						double Ffx_m = v00*u0 + v01*u1 + v02*u4 + v03*u5;
						double Ffy_m = v10*u0 + v11*u1 + v12*u4 + v13*u5;
						Fptemmx += (Frx + Ffx);
						Fptemmy += (Fry + Ffy);
						Fptemmx_m += (Frx_m + Ffx_m);
						Fptemmy_m += (Fry_m + Ffy_m);
					}
					ptemp->set_elastic_force(Fptempx, Fptempy);
					ptemm->set_elastic_force(Fptemmx, Fptemmy);
					///
					ptemp->set_elastic_force_migration(Fptempx_m, Fptempy_m);
					ptemm->set_elastic_force_migration(Fptemmx_m, Fptemmy_m);
					double ptempnx = ptemp->getab()->getnormal(0);
					double ptempny = ptemp->getab()->getnormal(1);
					double ptempnd = sqrt(ptempnx*ptempnx + ptempny*ptempny);
					double ptemp_normal_x = ptempnx/ptempnd;
					double ptemp_normal_y = ptempny/ptempnd;
					/***************************************
					s1    ptemp                t1
					o-->--X-->--o-->--o-->--o
					o--<--X--<--o--<--o--<--o
					t2    ptemm                s2
					***************************************/
					if ((ptemp->get_pop_out()>0 && ptemm->get_pop_out()>0) ||
						(ptemp->get_pop_out()>0 && ptemm->get_pop_out()==0) ||
						(ptemp->get_pop_out()==0 && ptemm->get_pop_out()>0))
						// ptemp-ptemm pair is still attached: both poppping outward or
						// one is static while the other is poping outward
					{
						ptemp->set_stored_apart_force(0);
						ptemm->set_stored_apart_force(0);
					}
					else if (ptemp->get_pop_out()==0 && ptemm->get_pop_out()==0)
						// ptemp-ptemm pair potentially breaks apart: both are static or shrinking backward or
						// one is static while the other is shrinking backward
					{
						double sF_ptemp = ptemp->get_stored_apart_force();
						double sF_ptemm = ptemm->get_stored_apart_force();
						double sF_ptemp_new = sF_ptemp;// + plot_ptemp;
						double sF_ptemm_new = sF_ptemm;// + plot_ptemm;
						int Dir_ptemp = -1;
						int Dir_ptemm = -1;
						double plot_ptemp = Fptempx*ptemp_normal_x + Fptempy*ptemp_normal_y;
						double plot_ptemm = -Fptemmx*ptemp_normal_x -Fptemmy*ptemp_normal_y;
						if (plot_ptemp>0) {Dir_ptemp = 1;}
						if (plot_ptemm>0) {Dir_ptemm = 1;}
						if (Dir_ptemp>0) {sF_ptemp_new += plot_ptemp;} else {sF_ptemp_new = 0;}
						if (Dir_ptemm>0) {sF_ptemm_new += plot_ptemm;} else {sF_ptemm_new = 0;}
						ptemp->set_stored_apart_force(sF_ptemp_new);
						ptemm->set_stored_apart_force(sF_ptemm_new);
						double s_f_p_x = sF_ptemp_new*ptemp_normal_x;
						double s_f_p_y = sF_ptemp_new*ptemp_normal_y;
						double s_f_m_x = -sF_ptemm_new*ptemp_normal_x;
						double s_f_m_y = -sF_ptemm_new*ptemp_normal_y;
						ptemp->set_stored_force_vis(s_f_p_x,s_f_p_y);
						ptemm->set_stored_force_vis(s_f_m_x,s_f_m_y);
						/////
						double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
						double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
						double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
						double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
						if (Dir_ptemp>0)
						{
							C1_sense_mech_vec_x += Fptempx;
							C1_sense_mech_vec_y += Fptempy;
							C2_sense_mech_vec_x += Fptemmx;
							C2_sense_mech_vec_y += Fptemmy;
							C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
							C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
							C11->set_sense_mech(true);
							C12->set_sense_mech(true);
						}
					}
					else if (ptemp->get_pop_out()<0 && ptemm->get_pop_out()>0)
					{
						double sF_ptemp = ptemp->get_stored_apart_force();
						double sF_ptemm = ptemm->get_stored_apart_force();
						double plot_ptemp = abs(Fptempx*ptemp_normal_x + Fptempy*ptemp_normal_y);
						double plot_ptemm = abs(-Fptemmx*ptemp_normal_x -Fptemmy*ptemp_normal_y);
						double sF_ptemp_new = sF_ptemp + plot_ptemp;
						double sF_ptemm_new = sF_ptemm + plot_ptemm;
						double diff_ptemp_ptemm = sF_ptemp_new - sF_ptemm_new;
						if (diff_ptemp_ptemm>0) {sF_ptemp_new = diff_ptemp_ptemm;sF_ptemm_new = diff_ptemp_ptemm;}
						else {sF_ptemp_new = 0;sF_ptemm_new = 0;}
						ptemp->set_stored_apart_force(sF_ptemp_new);
						ptemm->set_stored_apart_force(sF_ptemm_new);
						double s_f_p_x = -sF_ptemp_new*ptemp_normal_x;
						double s_f_p_y = -sF_ptemp_new*ptemp_normal_y;
						ptemp->set_stored_force_vis(s_f_p_x,s_f_p_y);
						ptemm->set_stored_force_vis(0,0);
						/////
						double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
						double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
						if (diff_ptemp_ptemm>0)
						{
							double Fptemp_d = sqrt(Fptempx*Fptempx + Fptempy*Fptempy);
							C2_sense_mech_vec_x += diff_ptemp_ptemm*Fptempx/Fptemp_d;
							C2_sense_mech_vec_y += diff_ptemp_ptemm*Fptempy/Fptemp_d;
							C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
							C12->set_sense_mech(true);
						}
					}
					else if (ptemp->get_pop_out()>0 && ptemm->get_pop_out()<0)
					{
						double sF_ptemp = ptemp->get_stored_apart_force();
						double sF_ptemm = ptemm->get_stored_apart_force();
						double plot_ptemp = abs(Fptempx*ptemp_normal_x + Fptempy*ptemp_normal_y);
						double plot_ptemm = abs(-Fptemmx*ptemp_normal_x -Fptemmy*ptemp_normal_y);
						double sF_ptemp_new = sF_ptemp + plot_ptemp;
						double sF_ptemm_new = sF_ptemm + plot_ptemm;
						double diff_ptemm_ptemp = sF_ptemm_new - sF_ptemp_new;
						if (diff_ptemm_ptemp>0) {sF_ptemp_new = diff_ptemm_ptemp;sF_ptemm_new = diff_ptemm_ptemp;}
						else {sF_ptemp_new = 0;sF_ptemm_new = 0;}
						ptemp->set_stored_apart_force(sF_ptemp_new);
						ptemm->set_stored_apart_force(sF_ptemm_new);
						double s_f_m_x = sF_ptemm_new*ptemp_normal_x;
						double s_f_m_y = sF_ptemm_new*ptemp_normal_y;
						ptemp->set_stored_force_vis(0,0);
						ptemm->set_stored_force_vis(s_f_m_x,s_f_m_y);
						/////
						double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
						double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
						if (diff_ptemm_ptemp>0)
						{
							double Fptemm_d = sqrt(Fptemmx*Fptemmx + Fptemmy*Fptemmy);
							C1_sense_mech_vec_x += diff_ptemm_ptemp*Fptemmx/Fptemm_d;
							C1_sense_mech_vec_y += diff_ptemm_ptemp*Fptemmy/Fptemm_d;
							C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
							C11->set_sense_mech(true);
						}
					}
					else if (ptemp->get_pop_out()<0 && ptemm->get_pop_out()==0)
					{
						double sF_ptemp = ptemp->get_stored_apart_force();
						double sF_ptemm = ptemm->get_stored_apart_force();
						double sF_ptemp_new = sF_ptemp;// + plot_ptemp;
						double sF_ptemm_new = sF_ptemm;// + plot_ptemm;
						int Dir_ptemm = -1;
						double plot_ptemp = abs(Fptempx*ptemp_normal_x + Fptempy*ptemp_normal_y);
						double plot_ptemm = -Fptemmx*ptemp_normal_x -Fptemmy*ptemp_normal_y;
						if (plot_ptemm>0) {Dir_ptemm = 1;}
						sF_ptemp_new += plot_ptemp;
						if (Dir_ptemm>0) {sF_ptemm_new += plot_ptemm;} else {sF_ptemm_new = 0;}
						ptemp->set_stored_apart_force(sF_ptemp_new);
						ptemm->set_stored_apart_force(sF_ptemm_new);
						double s_f_p_x = -sF_ptemp_new*ptemp_normal_x;
						double s_f_p_y = -sF_ptemp_new*ptemp_normal_y;
						double s_f_m_x = -sF_ptemm_new*ptemp_normal_x;
						double s_f_m_y = -sF_ptemm_new*ptemp_normal_y;
						ptemp->set_stored_force_vis(s_f_p_x,s_f_p_y);
						ptemm->set_stored_force_vis(s_f_m_x,s_f_m_y);
						/////
						double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
						double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
						if (Dir_ptemm>0)
						{
							C2_sense_mech_vec_x += Fptemmx;
							C2_sense_mech_vec_y += Fptemmy;
							C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);
							C12->set_sense_mech(true);
						}
					}
					else if (ptemp->get_pop_out()==0 && ptemm->get_pop_out()<0)
					{
						double sF_ptemp = ptemp->get_stored_apart_force();
						double sF_ptemm = ptemm->get_stored_apart_force();
						double sF_ptemp_new = sF_ptemp;// + plot_ptemp;
						double sF_ptemm_new = sF_ptemm;// + plot_ptemm;
						int Dir_ptemp = -1;
						double plot_ptemp = Fptempx*ptemp_normal_x + Fptempy*ptemp_normal_y;
						double plot_ptemm = abs(-Fptemmx*ptemp_normal_x -Fptemmy*ptemp_normal_y);
						if (plot_ptemp>0) {Dir_ptemp = 1;}
						if (Dir_ptemp>0) {sF_ptemp_new += plot_ptemp;} else {sF_ptemp_new = 0;}
						sF_ptemm_new += plot_ptemm;
						ptemp->set_stored_apart_force(sF_ptemp_new);
						ptemm->set_stored_apart_force(sF_ptemm_new);
						double s_f_p_x = sF_ptemp_new*ptemp_normal_x;
						double s_f_p_y = sF_ptemp_new*ptemp_normal_y;
						double s_f_m_x = sF_ptemm_new*ptemp_normal_x;
						double s_f_m_y = sF_ptemm_new*ptemp_normal_y;
						ptemp->set_stored_force_vis(s_f_p_x,s_f_p_y);
						ptemm->set_stored_force_vis(s_f_m_x,s_f_m_y);
						/////
						double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
						double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
						if (Dir_ptemp>0)
						{
							C1_sense_mech_vec_x += Fptempx;
							C1_sense_mech_vec_y += Fptempy;
							C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
							C11->set_sense_mech(true);
						}
					}
					else if (ptemp->get_pop_out()<0 && ptemm->get_pop_out()<0)
					{
						double sF_ptemp = ptemp->get_stored_apart_force();
						double sF_ptemm = ptemm->get_stored_apart_force();
						double plot_ptemp = abs(Fptempx*ptemp_normal_x + Fptempy*ptemp_normal_y);
						double plot_ptemm = abs(-Fptemmx*ptemp_normal_x -Fptemmy*ptemp_normal_y);
						double sF_ptemp_new = sF_ptemp + plot_ptemp;
						double sF_ptemm_new = sF_ptemm + plot_ptemm;
						ptemp->set_stored_apart_force(sF_ptemp_new);
						ptemm->set_stored_apart_force(sF_ptemm_new);
						double s_f_p_x = -sF_ptemp_new*ptemp_normal_x;
						double s_f_p_y = -sF_ptemp_new*ptemp_normal_y;
						double s_f_m_x = sF_ptemm_new*ptemp_normal_x;
						double s_f_m_y = sF_ptemm_new*ptemp_normal_y;
						ptemp->set_stored_force_vis(s_f_p_x,s_f_p_y);
						ptemm->set_stored_force_vis(s_f_m_x,s_f_m_y);
					}
					///
					if (!ptemp->get_migrate_mark() && !ptemm->get_migrate_mark())
					{
						if      (sense_migration_flag==0)
						{
							double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
							double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
							double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
							double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
							if (C11->Migrate())
							{
								double C11cx = C11->get_center(0);
								double C11cy = C11->get_center(1);
								double C12cx = C12->get_center(0);
								double C12cy = C12->get_center(1);
								double m21x = C11cx - C12cx;
								double m21y = C11cy - C12cy;
								double m21d = sqrt(m21x*m21x + m21y*m21y);
								C2_sense_migr_vec_x += m21x/m21d;
								C2_sense_migr_vec_y += m21y/m21d;
								C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C12->set_sense_migr(true);
							}
							if (C12->Migrate())
							{
								double C11cx = C11->get_center(0);
								double C11cy = C11->get_center(1);
								double C12cx = C12->get_center(0);
								double C12cy = C12->get_center(1);
								double m12x = C12cx - C11cx;
								double m12y = C12cy - C11cy;
								double m12d = sqrt(m12x*m12x + m12y*m12y);
								C1_sense_migr_vec_x += m12x/m12d;
								C1_sense_migr_vec_y += m12y/m12d;
								C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C11->set_sense_migr(true);
							}
						}
						else if (sense_migration_flag==1)
						{
							int Dir_ptemp = -1;
							int Dir_ptemm = -1;
							double plot_ptemp = Fptempx_m*ptemp_normal_x + Fptempy_m*ptemp_normal_y;
							double plot_ptemm = -Fptemmx_m*ptemp_normal_x -Fptemmy_m*ptemp_normal_y;
							if (plot_ptemp>0) {Dir_ptemp = 1;}
							if (plot_ptemm>0) {Dir_ptemm = 1;}
							double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
							double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
							double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
							double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
							if (Dir_ptemp>0)
							{
								/// strategy 2:
								if (C11->Migrate())
								{
									double C11_ma = C11->get_migrate_angle();
									double C11_ma_x = cos(C11_ma*PI/180);
									double C11_ma_y = sin(C11_ma*PI/180);
									double plot_C11_ma_p = -Fptempx_m*C11_ma_x -Fptempy_m*C11_ma_y;
									double plot_C11_ma_p_d = sqrt(Fptempx_m*Fptempx_m + Fptempy_m*Fptempy_m);
									if (plot_C11_ma_p>0.71*plot_C11_ma_p_d) // 45 degree
									{
										C2_sense_migr_vec_x -= Fptempx_m;
										C2_sense_migr_vec_y -= Fptempy_m;
										C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C12->set_sense_migr(true);
										if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> "<<C11->id()<<", "<<C12->id()<<" temp, temm: "<<Fptempx_m<<","<<Fptempy_m<<"; "<<Fptemmx_m<<","<<Fptemmy_m<<endl;}
									}
								}
							}
							if (Dir_ptemm>0)
							{
								/// strategy 2:
								if (C12->Migrate())
								{
									double C12_ma = C12->get_migrate_angle();
									double C12_ma_x = cos(C12_ma*PI/180);
									double C12_ma_y = sin(C12_ma*PI/180);
									double plot_C12_ma_m = -Fptemmx_m*C12_ma_x -Fptemmy_m*C12_ma_y;
									double plot_C12_ma_m_d = sqrt(Fptemmx_m*Fptemmx_m + Fptemmy_m*Fptemmy_m);
									if (plot_C12_ma_m>0.71*plot_C12_ma_m_d) // 45 degree
									{
										C1_sense_migr_vec_x -= Fptemmx_m;
										C1_sense_migr_vec_y -= Fptemmy_m;							
										C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
										C11->set_sense_migr(true);
										if (C11->id()==w_id1 || C12->id()==w_id1) {cout<<"  -> "<<C11->id()<<", "<<C12->id()<<" temp, temm: "<<Fptempx_m<<","<<Fptempy_m<<"; "<<Fptemmx_m<<","<<Fptemmy_m<<endl;}
									}
								}
							}
						}
					}
					///
					ptemp = ptemp->getfp();
					ptemm = ptemm->getrp();
				}
			}
		}
	}
	cout<<"  -> Two-endpoints done now for three endpoints check:\n";
	for (int i=0;i<nn;i++)
	{
		if (nodeList[i]->get_in_global())
		{
			int vertn = (int)nodeList[i]->getvertex()->size();
			if (vertn==3)
			{
				point3D *nv1 = nodeList[i]->getvertex()->at(0);
				point3D *nv2 = nodeList[i]->getvertex()->at(1);
				point3D *nv3 = nodeList[i]->getvertex()->at(2);
				cell *C21 = nv1->getab()->Cell();
				cell *C22 = nv2->getab()->Cell();
				cell *C23 = nv3->getab()->Cell();
				if (!C21->Migrate() && !C22->Migrate() && !C23->Migrate()) continue;
				edge *nv1_ab = nv1->getab();
				edge *nv2_ab = nv2->getab();
				edge *nv3_ab = nv3->getab();
				edge *nv1_ba = nv1->getba();
				edge *nv2_ba = nv2->getba();
				edge *nv3_ba = nv3->getba();
				bool nv1a = false;
				bool nv2a = false;
				bool nv3a = false;
				if (nv1_ab->get_attach() && nv1_ba->get_attach()) {nv1a = true;}
				if (nv2_ab->get_attach() && nv2_ba->get_attach()) {nv2a = true;}
				if (nv3_ab->get_attach() && nv3_ba->get_attach()) {nv3a = true;}
				if      (nv1a && nv2a && nv3a)
				{
					/*******************************
					o-<-oo      oo-<-o
						 \\ V1 //
						  \\  //
						   \\//
						 V2 oo V3
							||
							||
					    o->-oo->-o
					*******************************/
					point3D *V1 = nv1;
					point3D *V2 = nv1_ba->get_Neighbor()->p1();
					point3D *V3 = nv1_ab->get_Neighbor()->p2();
					cell *C11 = V1->getab()->Cell();
					cell *C12 = V2->getab()->Cell();
					cell *C13 = V3->getab()->Cell();
					interpair *Pl12 = nv1_ba->get_Pair();
					interpair *Pl13 = nv1_ab->get_Pair();
					interpair *Pl23 = V2->getba()->get_Pair();
					if (Pl12->MP()->size()==0 || Pl13->MP()->size()==0 || Pl23->MP()->size()==0)
					{
						cout<<"  -> error: 3 points joint: "<<V1->id()<<","<<V2->id()<<","<<V3->id()<<" of cell "<<V1->getab()->Cell()->id()<<","<<V2->getab()->Cell()->id()<<","<<V3->getab()->Cell()->id()<<endl;
						return V1->getab()->Cell()->id();
					}
					double break_force12 = Pl12->get_adhesion();
					double break_force13 = Pl13->get_adhesion();
					double break_force23 = Pl23->get_adhesion();
					double V1_fx = V1->get_elastic_force(0);
					double V1_fy = V1->get_elastic_force(1);
					double V2_fx = V2->get_elastic_force(0);
					double V2_fy = V2->get_elastic_force(1);
					double V3_fx = V3->get_elastic_force(0);
					double V3_fy = V3->get_elastic_force(1);
					///
					double V1_fx_m = V1->get_elastic_force_migration(0);
					double V1_fy_m = V1->get_elastic_force_migration(1);
					double V2_fx_m = V2->get_elastic_force_migration(0);
					double V2_fy_m = V2->get_elastic_force_migration(1);
					double V3_fx_m = V3->get_elastic_force_migration(0);
					double V3_fy_m = V3->get_elastic_force_migration(1);
					///
					double V1_ab_nx = V1->getab()->getnormal(0);
					double V1_ab_ny = V1->getab()->getnormal(1);
					double V2_ab_nx = V2->getab()->getnormal(0);
					double V2_ab_ny = V2->getab()->getnormal(1);
					double V3_ab_nx = V3->getab()->getnormal(0);
					double V3_ab_ny = V3->getab()->getnormal(1);
					double V1_ba_nx = V1->getba()->getnormal(0);
					double V1_ba_ny = V1->getba()->getnormal(1);
					double V2_ba_nx = V2->getba()->getnormal(0);
					double V2_ba_ny = V2->getba()->getnormal(1);
					double V3_ba_nx = V3->getba()->getnormal(0);
					double V3_ba_ny = V3->getba()->getnormal(1);
					double V1_ab_nd = sqrt(V1_ab_nx*V1_ab_nx + V1_ab_ny*V1_ab_ny);
					double V2_ab_nd = sqrt(V2_ab_nx*V2_ab_nx + V2_ab_ny*V2_ab_ny);
					double V3_ab_nd = sqrt(V3_ab_nx*V3_ab_nx + V3_ab_ny*V3_ab_ny);
					double V1_ba_nd = sqrt(V1_ba_nx*V1_ba_nx + V1_ba_ny*V1_ba_ny);
					double V2_ba_nd = sqrt(V2_ba_nx*V2_ba_nx + V2_ba_ny*V2_ba_ny);
					double V3_ba_nd = sqrt(V3_ba_nx*V3_ba_nx + V3_ba_ny*V3_ba_ny);
					double V1_ab_normal_x = V1_ab_nx/V1_ab_nd;
					double V1_ab_normal_y = V1_ab_ny/V1_ab_nd;
					double V2_ab_normal_x = V2_ab_nx/V2_ab_nd;
					double V2_ab_normal_y = V2_ab_ny/V2_ab_nd;
					double V3_ab_normal_x = V3_ab_nx/V3_ab_nd;
					double V3_ab_normal_y = V3_ab_ny/V3_ab_nd;
					double V1_ba_normal_x = V1_ba_nx/V1_ba_nd;
					double V1_ba_normal_y = V1_ba_ny/V1_ba_nd;
					double V2_ba_normal_x = V2_ba_nx/V2_ba_nd;
					double V2_ba_normal_y = V2_ba_ny/V2_ba_nd;
					double V3_ba_normal_x = V3_ba_nx/V3_ba_nd;
					double V3_ba_normal_y = V3_ba_ny/V3_ba_nd;
					bool V21_b = false; // V3 is protruding the boundary of V1 and V2
					bool V32_b = false; // V1 is protruding the boundary of V2 and V3
					bool V13_b = false; // V2 is protruding the boundary of V1 and V3
					double e32_x = V3->getfp()->x() - V3->x();
					double e32_y = V3->getfp()->y() - V3->y();
					double e32_d = sqrt(e32_x*e32_x + e32_y*e32_y);
					double e32_along_x = e32_x/e32_d;
					double e32_along_y = e32_y/e32_d;
					double plot_e32_along_V1 = V1_fx*e32_along_x + V1_fy*e32_along_y;
					double e13_x = V1->getfp()->x() - V1->x();
					double e13_y = V1->getfp()->y() - V1->y();
					double e13_d = sqrt(e13_x*e13_x + e13_y*e13_y);
					double e13_along_x = e13_x/e13_d;
					double e13_along_y = e13_y/e13_d;
					double plot_e13_along_V2 = V2_fx*e13_along_x + V2_fy*e13_along_y;
					double e21_x = V2->getfp()->x() - V2->x();
					double e21_y = V2->getfp()->y() - V2->y();
					double e21_d = sqrt(e21_x*e21_x + e21_y*e21_y);
					double e21_along_x = e21_x/e21_d;
					double e21_along_y = e21_y/e21_d;
					double plot_e21_along_V3 = V3_fx*e21_along_x + V3_fy*e21_along_y;
					double plot_V1_V1_ab = V1_fx*V1_ab_normal_x + V1_fy*V1_ab_normal_y;
					double plot_V1_V1_ba = V1_fx*V1_ba_normal_x + V1_fy*V1_ba_normal_y;
					double plot_V2_V2_ab = V2_fx*V2_ab_normal_x + V2_fy*V2_ab_normal_y;
					double plot_V2_V2_ba = V2_fx*V2_ba_normal_x + V2_fy*V2_ba_normal_y;
					double plot_V3_V3_ab = V3_fx*V3_ab_normal_x + V3_fy*V3_ab_normal_y;
					double plot_V3_V3_ba = V3_fx*V3_ba_normal_x + V3_fy*V3_ba_normal_y;
					double plot_V1_V1_ab_m = V1_fx_m*V1_ab_normal_x + V1_fy_m*V1_ab_normal_y;
					double plot_V1_V1_ba_m = V1_fx_m*V1_ba_normal_x + V1_fy_m*V1_ba_normal_y;
					double plot_V2_V2_ab_m = V2_fx_m*V2_ab_normal_x + V2_fy_m*V2_ab_normal_y;
					double plot_V2_V2_ba_m = V2_fx_m*V2_ba_normal_x + V2_fy_m*V2_ba_normal_y;
					double plot_V3_V3_ab_m = V3_fx_m*V3_ab_normal_x + V3_fy_m*V3_ab_normal_y;
					double plot_V3_V3_ba_m = V3_fx_m*V3_ba_normal_x + V3_fy_m*V3_ba_normal_y;
					bool V1_ab_b = false;
					bool V1_ba_b = false;
					bool V2_ab_b = false;
					bool V2_ba_b = false;
					bool V3_ab_b = false;
					bool V3_ba_b = false;
					double V1_ab_bf = 0; // the stored force to break edge of V1_ab
					double V1_ba_bf = 0; // the stored force to break edge of V1_ba
					double V2_ab_bf = 0; // the stored force to break edge of V2_ab
					double V2_ba_bf = 0; // the stored force to break edge of V2_ba
					double V3_ab_bf = 0; // the stored force to break edge of V3_ab
					double V3_ba_bf = 0; // the stored force to break edge of V3_ba
					double s_f_V1_x = 0;
					double s_f_V1_y = 0;
					double s_f_V2_x = 0;
					double s_f_V2_y = 0;
					double s_f_V3_x = 0;
					double s_f_V3_y = 0;
					/// for migration response ///
					if      (!V1->get_migrate_mark() && !V2->get_migrate_mark() && !V3->get_migrate_mark())
					{
						if (sense_migration_flag==0)
						{
							double C3_sense_migr_vec_x = C13->get_sense_migr_vec(0);
							double C3_sense_migr_vec_y = C13->get_sense_migr_vec(1);
							double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
							double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
							double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
							double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
							if (C11->Migrate())
							{
								double C11cx = C11->get_center(0);
								double C11cy = C11->get_center(1);
								double C12cx = C12->get_center(0);
								double C12cy = C12->get_center(1);
								double C13cx = C13->get_center(0);
								double C13cy = C13->get_center(1);
								double m21x = C11cx - C12cx;
								double m21y = C11cy - C12cy;
								double m21d = sqrt(m21x*m21x + m21y*m21y);
								double m31x = C11cx - C13cx;
								double m31y = C11cy - C13cy;
								double m31d = sqrt(m31x*m31x + m31y*m31y);
								C2_sense_migr_vec_x += m21x/m21d;
								C2_sense_migr_vec_y += m21y/m21d;
								C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C12->set_sense_migr(true);
								C3_sense_migr_vec_x += m31x/m31d;
								C3_sense_migr_vec_y += m31y/m31d;
								C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
								C13->set_sense_migr(true);
							}
							if (C12->Migrate())
							{
								double C11cx = C11->get_center(0);
								double C11cy = C11->get_center(1);
								double C12cx = C12->get_center(0);
								double C12cy = C12->get_center(1);
								double C13cx = C13->get_center(0);
								double C13cy = C13->get_center(1);
								double m12x = C12cx - C11cx;
								double m12y = C12cy - C11cy;
								double m12d = sqrt(m12x*m12x + m12y*m12y);
								double m32x = C12cx - C13cx;
								double m32y = C12cy - C13cy;
								double m32d = sqrt(m32x*m32x + m32y*m32y);
								C1_sense_migr_vec_x += m12x/m12d;
								C1_sense_migr_vec_y += m12y/m12d;
								C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C11->set_sense_migr(true);
								C3_sense_migr_vec_x += m32x/m32d;
								C3_sense_migr_vec_y += m32y/m32d;
								C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
								C13->set_sense_migr(true);
							}
							if (C13->Migrate())
							{
								double C11cx = C11->get_center(0);
								double C11cy = C11->get_center(1);
								double C12cx = C12->get_center(0);
								double C12cy = C12->get_center(1);
								double C13cx = C13->get_center(0);
								double C13cy = C13->get_center(1);
								double m13x = C13cx - C11cx;
								double m13y = C13cy - C11cy;
								double m13d = sqrt(m13x*m13x + m13y*m13y);
								double m23x = C13cx - C12cx;
								double m23y = C13cy - C12cy;
								double m23d = sqrt(m23x*m23x + m23y*m23y);
								C1_sense_migr_vec_x += m13x/m13d;
								C1_sense_migr_vec_y += m13y/m13d;
								C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C11->set_sense_migr(true);
								C2_sense_migr_vec_x += m23x/m23d;
								C2_sense_migr_vec_y += m23y/m23d;
								C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C12->set_sense_migr(true);
							}
						}
						else if (sense_migration_flag==1)
						{
							if (plot_V1_V1_ab_m>0 || plot_V1_V1_ba_m>0)
							{
								/// strategy 2:
								double C3_sense_migr_vec_x = C13->get_sense_migr_vec(0);
								double C3_sense_migr_vec_y = C13->get_sense_migr_vec(1);
								double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
								double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
								if (plot_V1_V1_ab_m>0 && plot_V3_V3_ba_m>0 && C11->Migrate())
								{
									double C11_ma = C11->get_migrate_angle();
									double C11_ma_x = cos(C11_ma*PI/180);
									double C11_ma_y = sin(C11_ma*PI/180);
									double C11_ma_1 = -V1_fx_m*C11_ma_x -V1_fy_m*C11_ma_y;
									double C11_ma_3 = V3_fx_m*C11_ma_x + V3_fy_m*C11_ma_y;
									double C11_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
									double C11_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
									if (C11_ma_1>0.71*C11_ma_1_d) // 45 degree
									{
										C3_sense_migr_vec_x -= V1_fx_m;
										C3_sense_migr_vec_y -= V1_fy_m;
										C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
										C13->set_sense_migr(true);
										if (C13->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V1: "<<-V1_fx_m<<","<<-V1_fy_m<<endl;}
									}
									else if (C11_ma_3>0.71*C11_ma_3_d) // 45 degree
									{
										C3_sense_migr_vec_x += V3_fx_m;
										C3_sense_migr_vec_y += V3_fy_m;
										C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
										C13->set_sense_migr(true);
										if (C13->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V1 but V3 itself: "<<V3_fx_m<<","<<V3_fy_m<<endl;}
									}
									if (C13->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
								if (plot_V1_V1_ba_m>0 && plot_V2_V2_ab_m>0 && C11->Migrate())
								{
									double C11_ma = C11->get_migrate_angle();
									double C11_ma_x = cos(C11_ma*PI/180);
									double C11_ma_y = sin(C11_ma*PI/180);
									double C11_ma_1 = -V1_fx_m*C11_ma_x -V1_fy_m*C11_ma_y;
									double C11_ma_2 = V2_fx_m*C11_ma_x + V2_fy_m*C11_ma_y;
									double C11_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
									double C11_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
									if (C11_ma_1>0.71*C11_ma_1_d) // 45 degree
									{
										C2_sense_migr_vec_x -= V1_fx_m;
										C2_sense_migr_vec_y -= V1_fy_m;
										C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C12->set_sense_migr(true);
										if (C12->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V1: "<<-V1_fx_m<<","<<-V1_fy_m<<endl;}
									}
									else if (C11_ma_2>0.71*C11_ma_2_d) // 45 degree
									{
										C2_sense_migr_vec_x += V2_fx_m;
										C2_sense_migr_vec_y += V2_fy_m;
										C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C12->set_sense_migr(true);
										if (C12->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V1 but V2 itself: "<<V2_fx_m<<","<<V2_fy_m<<endl;}
									}
									if (C12->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
							}
							if (plot_V2_V2_ab_m>0 || plot_V2_V2_ba_m>0)
							{
								/// strategy 2:
								double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
								double C3_sense_migr_vec_x = C13->get_sense_migr_vec(0);
								double C3_sense_migr_vec_y = C13->get_sense_migr_vec(1);
								if (plot_V2_V2_ab_m>0 && plot_V1_V1_ba_m>0 && C12->Migrate())
								{
									double C12_ma = C12->get_migrate_angle();
									double C12_ma_x = cos(C12_ma*PI/180);
									double C12_ma_y = sin(C12_ma*PI/180);
									double C12_ma_2 = -V2_fx_m*C12_ma_x -V2_fy_m*C12_ma_y;
									double C12_ma_1 = V1_fx_m*C12_ma_x + V1_fy_m*C12_ma_y;
									double C12_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
									double C12_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
									if (C12_ma_2>0.71*C12_ma_2_d) // 45 degree
									{
										C1_sense_migr_vec_x -= V2_fx_m;
										C1_sense_migr_vec_y -= V2_fy_m;
										C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
										C11->set_sense_migr(true);
										if (C11->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V2: "<<-V2_fx_m<<","<<-V2_fy_m<<endl;}
									}
									else if (C12_ma_1>0.71*C12_ma_1_d) // 45 degree
									{
										C1_sense_migr_vec_x += V1_fx_m;
										C1_sense_migr_vec_y += V1_fy_m;
										C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
										C11->set_sense_migr(true);
										if (C11->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V2 but V1 itself: "<<V1_fx_m<<","<<V1_fy_m<<endl;}
									}
									if (C11->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
								if (plot_V2_V2_ba_m>0 && plot_V3_V3_ab_m>0 && C12->Migrate())
								{
									double C12_ma = C12->get_migrate_angle();
									double C12_ma_x = cos(C12_ma*PI/180);
									double C12_ma_y = sin(C12_ma*PI/180);
									double C12_ma_2 = -V2_fx_m*C12_ma_x -V2_fy_m*C12_ma_y;
									double C12_ma_3 = V3_fx_m*C12_ma_x + V3_fy_m*C12_ma_y;
									double C12_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
									double C12_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
									if (C12_ma_2>0.71*C12_ma_2_d) // 45 degree
									{
										C3_sense_migr_vec_x -= V2_fx_m;
										C3_sense_migr_vec_y -= V2_fy_m;
										C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
										C13->set_sense_migr(true);
										if (C13->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V2: "<<-V2_fx_m<<","<<-V2_fy_m<<endl;}
									}
									else if (C12_ma_3>0.71*C12_ma_3_d) // 45 degree
									{
										C3_sense_migr_vec_x += V3_fx_m;
										C3_sense_migr_vec_y += V3_fy_m;
										C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
										C13->set_sense_migr(true);
										if (C13->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V2 but V3 itself: "<<V3_fx_m<<","<<V3_fy_m<<endl;}
									}								
									if (C13->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
							}
							if (plot_V3_V3_ab_m>0 || plot_V3_V3_ba_m>0)
							{
								/// strategy 2:
								double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
								double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
								double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
								if (plot_V3_V3_ab_m>0 && plot_V2_V2_ba_m>0 && C13->Migrate())
								{
									double C13_ma = C13->get_migrate_angle();
									double C13_ma_x = cos(C13_ma*PI/180);
									double C13_ma_y = sin(C13_ma*PI/180);
									double C13_ma_3 = -V3_fx_m*C13_ma_x -V3_fy_m*C13_ma_y;
									double C13_ma_2 = V2_fx_m*C13_ma_x + V2_fy_m*C13_ma_y;
									double C13_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
									double C13_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
									if (C13_ma_3>0.71*C13_ma_3_d) // 45 degree
									{
										C2_sense_migr_vec_x -= V3_fx_m;
										C2_sense_migr_vec_y -= V3_fy_m;
										C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C12->set_sense_migr(true);
										if (C12->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V3: "<<-V3_fx_m<<","<<-V3_fy_m<<endl;}
									}
									else if (C13_ma_2>0.71*C13_ma_2_d) // 45 degree
									{
										C2_sense_migr_vec_x += V2_fx_m;
										C2_sense_migr_vec_y += V2_fy_m;
										C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C12->set_sense_migr(true);
										if (C12->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V3 but V2 itself: "<<V2_fx_m<<","<<V2_fy_m<<endl;}
									}						
									if (C12->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
								if (plot_V3_V3_ba_m>0 && plot_V1_V1_ab_m>0 && C13->Migrate())
								{
									double C13_ma = C13->get_migrate_angle();
									double C13_ma_x = cos(C13_ma*PI/180);
									double C13_ma_y = sin(C13_ma*PI/180);
									double C13_ma_3 = -V3_fx_m*C13_ma_x -V3_fy_m*C13_ma_y;
									double C13_ma_1 = V1_fx_m*C13_ma_x + V1_fy_m*C13_ma_y;
									double C13_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
									double C13_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
									if (C13_ma_3>0.71*C13_ma_3_d) // 45 degree
									{
										C1_sense_migr_vec_x -= V3_fx_m;
										C1_sense_migr_vec_y -= V3_fy_m;
										C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
										C11->set_sense_migr(true);
										if (C11->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V3: "<<-V3_fx_m<<","<<-V3_fy_m<<endl;}
									}
									else if (C13_ma_1>0.71*C13_ma_1_d) // 45 degree
									{
										C1_sense_migr_vec_x += V1_fx_m;
										C1_sense_migr_vec_y += V1_fy_m;
										C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
										C11->set_sense_migr(true);
										if (C11->id()==w_id1) {cout<<"    -> cell "<<w_id1<<" from V3 but V1 itself: "<<V1_fx_m<<","<<V1_fy_m<<endl;}
									}
									if (C11->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
							}
						}
					}
					else if (V1->get_migrate_mark() && !V2->get_migrate_mark() && !V3->get_migrate_mark())
					{
						if (sense_migration_flag==1)
						{
							double plot_V1_V2_ab_m = V1_fx_m*V2_ab_normal_x + V1_fy_m*V2_ab_normal_y;
							double plot_V1_V2_ba_m = V1_fx_m*V2_ba_normal_x + V1_fy_m*V2_ba_normal_y;
							double plot_V1_V3_ab_m = V1_fx_m*V3_ab_normal_x + V1_fy_m*V3_ab_normal_y;
							double plot_V1_V3_ba_m = V1_fx_m*V3_ba_normal_x + V1_fy_m*V3_ba_normal_y;
							double C11_ma = C11->get_migrate_angle();
							double C11_ma_x = cos(C11_ma*PI/180);
							double C11_ma_y = sin(C11_ma*PI/180);
							double C11_ma_2 = V2_fx_m*C11_ma_x + V2_fy_m*C11_ma_y;
							double C11_ma_3 = V3_fx_m*C11_ma_x + V3_fy_m*C11_ma_y;
							if (C11_ma_2>0 && plot_V1_V2_ab_m>0 && plot_V1_V2_ba_m>0 && plot_V2_V2_ab_m>0 && plot_V2_V2_ba_m>0)
							{
								double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
								double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
								C2_sense_migr_vec_x += V2_fx_m;
								C2_sense_migr_vec_y += V2_fy_m;
								C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C12->set_sense_migr(true);
								if (C11->id()==w_id1 || C12->id()==w_id1 || C13->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
							if (C11_ma_3>0 && plot_V1_V3_ba_m>0 && plot_V1_V3_ab_m>0 && plot_V3_V3_ab_m>0 && plot_V3_V3_ba_m>0)
							{
								double C3_sense_migr_vec_x = C13->get_sense_migr_vec(0);
								double C3_sense_migr_vec_y = C13->get_sense_migr_vec(1);
								C3_sense_migr_vec_x += V3_fx_m;
								C3_sense_migr_vec_y += V3_fy_m;
								C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
								C13->set_sense_migr(true);
								if (C11->id()==w_id1 || C12->id()==w_id1 || C13->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
						}
					}
					else if (!V1->get_migrate_mark() && V2->get_migrate_mark() && !V3->get_migrate_mark())
					{
						if (sense_migration_flag==1)
						{
							double plot_V2_V1_ab_m = V2_fx_m*V1_ab_normal_x + V2_fy_m*V1_ab_normal_y;
							double plot_V2_V1_ba_m = V2_fx_m*V1_ba_normal_x + V2_fy_m*V1_ba_normal_y;
							double plot_V2_V3_ab_m = V2_fx_m*V3_ab_normal_x + V2_fy_m*V3_ab_normal_y;
							double plot_V2_V3_ba_m = V2_fx_m*V3_ba_normal_x + V2_fy_m*V3_ba_normal_y;
							double C12_ma = C12->get_migrate_angle();
							double C12_ma_x = cos(C12_ma*PI/180);
							double C12_ma_y = sin(C12_ma*PI/180);
							double C12_ma_1 = V1_fx_m*C12_ma_x + V1_fy_m*C12_ma_y;
							double C12_ma_3 = V3_fx_m*C12_ma_x + V3_fy_m*C12_ma_y;
							if (C12_ma_1>0 && plot_V2_V1_ab_m>0 && plot_V2_V1_ba_m>0 && plot_V1_V1_ab_m>0 && plot_V1_V1_ba_m>0)
							{
								double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
								C1_sense_migr_vec_x += V1_fx_m;
								C1_sense_migr_vec_y += V1_fy_m;
								C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C11->set_sense_migr(true);
								if (C11->id()==w_id1 || C12->id()==w_id1 || C13->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
							if (C12_ma_3>0 && plot_V2_V3_ab_m>0 && plot_V2_V3_ba_m>0 && plot_V3_V3_ab_m>0 && plot_V3_V3_ba_m>0)
							{
								double C3_sense_migr_vec_x = C13->get_sense_migr_vec(0);
								double C3_sense_migr_vec_y = C13->get_sense_migr_vec(1);
								C3_sense_migr_vec_x += V3_fx_m;
								C3_sense_migr_vec_y += V3_fy_m;
								C13->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
								C13->set_sense_migr(true);
								if (C11->id()==w_id1 || C12->id()==w_id1 || C13->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
						}
					}
					else if (!V1->get_migrate_mark() && !V2->get_migrate_mark() && V3->get_migrate_mark())
					{
						if (sense_migration_flag==1)
						{
							double plot_V3_V1_ab_m = V3_fx_m*V1_ab_normal_x + V3_fy_m*V1_ab_normal_y;
							double plot_V3_V1_ba_m = V3_fx_m*V1_ba_normal_x + V3_fy_m*V1_ba_normal_y;
							double plot_V3_V2_ab_m = V3_fx_m*V2_ab_normal_x + V3_fy_m*V2_ab_normal_y;
							double plot_V3_V2_ba_m = V3_fx_m*V2_ba_normal_x + V3_fy_m*V2_ba_normal_y;
							double C13_ma = C13->get_migrate_angle();
							double C13_ma_x = cos(C13_ma*PI/180);
							double C13_ma_y = sin(C13_ma*PI/180);
							double C13_ma_1 = V1_fx_m*C13_ma_x + V1_fy_m*C13_ma_y;
							double C13_ma_2 = V2_fx_m*C13_ma_x + V2_fy_m*C13_ma_y;
							if (C13_ma_1>0 && plot_V3_V1_ab_m>0 && plot_V3_V1_ba_m>0 && plot_V1_V1_ab_m>0 && plot_V1_V1_ba_m>0)
							{
								double C1_sense_migr_vec_x = C11->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C11->get_sense_migr_vec(1);
								C1_sense_migr_vec_x += V1_fx_m;
								C1_sense_migr_vec_y += V1_fy_m;
								C11->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C11->set_sense_migr(true);
								if (C11->id()==w_id1 || C12->id()==w_id1 || C13->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
							if (C13_ma_2>0 && plot_V3_V2_ab_m>0 && plot_V3_V2_ba_m>0 && plot_V2_V2_ab_m>0 && plot_V2_V2_ba_m>0)
							{
								double C2_sense_migr_vec_x = C12->get_sense_migr_vec(0);
								double C2_sense_migr_vec_y = C12->get_sense_migr_vec(1);
								C2_sense_migr_vec_x += V2_fx_m;
								C2_sense_migr_vec_y += V2_fy_m;
								C12->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C12->set_sense_migr(true);
								if (C11->id()==w_id1 || C12->id()==w_id1 || C13->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C11->id()<<","<<C12->id()<<","<<C13->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
						}
					}
					///
					if		(V1->get_pop_out()==0) 
					{
						V1_ab_bf += plot_V1_V1_ab;
						V1_ba_bf += plot_V1_V1_ba;
						/////
						double C1_sense_mech_vec_x = C11->get_sense_mech_vec(0);
						double C1_sense_mech_vec_y = C11->get_sense_mech_vec(1);
						C1_sense_mech_vec_x += V1_fx;
						C1_sense_mech_vec_y += V1_fy;
						C11->set_sense_mech_vec(C1_sense_mech_vec_x, C1_sense_mech_vec_y);
						C11->set_sense_mech(true);
					}
					else if (V1->get_pop_out()==1)
					{
						if (plot_V1_V1_ab>0) {V1_ab_bf -= plot_V1_V1_ab;}
						else                 {V1_ab_bf += plot_V1_V1_ab;}
						if (plot_V1_V1_ba>0) {V1_ba_bf -= plot_V1_V1_ba;}
						else                 {V1_ba_bf += plot_V1_V1_ba;}
						if (plot_e32_along_V1>0) {V3_ab_bf += plot_e32_along_V1;V2_ba_bf += plot_e32_along_V1;}
						else                     {V3_ab_bf -= plot_e32_along_V1;V2_ba_bf -= plot_e32_along_V1;}
					}
					else if (V1->get_pop_out()==-1)
					{
						if (plot_V1_V1_ab>0) {V1_ab_bf += plot_V1_V1_ab;}
						else                 {V1_ab_bf -= plot_V1_V1_ab;}
						if (plot_V1_V1_ba>0) {V1_ba_bf += plot_V1_V1_ba;}
						else                 {V1_ba_bf -= plot_V1_V1_ba;}
					}
					if		(V2->get_pop_out()==0)
					{
						V2_ab_bf += plot_V2_V2_ab;
						V2_ba_bf += plot_V2_V2_ba;
						/////
						double C2_sense_mech_vec_x = C12->get_sense_mech_vec(0);
						double C2_sense_mech_vec_y = C12->get_sense_mech_vec(1);
						C2_sense_mech_vec_x += V2_fx;
						C2_sense_mech_vec_y += V2_fy;
						C12->set_sense_mech_vec(C2_sense_mech_vec_x, C2_sense_mech_vec_y);\
						C12->set_sense_mech(true);
					}
					else if (V2->get_pop_out()==1)
					{
						if (plot_V2_V2_ab>0) {V2_ab_bf -= plot_V2_V2_ab;}
						else                 {V2_ab_bf += plot_V2_V2_ab;}
						if (plot_V2_V2_ba>0) {V2_ba_bf -= plot_V2_V2_ba;}
						else                 {V2_ba_bf += plot_V2_V2_ba;}
						if (plot_e13_along_V2>0) {V1_ab_bf += plot_e13_along_V2;V3_ba_bf += plot_e13_along_V2;}
						else                     {V1_ab_bf -= plot_e13_along_V2;V3_ba_bf -= plot_e13_along_V2;}
					}
					else if (V2->get_pop_out()==-1)
					{
						if (plot_V2_V2_ab>0) {V2_ab_bf += plot_V2_V2_ab;}
						else                 {V2_ab_bf -= plot_V2_V2_ab;}
						if (plot_V2_V2_ba>0) {V2_ba_bf += plot_V2_V2_ba;}
						else                 {V2_ba_bf -= plot_V2_V2_ba;}
					}
					if		(V3->get_pop_out()==0)
					{
						V3_ab_bf += plot_V3_V3_ab;
						V3_ba_bf += plot_V3_V3_ba;
						/////
						double C3_sense_mech_vec_x = C13->get_sense_mech_vec(0);
						double C3_sense_mech_vec_y = C13->get_sense_mech_vec(1);
						C3_sense_mech_vec_x += V3_fx;
						C3_sense_mech_vec_y += V3_fy;
						C13->set_sense_mech_vec(C3_sense_mech_vec_x, C3_sense_mech_vec_y);
						C13->set_sense_mech(true);
					}
					else if (V3->get_pop_out()==1)
					{
						if (plot_V3_V3_ab>0) {V3_ab_bf -= plot_V3_V3_ab;}
						else                 {V3_ab_bf += plot_V3_V3_ab;}
						if (plot_V3_V3_ba>0) {V3_ba_bf -= plot_V3_V3_ba;}
						else                 {V3_ba_bf += plot_V3_V3_ba;}
						if (plot_e21_along_V3>0) {V2_ab_bf += plot_e21_along_V3;V1_ba_bf += plot_e21_along_V3;}
						else                     {V2_ab_bf -= plot_e21_along_V3;V1_ba_bf -= plot_e21_along_V3;}
					}
					else if (V3->get_pop_out()==-1)
					{
						if (plot_V3_V3_ab>0) {V3_ab_bf += plot_V3_V3_ab;}
						else                 {V3_ab_bf -= plot_V3_V3_ab;}
						if (plot_V3_V3_ba>0) {V3_ba_bf += plot_V3_V3_ba;}
						else                 {V3_ba_bf -= plot_V3_V3_ba;}
					}
					if      (Pl12->MP()->at(0)->T1()==V1)
					{
						double Fb_t1 = Pl12->MP()->at(0)->get_F_t1();
						double Fb_s2 = Pl12->MP()->at(0)->get_F_s2();
						double Fb_t1_new = Fb_t1;
						double Fb_s2_new = Fb_s2;
						if (V1_ba_bf>0 && V2_ab_bf>0)
						{
							Fb_t1_new += V1_ba_bf;
							Fb_s2_new += V2_ab_bf;
						}
						else
						{
							double inc = V1_ba_bf + V2_ab_bf;
							if (inc>=0) {Fb_t1_new += inc;Fb_s2_new += inc;}
							else        {Fb_t1_new = 0;Fb_s2_new = 0;}
						}
						Pl12->MP()->at(0)->set_F_t1(Fb_t1_new);
						Pl12->MP()->at(0)->set_F_s2(Fb_s2_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ba_normal_x*Fb_t1_new;
							s_f_V1_y += V1_ba_normal_y*Fb_t1_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ba_normal_x*Fb_t1_new;
							s_f_V1_y += -V1_ba_normal_y*Fb_t1_new;
						}
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += V2_ab_normal_x*Fb_s2_new;
							s_f_V2_y += V2_ab_normal_y*Fb_s2_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += -V2_ab_normal_x*Fb_s2_new;
							s_f_V2_y += -V2_ab_normal_y*Fb_s2_new;
						}
						if (Fb_t1_new>break_force12 || Fb_s2_new>break_force12) {V21_b = true;}
					}
					else if (Pl12->MP()->at(0)->T2()==V1)
					{
						double Fb_t2 = Pl12->MP()->at(0)->get_F_t2();
						double Fb_s1 = Pl12->MP()->at(0)->get_F_s1();
						double Fb_t2_new = Fb_t2;
						double Fb_s1_new = Fb_s1;
						if (V1_ba_bf>0 && V2_ab_bf>0)
						{
							Fb_t2_new += V1_ba_bf;
							Fb_s1_new += V2_ab_bf;
						}
						else
						{
							double inc = V1_ba_bf + V2_ab_bf;
							if (inc>=0) {Fb_t2_new += inc;Fb_s1_new += inc;}
							else {Fb_t2_new = 0;Fb_s1_new = 0;}
						}
						Pl12->MP()->at(0)->set_F_t2(Fb_t2_new);
						Pl12->MP()->at(0)->set_F_s1(Fb_s1_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ba_normal_x*Fb_t2_new;
							s_f_V1_y += V1_ba_normal_y*Fb_t2_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ba_normal_x*Fb_t2_new;
							s_f_V1_y += -V1_ba_normal_y*Fb_t2_new;
						}
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += V2_ab_normal_x*Fb_s1_new;
							s_f_V2_y += V2_ab_normal_y*Fb_s1_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += -V2_ab_normal_x*Fb_s1_new;
							s_f_V2_y += -V2_ab_normal_y*Fb_s1_new;
						}
						if (Fb_t2_new>break_force12 || Fb_s1_new>break_force12) {V21_b = true;}
					}
					if      (Pl13->MP()->at(0)->S1()==V1)
					{
						double Fb_s1 = Pl13->MP()->at(0)->get_F_s1();
						double Fb_t2 = Pl13->MP()->at(0)->get_F_t2();
						double Fb_s1_new = Fb_s1;
						double Fb_t2_new = Fb_t2;
						if (V1_ab_bf>0 && V3_ba_bf>0)
						{
							Fb_s1_new += V1_ab_bf;
							Fb_t2_new += V3_ba_bf;
						}
						else
						{
							double inc = V1_ab_bf + V3_ba_bf;
							if (inc>=0) {Fb_s1_new += inc;Fb_t2_new += inc;}
							else {Fb_s1_new = 0;Fb_t2_new = 0;}
						}
						Pl13->MP()->at(0)->set_F_s1(Fb_s1_new);
						Pl13->MP()->at(0)->set_F_t2(Fb_t2_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ab_normal_x*Fb_s1_new;
							s_f_V1_y += V1_ab_normal_y*Fb_s1_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ab_normal_x*Fb_s1_new;
							s_f_V1_y += -V1_ab_normal_y*Fb_s1_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += V3_ba_normal_x*Fb_t2_new;
							s_f_V3_y += V3_ba_normal_y*Fb_t2_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += -V3_ba_normal_x*Fb_t2_new;
							s_f_V3_y += -V3_ba_normal_y*Fb_t2_new;
						}
						if (Fb_s1_new>break_force13 || Fb_t2_new>break_force13) {V13_b = true;}
					}
					else if (Pl13->MP()->at(0)->S2()==V1)
					{
						double Fb_s2 = Pl13->MP()->at(0)->get_F_s2();
						double Fb_t1 = Pl13->MP()->at(0)->get_F_t1();
						double Fb_s2_new = Fb_s2;
						double Fb_t1_new = Fb_t1;
						if (V1_ab_bf>0 && V3_ba_bf>0)
						{
							Fb_s2_new += V1_ab_bf;
							Fb_t1_new += V3_ba_bf;
						}
						else
						{
							double inc = V1_ab_bf + V3_ba_bf;
							if (inc>=0) {Fb_s2_new += inc;Fb_t1_new += inc;}
							else {Fb_s2_new = 0;Fb_t1_new = 0;}
						}
						Pl13->MP()->at(0)->set_F_s2(Fb_s2_new);
						Pl13->MP()->at(0)->set_F_t1(Fb_t1_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ab_normal_x*Fb_s2_new;
							s_f_V1_y += V1_ab_normal_y*Fb_s2_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ab_normal_x*Fb_s2_new;
							s_f_V1_y += -V1_ab_normal_y*Fb_s2_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += V3_ba_normal_x*Fb_t1_new;
							s_f_V3_y += V3_ba_normal_y*Fb_t1_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += -V3_ba_normal_x*Fb_t1_new;
							s_f_V3_y += -V3_ba_normal_y*Fb_t1_new;
						}
						if (Fb_s2_new>break_force13 || Fb_t1_new>break_force13) {V13_b = true;}
					}
					if      (Pl23->MP()->at(0)->T1()==V2)
					{
						double Fb_t1 = Pl23->MP()->at(0)->get_F_t1();
						double Fb_s2 = Pl23->MP()->at(0)->get_F_s2();
						double Fb_t1_new = Fb_t1;
						double Fb_s2_new = Fb_s2;
						if (V3_ab_bf>0 && V2_ba_bf>0)
						{
							Fb_t1_new += V2_ba_bf;
							Fb_s2_new += V3_ab_bf;
						}
						else
						{
							double inc = V3_ab_bf + V2_ba_bf;
							if (inc>=0) {Fb_t1_new += inc;Fb_s2_new += inc;}
							else {Fb_t1_new = 0;Fb_s2_new = 0;}
						}
						Pl23->MP()->at(0)->set_F_t1(Fb_t1_new);
						Pl23->MP()->at(0)->set_F_s2(Fb_s2_new);
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += V2_ba_normal_x*Fb_t1_new;
							s_f_V2_y += V2_ba_normal_y*Fb_t1_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += -V2_ba_normal_x*Fb_t1_new;
							s_f_V2_y += -V2_ba_normal_y*Fb_t1_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += V3_ab_normal_x*Fb_s2_new;
							s_f_V3_y += V3_ab_normal_y*Fb_s2_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += -V3_ab_normal_x*Fb_s2_new;
							s_f_V3_y += -V3_ab_normal_y*Fb_s2_new;
						}
						if (Fb_t1_new>break_force23 || Fb_s2_new>break_force23) {V32_b = true;}
					}
					else if (Pl23->MP()->at(0)->T2()==V2)
					{
						double Fb_t2 = Pl23->MP()->at(0)->get_F_t2();
						double Fb_s1 = Pl23->MP()->at(0)->get_F_s1();
						double Fb_t2_new = Fb_t2;
						double Fb_s1_new = Fb_s1;
						if (V3_ab_bf>0 && V2_ba_bf>0)
						{
							Fb_t2_new += V2_ba_bf;
							Fb_s1_new += V3_ab_bf;
						}
						else
						{
							double inc = V3_ab_bf + V2_ba_bf;
							if (inc>=0) {Fb_t2_new += inc;Fb_s1_new += inc;}
							else {Fb_t2_new = 0;Fb_s1_new = 0;}
						}
						Pl23->MP()->at(0)->set_F_t2(Fb_t2_new);
						Pl23->MP()->at(0)->set_F_s1(Fb_s1_new);
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += V2_ba_normal_x*Fb_t2_new;
							s_f_V2_y += V2_ba_normal_y*Fb_t2_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += -V2_ba_normal_x*Fb_t2_new;
							s_f_V2_y += -V2_ba_normal_y*Fb_t2_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += V3_ab_normal_x*Fb_s1_new;
							s_f_V3_y += V3_ab_normal_y*Fb_s1_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += -V3_ab_normal_x*Fb_s1_new;
							s_f_V3_y += -V3_ab_normal_y*Fb_s1_new;
						}
						if (Fb_t2_new>break_force23 || Fb_s1_new>break_force23) {V32_b = true;}
					}
					V1->set_stored_force_vis(s_f_V1_x,s_f_V1_y);
					V2->set_stored_force_vis(s_f_V2_x,s_f_V2_y);
					V3->set_stored_force_vis(s_f_V3_x,s_f_V3_y);
					if      (V21_b && !V32_b && !V13_b) 
					{
						if      (V3->get_pop_out()>0 && !V3->getab()->Cell()->Soften()) 
						{edgebreak1.push_back(V2->getab());}
						else if (V3->get_pop_out()<=0 && V1->getab()->Cell()->Migrate())
						{pointbreak2.push_back(V1);}
						else if (V3->get_pop_out()<=0 && V2->getab()->Cell()->Migrate())
						{
							pointbreak2.push_back(V2);
						}
					}
					else if (!V21_b && !V32_b && V13_b) 
					{
						if (V2->get_pop_out()>0 && !V2->getab()->Cell()->Soften())
						{edgebreak1.push_back(V1->getab());}
						else if (V2->get_pop_out()<=0 && V1->getab()->Cell()->Migrate())
						{pointbreak2.push_back(V1);}
						else if (V2->get_pop_out()<=0 && V3->getab()->Cell()->Migrate())
						{pointbreak2.push_back(V3);}
					}
					else if (!V21_b && V32_b && !V13_b) 
					{
						if (V1->get_pop_out()>0 && !V1->getab()->Cell()->Soften())
						{edgebreak1.push_back(V3->getab());}
						else if (V1->get_pop_out()<=0 && V2->getab()->Cell()->Migrate())
						{pointbreak2.push_back(V2);}
						else if (V1->get_pop_out()<=0 && V3->getab()->Cell()->Migrate())
						{pointbreak2.push_back(V3);}
					}
					else if (V21_b && !V32_b && V13_b)
					{
						if (V1->get_pop_out()<=0) {pointbreak2.push_back(V1);}
					}
					else if (V21_b && V32_b && !V13_b)
					{
						if (V2->get_pop_out()<=0) {pointbreak2.push_back(V2);}
					}
					else if (!V21_b && V32_b && V13_b)
					{
						if (V3->get_pop_out()<=0) {pointbreak2.push_back(V3);}
					}
				}
				else if ((nv1a && !nv2a && !nv3a) ||
					     (!nv1a && nv2a && !nv3a) ||
						 (!nv1a && !nv2a && nv3a))
				{
					/**********************************
					      o-<-o
						       \     o-<-o
							 V3 \   /     
							 o->-o /       
							 o-<-oo V2        
							   V1 \\        
								   \\        
								o->-oo->-o
					**********************************/
					point3D *V1 = nv1;
					if      (nv2a) {V1 = nv2;}
					else if (nv3a) {V1 = nv3;}
					point3D *V2 = V1->getba()->get_Neighbor()->p1();
					point3D *V3 = V1->getab()->get_Neighbor()->p2();
					//cout<<"3 points joint: "<<V1->id()<<","<<V2->id()<<","<<V3->id()<<" of cell "<<V1->getab()->Cell()->id()<<","<<V2->getab()->Cell()->id()<<","<<V3->getab()->Cell()->id()<<endl;
					cell *C1t = V1->getab()->Cell();
					cell *C2t = V2->getab()->Cell();
					cell *C3t = V3->getab()->Cell();
					interpair *Pl21 = V1->getba()->get_Pair();
					interpair *Pl13 = V1->getab()->get_Pair();
					interpair *Pl32 = Pl21;
					for (int j=0;j<(int)C2t->get_pairs()->size();j++)
					{
						interpair *Plt = C2t->get_pairs()->at(j);
						if ((Plt->get_I1()==C2t && Plt->get_I2()==C3t) ||
							(Plt->get_I1()==C2t && Plt->get_I2()==C3t)) {Pl32 = Plt;break;}
					}
					double break_force21 = Pl21->get_adhesion();
					double break_force13 = Pl13->get_adhesion();
					double break_force32 = Pl32->get_adhesion();
					double V1_fx = V1->get_elastic_force(0);
					double V1_fy = V1->get_elastic_force(1);
					double V2_fx = V2->get_elastic_force(0);
					double V2_fy = V2->get_elastic_force(1);
					double V3_fx = V3->get_elastic_force(0);
					double V3_fy = V3->get_elastic_force(1);
					///
					double V1_fx_m = V1->get_elastic_force_migration(0);
					double V1_fy_m = V1->get_elastic_force_migration(1);
					double V2_fx_m = V2->get_elastic_force_migration(0);
					double V2_fy_m = V2->get_elastic_force_migration(1);
					double V3_fx_m = V3->get_elastic_force_migration(0);
					double V3_fy_m = V3->get_elastic_force_migration(1);
					///
					double V1_ab_nx = V1->getab()->getnormal(0);
					double V1_ab_ny = V1->getab()->getnormal(1);
					double V1_ba_nx = V1->getba()->getnormal(0);
					double V1_ba_ny = V1->getba()->getnormal(1);
					double V2_ba_nx = V2->getba()->getnormal(0);
					double V2_ba_ny = V2->getba()->getnormal(1);
					double V3_ab_nx = V3->getab()->getnormal(0);
					double V3_ab_ny = V3->getab()->getnormal(1);
					double V1_ab_nd = sqrt(V1_ab_nx*V1_ab_nx + V1_ab_ny*V1_ab_ny);
					double V1_ba_nd = sqrt(V1_ba_nx*V1_ba_nx + V1_ba_ny*V1_ba_ny);
					double V3_ab_nd = sqrt(V3_ab_nx*V3_ab_nx + V3_ab_ny*V3_ab_ny);
					double V2_ba_nd = sqrt(V2_ba_nx*V2_ba_nx + V2_ba_ny*V2_ba_ny);
					double V1_ab_normal_x = V1_ab_nx/V1_ab_nd;
					double V1_ab_normal_y = V1_ab_ny/V1_ab_nd;
					double V1_ba_normal_x = V1_ba_nx/V1_ba_nd;
					double V1_ba_normal_y = V1_ba_ny/V1_ba_nd;
					double V3_ab_normal_x = V3_ab_nx/V3_ab_nd;
					double V3_ab_normal_y = V3_ab_ny/V3_ab_nd;
					double V2_ba_normal_x = V2_ba_nx/V2_ba_nd;
					double V2_ba_normal_y = V2_ba_ny/V2_ba_nd;
					double e32_along_x = (V2->getrp()->x() - V2->x()) + (V3->getfp()->x() - V3->x());
					double e32_along_y = (V2->getrp()->y() - V2->y()) + (V3->getfp()->y() - V3->y());
					double e32_along_d = sqrt(e32_along_x*e32_along_x + e32_along_y*e32_along_y);
					double e32_along_normal_x = e32_along_x/e32_along_d;
					double e32_along_normal_y = e32_along_y/e32_along_d;
					bool V21_b = false;
					bool V13_b = false;
					bool V32_b = false;
					double plot_V1_V1_ab = V1_fx*V1_ab_normal_x + V1_fy*V1_ab_normal_y;
					double plot_V3_V3_ba = -V3_fx*V1_ab_normal_x -V3_fy*V1_ab_normal_y;
					double plot_V1_V1_ba = V1_fx*V1_ba_normal_x + V1_fy*V1_ba_normal_y;
					double plot_V2_V2_ab = -V2_fx*V1_ba_normal_x -V2_fy*V1_ba_normal_y;
					double plot_V3_V3_ab = V3_fx*V3_ab_normal_x + V3_fy*V3_ab_normal_y;
					double plot_V2_V2_ba = V2_fx*V2_ba_normal_x + V2_fy*V2_ba_normal_y;
					///
					double plot_V1_V1_ab_m = V1_fx_m*V1_ab_normal_x + V1_fy_m*V1_ab_normal_y;
					double plot_V3_V3_ba_m = -V3_fx_m*V1_ab_normal_x -V3_fy_m*V1_ab_normal_y;
					double plot_V1_V1_ba_m = V1_fx_m*V1_ba_normal_x + V1_fy_m*V1_ba_normal_y;
					double plot_V2_V2_ab_m = -V2_fx_m*V1_ba_normal_x -V2_fy_m*V1_ba_normal_y;
					///
					double plot_e32_along_V1 = V1_fx*e32_along_normal_x + V1_fy*e32_along_normal_y;
					double V1_ab_bf = 0;
					double V1_ba_bf = 0;
					double V2_ab_bf = 0;
					double V2_ba_bf = 0;
					double V3_ab_bf = 0;
					double V3_ba_bf = 0;
					double s_f_V1_x = 0;
					double s_f_V1_y = 0;
					double s_f_V2_x = 0;
					double s_f_V2_y = 0;
					double s_f_V3_x = 0;
					double s_f_V3_y = 0;
					/// in response to migration
					if (!V1->get_migrate_mark() && !V2->get_migrate_mark() && !V3->get_migrate_mark())
					{
						if (sense_migration_flag==0)
						{
							double C3_sense_migr_vec_x = C3t->get_sense_migr_vec(0);
							double C3_sense_migr_vec_y = C3t->get_sense_migr_vec(1);
							double C2_sense_migr_vec_x = C2t->get_sense_migr_vec(0);
							double C2_sense_migr_vec_y = C2t->get_sense_migr_vec(1);
							double C1_sense_migr_vec_x = C1t->get_sense_migr_vec(0);
							double C1_sense_migr_vec_y = C1t->get_sense_migr_vec(1);
							if (C1t->Migrate())
							{
								double C11cx = C1t->get_center(0);
								double C11cy = C1t->get_center(1);
								double C12cx = C2t->get_center(0);
								double C12cy = C2t->get_center(1);
								double C13cx = C3t->get_center(0);
								double C13cy = C3t->get_center(1);
								double m21x = C11cx - C12cx;
								double m21y = C11cy - C12cy;
								double m21d = sqrt(m21x*m21x + m21y*m21y);
								double m31x = C11cx - C13cx;
								double m31y = C11cy - C13cy;
								double m31d = sqrt(m31x*m31x + m31y*m31y);
								C2_sense_migr_vec_x += m21x/m21d;
								C2_sense_migr_vec_y += m21y/m21d;
								C2t->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C2t->set_sense_migr(true);
								C3_sense_migr_vec_x += m31x/m31d;
								C3_sense_migr_vec_y += m31y/m31d;
								C3t->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
								C3t->set_sense_migr(true);
							}
							if (C2t->Migrate())
							{
								double C11cx = C1t->get_center(0);
								double C11cy = C1t->get_center(1);
								double C12cx = C2t->get_center(0);
								double C12cy = C2t->get_center(1);
								double C13cx = C3t->get_center(0);
								double C13cy = C3t->get_center(1);
								double m12x = C12cx - C11cx;
								double m12y = C12cy - C11cy;
								double m12d = sqrt(m12x*m12x + m12y*m12y);
								double m32x = C12cx - C13cx;
								double m32y = C12cy - C13cy;
								double m32d = sqrt(m32x*m32x + m32y*m32y);
								C1_sense_migr_vec_x += m12x/m12d;
								C1_sense_migr_vec_y += m12y/m12d;
								C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C1t->set_sense_migr(true);
								C3_sense_migr_vec_x += m32x/m32d;
								C3_sense_migr_vec_y += m32y/m32d;
								C3t->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
								C3t->set_sense_migr(true);
							}
							if (C3t->Migrate())
							{
								double C11cx = C1t->get_center(0);
								double C11cy = C1t->get_center(1);
								double C12cx = C2t->get_center(0);
								double C12cy = C2t->get_center(1);
								double C13cx = C3t->get_center(0);
								double C13cy = C3t->get_center(1);
								double m13x = C13cx - C11cx;
								double m13y = C13cy - C11cy;
								double m13d = sqrt(m13x*m13x + m13y*m13y);
								double m23x = C13cx - C12cx;
								double m23y = C13cy - C12cy;
								double m23d = sqrt(m23x*m23x + m23y*m23y);
								C1_sense_migr_vec_x += m13x/m13d;
								C1_sense_migr_vec_y += m13y/m13d;
								C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
								C1t->set_sense_migr(true);
								C2_sense_migr_vec_x += m23x/m23d;
								C2_sense_migr_vec_y += m23y/m23d;
								C2t->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
								C2t->set_sense_migr(true);
							}
						}
						else if (sense_migration_flag==1)
						{
							if (plot_V1_V1_ab_m>0 || plot_V1_V1_ba_m>0)
							{
								/// strategy 2:
								double C2_sense_migr_vec_x = C2t->get_sense_migr_vec(0);
								double C2_sense_migr_vec_y = C2t->get_sense_migr_vec(1);
								double C3_sense_migr_vec_x = C3t->get_sense_migr_vec(0);
								double C3_sense_migr_vec_y = C3t->get_sense_migr_vec(1);
								if (plot_V1_V1_ab_m>0 && plot_V3_V3_ba_m>0 && C1t->Migrate())
								{
									double C1t_ma = C1t->get_migrate_angle();
									double C1t_ma_x = cos(C1t_ma*PI/180);
									double C1t_ma_y = sin(C1t_ma*PI/180);
									double C1t_ma_1 = -V1_fx_m*C1t_ma_x -V1_fy_m*C1t_ma_y;
									double C1t_ma_3 = V3_fx_m*C1t_ma_x + V3_fy_m*C1t_ma_y;
									double C1t_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
									double C1t_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
									if (C1t_ma_1>0.71*C1t_ma_1_d) // 45 degree
									{
										C3_sense_migr_vec_x -= V1_fx_m;
										C3_sense_migr_vec_y -= V1_fy_m;
										C3t->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
										C3t->set_sense_migr(true);
									}
									else if (C1t_ma_3>0.71*C1t_ma_3_d) //  45 degree
									{
										C3_sense_migr_vec_x += V3_fx_m;
										C3_sense_migr_vec_y += V3_fy_m;
										C3t->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
										C3t->set_sense_migr(true);
									}
									if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
								if (plot_V1_V1_ba_m>0 && plot_V2_V2_ab_m>0 && C1t->Migrate())
								{
									double C1t_ma = C1t->get_migrate_angle();
									double C1t_ma_x = cos(C1t_ma*PI/180);
									double C1t_ma_y = sin(C1t_ma*PI/180);
									double C1t_ma_1 = -V1_fx_m*C1t_ma_x -V1_fy_m*C1t_ma_y;
									double C1t_ma_2 = V2_fx_m*C1t_ma_x + V2_fy_m*C1t_ma_y;
									double C1t_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
									double C1t_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
									if (C1t_ma_1>0.71*C1t_ma_1_d) //  45 degree
									{
										C2_sense_migr_vec_x -= V1_fx_m;
										C2_sense_migr_vec_y -= V1_fy_m;
										C2t->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C2t->set_sense_migr(true);
									}
									else if (C1t_ma_2>0.71*C1t_ma_2_d) //  45 degree
									{
										C2_sense_migr_vec_x += V2_fx_m;
										C2_sense_migr_vec_y += V2_fy_m;
										C2t->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
										C2t->set_sense_migr(true);
									}
									if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
									{
										cout<<"    -> "<<w_id1<<" sense mech: cells ";
										cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
										cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
										cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
										cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
									}
								}
							}
							if (plot_V2_V2_ab_m>0 && plot_V1_V1_ba_m>0 && C2t->Migrate())
							{
								/// strategy 2:
								double C1_sense_migr_vec_x = C1t->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C1t->get_sense_migr_vec(1);
								double C2t_ma = C2t->get_migrate_angle();
								double C2t_ma_x = cos(C2t_ma*PI/180);
								double C2t_ma_y = sin(C2t_ma*PI/180);
								double C2t_ma_2 = -V2_fx_m*C2t_ma_x -V2_fy_m*C2t_ma_y;
								double C2t_ma_1 = V1_fx_m*C2t_ma_x + V1_fy_m*C2t_ma_y;
								double C2t_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
								double C2t_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
								if (C2t_ma_2>0.71*C2t_ma_2_d) // 45 degree
								{
									C1_sense_migr_vec_x -= V2_fx_m;
									C1_sense_migr_vec_y -= V2_fy_m;
									C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
									C1t->set_sense_migr(true);
								}
								else if (C2t_ma_1>0.71*C2t_ma_1_d) // 45 degree
								{
									C1_sense_migr_vec_x += V1_fx_m;
									C1_sense_migr_vec_y += V1_fy_m;
									C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
									C1t->set_sense_migr(true);
								}
								if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
							if (plot_V3_V3_ba_m>0 && plot_V1_V1_ab_m>0 && C3t->Migrate())
							{
								/// strategy 2:
								double C1_sense_migr_vec_x = C1t->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C1t->get_sense_migr_vec(1);
								double C3t_ma = C3t->get_migrate_angle();
								double C3t_ma_x = cos(C3t_ma*PI/180);
								double C3t_ma_y = sin(C3t_ma*PI/180);
								double C3t_ma_3 = -V3_fx_m*C3t_ma_x -V3_fy_m*C3t_ma_y;
								double C3t_ma_1 = V1_fx_m*C3t_ma_x + V1_fy_m*C3t_ma_y;
								double C3t_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
								double C3t_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
								if (C3t_ma_3>0.71*C3t_ma_3_d) // 45 degree
								{
									C1_sense_migr_vec_x -= V3_fx_m;
									C1_sense_migr_vec_y -= V3_fy_m;
									C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
									C1t->set_sense_migr(true);
								}
								else if (C3t_ma_1>0.71*C3t_ma_1_d) // 45 degree
								{
									C1_sense_migr_vec_x += V1_fx_m;
									C1_sense_migr_vec_y += V1_fy_m;
									C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
									C1t->set_sense_migr(true);
								}
								if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
						}
					}
					else if (!V1->get_migrate_mark() && V2->get_migrate_mark() && !V3->get_migrate_mark())
					{
						if (sense_migration_flag==1)
						{
							if (plot_V1_V1_ab_m>0 && plot_V1_V1_ba_m>0)
							{
								double C1_sense_migr_vec_x = C1t->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C1t->get_sense_migr_vec(1);
								double C2t_ma = C2t->get_migrate_angle();
								double C2t_ma_x = cos(C2t_ma*PI/180);
								double C2t_ma_y = sin(C2t_ma*PI/180);
								double C2t_ma_1 = V1_fx_m*C2t_ma_x + V1_fy_m*C2t_ma_y;
								double C2t_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
								if (C2t_ma_1>0.71*C2t_ma_1_d) // 45 degree
								{
									C1_sense_migr_vec_x += V1_fx_m;
									C1_sense_migr_vec_y += V1_fy_m;
									C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
									C1t->set_sense_migr(true);
								}
								if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
							if (plot_V3_V3_ba_m>0)
							{
								double C3_sense_migr_vec_x = C3t->get_sense_migr_vec(0);
								double C3_sense_migr_vec_y = C3t->get_sense_migr_vec(1);
								double C2t_ma = C2t->get_migrate_angle();
								double C2t_ma_x = cos(C2t_ma*PI/180);
								double C2t_ma_y = sin(C2t_ma*PI/180);
								double C2t_ma_3 = V3_fx_m*C2t_ma_x + V3_fy_m*C2t_ma_y;
								double C2t_ma_3_d = sqrt(V3_fx_m*V3_fx_m + V3_fy_m*V3_fy_m);
								if (C2t_ma_3>0.71*C2t_ma_3_d) // 45 degree
								{
									C3_sense_migr_vec_x += V3_fx_m;
									C3_sense_migr_vec_y += V3_fy_m;
									C3t->set_sense_migr_vec(C3_sense_migr_vec_x, C3_sense_migr_vec_y);
									C3t->set_sense_migr(true);
								}
								if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
						}
					}
					else if (!V1->get_migrate_mark() && !V2->get_migrate_mark() && V3->get_migrate_mark())
					{
						if (sense_migration_flag==1)
						{
							if (plot_V1_V1_ab_m>0 && plot_V1_V1_ba_m>0)
							{
								double C1_sense_migr_vec_x = C1t->get_sense_migr_vec(0);
								double C1_sense_migr_vec_y = C1t->get_sense_migr_vec(1);
								double C3t_ma = C3t->get_migrate_angle();
								double C3t_ma_x = cos(C3t_ma*PI/180);
								double C3t_ma_y = sin(C3t_ma*PI/180);
								double C3t_ma_1 = V1_fx_m*C3t_ma_x + V1_fy_m*C3t_ma_y;
								double C3t_ma_1_d = sqrt(V1_fx_m*V1_fx_m + V1_fy_m*V1_fy_m);
								if (C3t_ma_1>0)
								{
									C1_sense_migr_vec_x += V1_fx_m;
									C1_sense_migr_vec_y += V1_fy_m;
									C1t->set_sense_migr_vec(C1_sense_migr_vec_x, C1_sense_migr_vec_y);
									C1t->set_sense_migr(true);
								}
								if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
							if (plot_V2_V2_ab_m>0)
							{
								double C2_sense_migr_vec_x = C2t->get_sense_migr_vec(0);
								double C2_sense_migr_vec_y = C2t->get_sense_migr_vec(1);
								double C3t_ma = C3t->get_migrate_angle();
								double C3t_ma_x = cos(C3t_ma*PI/180);
								double C3t_ma_y = sin(C3t_ma*PI/180);
								double C3t_ma_2 = V2_fx_m*C3t_ma_x + V2_fy_m*C3t_ma_y;
								double C3t_ma_2_d = sqrt(V2_fx_m*V2_fx_m + V2_fy_m*V2_fy_m);
								if (C3t_ma_2>0.71*C3t_ma_2_d) //  45 degree
								{
									C2_sense_migr_vec_x += V2_fx_m;
									C2_sense_migr_vec_y += V2_fy_m;
									C2t->set_sense_migr_vec(C2_sense_migr_vec_x, C2_sense_migr_vec_y);
									C2t->set_sense_migr(true);
								}
								if (C1t->id()==w_id1 || C2t->id()==w_id1 || C3t->id()==w_id1)
								{
									cout<<"    -> "<<w_id1<<" sense mech: cells ";
									cout<<" "<<C1t->id()<<","<<C2t->id()<<","<<C3t->id()<<": "<<V1->get_migrate_mark()<<","<<V2->get_migrate_mark()<<","<<V3->get_migrate_mark()<<endl;
									cout<<"    -> "<<V1_fx_m<<","<<V1_fy_m<<endl;
									cout<<"    -> "<<V2_fx_m<<","<<V2_fy_m<<endl;
									cout<<"    -> "<<V3_fx_m<<","<<V3_fy_m<<endl;
								}
							}
						}
					}
					////////// V1 //////////
					if      (V1->get_pop_out()==0)
					{
						V1_ab_bf += plot_V1_V1_ab;
						V1_ba_bf += plot_V1_V1_ba;
					}
					else if (V1->get_pop_out()==1)
					{
						if (plot_V1_V1_ab>0) {V1_ab_bf -= plot_V1_V1_ab;}
						else                 {V1_ab_bf += plot_V1_V1_ab;}
						if (plot_V1_V1_ba>0) {V1_ba_bf -= plot_V1_V1_ba;}
						else                 {V1_ba_bf += plot_V1_V1_ba;}
						if (plot_e32_along_V1>0) {V3_ab_bf += plot_e32_along_V1;V2_ba_bf += plot_e32_along_V1;}
						else                     {V3_ab_bf -= plot_e32_along_V1;V2_ba_bf -= plot_e32_along_V1;}
					}
					else if (V1->get_pop_out()==-1)
					{
						if (plot_V1_V1_ab>0) {V1_ab_bf += plot_V1_V1_ab;}
						else                 {V1_ab_bf -= plot_V1_V1_ab;}
						if (plot_V1_V1_ba>0) {V1_ba_bf += plot_V1_V1_ba;}
						else                 {V1_ba_bf -= plot_V1_V1_ba;}
					}
					///////////// V2 /////////////
					if      (V2->get_pop_out()==0)
					{
						V2_ab_bf += plot_V2_V2_ab;
						V2_ba_bf += plot_V2_V2_ba;
					}
					else if (V2->get_pop_out()==1)
					{
						if (plot_V2_V2_ab>0) {V2_ab_bf -= plot_V2_V2_ab;}
						else                 {V2_ab_bf += plot_V2_V2_ab;}
						if (plot_V2_V2_ba>0) {V2_ba_bf -= plot_V2_V2_ba;}
						else                 {V2_ba_bf += plot_V2_V2_ba;}
					}
					else if (V2->get_pop_out()==-1)
					{
						if (plot_V2_V2_ab>0) {V2_ab_bf += plot_V2_V2_ab;}
						else                 {V2_ab_bf -= plot_V2_V2_ab;}
						if (plot_V2_V2_ba>0) {V2_ba_bf += plot_V2_V2_ba;}
						else                 {V2_ba_bf -= plot_V2_V2_ba;}
					}
					////////////// V3 //////////////
					if      (V3->get_pop_out()==0)
					{
						V3_ba_bf += plot_V3_V3_ba;
						V3_ab_bf += plot_V3_V3_ab;
					}
					else if (V3->get_pop_out()==1)
					{
						if (plot_V3_V3_ab>0) {V3_ab_bf -= plot_V3_V3_ab;}
						else                 {V3_ab_bf += plot_V3_V3_ab;}
						if (plot_V3_V3_ba>0) {V3_ba_bf -= plot_V3_V3_ba;}
						else                 {V3_ba_bf += plot_V3_V3_ba;}
					}
					else if (V3->get_pop_out()==-1)
					{
						if (plot_V3_V3_ab>0) {V3_ab_bf += plot_V3_V3_ab;}
						else                 {V3_ab_bf -= plot_V3_V3_ab;}
						if (plot_V3_V3_ba>0) {V3_ba_bf += plot_V3_V3_ba;}
						else                 {V3_ba_bf -= plot_V3_V3_ba;}
					}
					if (Pl21->MP()->size()==0)
					{
						cout<<"  -> error: 3 points joint, V1-V2 joint apart: "<<V1->id()<<","<<V2->id()<<","<<V3->id()<<" of cell "<<V1->getab()->Cell()->id()<<","<<V2->getab()->Cell()->id()<<","<<V3->getab()->Cell()->id()<<endl;
						return V1->getab()->Cell()->id();
					}
					if      (Pl21->MP()->at(0)->T1()==V1)
					{
						double Fb_t1 = Pl21->MP()->at(0)->get_F_t1();
						double Fb_s2 = Pl21->MP()->at(0)->get_F_s2();
						double Fb_t1_new = Fb_t1;
						double Fb_s2_new = Fb_s2;
						if (V1_ba_bf>0 && V2_ab_bf>0)
						{
							Fb_t1_new += V1_ba_bf;
							Fb_s2_new += V2_ab_bf;
						}
						else 
						{
							double inc = V1_ba_bf + V2_ab_bf;
							if (inc>=0) {Fb_t1_new += inc;Fb_s2_new += inc;}
							else {Fb_t1_new = 0;Fb_s2_new = 0;}
						}
						Pl21->MP()->at(0)->set_F_t1(Fb_t1_new);
						Pl21->MP()->at(0)->set_F_s2(Fb_s2_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ba_normal_x*Fb_t1_new;
							s_f_V1_y += V1_ba_normal_y*Fb_t1_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ba_normal_x*Fb_t1_new;
							s_f_V1_y += -V1_ba_normal_y*Fb_t1_new;
						}
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += -V1_ba_normal_x*Fb_s2_new;
							s_f_V2_y += -V1_ba_normal_y*Fb_s2_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += V1_ba_normal_x*Fb_s2_new;
							s_f_V2_y += V1_ba_normal_y*Fb_s2_new;
						}
						if (Fb_t1_new>break_force21 || Fb_s2_new>break_force21)
						{
							V21_b = true;
						}
					}
					else if (Pl21->MP()->at(0)->T2()==V1)
					{
						double Fb_t2 = Pl21->MP()->at(0)->get_F_t2();
						double Fb_s1 = Pl21->MP()->at(0)->get_F_s1();
						double Fb_t2_new = Fb_t2;
						double Fb_s1_new = Fb_s1;
						if (V1_ba_bf>0 && V2_ab_bf>0)
						{
							Fb_t2_new += V1_ba_bf;
							Fb_s1_new += V2_ab_bf;
						}
						else
						{
							double inc = V1_ba_bf + V2_ab_bf;
							if (inc>=0) {Fb_t2_new += inc;Fb_s1_new += inc;}
							else {Fb_t2_new = 0;Fb_s1_new = 0;}
						}
						Pl21->MP()->at(0)->set_F_t2(Fb_t2_new);
						Pl21->MP()->at(0)->set_F_s1(Fb_s1_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ba_normal_x*Fb_t2_new;
							s_f_V1_y += V1_ba_normal_y*Fb_t2_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ba_normal_x*Fb_t2_new;
							s_f_V1_y += -V1_ba_normal_y*Fb_t2_new;
						}
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += -V1_ba_normal_x*Fb_s1_new;
							s_f_V2_y += -V1_ba_normal_y*Fb_s1_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += V1_ba_normal_x*Fb_s1_new;
							s_f_V2_y += V1_ba_normal_y*Fb_s1_new;
						}
						if (Fb_t2_new>break_force21 || Fb_s1_new>break_force21)
						{
							V21_b = true;
						}
					}
					if (Pl13->MP()->size()==0)
					{
						cout<<"  -> error: 3 points joint, V1-V3 joint apart: "<<V1->id()<<","<<V2->id()<<","<<V3->id()<<" of cell "<<V1->getab()->Cell()->id()<<","<<V2->getab()->Cell()->id()<<","<<V3->getab()->Cell()->id()<<endl;
						return V1->getab()->Cell()->id();
					}
					if      (Pl13->MP()->at(0)->S1()==V1)
					{
						double Fb_s1 = Pl13->MP()->at(0)->get_F_s1();
						double Fb_t2 = Pl13->MP()->at(0)->get_F_t2();
						double Fb_s1_new = Fb_s1;
						double Fb_t2_new = Fb_t2;
						if (V1_ab_bf>0 && V3_ba_bf>0)
						{
							Fb_s1_new += V1_ab_bf;
							Fb_t2_new += V3_ba_bf;
						}
						else
						{
							double inc = V1_ab_bf + V3_ba_bf;
							if (inc>=0) {Fb_s1_new += inc;Fb_t2_new += inc;}
							else {Fb_s1_new = 0;Fb_t2_new = 0;}
						}
						Pl13->MP()->at(0)->set_F_s1(Fb_s1_new);
						Pl13->MP()->at(0)->set_F_t2(Fb_t2_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ab_normal_x*Fb_s1_new;
							s_f_V1_y += V1_ab_normal_y*Fb_s1_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ab_normal_x*Fb_s1_new;
							s_f_V1_y += -V1_ab_normal_y*Fb_s1_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += -V1_ab_normal_x*Fb_t2_new;
							s_f_V3_y += -V1_ab_normal_y*Fb_t2_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += V1_ab_normal_x*Fb_t2_new;
							s_f_V3_y += V1_ab_normal_y*Fb_t2_new;
						}
						if (Fb_s1_new>break_force13 || Fb_t2_new>break_force13) {V13_b = true;}
					}
					else if (Pl13->MP()->at(0)->S2()==V1)
					{
						double Fb_s2 = Pl13->MP()->at(0)->get_F_s2();
						double Fb_t1 = Pl13->MP()->at(0)->get_F_t1();
						double Fb_s2_new = Fb_s2;
						double Fb_t1_new = Fb_t1;
						if (V1_ab_bf>0 && V3_ba_bf>0)
						{
							Fb_s2_new += V1_ab_bf;
							Fb_t1_new += V3_ba_bf;
						}
						else
						{
							double inc = V1_ab_bf + V3_ba_bf;
							if (inc>=0) {Fb_s2_new += inc;Fb_t1_new += inc;}
							else {Fb_s2_new = 0;Fb_t1_new = 0;}
						}
						Pl13->MP()->at(0)->set_F_s2(Fb_s2_new);
						Pl13->MP()->at(0)->set_F_t1(Fb_t1_new);
						if      (V1->get_pop_out()==0)
						{
							s_f_V1_x += V1_ab_normal_x*Fb_s2_new;
							s_f_V1_y += V1_ab_normal_y*Fb_s2_new;
						}
						else if (V1->get_pop_out()==-1)
						{
							s_f_V1_x += -V1_ab_normal_x*Fb_s2_new;
							s_f_V1_y += -V1_ab_normal_y*Fb_s2_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += -V1_ab_normal_x*Fb_t1_new;
							s_f_V3_y += -V1_ab_normal_y*Fb_t1_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += V1_ab_normal_x*Fb_t1_new;
							s_f_V3_y += V1_ab_normal_y*Fb_t1_new;
						}
						if (Fb_s2_new>break_force13 || Fb_t1_new>break_force13) {V13_b = true;}
					}
					if (Pl32->MP()->size()==0)
					{
						cout<<"  -> error: 3 points joint, V2-V3 joint apart: "<<V1->id()<<","<<V2->id()<<","<<V3->id()<<" of cell "<<V1->getab()->Cell()->id()<<","<<V2->getab()->Cell()->id()<<","<<V3->getab()->Cell()->id()<<endl;
						return V1->getab()->Cell()->id();
					}
					if      (Pl32->MP()->at(0)->S1()==V2)
					{
						double Fb_t1 = Pl32->MP()->at(0)->get_F_t1();
						double Fb_s2 = Pl32->MP()->at(0)->get_F_s2();
						double Fb_s1 = Pl32->MP()->at(0)->get_F_s1();
						double Fb_t2 = Pl32->MP()->at(0)->get_F_t2();
						double Fb_t1_new = Fb_t1;
						double Fb_s2_new = Fb_s2;
						double Fb_s1_new = Fb_s1;
						double Fb_t2_new = Fb_t2;
						if (V3_ab_bf>0 && V2_ba_bf>0)
						{
							Fb_t1_new += V2_ba_bf;
							Fb_s2_new += V3_ab_bf;
							Fb_s1_new += V2_ba_bf;
							Fb_t2_new += V3_ab_bf;
						}
						else 
						{
							double inc = V3_ab_bf + V2_ba_bf;
							if (inc>=0)
							{
								Fb_t1_new += inc;Fb_s2_new += inc;
								Fb_s1_new += inc;Fb_t2_new += inc;
							}
							else 
							{
								Fb_t1_new = 0;Fb_s2_new = 0;
								Fb_s1_new = 0;Fb_t2_new = 0;
							}
						}
						Pl32->MP()->at(0)->set_F_t1(Fb_t1_new);
						Pl32->MP()->at(0)->set_F_s1(Fb_s1_new);
						Pl32->MP()->at(0)->set_F_s2(Fb_s2_new);
						Pl32->MP()->at(0)->set_F_t2(Fb_t2_new);
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += V2_ba_normal_x*Fb_s1_new;
							s_f_V2_y += V2_ba_normal_y*Fb_s1_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += -V2_ba_normal_x*Fb_s1_new;
							s_f_V2_y += -V2_ba_normal_y*Fb_s1_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += V3_ab_normal_x*Fb_t2_new;
							s_f_V3_y += V3_ab_normal_y*Fb_t2_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += -V3_ab_normal_x*Fb_t2_new;
							s_f_V3_y += -V3_ab_normal_y*Fb_t2_new;
						}
						if (Fb_t1_new>break_force32 ||
							Fb_s2_new>break_force32 ||
							Fb_s1_new>break_force32 ||
							Fb_t2_new>break_force32) 
						{V32_b = true;}
					}
					else if (Pl32->MP()->at(0)->S2()==V2)
					{
						double Fb_t2 = Pl32->MP()->at(0)->get_F_t2();
						double Fb_s1 = Pl32->MP()->at(0)->get_F_s1();
						double Fb_s2 = Pl32->MP()->at(0)->get_F_s2();
						double Fb_t1 = Pl32->MP()->at(0)->get_F_t1();
						double Fb_t2_new = Fb_t2;
						double Fb_s1_new = Fb_s1;
						double Fb_s2_new = Fb_s2;
						double Fb_t1_new = Fb_t1;
						if (V3_ab_bf>0 && V2_ba_bf>0)
						{
							Fb_t2_new += V2_ba_bf;
							Fb_s1_new += V3_ab_bf;
							Fb_s2_new += V2_ba_bf;
							Fb_t1_new += V3_ab_bf;
						}
						else
						{
							double inc = V3_ab_bf + V2_ba_bf;
							if (inc>=0)
							{
								Fb_t2_new += inc;Fb_s1_new += inc;
								Fb_s2_new += inc;Fb_t1_new += inc;
							}
							else
							{
								Fb_t2_new = 0;Fb_s1_new = 0;
								Fb_s2_new = 0;Fb_t1_new = 0;
							}
						}
						Pl32->MP()->at(0)->set_F_t2(Fb_t2_new);
						Pl32->MP()->at(0)->set_F_s2(Fb_s2_new);
						Pl32->MP()->at(0)->set_F_s1(Fb_s1_new);
						Pl32->MP()->at(0)->set_F_t1(Fb_t1_new);
						if      (V2->get_pop_out()==0)
						{
							s_f_V2_x += V2_ba_normal_x*Fb_s2_new;
							s_f_V2_y += V2_ba_normal_y*Fb_s2_new;
						}
						else if (V2->get_pop_out()==-1)
						{
							s_f_V2_x += -V2_ba_normal_x*Fb_s2_new;
							s_f_V2_y += -V2_ba_normal_y*Fb_s2_new;
						}
						if      (V3->get_pop_out()==0)
						{
							s_f_V3_x += V3_ab_normal_x*Fb_t1_new;
							s_f_V3_y += V3_ab_normal_y*Fb_t1_new;
						}
						else if (V3->get_pop_out()==-1)
						{
							s_f_V3_x += -V3_ab_normal_x*Fb_t1_new;
							s_f_V3_y += -V3_ab_normal_y*Fb_t1_new;
						}
						if (Fb_t2_new>break_force32 ||
							Fb_s1_new>break_force32 ||
							Fb_s2_new>break_force32 ||
							Fb_t1_new>break_force32)
						{V32_b = true;}
					}
					V1->set_stored_force_vis(s_f_V1_x,s_f_V1_y);
					V2->set_stored_force_vis(s_f_V2_x,s_f_V2_y);
					V3->set_stored_force_vis(s_f_V3_x,s_f_V3_y);
					if      (V21_b && !V13_b && !V32_b)
					{
						edgebreak0.push_back(V2->getab());
					} // edge_break0
					else if (!V21_b && V13_b && !V32_b)
					{
						edgebreak0.push_back(V1->getab());
					} // edge_break0
					else if (V21_b && V13_b && !V32_b)  {pointbreak1.push_back(V1);} // point_break1
					else if (V32_b)                     {pointbreak1.push_back(V1);} // point_break1
				}
			}
		}
	}
	t_e_re = clock();
	cout<<"Elastic force setup done: "<<(t_e_re-t_b_re)*0.001<<" seconds!"<<endl;
	/////// free up memory ////////
	time_t t_b_free1, t_e_free1, t_b_free2, t_e_free2;
	t_b_free1 = clock();
	delete[] SOR_II_V;
	delete[] SOR_N;
	delete[] JA;
	delete[] VA;
	delete[] X;
	delete[] Xn;
	delete[] Fn;
	delete[] F_relax;
	delete[] X_relax;
	delete[] X_relax_n;
	t_e_free1 = clock();
	t_b_free2 = clock();
	int smn = (int)sms.size();
	for (int i=0;i<smn;i++)
	{
		delete sms[i];
		sms[i] = NULL;
	}
	sms.clear();
	nodepoints.clear();
	t_e_free2 = clock();
	cout<<"Free up memory done: "<<(t_e_free2-t_b_free1)*0.001<<" seconds!"<<endl;
	time_t t_b_ad, t_e_ad;
	t_b_ad = clock();

	for (int i=0;i<(int)edgebreak0.size();i++) 
	{
		//cout<<"  -> ";
		edge_response_b0(edgebreak0[i]);
	}
	for (int i=0;i<(int)edgebreak1.size();i++) 
	{
		//cout<<"  -> ";
		edge_response_b1(edgebreak1[i]);
	}
	for (int i=0;i<(int)pointbreak0.size();i++) 
	{
		//cout<<"  -> ";
		point_response_b0(pointbreak0[i]);
	}
	for (int i=0;i<(int)pointbreak1.size();i++) 
	{
		//cout<<"  -> ";
		point_response_b1(pointbreak1[i]);
	}
	for (int i=0;i<(int)pointbreak2.size();i++) 
	{
		//cout<<"  -> ";
		point_response_b2(pointbreak2[i]);
	}
	t_e_ad = clock();
	cout<<"  -> Adhesion force part done in "<<(t_e_ad-t_b_ad)*0.001<<" seconds!"<<endl;
	///////// friction force update for migrating cells /////////
	time_t t_b_fr, t_e_fr;
	t_b_fr = clock();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (cellList[i]->Migrate())
		{
			//cout<<"    -> friction response for cell "<<i<<endl;
			cell_friction_response(cellList[i],number);
			cellList[i]->set_sense_migr(false); // all migrating cells don't respond with adjacent migrating signals
			cellList[i]->set_sense_migr_vec(0, 0);
		}
	}
	t_e_fr = clock();
	cout<<"  -> Friction force part done in "<<(t_e_fr-t_b_fr)*0.001<<" seconds!"<<endl;
	//////////// relax the migrating and soften cell ////////////
	time_t t_b_msr, t_e_msr;
	t_b_msr = clock();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (cellList[i]->get_migrate_relax())
		{
			//cout<<"    -> simple relax for cell "<<i<<endl;
			cell_simple_relaxation(cellList[i],number);
		}
		else if (cellList[i]->Soften())
		{
			cell_simple_relaxation(cellList[i],number);
		}
	}
	t_e_msr = clock();
	cout<<"  -> Migrating and soften cell relaxation in "<<(t_e_msr-t_b_msr)*0.001<<" seconds!"<<endl;
	///////// self collision check for soften cells /////////////
	time_t t_b_csc, t_e_csc;
	t_b_csc = clock();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (cellList[i]->Soften())
		{
			bool scsc = cell_self_collision(cellList[i],number);
			if (scsc) {cell_death_response(cellList[i]);}
		}
	}
	t_e_csc = clock();
	cout<<"  -> Soften cell self-collision detection done in "<<(t_e_csc-t_b_csc)*0.001<<" seconds!"<<endl;
	/////// new cell born at rear side of migrating cell ///////
	time_t t_b_gr, t_e_gr;
	t_b_gr = clock();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		if (cellList[i]->Migrate())
		{
			int flag_mcr = migrating_cell_rear(cellList[i],0,number); // hole filled up by ECM: 0
			if (flag_mcr>=0) {return flag_mcr;}
		}
	}
	/////// update migrating status ///////
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue;
		//////// clear the old Growth rate ////////
		if (cellList[i]->get_behavior()!=2) {cellList[i]->set_GR(0);}
		////////
		if (cellList[i]->Migrate())
		{
			cellList[i]->set_migrate(0);
			cellList[i]->set_migrate_relax(0);
			cellList[i]->set_migrate_p(NULL,NULL);
		}
	}
	t_e_gr = clock();
	cout<<"  -> New cells born at rear side of migrating cells done in "<<(t_e_gr-t_b_gr)*0.001<<" seconds!"<<endl;
	t_e_re = clock();
	cout<<"Internal force recovery done: "<<(t_e_re-t_b_re)*0.001<<" seconds!"<<endl;
	edgebreak0.clear();
	edgebreak1.clear();
	pointbreak0.clear();
	pointbreak1.clear();
	pointbreak2.clear();

	return -1;
}

#endif