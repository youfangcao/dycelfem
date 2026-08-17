// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     point3D.cpp         ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef POINT3D_CPP
#define POINT3D_CPP

#include <cmath>
#include "point3D.h"
#include "biology.h"

using namespace std;

point3D::point3D(int id, double x, double y)
{
	index = id;
	u = x;
	v = y;
	u_3[0] = x;u_3[1] = x;u_3[2] = x;u_3[3] = x;
	v_3[0] = y;v_3[1] = y;v_3[2] = y;v_3[3] = y;
	merge = false;
	in_cell = true;
	node_check = false;
	cell_boundary = true;
	burry = false;
	end_point = false;
	migrate_mark = false;
	node_sm_n = -10000;
	pop_out = 1; // default growing status: grow
	adhesion_break[0] = 0;adhesion_break[1] = 0;
	elastic_force[0] = 0;elastic_force[1] = 0;
	migrate_response_force[0] = 0;migrate_response_force[1] = 0;
	single_Force[0] = 0;single_Force[1] = 0;
	P_in[0] = 0;P_in[1] = 0;
	A_ij[0] = 0;A_ij[1] = 0;A_ij[2] = 0;A_ij[3] = 0;
	A_ij_single[0] = 0;A_ij_single[1] = 0;A_ij_single[2] = 0;A_ij_single[3] = 0;
	Dp[0] = 0;Dp[1] = 0;
	stored_apart_force = 0;
	stored_force_vis[0] = 0;stored_force_vis[1] = 0;
	stored_friction_force = 0;
	migration_rate = 0;
	single_growth_rate = 0;
}

point3D::~point3D()
{
	C = NULL;
	fp = NULL;
	rp = NULL;
	int mpn = (int)mp.size();
	if (mpn>0)
	{
		for (int i=0;i<mpn;i++)
		{
			mp.at(i) = NULL;
		}
		mp.clear();
	}
	int iatn = (int)inner_at.size();
	if (iatn>0)
	{
		for (int i=0;i<iatn;i++)
		{
			inner_at.at(i) = NULL;
		}
		inner_at.clear();
	}
	int smn = (int)node_sm.size();
	if (smn>0)
	{
		for (int i=0;i<smn;i++)
		{
			node_sm.at(i)= NULL;
		}
		node_sm.clear();
	}
}

void point3D::setPressure_increment()
{
	/******************************************************************** 
	                       setup incremental area 
	  =================================================================
                                
	           o----o----o----o  
              /                \  
             o                  o
			 |                  |
			 o                  o
			  \  ba V ab       / 
			   o->--o->--o----o  
                 |__| => incremental area due to pressure
	********************************************************************/
	double area = GAR/ab->Cell()->get_sides()->size();
	double ab_l = ab->getlength();
	double ba_l = ba->getlength();
	double ab_m = ab_l/BD*ab->Cell()->get_mass();
	double ba_m = ba_l/BD*ba->Cell()->get_mass();
	double mass = (ab_m + ba_m)*0.5;
	
	double ab_x = (fp->x() - u)*0.5;
	double ab_y = (fp->y() - v)*0.5;
	double ba_x = (u - rp->x())*0.5;
	double ba_y = (v - rp->y())*0.5;

	double normal_ab_x = ab->getnormal(0);
	double normal_ab_y = ab->getnormal(1);
	double normal_ba_x = ba->getnormal(0);
	double normal_ba_y = ba->getnormal(1);
	double d_ab = sqrt(normal_ab_x*normal_ab_x + normal_ab_y*normal_ab_y);
	double d_ba = sqrt(normal_ba_x*normal_ba_x + normal_ba_y*normal_ba_y);

	double T_ab[2],T_ba[2]; // Tension
	T_ab[0] = (ab->p2()->x() - ab->p1()->x())*ab->Cell()->get_tcoef();
	T_ab[1] = (ab->p2()->y() - ab->p1()->y())*ab->Cell()->get_tcoef();
	T_ba[0] = (ba->p1()->x() - ba->p2()->x())*ba->Cell()->get_tcoef();
	T_ba[1] = (ba->p1()->y() - ba->p2()->y())*ba->Cell()->get_tcoef();

	double P_ab[2],P_ba[2]; // Pressure
	P_ab[0] = 0.5*ab->Cell()->get_pcoef()*ab_l*normal_ab_x/d_ab;
	P_ab[1] = 0.5*ab->Cell()->get_pcoef()*ab_l*normal_ab_y/d_ab;
	P_ba[0] = 0.5*ba->Cell()->get_pcoef()*ba_l*normal_ba_x/d_ba;
	P_ba[1] = 0.5*ba->Cell()->get_pcoef()*ba_l*normal_ba_y/d_ba;

	double P_increment_ab[2];	// Pressure increment
	double P_increment_ba[2];	// Pressure increment
	double F_ab_volume = area*4*ab_m/((ab_l+ba_l)*TM*TM);
	double F_ba_volume = area*4*ba_m/((ba_l+ab_l)*TM*TM);
	P_increment_ab[0] = normal_ab_x/d_ab*F_ab_volume;
	P_increment_ab[1] = normal_ab_y/d_ab*F_ab_volume;
	P_increment_ba[0] = normal_ba_x/d_ba*F_ba_volume;
	P_increment_ba[1] = normal_ba_y/d_ba*F_ba_volume;

	double F[2]; // Composition Force
	F[0] = T_ab[0] + T_ba[0] + P_ab[0] + P_ba[0];
	F[1] = T_ab[1] + T_ba[1] + P_ab[1] + P_ba[1];

	setForce(F[0],F[1]);
	setP_in_ab(P_increment_ab[0],P_increment_ab[1]);
	setP_in_ba(P_increment_ba[0],P_increment_ba[1]);
}

void point3D::pushmp(point3D* P)
{
	if (P->id()!=index)
	{
		bool exist = false;
		for (int i=0;i<(int)mp.size();i++)
		{
			if (P==mp.at(i)) {exist = true;break;}
		}
		if (!exist) {mp.push_back(P);}
	}
}

bool point3D::checkmp(point3D *P)
{
	bool valid = false;
	if (P->id()!=index)
	{
		for (int i=0;i<(int)mp.size();i++)
		{
			if (P==mp.at(i)) {valid = true;break;}
		}
	}
	return valid;
}
 
void point3D::removemp(point3D* P)
{
	int flag = -1;
	int mn = (int)mp.size();
	for (int i=0;i<(int)mp.size();i++)
	{
		if (P==mp.at(i))
		{
			flag = i;
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==mn-1)
		{
			mp.pop_back();
		}
		else
		{
			for (int i=flag;i<mn-1;i++)
			{
				mp.at(i) = mp.at(i+1);
			}
			mp.pop_back();
		}
	}
}

void point3D::clearmp()
{
	while(mp.size())
	{
		mp.pop_back();
	}
}

void point3D::clear_inner_at()
{
	while (inner_at.size())
	{
		inner_at.pop_back();
	}
}

double point3D::point_angle()
{
	double angle = 0;
	double x1 = rp->x() - u;
	double y1 = rp->y() - v;
	double x2 = fp->x() - u;
	double y2 = fp->y() - v;
	double D1 = sqrt(x1*x1 + y1*y1);
	double D2 = sqrt(x2*x2 + y2*y2);
	double angle1 = 0,angle2 = 0;
	
	/**************************
	           /| y
				|     
				|
      __________|___________
	            |\         / x
				| \
				|  \
				|   o
	**************************/
	if      (x1>=0 && y1>=0)
	{angle1 = asin(y1/D1)*180.0/PI;}
	else if (x1>=0 && y1<0)
	{angle1 = asin(y1/D1)*180.0/PI + 360.0;}
	else if (x1<0 && y1>=0)
	{angle1 = 180.0 - asin(y1/D1)*180.0/PI;}
	else if (x1<0 && y1<0)
	{angle1 = 180.0 - asin(y1/D1)*180.0/PI;}

	if      (x2>=0 && y2>=0)
	{angle2 = asin(y2/D2)*180.0/PI;}
	else if (x2>=0 && y2<0)
	{angle2 = asin(y2/D2)*180.0/PI + 360.0;}
	else if (x2<0 && y2>=0)
	{angle2 = 180.0 - asin(y2/D2)*180.0/PI;}
	else if (x2<0 && y2<0)
	{angle2 = 180.0 - asin(y2/D2)*180.0/PI;}

	angle = angle1 - angle2;
	if (angle<0) {angle = 360.0 + angle;}
	
	return angle;
}

void point3D::clear_node_sm()
{
	while(node_sm.size())
	{
		int nmn = (int)node_sm.size();
		node_sm.at(nmn-1) = NULL;
		node_sm.pop_back();
	}
}

void point3D::clear_node_sm_single()
{
	while (node_sm_single.size())
	{
		int nmn = (int)node_sm_single.size();
		node_sm_single.at(nmn-1) = NULL;
		node_sm_single.pop_back();
	}
}

edge::edge(int id, point3D *A, point3D *B, cell *Cell)
{
	index = id;
	P1 = A;
	P2 = B;
	C = Cell;
	attach = false;
	in_cell = true;
	virtual_stiffness_matrix[0] = 0;
	virtual_stiffness_matrix[1] = 0;
	virtual_stiffness_matrix[2] = 0;
	virtual_stiffness_matrix[3] = 0;
}

edge::~edge()
{
	P1 = NULL;
	P2 = NULL;
	C = NULL;
	Neighbor = NULL;
}

double edge::getnormal(int i)
{
	double v = 0;
	if (i==0)
	{
		v = P2->y() - P1->y();
	}
	else if (i==1)
	{
		v = P1->x() - P2->x();
	}
	return v;
}

double edge::getlength()
{
	double l = 0;
	double x1 = P1->x();
	double y1 = P1->y();
	double x2 = P2->x();
	double y2 = P2->y();
	l = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	return l;
}

bool edge::edge_twist()
{
	/************************

     P3 o---->------o P1
	   /           /
      /           /
     /           /
	o		    o---->---o P4
	           P2         \
						   \
							o 
	************************/
	bool valid = false;
	double x1 = P1->x();
	double y1 = P1->y();
	double x2 = P2->x();
	double y2 = P2->y();
	double x3 = P1->getrp()->x();
	double y3 = P1->getrp()->y();
	double x4 = P2->getfp()->x();
	double y4 = P2->getfp()->y();
	double cross1 = (x1-x3)*(x2-x1) + (y1-y3)*(y2-y1);
	double cross2 = (x2-x1)*(x4-x2) + (y2-y1)*(y4-y2);
	double angle2 = P2->point_angle();
	if ((cross1<0 && cross2<0) || angle2>=300) {valid = true;}
	return valid;
}

bool edge::judge_vert_out_side(double x, double y)
{
	/****************************************************
                   N  P
				   | /
				   |/
	  P2 o----<<---o P1
	****************************************************/
	bool valid = false;
	double norm_x = P2->y() - P1->y();
	double norm_y = P1->x() - P2->x();
	double v_x = x - P1->x();
	double v_y = y - P1->y();
	double prod = v_x*norm_x + v_y*norm_y;
	if (abs(prod)<0.00000000001) {prod=0;}
	if (prod>0) {valid = true;}
	return valid;
}

void edge::set_virtual_stiffness_matrix()
{
	/************************************
                  A
	              o    
		         / \    normal
		        /   \  |\ 
		       /     \ |
		      /       \|
		 P2  o---<<----o  P1  
	************************************/
	point3D *p1 = P1->get_node_id()->getvertex()->at(0);
	point3D *p2 = P2->get_node_id()->getvertex()->at(0);
	double x1 = p1->x();
	double y1 = p1->y();
	double x2 = p2->x();
	double y2 = p2->y();
	setnormal();
	double d = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
	double len = getlength();
	double xA = (x1 + x2)/2 + normal[0]/d*0.866*len;
	double yA = (y1 + y2)/2 + normal[1]/d*0.866*len;
	////////// triangle 12A //////////
	double b1t = y2 - yA;
	double c1t = xA - x2;
	double b2t = yA - y1;
	double c2t = x1 - xA;
	double b3t = y1 - y2;
	double c3t = x2 - x1;
	double A2t = c1t*y1 + c2t*y2 + c3t*yA;if (A2t<0) {A2t = -A2t;}
	double b1_b1 = b1t*b1t/A2t;double b1_c1 = b1t*c1t/A2t;double b1_b2 = b1t*b2t/A2t;double b1_c2 = b1t*c2t/A2t;double b1_b3 = b1t*b3t/A2t;double b1_c3 = b1t*c3t/A2t;
						       double c1_c1 = c1t*c1t/A2t;double c1_b2 = c1t*b2t/A2t;double c1_c2 = c1t*c2t/A2t;double c1_b3 = c1t*b3t/A2t;double c1_c3 = c1t*c3t/A2t;
											              double b2_b2 = b2t*b2t/A2t;double b2_c2 = b2t*c2t/A2t;double b2_b3 = b2t*b3t/A2t;double b2_c3 = b2t*c3t/A2t;
													    					         double c2_c2 = c2t*c2t/A2t;double c2_b3 = c2t*b3t/A2t;double c2_c3 = c2t*c3t/A2t;
																												double b3_b3 = b3t*b3t/A2t;double b3_c3 = b3t*c3t/A2t;
																																	       double c3_c3 = c3t*c3t/A2t;
    /////////////
	double K00 = (MU_virtual2_LAMBDA_virtual*b1_b1 + MU_virtual*c1_c1)/2;
	double K01 = (MU_virtual1_LAMBDA_virtual*b1_c1)/2;
	double K02 = (MU_virtual2_LAMBDA_virtual*b1_b2 + MU_virtual*c1_c2)/2;
	double K03 = (LAMBDA_virtual*b1_c2             + MU_virtual*c1_b2)/2;
	double K04 = (MU_virtual2_LAMBDA_virtual*b1_b3 + MU_virtual*c1_c3)/2;
	double K05 = (LAMBDA_virtual*b1_c3             + MU_virtual*c1_b3)/2;
	/////////////
	double K11 = (MU_virtual2_LAMBDA_virtual*c1_c1 + MU_virtual*b1_b1)/2;
	double K12 = (LAMBDA_virtual*c1_b2			   + MU_virtual*b1_c2)/2;
	double K13 = (MU_virtual2_LAMBDA_virtual*c1_c2 + MU_virtual*b1_b2)/2;
	double K14 = (LAMBDA_virtual*c1_b3             + MU_virtual*b1_c3)/2;
	double K15 = (MU_virtual2_LAMBDA_virtual*c1_c3 + MU_virtual*b1_b3)/2;
	/////////////
	double K22 = (MU_virtual2_LAMBDA_virtual*b2_b2 + MU_virtual*c2_c2)/2;
	double K23 = (MU_virtual1_LAMBDA_virtual*b2_c2)/2;
	double K24 = (MU_virtual2_LAMBDA_virtual*b2_b3 + MU_virtual*c2_c3)/2;
	double K25 = (LAMBDA_virtual*b2_c3             + MU_virtual*c2_b3)/2;
	/////////////
	double K33 = (MU_virtual2_LAMBDA_virtual*c2_c2 + MU_virtual*b2_b2)/2;
	double K34 = (LAMBDA_virtual*c2_b3             + MU_virtual*b2_c3)/2;
	double K35 = (MU_virtual2_LAMBDA_virtual*c2_c3 + MU_virtual*b2_b3)/2;
	/////////////
	double K44 = (MU_virtual2_LAMBDA_virtual*b3_b3 + MU_virtual*c3_c3)/2;
	double K45 = (MU_virtual1_LAMBDA_virtual*b3_c3)/2;
	/////////////
	double K55 = (MU_virtual2_LAMBDA_virtual*c3_c3 + MU_virtual*b3_b3)/2;
	//////////////////////////////////////////////////////////////////////
	virtual_triangle_stiffness_matrix[0][0] = K00;
	virtual_triangle_stiffness_matrix[0][1] = K01;
	virtual_triangle_stiffness_matrix[0][2] = K02;
	virtual_triangle_stiffness_matrix[0][3] = K03;
	virtual_triangle_stiffness_matrix[0][4] = K04;
	virtual_triangle_stiffness_matrix[0][5] = K05;
	///////
	virtual_triangle_stiffness_matrix[1][0] = K01;
	virtual_triangle_stiffness_matrix[1][1] = K11;
	virtual_triangle_stiffness_matrix[1][2] = K12;
	virtual_triangle_stiffness_matrix[1][3] = K13;
	virtual_triangle_stiffness_matrix[1][4] = K14;
	virtual_triangle_stiffness_matrix[1][5] = K15;
	///////
	virtual_triangle_stiffness_matrix[2][0] = K02;
	virtual_triangle_stiffness_matrix[2][1] = K12;
	virtual_triangle_stiffness_matrix[2][2] = K22;
	virtual_triangle_stiffness_matrix[2][3] = K23;
	virtual_triangle_stiffness_matrix[2][4] = K24;
	virtual_triangle_stiffness_matrix[2][5] = K25;
	///////
	virtual_triangle_stiffness_matrix[3][0] = K03;
	virtual_triangle_stiffness_matrix[3][1] = K13;
	virtual_triangle_stiffness_matrix[3][2] = K23;
	virtual_triangle_stiffness_matrix[3][3] = K33;
	virtual_triangle_stiffness_matrix[3][4] = K34;
	virtual_triangle_stiffness_matrix[3][5] = K35;
	///////
	virtual_triangle_stiffness_matrix[4][0] = K04;
	virtual_triangle_stiffness_matrix[4][1] = K14;
	virtual_triangle_stiffness_matrix[4][2] = K24;
	virtual_triangle_stiffness_matrix[4][3] = K34;
	virtual_triangle_stiffness_matrix[4][4] = K44;
	virtual_triangle_stiffness_matrix[4][5] = K45;
	///////
	virtual_triangle_stiffness_matrix[5][0] = K05;
	virtual_triangle_stiffness_matrix[5][1] = K15;
	virtual_triangle_stiffness_matrix[5][2] = K25;
	virtual_triangle_stiffness_matrix[5][3] = K35;
	virtual_triangle_stiffness_matrix[5][4] = K45;
	virtual_triangle_stiffness_matrix[5][5] = K55;
	///////// add the virtual element matrix into the global matrix ////////
	for (int i=0;i<(int)p1->get_node_sm()->size();i++)
	{
		if (p1->get_node_sm()->at(i)->getN1()==p1 &&
			p1->get_node_sm()->at(i)->getN2()==p2)
		{
			double a11 = K02 + p1->get_node_sm()->at(i)->get_stiffness_matrix_12(0);
			double a12 = K03 + p1->get_node_sm()->at(i)->get_stiffness_matrix_12(1);
			double a21 = K12 + p1->get_node_sm()->at(i)->get_stiffness_matrix_12(2);
			double a22 = K13 + p1->get_node_sm()->at(i)->get_stiffness_matrix_12(3);
			p1->get_node_sm()->at(i)->set_sm_12(a11,a12,a21,a22);
			p1->get_node_sm()->at(i)->set_sm_21(a11,a21,a12,a22);
			virtual_stiffness_matrix[0] = a11;
			virtual_stiffness_matrix[1] = a12;
			virtual_stiffness_matrix[2] = a21;
			virtual_stiffness_matrix[3] = a22;
			break;
		}
		else if (p1->get_node_sm()->at(i)->getN1()==p2 &&
				 p1->get_node_sm()->at(i)->getN2()==p1)
		{
			double a11 = K02 + p1->get_node_sm()->at(i)->get_stiffness_matrix_21(0);
			double a12 = K03 + p1->get_node_sm()->at(i)->get_stiffness_matrix_21(1);
			double a21 = K12 + p1->get_node_sm()->at(i)->get_stiffness_matrix_21(2);
			double a22 = K13 + p1->get_node_sm()->at(i)->get_stiffness_matrix_21(3);
			p1->get_node_sm()->at(i)->set_sm_21(a11,a12,a21,a22);
			p1->get_node_sm()->at(i)->set_sm_12(a11,a21,a12,a22);
			virtual_stiffness_matrix[0] = a11;
			virtual_stiffness_matrix[1] = a12;
			virtual_stiffness_matrix[2] = a21;
			virtual_stiffness_matrix[3] = a22;
			break;
		}
	}
	double a11_1 = K00 + p1->getA_ij(0);
	double a12_1 = K01 + p1->getA_ij(1);
	double a21_1 = K01 + p1->getA_ij(2);
	double a22_1 = K11 + p1->getA_ij(3);
	double a11_2 = K22 + p2->getA_ij(0);
	double a12_2 = K23 + p2->getA_ij(1);
	double a21_2 = K23 + p2->getA_ij(2);
	double a22_2 = K33 + p2->getA_ij(3);
	p1->setA_ij(a11_1,a12_1,a21_1,a22_1);
	p2->setA_ij(a11_2,a12_2,a21_2,a22_2);
}

void edge::set_virtual_stiffness_matrix_value_only()
{
	/************************************
                  A
	              o    
		         / \    normal
		        /   \  |\ 
		       /     \ |
		      /       \|
		 P2  o---<<----o  P1  
	************************************/
	point3D *p1 = P1->get_node_id()->getvertex()->at(0);
	point3D *p2 = P2->get_node_id()->getvertex()->at(0);
	double x1 = p1->x();
	double y1 = p1->y();
	double x2 = p2->x();
	double y2 = p2->y();
	setnormal();
	double d = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
	double len = getlength();
	double xA = (x1 + x2)/2 + normal[0]/d*0.866*len;
	double yA = (y1 + y2)/2 + normal[1]/d*0.866*len;
	////////// triangle 12A //////////
	double b1t = y2 - yA;
	double c1t = xA - x2;
	double b2t = yA - y1;
	double c2t = x1 - xA;
	double b3t = y1 - y2;
	double c3t = x2 - x1;
	double A2t = c1t*y1 + c2t*y2 + c3t*yA;if (A2t<0) {A2t = -A2t;}
	double b1_b1 = b1t*b1t/A2t;double b1_c1 = b1t*c1t/A2t;double b1_b2 = b1t*b2t/A2t;double b1_c2 = b1t*c2t/A2t;double b1_b3 = b1t*b3t/A2t;double b1_c3 = b1t*c3t/A2t;
						       double c1_c1 = c1t*c1t/A2t;double c1_b2 = c1t*b2t/A2t;double c1_c2 = c1t*c2t/A2t;double c1_b3 = c1t*b3t/A2t;double c1_c3 = c1t*c3t/A2t;
											              double b2_b2 = b2t*b2t/A2t;double b2_c2 = b2t*c2t/A2t;double b2_b3 = b2t*b3t/A2t;double b2_c3 = b2t*c3t/A2t;
													    					         double c2_c2 = c2t*c2t/A2t;double c2_b3 = c2t*b3t/A2t;double c2_c3 = c2t*c3t/A2t;
																												double b3_b3 = b3t*b3t/A2t;double b3_c3 = b3t*c3t/A2t;
																																	       double c3_c3 = c3t*c3t/A2t;
    /////////////
	double K00 = (MU_virtual2_LAMBDA_virtual*b1_b1 + MU_virtual*c1_c1)/2;
	double K01 = (MU_virtual1_LAMBDA_virtual*b1_c1)/2;
	double K02 = (MU_virtual2_LAMBDA_virtual*b1_b2 + MU_virtual*c1_c2)/2;
	double K03 = (LAMBDA_virtual*b1_c2             + MU_virtual*c1_b2)/2;
	double K04 = (MU_virtual2_LAMBDA_virtual*b1_b3 + MU_virtual*c1_c3)/2;
	double K05 = (LAMBDA_virtual*b1_c3             + MU_virtual*c1_b3)/2;
	/////////////
	double K11 = (MU_virtual2_LAMBDA_virtual*c1_c1 + MU_virtual*b1_b1)/2;
	double K12 = (LAMBDA_virtual*c1_b2			   + MU_virtual*b1_c2)/2;
	double K13 = (MU_virtual2_LAMBDA_virtual*c1_c2 + MU_virtual*b1_b2)/2;
	double K14 = (LAMBDA_virtual*c1_b3             + MU_virtual*b1_c3)/2;
	double K15 = (MU_virtual2_LAMBDA_virtual*c1_c3 + MU_virtual*b1_b3)/2;
	/////////////
	double K22 = (MU_virtual2_LAMBDA_virtual*b2_b2 + MU_virtual*c2_c2)/2;
	double K23 = (MU_virtual1_LAMBDA_virtual*b2_c2)/2;
	double K24 = (MU_virtual2_LAMBDA_virtual*b2_b3 + MU_virtual*c2_c3)/2;
	double K25 = (LAMBDA_virtual*b2_c3             + MU_virtual*c2_b3)/2;
	/////////////
	double K33 = (MU_virtual2_LAMBDA_virtual*c2_c2 + MU_virtual*b2_b2)/2;
	double K34 = (LAMBDA_virtual*c2_b3             + MU_virtual*b2_c3)/2;
	double K35 = (MU_virtual2_LAMBDA_virtual*c2_c3 + MU_virtual*b2_b3)/2;
	/////////////
	double K44 = (MU_virtual2_LAMBDA_virtual*b3_b3 + MU_virtual*c3_c3)/2;
	double K45 = (MU_virtual1_LAMBDA_virtual*b3_c3)/2;
	/////////////
	double K55 = (MU_virtual2_LAMBDA_virtual*c3_c3 + MU_virtual*b3_b3)/2;
	//////////////////////////////////////////////////////////////////////
	virtual_triangle_stiffness_matrix[0][0] = K00;
	virtual_triangle_stiffness_matrix[0][1] = K01;
	virtual_triangle_stiffness_matrix[0][2] = K02;
	virtual_triangle_stiffness_matrix[0][3] = K03;
	virtual_triangle_stiffness_matrix[0][4] = K04;
	virtual_triangle_stiffness_matrix[0][5] = K05;
	///////
	virtual_triangle_stiffness_matrix[1][0] = K01;
	virtual_triangle_stiffness_matrix[1][1] = K11;
	virtual_triangle_stiffness_matrix[1][2] = K12;
	virtual_triangle_stiffness_matrix[1][3] = K13;
	virtual_triangle_stiffness_matrix[1][4] = K14;
	virtual_triangle_stiffness_matrix[1][5] = K15;
	///////
	virtual_triangle_stiffness_matrix[2][0] = K02;
	virtual_triangle_stiffness_matrix[2][1] = K12;
	virtual_triangle_stiffness_matrix[2][2] = K22;
	virtual_triangle_stiffness_matrix[2][3] = K23;
	virtual_triangle_stiffness_matrix[2][4] = K24;
	virtual_triangle_stiffness_matrix[2][5] = K25;
	///////
	virtual_triangle_stiffness_matrix[3][0] = K03;
	virtual_triangle_stiffness_matrix[3][1] = K13;
	virtual_triangle_stiffness_matrix[3][2] = K23;
	virtual_triangle_stiffness_matrix[3][3] = K33;
	virtual_triangle_stiffness_matrix[3][4] = K34;
	virtual_triangle_stiffness_matrix[3][5] = K35;
	///////
	virtual_triangle_stiffness_matrix[4][0] = K04;
	virtual_triangle_stiffness_matrix[4][1] = K14;
	virtual_triangle_stiffness_matrix[4][2] = K24;
	virtual_triangle_stiffness_matrix[4][3] = K34;
	virtual_triangle_stiffness_matrix[4][4] = K44;
	virtual_triangle_stiffness_matrix[4][5] = K45;
	///////
	virtual_triangle_stiffness_matrix[5][0] = K05;
	virtual_triangle_stiffness_matrix[5][1] = K15;
	virtual_triangle_stiffness_matrix[5][2] = K25;
	virtual_triangle_stiffness_matrix[5][3] = K35;
	virtual_triangle_stiffness_matrix[5][4] = K45;
	virtual_triangle_stiffness_matrix[5][5] = K55;
}

void edge::set_soften_stiffness_matrix()
{
	/************************************
	                 |\
	                 |
	                 |
	    P2 o---<<----o P1
		    \       /
			 \     /
			  \   /
			   \ /
			    o  A
	************************************/
	point3D *p1 = P1->get_node_id()->getvertex()->at(0);
	point3D *p2 = P2->get_node_id()->getvertex()->at(0);
	double x1 = p1->x();
	double y1 = p1->y();
	double x2 = p2->x();
	double y2 = p2->y();
	this->setnormal();
	double d = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
	double len = getlength();
	double xA = (x1 + x2)/2 - normal[0]/d*0.866*len;
	double yA = (y1 + y2)/2 - normal[1]/d*0.866*len;
	double b1t = y2 - yA;
	double c1t = xA - x2;
	double b2t = yA - y1;
	double c2t = x1 - xA;
	double b3t = y1 - y2;
	double c3t = x2 - x1;
	double A2t = c1t*y1 + c2t*y2 + c3t*yA;if (A2t<0) {A2t = -A2t;}
	double b1_b1 = b1t*b1t/A2t;double b1_c1 = b1t*c1t/A2t;double b1_b2 = b1t*b2t/A2t;double b1_c2 = b1t*c2t/A2t;double b1_b3 = b1t*b3t/A2t;double b1_c3 = b1t*c3t/A2t;
						       double c1_c1 = c1t*c1t/A2t;double c1_b2 = c1t*b2t/A2t;double c1_c2 = c1t*c2t/A2t;double c1_b3 = c1t*b3t/A2t;double c1_c3 = c1t*c3t/A2t;
											              double b2_b2 = b2t*b2t/A2t;double b2_c2 = b2t*c2t/A2t;double b2_b3 = b2t*b3t/A2t;double b2_c3 = b2t*c3t/A2t;
													    					         double c2_c2 = c2t*c2t/A2t;double c2_b3 = c2t*b3t/A2t;double c2_c3 = c2t*c3t/A2t;
																												double b3_b3 = b3t*b3t/A2t;double b3_c3 = b3t*c3t/A2t;
																																	       double c3_c3 = c3t*c3t/A2t;
	/////////////
	double K00 = (MU_soft2_LAMBDA_soft*b1_b1 + MU_soft*c1_c1)/2;
	double K01 = (MU_soft1_LAMBDA_soft*b1_c1)/2;
	double K02 = (MU_soft2_LAMBDA_soft*b1_b2 + MU_soft*c1_c2)/2;
	double K03 = (LAMBDA_soft*b1_c2          + MU_soft*c1_b2)/2;
	double K04 = (MU_soft2_LAMBDA_soft*b1_b3 + MU_soft*c1_c3)/2;
	double K05 = (LAMBDA_soft*b1_c3          + MU_soft*c1_b3)/2;
	/////////////
	double K11 = (MU_soft2_LAMBDA_soft*c1_c1 + MU_soft*b1_b1)/2;
	double K12 = (LAMBDA_soft*c1_b2          + MU_soft*b1_c2)/2;
	double K13 = (MU_soft2_LAMBDA_soft*c1_c2 + MU_soft*b1_b2)/2;
	double K14 = (LAMBDA_soft*c1_b3          + MU_soft*b1_c3)/2;
	double K15 = (MU_soft2_LAMBDA_soft*c1_c3 + MU_soft*b1_b3)/2;
	/////////////
	double K22 = (MU_soft2_LAMBDA_soft*b2_b2 + MU_soft*c2_c2)/2;
	double K23 = (MU_soft1_LAMBDA_soft*b2_c2)/2;
	double K24 = (MU_soft2_LAMBDA_soft*b2_b3 + MU_soft*c2_c3)/2;
	double K25 = (LAMBDA_soft*b2_c3          + MU_soft*c2_b3)/2;
	/////////////
	double K33 = (MU_soft2_LAMBDA_soft*c2_c2 + MU_soft*b2_b2)/2;
	double K34 = (LAMBDA_soft*c2_b3          + MU_soft*b2_c3)/2;
	double K35 = (MU_soft2_LAMBDA_soft*c2_c3 + MU_soft*b2_b3)/2;
	/////////////
	double K44 = (MU_soft2_LAMBDA_soft*b3_b3 + MU_soft*c3_c3)/2;
	double K45 = (MU_soft1_LAMBDA_soft*b3_c3)/2;
	/////////////
	double K55 = (MU_soft2_LAMBDA_soft*c3_c3 + MU_soft*b3_b3)/2;
	/////////////
	soften_triangle_stiffness_matrix[0][0] = K00;
	soften_triangle_stiffness_matrix[0][1] = K01;
	soften_triangle_stiffness_matrix[0][2] = K02;
	soften_triangle_stiffness_matrix[0][3] = K03;
	soften_triangle_stiffness_matrix[0][4] = K04;
	soften_triangle_stiffness_matrix[0][5] = K05;
	//////
	soften_triangle_stiffness_matrix[1][0] = K01;
	soften_triangle_stiffness_matrix[1][1] = K11;
	soften_triangle_stiffness_matrix[1][2] = K12;
	soften_triangle_stiffness_matrix[1][3] = K13;
	soften_triangle_stiffness_matrix[1][4] = K14;
	soften_triangle_stiffness_matrix[1][5] = K15;
	//////
	soften_triangle_stiffness_matrix[2][0] = K02;
	soften_triangle_stiffness_matrix[2][1] = K12;
	soften_triangle_stiffness_matrix[2][2] = K22;
	soften_triangle_stiffness_matrix[2][3] = K23;
	soften_triangle_stiffness_matrix[2][4] = K24;
	soften_triangle_stiffness_matrix[2][5] = K25;
	//////
	soften_triangle_stiffness_matrix[3][0] = K03;
	soften_triangle_stiffness_matrix[3][1] = K13;
	soften_triangle_stiffness_matrix[3][2] = K23;
	soften_triangle_stiffness_matrix[3][3] = K33;
	soften_triangle_stiffness_matrix[3][4] = K34;
	soften_triangle_stiffness_matrix[3][5] = K35;
	//////
	soften_triangle_stiffness_matrix[4][0] = K04;
	soften_triangle_stiffness_matrix[4][1] = K14;
	soften_triangle_stiffness_matrix[4][2] = K24;
	soften_triangle_stiffness_matrix[4][3] = K34;
	soften_triangle_stiffness_matrix[4][4] = K44;
	soften_triangle_stiffness_matrix[4][5] = K45;
	//////
	soften_triangle_stiffness_matrix[5][0] = K05;
	soften_triangle_stiffness_matrix[5][1] = K15;
	soften_triangle_stiffness_matrix[5][2] = K25;
	soften_triangle_stiffness_matrix[5][3] = K35;
	soften_triangle_stiffness_matrix[5][4] = K45;
	soften_triangle_stiffness_matrix[5][5] = K55;
}

cell::cell(int id) 
{
	index = id;
	burried = false;
	dead = false;
	set_dead = false;
	migrate = false;
	migrate_relax = false;
	migrate_angle = -1;
	mark_number = 0; // default
	mAngle = 0;
	soften = false;
	migrate_response = false;
	circle = true;
	interior_refresh = false;
	cell_status = 0;
	behavior = 0; // default set: static
	bhsteps[0] = bhsteps[1] = bhsteps[2] = bhsteps[3] = bhsteps[4] = 0;
	AABB[0] = NULL;AABB[1] = NULL;AABB[2] = NULL;AABB[3] = NULL;
	migrate_p[0] = NULL;migrate_p[0] = NULL;
	int N_SPECIES = 20; // YC: NOT GOOD, SHOULD GET THE VALUE FROM BIOLOGY CLASS.
	for (int i=0; i<N_SPECIES; i++)
	{
		state.push_back(0);
	}
}

cell::~cell()
{
	if (sides.size()>0) {sides.clear();}
	if (pairs.size()>0) {pairs.clear();}
	if (inner_t.size()>0) {inner_t.clear();}
	if (inner_p.size()>0) {inner_p.clear();}
	if (inner_p_slip.size()>0) {inner_p_slip.clear();}
}

void cell::setup_Circle()
{
	bool valid = true;
	if ((int)pairs.size()>0)
	{
		for (int i=0;i<(int)pairs.size();i++)
		{
			if ((int)pairs.at(i)->MP()->size()>0) 
			{
				valid = false;
				break;
			}
		}
	}
	circle = valid;
}

void cell::set_cell_type(int i)
{
	this->cell_type = i;
	double ri = sqrt(this->get_area()/2/PI);
	this->set_r(ri);
	double color_r = color_class_R[i];
	double color_g = color_class_G[i];
	double color_b = color_class_B[i];
	this->set_initial_area(this->get_initial_area());
	this->setup_tcoef(tension_edge_class[i]);
	this->setup_pcoef(pressure_edge_class[i]);
	this->setup_mcoef(mass_edge_class[i]);
	this->set_lame(0,elastic_MU[i]);
	this->set_lame(1,elastic_LAMBDA[i]);
	this->set_cell_color(color_r,color_g,color_b);
	int pn = (int)this->get_pairs()->size();
	for (int j=0;j<pn;j++)
	{
		if (this->get_pairs()->at(j)->get_redundant()) continue;
		cell *C1 = this;
		cell *C2 = NULL;
		if (this->get_pairs()->at(j)->get_I1()==this) {C2 = this->get_pairs()->at(j)->get_I2();}
		else                                          {C2 = this->get_pairs()->at(j)->get_I1();}
		//int iptype = return_type(C1->get_cell_type(),C2->get_cell_type());
		//this->get_pairs()->at(j)->set_adhesion(adhesion_pair_class[iptype]);
		//this->get_pairs()->at(j)->set_friction(friction_pair_class[iptype]);
		double adhesion_force = adhesion_pair_array[C1->get_cell_type()][C2->get_cell_type()];
		this->get_pairs()->at(j)->set_adhesion(adhesion_force);
		this->get_pairs()->at(j)->set_friction(adhesion_force);
	}
}

double cell::shape_based_division_angle()
{
	double angle = 0;
	double D02 = (AABB[0]->x() - AABB[2]->x())*(AABB[0]->x() - AABB[2]->x()) + (AABB[0]->y() - AABB[2]->y())*(AABB[0]->y() - AABB[2]->y());
	double D13 = (AABB[1]->x() - AABB[3]->x())*(AABB[1]->x() - AABB[3]->x()) + (AABB[1]->y() - AABB[3]->y())*(AABB[1]->y() - AABB[3]->y());
	if (D02>D13) // D02
	{
		double vector_x = AABB[0]->x() - AABB[2]->x();
		double vector_y = AABB[0]->y() - AABB[2]->y();
		double dist = sqrt(vector_x*vector_x +vector_y*vector_y);
		double x1 = vector_x/dist;
		double y1 = vector_y/dist;
		if      (x1>=0 && y1>=0)
		{
			angle = asin(y1)*180.0/PI;
		}
		else if (x1<0 && y1>=0)
		{
			angle = 180 - asin(y1)*180.0/PI;
		}
		else if (x1<0 && y1<0)
		{
			angle = 180 + asin(-y1)*180.0/PI;
		}
		else if (x1>=0 && y1<0)
		{
			angle = 360 - asin(-y1)*180.0/PI;
		}
		angle-=90;
		if (angle<0) {angle+=360;}
	}
	else // D13
	{
		double vector_x = AABB[3]->x() - AABB[1]->x();
		double vector_y = AABB[3]->y() - AABB[1]->y();
		double dist = sqrt(vector_x*vector_x +vector_y*vector_y);
		double x1 = vector_x/dist;
		double y1 = vector_y/dist;
		if      (x1>=0 && y1>=0)
		{
			angle = asin(y1)*180.0/PI;
		}
		else if (x1<0 && y1>=0)
		{
			angle = 180 - asin(y1)*180.0/PI;
		}
		else if (x1<0 && y1<0)
		{
			angle = 180 + asin(-y1)*180.0/PI;
		}
		else if (x1>=0 && y1<0)
		{
			angle = 360 - asin(-y1)*180.0/PI;
		}
		angle-=90;
		if (angle<0) {angle+=360;}
	}
	return angle;
}

void cell::set_area()
{
	double A = 0;
	for (int i=0;i<(int)sides.size();i++)
	{
		double x1 = sides.at(i)->p1()->x();
		double y1 = sides.at(i)->p1()->y();
		double x2 = sides.at(i)->p2()->x();
		double y2 = sides.at(i)->p2()->y();
		A+=0.5*(x1*y2 - x2*y1);
	}
	A = abs(A);
	area = A;
}

void cell::refresh_AABB()
{   
	// 0:xmin, 1:ymin, 2:xmax, 3:ymax
	double xmin = 100000;
	double ymin = 100000;
	double xmax =-100000;
	double ymax =-100000;

	for (int i = 0;i<(int)sides.size();i++)
	{
		double x = sides.at(i)->p1()->x();
		double y = sides.at(i)->p1()->y();
		if      (x<xmin) {xmin = x;AABB[0] = sides.at(i)->p1();}
		if      (x>xmax) {xmax = x;AABB[2] = sides.at(i)->p1();}
		if      (y<ymin) {ymin = y;AABB[1] = sides.at(i)->p1();}
		if      (y>ymax) {ymax = y;AABB[3] = sides.at(i)->p1();}
	}
}

void cell::center_refresh()
{
	double x = 0;
	double y = 0;
	int ns = (int)sides.size();
	for (int i=0;i<ns;i++)
	{
		x += sides.at(i)->p1()->x();
		y += sides.at(i)->p1()->y();
	}
	x = x/ns;
	y = y/ns;
	center[0] = x;
	center[1] = y;
}

void cell::setup_Burry()
{
	bool valid = true;
	int szn = sides.size();
	for (int i=0;i<szn;i++)
	{
		if (!sides[i]->get_attach())
		{
			valid = false;
			break;
		}
	}
	burried = valid;
}

void cell::removeside(edge *Edge)
{
	int flag = -1;
	int sn = (int)sides.size();
	for (int i=0;i<(int)sides.size();i++)
	{
		if (Edge==sides.at(i)) 
		{
			flag = i;
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==sn-1)
		{
			sides.pop_back();
		}
		else
		{
			for (int i=flag;i<sn-1;i++)
			{
				sides.at(i) = sides.at(i+1);
			}
			sides.pop_back();
		}
	}
}

void cell::removepair(interpair *Pair)
{
	int flag = -1;
	int pn = (int)pairs.size();
	for (int i=0;i<(int)pairs.size();i++)
	{
		if (Pair==pairs.at(i))
		{
			flag = i;
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==pn-1)
		{
			pairs.pop_back();
		}
		else
		{
			for (int i=flag;i<pn-1;i++)
			{
				pairs.at(i) = pairs.at(i+1);
			}
			pairs.pop_back();
		}
	}
}

bool cell::test_point_inside(double x, double y)
{
	bool valid = false;
	int sn = (int)sides.size();
	int pass_before = 0;
	int pass_after = 0;
	for (int i=0;i<sn;i++)
	{
		double x1 = sides.at(i)->p1()->x();
		double y1 = sides.at(i)->p1()->y();
		double x2 = sides.at(i)->p2()->x();
		double y2 = sides.at(i)->p2()->y();
		double t = 0;
		if (y1==y2) continue;
		t = (y-y2)/(y1-y2);
		if (t<=0 || t>1) continue;
		double xt = t*x1 + (1-t)*x2;
		if (xt<x) {pass_before++;}
		if (xt>x) {pass_after++;}
	}
	if (pass_before%2==1 && pass_after%2==1) {valid = true;}
	return valid;
}

void cell::clear_inner_p()
{
	int pn = (int)inner_p.size();
	for (int i=pn-1;i>=0;i--) 
	{
		delete inner_p.at(i);
		inner_p.at(i) = NULL;
		inner_p.pop_back();
	} 
}

void cell::clear_inner_t()
{
	int tn = (int)inner_t.size();
	for (int i=tn-1;i>=0;i--) 
	{
		delete inner_t.at(i);
		inner_t.at(i) = NULL;
		inner_t.pop_back();
	} 
}

void cell::clear_inner_p_slip()
{
	int pn = (int)inner_p_slip.size();
	for (int i=pn-1;i>=0;i--) 
	{
		delete inner_p_slip.at(i);
		inner_p_slip.at(i) = NULL;
		inner_p_slip.pop_back();
	} 
}

void cell::remove_filling_cell(cell *DC)
{
	int flag = -1;
	int mn = (int)migrate_dead_filling_cells.size();
	rear_fc *FC = NULL;
	for (int i=0;i<(int)migrate_dead_filling_cells.size();i++)
	{
		if (DC==migrate_dead_filling_cells.at(i)->getC()) 
		{
			flag = i;
			FC = migrate_dead_filling_cells.at(i);
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==mn-1)
		{
			migrate_dead_filling_cells.pop_back();
		}
		else
		{
			for (int i=flag;i<mn-1;i++)
			{
				migrate_dead_filling_cells.at(i) = migrate_dead_filling_cells.at(i+1);
			}
			migrate_dead_filling_cells.pop_back();
		}
		delete FC;
		FC = NULL;
	}
}

void cell::clear_filling_cells()
{
	int sn = (int)migrate_dead_filling_cells.size();
	for (int i=sn-1;i>=0;i--)
	{
		delete migrate_dead_filling_cells.at(i);
		migrate_dead_filling_cells.at(i) = NULL;
		migrate_dead_filling_cells.pop_back();
	}
}

interpair::interpair(int n, cell* A, cell* B) 
{
	index = n;
	I1 = A;
	I2 = B;
	attach = false;
	redundant = false;
}

interpair::~interpair()
{
	int pn = (int)P.size();
	for (int i=0;i<pn;i++)
	{
		P.at(i) = NULL;
	}
	P.clear();
}

void interpair::remove_MP(mergepair *MP)
{
	int flag = -1;
	int mn = (int)P.size();
	for (int i=0;i<(int)P.size();i++)
	{
		if (MP==P.at(i)) 
		{
			flag = i;
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==mn-1)
		{
			P.pop_back();
		}
		else
		{
			for (int i=flag;i<mn-1;i++)
			{
				P.at(i) = P.at(i+1);
			}
			P.pop_back();
		}
	}
	delete MP;
	MP = NULL;
}

bool interpair::check_detection()
{
	bool valid = false;
	int pn = P.size();
	if (pn>0)
	{
		point3D *s1 = P[0]->S1();
		point3D *t1 = P[0]->T1();
		point3D *s2 = P[0]->S2();
		point3D *t2 = P[0]->T2();
		if (s1->getba()->get_attach() &&
			t1->getab()->get_attach() &&
			s2->getba()->get_attach() &&
			t2->getab()->get_attach())
		{
			valid = true;
		}
	}
	return valid;
}

mergepair::mergepair(int n, point3D* A, point3D* B, point3D* C, point3D* D)
{
	index = n;
	s1 = A;
	t1 = B;
	s2 = C;
	t2 = D;
	F_s1 = 0;
	F_t1 = 0;
	F_s2 = 0;
	F_t2 = 0;
	Fr_s1 = 0;
	Fr_t1 = 0;
	Fr_dir = -1; // null
}

mergepair::~mergepair()
{
	/*s1 = NULL;
	t1 = NULL;
	s2 = NULL;
	t2 = NULL;*/
}

node::node(int n)
{
	index = n;
	in_global = true;
	joint = false;
};

node::~node()
{
	for (int i=0;i<(int)P.size();i++)
	{
		P.at(i) = NULL;
	}
	P.clear();
}

void node::removevertex(point3D* A)
{
	int flag = -1;
	int sn = (int)P.size();
	for (int i=0;i<(int)P.size();i++)
	{
		if (A==P.at(i)) 
		{
			flag = i;
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==sn-1)
		{
			P.pop_back();
		}
		else
		{
			for (int i=flag;i<sn-1;i++)
			{
				P.at(i) = P.at(i+1);
			}
			P.pop_back();
		}
	}
}

void node::pushvertex(point3D* A)
{
	bool exist = false;
	for (int i=0;i<(int)P.size();i++)
	{
		if (A==P.at(i)) {exist = true;break;}
	}
	if (!exist) {P.push_back(A);}
}

triangle::triangle(int id, point3D *a, point3D *b, point3D *c) 
{
	index = id;
	A = a;
	B = b;
	C = c;
	pick = false;
	Stiffness_matrix[0][0] = 0;Stiffness_matrix[0][1] = 0;Stiffness_matrix[0][2] = 0;Stiffness_matrix[0][3] = 0;Stiffness_matrix[0][4] = 0;Stiffness_matrix[0][5] = 0;
	Stiffness_matrix[1][0] = 0;Stiffness_matrix[1][1] = 0;Stiffness_matrix[1][2] = 0;Stiffness_matrix[1][3] = 0;Stiffness_matrix[1][4] = 0;Stiffness_matrix[1][5] = 0;
	Stiffness_matrix[2][0] = 0;Stiffness_matrix[2][1] = 0;Stiffness_matrix[2][2] = 0;Stiffness_matrix[2][3] = 0;Stiffness_matrix[2][4] = 0;Stiffness_matrix[2][5] = 0;
	Stiffness_matrix[3][0] = 0;Stiffness_matrix[3][1] = 0;Stiffness_matrix[3][2] = 0;Stiffness_matrix[3][3] = 0;Stiffness_matrix[3][4] = 0;Stiffness_matrix[3][5] = 0;
	Stiffness_matrix[4][0] = 0;Stiffness_matrix[4][1] = 0;Stiffness_matrix[4][2] = 0;Stiffness_matrix[4][3] = 0;Stiffness_matrix[4][4] = 0;Stiffness_matrix[4][5] = 0;
	Stiffness_matrix[5][0] = 0;Stiffness_matrix[5][1] = 0;Stiffness_matrix[5][2] = 0;Stiffness_matrix[5][3] = 0;Stiffness_matrix[5][4] = 0;Stiffness_matrix[5][5] = 0;
}

triangle::~triangle()
{
	/*A = NULL;
	B = NULL;
	C = NULL;*/
}

void triangle::set_Stiffness_matrix()
{
	double x1 = A->x();
	double y1 = A->y();
	double x2 = B->x();
	double y2 = B->y();
	double x3 = C->x();
	double y3 = C->y();
	double T_MU = this->get_lame(0);
	double T_LAMBDA = this->get_lame(1);
	double T_MU2_LAMBDA = T_MU*2 + T_LAMBDA;
	double T_MU1_LAMBDA = T_MU + T_LAMBDA;
	/*************************************************
	                 
	 -1  | 00 01 02 | -1                | a1 b1 c1 |
	A  = | 10 11 12 |   ==inv=> 1/det(A)| a2 b2 c2 |
	     | 20 21 22 |                   | a3 b3 c3 |
	
	a1 = 11X22 - 12X21 b1 = 02X21 - 01X22 c1 = 01X12 - 02X11
	a2 = 12X20 - 10X22 b2 = 00X22 - 02X20 c2 = 02X10 - 00X12
	a3 = 10X21 - 11X20 b3 = 20X01 - 00X21 c3 = 00X11 - 01X10
	-------------------------------------------------
	a1 = | 11 12 | b1 = -| 01 02 | c1 = | 01 02 |
	     | 21 22 |		 | 21 22 |		| 11 12 |

    a2 = -| 10 12 | b2 = | 00 02 | c2 = -| 00 02 |
		  | 20 22 |		 | 20 22 |		 | 10 12 |
    
	a3 = | 10 11 | b3 = -| 00 01 | c3 = | 00 01 | 
		 | 20 21 |		 | 20 21 |		| 10 11 |
	*************************************************/
	double b1t = y2 - y3;
	double c1t = x3 - x2;
	double b2t = y3 - y1;
	double c2t = x1 - x3;
	double b3t = y1 - y2;
	double c3t = x2 - x1; // 10E-6 m
	double A2t = c1t*y1 + c2t*y2 + c3t*y3;if (A2t<0) {A2t = -A2t;} // 10E-12 m^2
	/////////////////////////////////
	double b1_b1 = b1t*b1t/A2t;double b1_c1 = b1t*c1t/A2t;double b1_b2 = b1t*b2t/A2t;double b1_c2 = b1t*c2t/A2t;double b1_b3 = b1t*b3t/A2t;double b1_c3 = b1t*c3t/A2t;
						       double c1_c1 = c1t*c1t/A2t;double c1_b2 = c1t*b2t/A2t;double c1_c2 = c1t*c2t/A2t;double c1_b3 = c1t*b3t/A2t;double c1_c3 = c1t*c3t/A2t;
														  double b2_b2 = b2t*b2t/A2t;double b2_c2 = b2t*c2t/A2t;double b2_b3 = b2t*b3t/A2t;double b2_c3 = b2t*c3t/A2t;
													    						     double c2_c2 = c2t*c2t/A2t;double c2_b3 = c2t*b3t/A2t;double c2_c3 = c2t*c3t/A2t;
																											    double b3_b3 = b3t*b3t/A2t;double b3_c3 = b3t*c3t/A2t;
																																	       double c3_c3 = c3t*c3t/A2t;
	/////////////
	double K00 = T_MU2_LAMBDA*b1_b1 + T_MU*c1_c1;//MU2_LAMBDA*b1_b1 + MU*c1_c1;
	double K01 = T_MU1_LAMBDA*b1_c1;             //MU1_LAMBDA*b1_c1;
	double K02 = T_MU2_LAMBDA*b1_b2 + T_MU*c1_c2;//MU2_LAMBDA*b1_b2 + MU*c1_c2;
	double K03 = T_LAMBDA*b1_c2     + T_MU*c1_b2;//LAMBDA*b1_c2     + MU*c1_b2;
	double K04 = T_MU2_LAMBDA*b1_b3 + T_MU*c1_c3;//MU2_LAMBDA*b1_b3 + MU*c1_c3;
	double K05 = T_LAMBDA*b1_c3     + T_MU*c1_b3;//LAMBDA*b1_c3     + MU*c1_b3;
	/////////////
	double K11 = T_MU2_LAMBDA*c1_c1 + T_MU*b1_b1;//MU2_LAMBDA*c1_c1 + MU*b1_b1;
	double K12 = T_LAMBDA*c1_b2     + T_MU*b1_c2;//LAMBDA*c1_b2     + MU*b1_c2;
	double K13 = T_MU2_LAMBDA*c1_c2 + T_MU*b1_b2;//MU2_LAMBDA*c1_c2 + MU*b1_b2;
	double K14 = T_LAMBDA*c1_b3     + T_MU*b1_c3;//LAMBDA*c1_b3     + MU*b1_c3;
	double K15 = T_MU2_LAMBDA*c1_c3 + T_MU*b1_b3;//MU2_LAMBDA*c1_c3 + MU*b1_b3;
	/////////////
	double K22 = T_MU2_LAMBDA*b2_b2 + T_MU*c2_c2;//MU2_LAMBDA*b2_b2 + MU*c2_c2;
	double K23 = T_MU1_LAMBDA*b2_c2;             //MU1_LAMBDA*b2_c2;
	double K24 = T_MU2_LAMBDA*b2_b3 + T_MU*c2_c3;//MU2_LAMBDA*b2_b3 + MU*c2_c3;
	double K25 = T_LAMBDA*b2_c3     + T_MU*c2_b3;//LAMBDA*b2_c3     + MU*c2_b3;
	/////////////
	double K33 = T_MU2_LAMBDA*c2_c2 + T_MU*b2_b2;//MU2_LAMBDA*c2_c2 + MU*b2_b2;
	double K34 = T_LAMBDA*c2_b3     + T_MU*b2_c3;//LAMBDA*c2_b3     + MU*b2_c3;
	double K35 = T_MU2_LAMBDA*c2_c3 + T_MU*b2_b3;//MU2_LAMBDA*c2_c3 + MU*b2_b3;
	/////////////
	double K44 = T_MU2_LAMBDA*b3_b3 + T_MU*c3_c3;//MU2_LAMBDA*b3_b3 + MU*c3_c3;
	double K45 = T_MU1_LAMBDA*b3_c3;             //MU1_LAMBDA*b3_c3;
	/////////////
	double K55 = T_MU2_LAMBDA*c3_c3 + T_MU*b3_b3;//MU2_LAMBDA*c3_c3 + MU*b3_b3;
	/////////////
	Stiffness_matrix[0][0] = K00/2;Stiffness_matrix[0][1] = K01/2;Stiffness_matrix[0][2] = K02/2;Stiffness_matrix[0][3] = K03/2;Stiffness_matrix[0][4] = K04/2;Stiffness_matrix[0][5] = K05/2;
	Stiffness_matrix[1][0] = K01/2;Stiffness_matrix[1][1] = K11/2;Stiffness_matrix[1][2] = K12/2;Stiffness_matrix[1][3] = K13/2;Stiffness_matrix[1][4] = K14/2;Stiffness_matrix[1][5] = K15/2;
	Stiffness_matrix[2][0] = K02/2;Stiffness_matrix[2][1] = K12/2;Stiffness_matrix[2][2] = K22/2;Stiffness_matrix[2][3] = K23/2;Stiffness_matrix[2][4] = K24/2;Stiffness_matrix[2][5] = K25/2;
	Stiffness_matrix[3][0] = K03/2;Stiffness_matrix[3][1] = K13/2;Stiffness_matrix[3][2] = K23/2;Stiffness_matrix[3][3] = K33/2;Stiffness_matrix[3][4] = K34/2;Stiffness_matrix[3][5] = K35/2;
	Stiffness_matrix[4][0] = K04/2;Stiffness_matrix[4][1] = K14/2;Stiffness_matrix[4][2] = K24/2;Stiffness_matrix[4][3] = K34/2;Stiffness_matrix[4][4] = K44/2;Stiffness_matrix[4][5] = K45/2;
	Stiffness_matrix[5][0] = K05/2;Stiffness_matrix[5][1] = K15/2;Stiffness_matrix[5][2] = K25/2;Stiffness_matrix[5][3] = K35/2;Stiffness_matrix[5][4] = K45/2;Stiffness_matrix[5][5] = K55/2;
}

le_sm::le_sm(int n)
{
	 index = n;
	 stiffness_matrix_12[0] = 0;
	 stiffness_matrix_12[1] = 0;
	 stiffness_matrix_12[2] = 0;
	 stiffness_matrix_12[2] = 0;
	 stiffness_matrix_21[0] = 0;
	 stiffness_matrix_21[1] = 0;
	 stiffness_matrix_21[2] = 0;
	 stiffness_matrix_21[2] = 0;
}

le_sm::~le_sm()
{
	N1 = NULL;
	N2 = NULL;
}

rear_fc::rear_fc(int n)
{
	index = n;
	initial = false;
}

rear_fc::~rear_fc()
{
	C = NULL;
	int ncn = (int)NCs.size();
	for (int i=0;i<ncn;i++)
	{
		NCs.at(i) = NULL;
	}
	NCs.clear();
}

void rear_fc::push_NCs(cell *C1)
{
	bool exist = false;
	for (int i=0;i<(int)NCs.size();i++)
	{
		if (C1==NCs.at(i)) {exist = true;break;}
	}
	if (!exist) {NCs.push_back(C1);}
}

void rear_fc::remove_NCs(cell *C1)
{
	int flag = -1;
	int mn = (int)NCs.size();
	for (int i=0;i<(int)NCs.size();i++)
	{
		if (C1==NCs.at(i)) 
		{
			flag = i;
			break;
		}
	}
	if (flag>=0)
	{
		if (flag==mn-1)
		{
			NCs.pop_back();
		}
		else
		{
			for (int i=flag;i<mn-1;i++)
			{
				NCs.at(i) = NCs.at(i+1);
			}
			NCs.pop_back();
		}
	}
}

void rear_fc::clear_NCs()
{
	while (NCs.size())
	{
		NCs.pop_back();
	}
}

#endif

