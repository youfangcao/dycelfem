// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/***************************************
***    Project:  Cell Growth         ***
***    File:     dbReader.cpp        ***
***    Author:   Jieling Zhao        ***
***                                  ***
***    Created on March 19, 2012     ***
***************************************/

#ifndef DBREADER_CPP
#define DBREADER_CPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstring>
#include "dbReader.h"
#include "biology.h"
#include <common/dycelfem_config.h>

using namespace std;

dbReader::dbReader()
{
	viewVolume *V = new viewVolume();
	vv = V;
	initial = 0;
	Complex = NULL;
	biology = new Biology();
}

viewVolume* dbReader::getViewVolume()
{
	return vv;
}

void dbReader::read(char* fN, int number) 
{
	construct_map_of_para();
	UserInitialParas ();

	ifstream fileIn;
	fileIn.open(fN);

	if(!fileIn.is_open())
	{
		printf("Error: %s does not exist!\n", fN);
		exit(1);
	}

	// boolean flag
	bool point_based_flag = false;
	bool center_based_flag = false;
	bool vert_list_flag = false;
	bool edge_list_flag = false;
	bool cell_list_flag = false;
	bool pair_list_flag = false;
	bool node_list_flag = false;
	bool dead_vert_list_flag = false;
	bool edge_pair_list_flag = false;
	bool migr_cell_list_flag = false;
	int cn=0;
	double range[4]; // 0:xmin 1:ymin 2:xmax 3:ymax
	range[0] = 1000000000000;
	range[1] = 1000000000000;
	range[2] =-1000000000000;
	range[3] =-1000000000000;

	string ln;

	// read database file line-by-line
	while(getline(fileIn, ln))
	{
		if (cell_list_flag && center_based_flag)
		{
			string strtmp_END_VERTEX_LIST = ln.substr(0,13);
			if (strcmp(strtmp_END_VERTEX_LIST.c_str(),"END_CELL_LIST"))
			{
				vector<string> array1;
				vector<string> array2;
				str_split(ln,array1,"#");
				str_split(array1[0],array2,",");
				double vx = atof(array2[0].c_str());
				double vy = atof(array2[1].c_str());
				int ct = atoi(array2[2].c_str());
				double ri = CR;
				cell* A = new cell(cn);
				A->set_center(vx,vy);
				A->set_r(ri);
				A->set_cell_type(ct);
				A->set_life(number);
				A->setup_tcoef(tension_edge_class[ct]);
				A->setup_pcoef(pressure_edge_class[ct]);
				A->setup_mcoef(mass_edge_class[ct]);
				cellList.push_back(A);
				if (range[0]>vx) {range[0] = vx;}
				if (range[1]>vy) {range[1] = vy;}
				if (range[2]<vx) {range[2] = vx;}
				if (range[3]<vy) {range[3] = vy;}
				double bin = PI*2.0/SN;
				for (int i=0;i<=SN;i++)
				{
					if (i==0)
					{ 
						double xi = vx + cos(i*bin)*ri;
						double yi = vy + sin(i*bin)*ri;
						point3D *B = new point3D(cn*SN+i,xi,yi);
						vertexList.push_back(B);
						node *B1 = new node(cn*SN+i);
						nodeList.push_back(B1);
						nodeList[cn*SN+i]->getvertex()->push_back(vertexList[cn*SN+i]);
						vertexList[cn*SN+i]->set_node_check(true);
						vertexList[cn*SN+i]->set_node_id(nodeList[cn*SN+i]);
					}
					else if (i>0 && i<SN)
					{
						double xi = vx + cos(i*bin)*ri;
						double yi = vy + sin(i*bin)*ri;
						point3D *B = new point3D(cn*SN+i,xi,yi);
						vertexList.push_back(B);
						node *B1 = new node(cn*SN+i);
						nodeList.push_back(B1);
						nodeList[cn*SN+i]->pushvertex(vertexList[cn*SN+i]);
						vertexList[cn*SN+i]->set_node_check(true);
						vertexList[cn*SN+i]->set_node_id(nodeList[cn*SN+i]);
						edge *C = new edge(cn*SN+i-1,vertexList[cn*SN+i-1],vertexList[cn*SN+i],cellList[cn]);
						C->setnormal();
						C->set_in_cell(1);
						edgeList.push_back(C);
						cellList[cn]->pushside(C);
						vertexList[cn*SN+i-1]->setab(C);
						vertexList[cn*SN+i]->setba(C);
						vertexList[cn*SN+i-1]->setfp(vertexList[cn*SN+i]);
						vertexList[cn*SN+i]->setrp(vertexList[cn*SN+i-1]);
					}
					else if (i==SN)
					{
						edge *C = new edge(cn*SN+i-1,vertexList[cn*SN+i-1],vertexList[cn*SN],cellList[cn]);
						C->setnormal();
						C->set_in_cell(1);
						edgeList.push_back(C);
						cellList[cn]->pushside(C);
						vertexList[cn*SN+i-1]->setab(C);
						vertexList[cn*SN]->setba(C);
						vertexList[cn*SN+i-1]->setfp(vertexList[cn*SN]);
						vertexList[cn*SN]->setrp(vertexList[cn*SN+i-1]);
					}
				}
				cellList[cn]->set_AABB(vertexList[cn*SN+SN/2],vertexList[cn*SN+SN/4*3],vertexList[cn*SN],vertexList[cn*SN+SN/4]);
				double color_r = color_class_R[ct];
				double color_g = color_class_G[ct];
				double color_b = color_class_B[ct];
				cellList[cn]->set_cell_color(color_r,color_g,color_b);
				cellList[cn]->set_cell_diffusion_color(1,1,1);
				cellList[cn]->set_GR(growth_rate_class[ct]);
				double lame_MU = elastic_MU[ct];
				double lame_LAMBDA = elastic_LAMBDA[ct];
				cellList[cn]->set_lame(0,lame_MU);
				cellList[cn]->set_lame(1,lame_LAMBDA);
				if      (cellList[cn]->get_GR()>0) {cellList[cn]->set_cell_status(1);}
				else if (cellList[cn]->get_GR()<0) {cellList[cn]->set_cell_status(-1);}
				for (int i=0;i<biology->getNumSpecies();i++) {cellList[cn]->set_cell_state(i,0);}
				cn++;
			}
			if (!strcmp(strtmp_END_VERTEX_LIST.c_str(),"END_CELL_LIST"))
			{cell_list_flag = false;}
		}
		if (vert_list_flag && point_based_flag)
		{
			string strtmp_END_VERT_LIST = ln.substr(0,13);
			if (strcmp(strtmp_END_VERT_LIST.c_str(),"END_VERT_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				double vx = atof(array1[0].c_str());
				double vy = atof(array1[1].c_str());
				int vid = atoi(array1[2].c_str());
				int in_cell = atoi(array1[3].c_str());
				double sf = atof(array1[4].c_str());
				double sff = atof(array1[5].c_str());
				double efx = atof(array1[6].c_str());
				double efy = atof(array1[7].c_str());
				double bfx = atof(array1[8].c_str());
				double bfy = atof(array1[9].c_str());
				point3D *V = new point3D(vid,vx,vy);
				vertexList.push_back(V);
				if      (in_cell==0) {vertexList[vid]->set_in_cell(false);}
				else if (in_cell==1) {vertexList[vid]->set_in_cell(true);}
				vertexList[vid]->set_stored_apart_force(sf);
				vertexList[vid]->set_stored_friction_force(sff);
				vertexList[vid]->set_elastic_force(efx,efy);
				vertexList[vid]->set_adhesion_break(bfx,bfy);
				vertexList[vid]->set_node_id(NULL);
			}
			else if (!strcmp(strtmp_END_VERT_LIST.c_str(),"END_VERT_LIST"))
			{vert_list_flag = false;}
		}
		if (edge_list_flag && point_based_flag)
		{
			string strtmp_END_EDGE_LIST = ln.substr(0,13);
			if (strcmp(strtmp_END_EDGE_LIST.c_str(),"END_EDGE_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				int eid = atoi(array1[0].c_str());
				int p1id = atoi(array1[1].c_str());
				int p2id = atoi(array1[2].c_str());
				edge *E = new edge(eid,vertexList[p1id],vertexList[p2id],NULL);
				edgeList.push_back(E);
				point3D *p1 = vertexList[p1id];
				point3D *p2 = vertexList[p2id];
				bool ing1 = p1->get_in_cell();
				bool ing2 = p2->get_in_cell();
				if (ing1 && ing2)
				{
					vertexList[p1id]->setfp(vertexList[p2id]);
					vertexList[p2id]->setrp(vertexList[p1id]);
					vertexList[p1id]->setab(edgeList[eid]);
					vertexList[p2id]->setba(edgeList[eid]);
				}
			}
			else if (!strcmp(strtmp_END_EDGE_LIST.c_str(),"END_EDGE_LIST"))
			{edge_list_flag = false;}
		}
		if (cell_list_flag && point_based_flag)
		{
			string strtmp_END_CELL_LIST = ln.substr(0,13);
			if (strcmp(strtmp_END_CELL_LIST.c_str(),"END_CELL_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				int cid = atoi(array1[0].c_str());
				int dead = atoi(array1[1].c_str());
				int set_dead = atoi(array1[2].c_str());
				int migrate = atoi(array1[3].c_str());
				int mAngle = atoi(array1[4].c_str());
				int type = atoi(array1[5].c_str());
				double iarea = atof(array1[8].c_str());
				int nspecies = atoi(array1[9].c_str()); // YFC: READ STATE VECTOR LENGTH
				cell *C = new cell(cid);
				cellList.push_back(C);
				if      (dead==0) {cellList[cid]->setup_Dead(false);}
				else if (dead==1) {cellList[cid]->setup_Dead(true);}
				if      (set_dead==0) {cellList[cid]->set_set_dead(false);}
				else if (set_dead==1) {cellList[cid]->set_set_dead(true);}
				if      (migrate==0) {cellList[cid]->set_migrate(false);}
				else if (migrate==1) {cellList[cid]->set_migrate(true);}
				//if      (soften==0) {cellList[cid]->set_soften(false);}
				//else if (soften==1) {cellList[cid]->set_soften(true);}
				cellList[cid]->set_cell_type(type);
				cellList[cid]->set_r(CR);
				cellList[cid]->set_life(number);
				cellList[cid]->set_initial_area(AR/2);
				cellList[cid]->setup_tcoef(tension_edge_class[type]);
				cellList[cid]->setup_pcoef(pressure_edge_class[type]);
				cellList[cid]->setup_mcoef(mass_edge_class[type]);
                cellList[cid]->set_center(atof(array1[6].c_str()), atof(array1[7].c_str()));
				// YFC: READ CELL STATE
				for (int j=0; j<nspecies; j++)
				{
					cellList[cid]->set_cell_state(j,(int)atof(array1[10+j].c_str()));
				}
				double color_r = color_class_R[type];
				double color_g = color_class_G[type];
				double color_b = color_class_B[type];
				cellList[cid]->set_cell_color(color_r,color_g,color_b);
				cellList[cid]->set_cell_diffusion_color(1,1,1);
				double lame_MU = elastic_MU[type];
				double lame_LAMBDA = elastic_LAMBDA[type];
				cellList[cid]->set_lame(0,lame_MU);
				cellList[cid]->set_lame(1,lame_LAMBDA);
				double gr = growth_rate_class[type];
				cellList[cid]->set_GR(gr);
				int cstatus = 0;
				if      (gr>0) {cstatus = 1;}
				else if (gr<0) {cstatus = -1;}
				cellList[cid]->set_cell_status(cstatus);
				int sdn = (int)array1.size();
				for (int ii=10+nspecies;ii<sdn;ii++)
				{
					int eid = atoi(array1[ii].c_str());
					edgeList[eid]->reset_Cell(cellList[cid]);
					cellList[cid]->pushside(edgeList[eid]);
				}
			}
			else if (!strcmp(strtmp_END_CELL_LIST.c_str(),"END_CELL_LIST"))
			{cell_list_flag = false;}
		}
		if (pair_list_flag && point_based_flag)
		{
			string strtmp_END_PAIR_LIST = ln.substr(0,13);
			if (strcmp(strtmp_END_PAIR_LIST.c_str(),"END_PAIR_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				int pid = atoi(array1[0].c_str());
				int reid = atoi(array1[1].c_str());
				int c1id = atoi(array1[2].c_str());
				int c2id = atoi(array1[3].c_str());
				interpair *Pl = new interpair(pid,cellList[c1id],cellList[c2id]);
				collisionpairList.push_back(Pl);
				//int iptype = interpair_type(cellList[c1id],cellList[c2id]);
				//collisionpairList[pid]->set_adhesion(adhesion_pair_class[iptype]);
				//collisionpairList[pid]->set_friction(friction_pair_class[iptype]);
				double adhesion_force = adhesion_pair_array[cellList[c1id]->get_cell_type()][cellList[c2id]->get_cell_type()];
				collisionpairList[pid]->set_adhesion(adhesion_force);
				collisionpairList[pid]->set_friction(adhesion_force);
				if      (reid==0) {collisionpairList[pid]->set_redundant(false);}
				else if (reid==1) {collisionpairList[pid]->set_redundant(true);}	
				int an = (int)array1.size();
				if (an==16)
				{
					int s1id = atoi(array1[5].c_str());
					int t1id = atoi(array1[6].c_str());
					int s2id = atoi(array1[7].c_str());
					int t2id = atoi(array1[8].c_str());
					double fs1 = atof(array1[9].c_str());
					double ft1 = atof(array1[10].c_str());
					double fs2 = atof(array1[11].c_str());
					double ft2 = atof(array1[12].c_str());
					double ffs1 = atof(array1[13].c_str());
					double fft1 = atof(array1[14].c_str());
					int ffd = atoi(array1[15].c_str());
					mergepair *Mp = new mergepair(0,vertexList[s1id],vertexList[t1id],vertexList[s2id],vertexList[t2id]);
					collisionpairList[pid]->push_MP(Mp);
					Mp->set_F_s1(fs1);
					Mp->set_F_t1(ft1);
					Mp->set_F_s2(fs2);
					Mp->set_F_t2(ft2);
					Mp->set_Fr_s1(ffs1);
					Mp->set_Fr_t1(fft1);
					Mp->set_Fr_dir(ffd);
					point3D *s1 = vertexList[s1id];
					point3D *t1 = vertexList[t1id];
					point3D *s2 = vertexList[s2id];
					point3D *t2 = vertexList[t2id];
				}
				if (reid==0)
				{
					cellList[c1id]->pushpair(collisionpairList[pid]);
					cellList[c2id]->pushpair(collisionpairList[pid]);
				}
			}
			else if (!strcmp(strtmp_END_PAIR_LIST.c_str(),"END_PAIR_LIST"))
			{pair_list_flag = false;}
		}
		if (node_list_flag && point_based_flag)
		{
			string strtmp_END_NODE_LIST = ln.substr(0,13);
			if (strcmp(strtmp_END_NODE_LIST.c_str(),"END_NODE_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				int nid = atoi(array1[0].c_str());
				int in_global = atoi(array1[1].c_str());
				node *Nd = new node(nid);
				nodeList.push_back(Nd);
				if     (in_global==0) {nodeList[nid]->set_in_global(false);}
				else if(in_global==1) {nodeList[nid]->set_in_global(true);}
				if (in_global)
				{
					int vpn = (int)array1.size();
					vector<point3D*> vps;
					for (int ii=2;ii<vpn;ii++)
					{
						int vid = atoi(array1[ii].c_str());
						nodeList[nid]->pushvertex(vertexList[vid]);
						vps.push_back(vertexList[vid]);
						vertexList[vid]->set_node_id(nodeList[nid]);
					}
					int vpsn = (int)vps.size();
					if (vpsn>1)
					{
						for (int ii=0;ii<vpsn-1;ii++)
						{
							point3D *vpii = vps[ii];
							for (int jj=ii+1;jj<vpsn;jj++)
							{
								point3D *vpjj = vps[jj];
								vpii->pushmp(vpjj);
								vpjj->pushmp(vpii);
							}
						}
					}
					vps.clear();
				}
			}
			else if (!strcmp(strtmp_END_NODE_LIST.c_str(),"END_NODE_LIST"))
			{node_list_flag = false;}
		}
		if (dead_vert_list_flag && point_based_flag)
		{
			string strtmp_END_DEAD_VERT_LIST = ln.substr(0,18);
			if (strcmp(strtmp_END_DEAD_VERT_LIST.c_str(),"END_DEAD_VERT_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				int dvid = atoi(array1[0].c_str());
				int fpid = atoi(array1[1].c_str());
				int rpid = atoi(array1[2].c_str());
				int abid = atoi(array1[3].c_str());
				int baid = atoi(array1[4].c_str());
				vertexList[dvid]->setfp(vertexList[fpid]);
				vertexList[dvid]->setrp(vertexList[rpid]);
				vertexList[dvid]->setab(edgeList[abid]);
				vertexList[dvid]->setba(edgeList[baid]);
			}
			else if (!strcmp(strtmp_END_DEAD_VERT_LIST.c_str(),"END_DEAD_VERT_LIST"))
			{dead_vert_list_flag = false;}
		}
		if (edge_pair_list_flag && point_based_flag)
		{
			string strtmp_END_EDGE_PAIR_LIST = ln.substr(0,18);
			if (strcmp(strtmp_END_EDGE_PAIR_LIST.c_str(),"END_EDGE_PAIR_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,",");
				int eid = atoi(array1[0].c_str());
				int enid = atoi(array1[1].c_str());
				int pid = atoi(array1[2].c_str());
				if (enid>=0) 
				{
					edgeList[eid]->set_attach(1);
					edgeList[eid]->set_Neighbor(edgeList[enid]);
					edgeList[eid]->set_Pair(collisionpairList[pid]);
				}
				else
				{
					edgeList[eid]->set_attach(0);
				}
			}
			else if (!strcmp(strtmp_END_EDGE_PAIR_LIST.c_str(),"END_EDGE_PAIR_LIST"))
			{edge_pair_list_flag = false;}
		}
		if (migr_cell_list_flag && point_based_flag)
		{
			string strtmp_END_MIGR_CELL_LIST = ln.substr(0,18);
			if (strcmp(strtmp_END_MIGR_CELL_LIST.c_str(),"END_MIGR_CELL_LIST"))
			{
				vector<string> array1;
				str_split(ln,array1,";");
				vector<string> array2;
				str_split(array1[0],array2,",");
				int cid = atoi(array2[0].c_str());
				int pid = atoi(array2[1].c_str());
				if (pid>0)
				{
					for (int i=0;i<pid;i++)
					{
						vector<string> array3;
						str_split(array1[i+1],array3,",");
						int cfid = atoi(array3[0].c_str());
						int initial_c = atoi(array3[1].c_str());
						double fccx = atof(array3[2].c_str());
						double fccy = atof(array3[3].c_str());
						double angle = atof(array3[4].c_str());
						int fcn = atoi(array3[5].c_str());
						rear_fc *FC = new rear_fc(i);
						FC->set_C(cellList[cfid]);
						if      (initial_c==0) {FC->set_initial(false);}
						else if (initial_c==1) {FC->set_initial(true);}	
						FC->set_center(fccx,fccy);
						FC->set_angle(angle);
						for (int j=0;j<fcn;j++)
						{
							int cncid = atoi(array3[6+j].c_str());
							FC->push_NCs(cellList[cncid]);
						}
						cellList[cid]->push_filling_cell(FC);
					}
				}
			}
			else if (!strcmp(strtmp_END_MIGR_CELL_LIST.c_str(),"END_MIGR_CELL_LIST"))
			{migr_cell_list_flag = false;}
		}
		string strtmp_center_based = ln.substr(0,12);
		if (!strcmp(strtmp_center_based.c_str(),"Center-based"))
		{center_based_flag = true;}
		string strtmp_point_based = ln.substr(0,11);
		if (!strcmp(strtmp_point_based.c_str(),"Point-based"))
		{point_based_flag = true;}
		string strtmp_u = ln.substr(0,1);
		if (!strcmp(strtmp_u.c_str(),"u")) 
		{
			vector<string> array1;
			vector<string> array2;
			str_split(ln,array1," ");
			str_split(array1[1],array2,",");
			double umin = atof(array2[0].c_str());
			double umax = atof(array2[1].c_str());
			vv->setUMin(umin);
			vv->setUMax(umax);
			vv->setUMin_default(umin);
			vv->setUMax_default(umax);
		}
		string strtmp_v = ln.substr(0,1);
		if (!strcmp(strtmp_v.c_str(),"v")) 
		{
			vector<string> array1;
			vector<string> array2;
			str_split(ln,array1," ");
			str_split(array1[1],array2,",");
			double vmin = atof(array2[0].c_str());
			double vmax = atof(array2[1].c_str());
			vv->setVMin(vmin);
			vv->setVMax(vmax);
			vv->setVMin_default(vmin);
			vv->setVMax_default(vmax);
		}
		string strtmp_VPN = ln.substr(0,3);
		if (!strcmp(strtmp_VPN.c_str(),"VPN")) 
		{
			vector<string> array1;
			vector<string> array2;
			str_split(ln,array1," ");
			str_split(array1[1],array2,",");
			double vpn1 = atof(array2[0].c_str());
			double vpn2 = atof(array2[1].c_str());
			point2D *A=new point2D(vpn1,vpn2);
			vv->setVPN(A);
		}
		string strtmp_VUP = ln.substr(0,3);
		if (!strcmp(strtmp_VUP.c_str(),"VUP")) 
		{
			vector<string> array1;
			vector<string> array2;
			str_split(ln,array1," ");
			str_split(array1[1],array2,",");
			double vup1 = atof(array2[0].c_str());
			double vup2 = atof(array2[1].c_str());
			point2D *A=new point2D(vup1,vup2);
			vv->setVUP(A);
		}
		string strtmp_VRP = ln.substr(0,3);
		if (!strcmp(strtmp_VRP.c_str(),"VRP")) 
		{
			vector<string> array1;
			vector<string> array2;
			str_split(ln,array1," ");
			str_split(array1[1],array2,",");
			double vrp1 = atof(array2[0].c_str());
			double vrp2 = atof(array2[1].c_str());
			point2D *A=new point2D(vrp1,vrp2);
			vv->setVRP(A);
		}
		string strtmp_PRP = ln.substr(0,3);
		if (!strcmp(strtmp_PRP.c_str(),"PRP")) 
		{
			vector<string> array1;
			vector<string> array2;
			str_split(ln,array1," ");
			str_split(array1[1],array2,",");
			double prp1 = atof(array2[0].c_str());
			double prp2 = atof(array2[1].c_str());
			point2D *A=new point2D(prp1,prp2);
			vv->setPRP(A);
		}
		//
		string strtmp_VERT_LIST = ln.substr(0,9);
		if (!strcmp(strtmp_VERT_LIST.c_str(),"VERT_LIST")) {vert_list_flag = true;}
		string strtmp_EDGE_LIST = ln.substr(0,9);
		if (!strcmp(strtmp_EDGE_LIST.c_str(),"EDGE_LIST")) {edge_list_flag = true;}
		string strtmp_CELL_LIST = ln.substr(0,9);
		if (!strcmp(strtmp_CELL_LIST.c_str(),"CELL_LIST")) {cell_list_flag = true;}
		string strtmp_PAIR_LIST = ln.substr(0,9);
		if (!strcmp(strtmp_PAIR_LIST.c_str(),"PAIR_LIST")) {pair_list_flag = true;}
		string strtmp_NODE_LIST = ln.substr(0,9);
		if (!strcmp(strtmp_NODE_LIST.c_str(),"NODE_LIST")) {node_list_flag = true;}
		string strtmp_DEAD_VERT_LIST = ln.substr(0,14);
		if (!strcmp(strtmp_DEAD_VERT_LIST.c_str(),"DEAD_VERT_LIST")) {dead_vert_list_flag = true;}
		string strtmp_EDGE_PAIR_LIST = ln.substr(0,14);
		if (!strcmp(strtmp_EDGE_PAIR_LIST.c_str(),"EDGE_PAIR_LIST")) {edge_pair_list_flag = true;}
		string strtmp_MIGR_CELL_LIST = ln.substr(0,14);
		if (!strcmp(strtmp_MIGR_CELL_LIST.c_str(),"MIGR_CELL_LIST")) {migr_cell_list_flag = true;}
	}
	fileIn.close();

	range[0] -= 1000;
	range[1] -= 1000;
	range[2] += 1000;
	range[3] += 1000;

	/************************************
	 Delaunay triangulaiton construction 
	 Collision detection
	 Mesh generation
	************************************/
	if (center_based_flag)
	{
		delaunay_triangulation(range);
		detect_collision(number);
		cell_mesh_generation(number);
		node_test();
		Complex->~CComplex();
	}
	if (point_based_flag)
	{
		cell_mesh_generation(number);
		node_test();
		initial = 1;
	}
	int cell_n = (int)cellList.size();
	for (int i=0;i<cell_n;i++)
	{
		cellList[i]->refresh_AABB();
		cellList[i]->set_area();
		cellList[i]->set_initial_area(AR/2);
		for (int j=0;j<10;j++)
		{
			cellList[i]->set_area_p(j,cellList[i]->get_area());
		}
	}
}

void dbReader::node_test()
{
	int n_size = 0;
	int node_v_size[10];
	for (int i=0;i<10;i++) {node_v_size[i] = 0;}
	for (int i=0;i<(int)cellList.size();i++)
	{
		cellList[i]->center_refresh();
	}
	for (int i=0;i<(int)nodeList.size();i++)
	{
		if (nodeList[i]->get_in_global()) 
		{
			n_size++;
			int vpn = (int)nodeList[i]->getvertex()->size();
			if      (vpn==1) {node_v_size[0]++;}
			else if (vpn==2) {node_v_size[1]++;}
			else if (vpn==3) {node_v_size[2]++;}
			else if (vpn==4) {node_v_size[3]++;}
			else if (vpn==5) {node_v_size[4]++;}
			else if (vpn==6) {node_v_size[5]++;}
			else if (vpn==7) {node_v_size[6]++;}
			else if (vpn==8) {node_v_size[7]++;}
			else if (vpn==9) {node_v_size[8]++;}
		}
	}
	int n_size_all = (int)vertexList.size();
	int n_size_r[10];
	int n_remove_p = 0;
	for (int i=0;i<10;i++) {n_size_r[i] = 0;}
	for (int i=0;i<(int)vertexList.size();i++)
	{
		if (vertexList[i]->get_in_cell()) 
		{
			int nei_n = (int)vertexList[i]->getmp()->size();
			n_size_r[nei_n]++;
		}
		else {n_remove_p++;}
	}
	n_size_all -= n_remove_p;
	for (int i=0;i<10;i++)
	{
		n_size_all -= i*n_size_r[i]/(i+1);
	}
	cout<<"The node number: "<<n_size<<" "<<n_size_all<<endl;
}

void dbReader::print_matlab_file()
{
	/***********************************
	    The points only for cell 0
	***********************************/
	string mname = "fem1.m";
	ofstream OUT (mname.c_str());
	vector<double> xc;
	vector<double> yc;
	vector< vector <int> > rel;
	int sm_id = 0;
	for (int i=0;i<(int)cellList[0]->get_sides()->size();i++)
	{
		xc.push_back(cellList[0]->get_sides()->at(i)->p1()->x());
		yc.push_back(cellList[0]->get_sides()->at(i)->p1()->y());
		cellList[0]->get_sides()->at(i)->p1()->get_node_id()->set_SM_ID(sm_id);
		sm_id++;
	}
	for (int i=0;i<(int)cellList[0]->get_inner_p()->size();i++)
	{
		cellList[0]->get_inner_p()->at(i)->setid(sm_id);
		xc.push_back(cellList[0]->get_inner_p()->at(i)->x());
		yc.push_back(cellList[0]->get_inner_p()->at(i)->y());
		sm_id++;
	}
	for (int i=0;i<(int)cellList[0]->get_sides()->size();i++)
	{
		double x1 = cellList[0]->get_sides()->at(i)->p1()->x();
		double y1 = cellList[0]->get_sides()->at(i)->p1()->y();
		double x2 = cellList[0]->get_sides()->at(i)->p2()->x();
		double y2 = cellList[0]->get_sides()->at(i)->p2()->y();
		cellList[0]->get_sides()->at(i)->setnormal();
		double n_x = cellList[0]->get_sides()->at(i)->getnormal(0);
		double n_y = cellList[0]->get_sides()->at(i)->getnormal(1);
		double d = sqrt(n_x*n_x + n_y*n_y);
		double len = cellList[0]->get_sides()->at(i)->getlength();
		double xm = (x1+x2)/2 + n_x/d*len*0.866;
		double ym = (y1+y2)/2 + n_y/d*len*0.866;
		xc.push_back(xm);
		yc.push_back(ym);
	}
	for (int i=0;i<(int)cellList[0]->get_inner_t()->size();i++)
	{
		point3D* A = cellList[0]->get_inner_t()->at(i)->getA();
		point3D* B = cellList[0]->get_inner_t()->at(i)->getB();
		point3D* C = cellList[0]->get_inner_t()->at(i)->getC();
		int Aid = 0;
		int Bid = 0;
		int Cid = 0;
		if (A->get_cell_boundary()) {Aid = A->get_node_id()->SM_ID();}
		else {Aid = A->id();}
		if (B->get_cell_boundary()) {Bid = B->get_node_id()->SM_ID();}
		else {Bid = B->id();}
		if (C->get_cell_boundary()) {Cid = C->get_node_id()->SM_ID();}
		else {Cid = C->id();}
		vector<int> relrow;
		Aid += 1;Bid += 1;Cid += 1;
		relrow.push_back(Aid);
		relrow.push_back(Bid);
		relrow.push_back(Cid);
		rel.push_back(relrow);
	}
	for (int i=0;i<(int)cellList[0]->get_sides()->size();i++)
	{
		int i1 = i+1;
		int i2 = i+2;if (i==(int)cellList[0]->get_sides()->size()-1) {i2 = 1;}
		int i3 = sm_id + 1;
		sm_id++;
		vector<int> relrow;
		relrow.push_back(i1);
		relrow.push_back(i2);
		relrow.push_back(i3);
		rel.push_back(relrow);
	}
	OUT<<"x=[";
	for (int i=0;i<(int)xc.size();i++) {OUT<<xc[i]<<" ";}
	OUT<<"];\n";
	OUT<<"y=[";
	for (int i=0;i<(int)yc.size();i++) {OUT<<yc[i]<<" ";}
	OUT<<"];\n";
	OUT<<"rel = [";
	for (int i=0;i<(int)rel.size();i++)
	{
		for (int j=0;j<(int)rel[i].size();j++) {OUT<<rel[i][j]<<" ";}
		OUT<<";\n";
	}
	OUT<<"];\n";
	OUT<<"reln = size(rel);\n";
	for (int i=0;i<(int)rel.size();i++)
	{
		int id1 = rel[i][0];
		int id2 = rel[i][1];
		int id3 = rel[i][2];
		OUT<<"b"<<i<<" = [1 1 1;\n";
		OUT<<"x("<<id1<<") x("<<id2<<") x("<<id3<<");\n";
		OUT<<"y("<<id1<<") y("<<id2<<") y("<<id3<<")];\n";
	}
	for (int i=0;i<(int)rel.size();i++)
	{
		OUT<<"det"<<i<<" = det(b"<<i<<");if (det"<<i<<"<0) det"<<i<<" = det"<<i<<"*-1;end\n";
	}
	for (int i=0;i<(int)rel.size();i++)
	{
		OUT<<"B"<<i<<" = det"<<i<<"*inv(b"<<i<<");\n";
	}
	for (int i=0;i<(int)rel.size();i++)
	{
		OUT<<"Be"<<i<<" = [B"<<i<<"(1,2) 0 B"<<i<<"(2,2) 0 B"<<i<<"(3,2) 0;\n";
		OUT<<"0 B"<<i<<"(1,3) 0 B"<<i<<"(2,3) 0 B"<<i<<"(3,3);\n";
		OUT<<"B"<<i<<"(1,3) B"<<i<<"(1,2) B"<<i<<"(2,3) B"<<i<<"(2,2) B"<<i<<"(3,3) B"<<i<<"(3,2)];\n";
	}
	for (int i=0;i<(int)rel.size();i++)
	{
		OUT<<"Be"<<i<<" = Be"<<i<<"/det"<<i<<";\n";
	}
	OUT<<"mu = 1.91;lambda = 2.87;mu1 = 0.034;lambda1 = 1.63;\n";
	OUT<<"D = [lambda+2*mu lambda 0;\nlambda lambda+2*mu 0;\n0 0 mu];\n";
	OUT<<"D1 = [lambda1+2*mu1 lambda1 0;\nlambda1 lambda1+2*mu1 0;\n0 0 mu1];\n";
	for (int i=0;i<(int)rel.size();i++)
	{
		if (i<(int)cellList[0]->get_inner_t()->size()) 
		{OUT<<"Ke"<<i<<" = Be"<<i<<"'*D*Be"<<i<<"*det"<<i<<"/2;\n";}
		else
		{OUT<<"Ke"<<i<<" = Be"<<i<<"'*D1*Be"<<i<<"*det"<<i<<"/2;\n";}
	}
	int n_t = (int)cellList[0]->get_sides()->size()*2 + (int)cellList[0]->get_inner_p()->size();
	n_t = n_t*2;
	OUT<<"K = zeros("<<n_t<<","<<n_t<<");\n";
	OUT<<"for i=1:reln(1,1)\np1 = rel(i,1);\np2 = rel(i,2);\np3 = rel(i,3);\nKei = Ke0;\n";
	OUT<<"if i==1\nKei = Ke0;\n";
	for (int i=1;i<(int)rel.size();i++)
	{
		OUT<<"else if i=="<<i+1<<"\n";
		OUT<<"Kei = Ke"<<i<<";\n";
	}
	for (int i=0;i<(int)rel.size();i++)
	{
		OUT<<"end\n";
	}
	OUT<<"i1 = p1*2-1;\nj1 = p1*2;\ni2 = p2*2-1;\nj2 = p2*2;\ni3 = p3*2-1;\nj3 = p3*2;\n";
	OUT<<"K(i1,i1) = K(i1,i1) + Kei(1,1);\n";
	OUT<<"K(i1,j1) = K(i1,j1) + Kei(1,2);\n";
	OUT<<"K(j1,i1) = K(j1,i1) + Kei(2,1);\n";
	OUT<<"K(j1,j1) = K(j1,j1) + Kei(2,2);\n";
	OUT<<"K(i1,i2) = K(i1,i2) + Kei(1,3);\n";
	OUT<<"K(i1,j2) = K(i1,j2) + Kei(1,4);\n";
	OUT<<"K(j1,i2) = K(j1,i2) + Kei(2,3);\n";
	OUT<<"K(j1,j2) = K(j1,j2) + Kei(2,4);\n";
	OUT<<"K(i1,i3) = K(i1,i3) + Kei(1,5);\n";
	OUT<<"K(i1,j3) = K(i1,j3) + Kei(1,6);\n";
	OUT<<"K(j1,i3) = K(j1,i3) + Kei(2,5);\n";
	OUT<<"K(j1,j3) = K(j1,j3) + Kei(2,6);\n";
	OUT<<"K(i2,i1) = K(i2,i1) + Kei(3,1);\n";
	OUT<<"K(i2,j1) = K(i2,j1) + Kei(3,2);\n";
	OUT<<"K(j2,i1) = K(j2,i1) + Kei(4,1);\n";
	OUT<<"K(j2,j1) = K(j2,j1) + Kei(4,2);\n";
	OUT<<"K(i2,i2) = K(i2,i2) + Kei(3,3);\n";
	OUT<<"K(i2,j2) = K(i2,j2) + Kei(3,4);\n";
	OUT<<"K(j2,i2) = K(j2,i2) + Kei(4,3);\n";
	OUT<<"K(j2,j2) = K(j2,j2) + Kei(4,4);\n";
	OUT<<"K(i2,i3) = K(i2,i3) + Kei(3,5);\n";
	OUT<<"K(i2,j3) = K(i2,j3) + Kei(3,6);\n";
	OUT<<"K(j2,i3) = K(j2,i3) + Kei(4,5);\n";
	OUT<<"K(j2,j3) = K(j2,j3) + Kei(4,6);\n";
	OUT<<"K(i3,i1) = K(i3,i1) + Kei(5,1);\n";
	OUT<<"K(i3,j1) = K(i3,j1) + Kei(5,2);\n";
	OUT<<"K(j3,i1) = K(j3,i1) + Kei(6,1);\n";
	OUT<<"K(j3,j1) = K(j3,j1) + Kei(6,2);\n";
	OUT<<"K(i3,i2) = K(i3,i2) + Kei(5,3);\n";
	OUT<<"K(i3,j2) = K(i3,j2) + Kei(5,4);\n";
	OUT<<"K(j3,i2) = K(j3,i2) + Kei(6,3);\n";
	OUT<<"K(j3,j2) = K(j3,j2) + Kei(6,4);\n";
	OUT<<"K(i3,i3) = K(i3,i3) + Kei(5,5);\n";
	OUT<<"K(i3,j3) = K(i3,j3) + Kei(5,6);\n";
	OUT<<"K(j3,i3) = K(j3,i3) + Kei(6,5);\n";
	OUT<<"K(j3,j3) = K(j3,j3) + Kei(6,6);\n";
	OUT<<"end\n";
	OUT<<"eig(K)\n";
	OUT<<"for i=1:"<<n_t<<"\n";
	for (int i=0;i<(int)cellList[0]->get_sides()->size();i++)
	{
		int i1 = n_t - 2*i;
		int i2 = n_t - 2*i - 1;
		OUT<<"K("<<i1<<",i)=0;K(i,"<<i1<<")=0;\n";
		OUT<<"K("<<i2<<",i)=0;K(i,"<<i2<<")=0;\n";
	}
	OUT<<"end\n";
	for (int i=0;i<(int)cellList[0]->get_sides()->size();i++)
	{
		int i1 = n_t - 2*i;
		int i2 = n_t - 2*i - 1;
		OUT<<"K("<<i1<<","<<i1<<")=1;\n";
		OUT<<"K("<<i2<<","<<i2<<")=1;\n";
	}
	cellList[0]->center_refresh();
	OUT<<"F=zeros("<<n_t<<",1);\n";
	for (int i=0;i<(int)cellList[0]->get_sides()->size();i++)
	{
		double x1 = cellList[0]->get_sides()->at(i)->p1()->x() - cellList[0]->get_center(0);
		double y1 = cellList[0]->get_sides()->at(i)->p1()->y() - cellList[0]->get_center(1);
		double d = sqrt(x1*x1 + y1*y1);
		OUT<<"F("<<i*2+1<<","<<1<<")="<<x1/d<<";\n";
		OUT<<"F("<<i*2+2<<","<<1<<")="<<y1/d<<";\n";
	}
	OUT<<"u = inv(K)*F;\n";
	OUT<<"x1 = x;\ny1 = y;\n";
	OUT<<"x2=zeros("<<n_t/2-cellList[0]->get_sides()->size()<<",1);\n";
	OUT<<"y2=zeros("<<n_t/2-cellList[0]->get_sides()->size()<<",1);\n";
	OUT<<"for i=1:"<<n_t/2-cellList[0]->get_sides()->size()<<"\n";
	OUT<<"x1(i) = x(i)+u(2*i-1);\n";
	OUT<<"y1(i) = y(i)+u(2*i);\n";
	OUT<<"x2(i) = x(i)+u(2*i-1);\n";
	OUT<<"y2(i) = y(i)+u(2*i);\n";
	OUT<<"end\nplot(x2,y2,'o');\n";
	OUT.close();
}

void dbReader::print_cell_file(int number) // INTERFACE: output format
{
	time_t t_b_c, t_e_c;
	t_b_c = clock();
	int vn = (int)vertexList.size();
	int en = (int)edgeList.size();
	int cn = (int)cellList.size();
	int pn = (int)collisionpairList.size();
	int nn = (int)nodeList.size();
	ostringstream temp;
	temp<<number;
	string tmps = temp.str();
	string cname = "printout"+tmps+".txt";
	/* output.dir also holds the per-step cell dumps when it is set. */
	{
		string od = dycfg::getStr("output.dir", "");
		if (!od.empty() && od != ".") cname = od + "/" + cname;
	}
	ofstream OUT (cname.c_str());
	if (!OUT.good())
	{
		cerr<<"ERROR: cannot write "<<cname<<" - does the output directory exist "
		    <<"and is it writable?"<<endl;
		return;
	}
	/* output.sections:
	     full   - VERT/EDGE/NODE/MIGR lists included. The file can be fed back in
	              as a restart input. ~15 MB/step for a 7000-cell tissue.
	     simple - CELL_LIST and PAIR_LIST only, which is everything the analysis
	              scripts read. ~2 MB/step. This is what the archived
	              branch_dsq9_simpleoutput tree hardcoded, and what produced the
	              published results.
	   Both write identical CELL_LIST and PAIR_LIST blocks. */
	bool full = (dycfg::getStr("output.sections", "full") != "simple");
	if (full) {
	OUT<<"#######################\n# Cell Growth Project #\n# Sample Database     #\n#######################"<<endl;
	OUT<<"Point-based\n"<<endl;
	OUT<<"# This is the bounds of your viewing window plane in VRC coordinates"<<endl;
	OUT<<"u -100.0,100.0\nv -75.0,75.0"<<endl;
	OUT<<"# VRC specifying vectors\nVPN 0.0,1.0\nVUP 1.0,0.0"<<endl;
	OUT<<"# World coordinate point specifying the VRC origin\nVRP 200.0,100.0"<<endl;
	OUT<<"# Specify the center of project in reference to the VRC\nPRP 0.0,0.0"<<endl;
	OUT<<"# List of vertices of objects in world coordinates"<<endl;
	OUT<<"VERT_LIST"<<endl;
	for (int i=0;i<vn;i++)
	{
		double x = vertexList[i]->x();
		double y = vertexList[i]->y();
		int index = vertexList[i]->id();
		int in_cell = vertexList[i]->get_in_cell();
		double store_force = vertexList[i]->get_stored_apart_force();
		double store_fr_force = vertexList[i]->get_stored_friction_force();
		double elastic_x = vertexList[i]->get_elastic_force(0);
		double elastic_y = vertexList[i]->get_elastic_force(1);
		double break_f_x = vertexList[i]->get_adhesion_break(0);
		double break_f_y = vertexList[i]->get_adhesion_break(1);
		OUT<<x<<","<<y<<","<<index<<","<<in_cell<<","<<store_force<<","<<store_fr_force<<","<<elastic_x<<","<<elastic_y<<","<<break_f_x<<","<<break_f_y<<endl;
	}
	OUT<<"END_VERT_LIST"<<endl;
	OUT<<"EDGE_LIST"<<endl;
	for (int i=0;i<en;i++)
	{
		int index = edgeList[i]->id();
		int p1_id = edgeList[i]->p1()->id();
		int p2_id = edgeList[i]->p2()->id();
		OUT<<index<<","<<p1_id<<","<<p2_id<<endl;
	}
	OUT<<"END_EDGE_LIST"<<endl;
	} /* end if (full) */
	OUT<<"#ID,Dead,S_Dead,Mig,mAng,Type,X,Y,Area,N_Species,Species0,Species1,..."<<endl;
	OUT<<"CELL_LIST"<<endl;
	for (int i=0;i<cn;i++)
	{
		int index = cellList[i]->id();
		int dead = cellList[i]->Dead();
		int set_dead = cellList[i]->Set_Dead();
		int migrate = cellList[i]->Mark_number(); // 0: static; 1: migrating
		//int soften = cellList[i]->Soften();
		double mAngle = cellList[i]->get_migrate_angle();
		int type = cellList[i]->get_cell_type();
		double cx = cellList[i]->get_center(0);
		double cy = cellList[i]->get_center(1);
		double iarea = cellList[i]->get_area();
		OUT<<index<<","<<dead<<","<<set_dead<<","<<migrate<<","<<mAngle/*soften*/<<","<<type<<","<<cx<<","<<cy<<","<<iarea;
		// YFC: ADD OUTPUT FOR CELL STATE
		int maxlen = max(cellList[i]->get_cell_statelen(), biology->getNumSpecies());
		int len = cellList[i]->get_cell_statelen();
		OUT<<","<<maxlen;
		for (int j=0;j<maxlen;j++)
		{
			if (j<len)
			{
				OUT<<","<<cellList[i]->get_cell_state(j);
			}
			else
			{
				OUT<<",0";
			}
		}
		/*for (int j=0;j<(int)cellList[i]->get_sides()->size();j++)
		{
			int eid = cellList[i]->get_sides()->at(j)->id();
			OUT<<","<<eid;
		}*/
		OUT<<endl;
	}
	OUT<<"END_CELL_LIST"<<endl;
	OUT<<"PAIR_LIST"<<endl;
	for (int i=0;i<pn;i++)
	{
		int index = collisionpairList[i]->id();
		int redundant = collisionpairList[i]->get_redundant();
		int c1id = collisionpairList[i]->get_I1()->id();
		int c2id = collisionpairList[i]->get_I2()->id();
		int mpn = (int)collisionpairList[i]->MP()->size();
		OUT<<index<<","<<redundant<<","<<c1id<<","<<c2id<<","<<mpn;
		if (mpn>0)
		{
			int s1id = collisionpairList[i]->MP()->at(0)->S1()->id();
			int t1id = collisionpairList[i]->MP()->at(0)->T1()->id();
			int s2id = collisionpairList[i]->MP()->at(0)->S2()->id();
			int t2id = collisionpairList[i]->MP()->at(0)->T2()->id();
			double fs1 = collisionpairList[i]->MP()->at(0)->get_F_s1();
			double ft1 = collisionpairList[i]->MP()->at(0)->get_F_t1();
			double fs2 = collisionpairList[i]->MP()->at(0)->get_F_s2();
			double ft2 = collisionpairList[i]->MP()->at(0)->get_F_t2();
			double frs1 = collisionpairList[i]->MP()->at(0)->get_Fr_s1();
			double frt1 = collisionpairList[i]->MP()->at(0)->get_Fr_t1();
			int frd = collisionpairList[i]->MP()->at(0)->get_Fr_dir();
			OUT<<","<<s1id<<","<<t1id<<","<<s2id<<","<<t2id<<","<<fs1<<","<<ft1<<","<<fs2<<","<<ft2<<","<<frs1<<","<<frt1<<","<<frd;
		}
		OUT<<endl;
	}
	OUT<<"END_PAIR_LIST"<<endl;
	if (full) {
	OUT<<"NODE_LIST"<<endl;
	for (int i=0;i<nn;i++)
	{
		int index = nodeList[i]->id();
		int in_global = nodeList[i]->get_in_global();
		OUT<<index<<","<<in_global;
		int vpn = (int)nodeList[i]->getvertex()->size();
		for (int j=0;j<vpn;j++)
		{
			int vpid = nodeList[i]->getvertex()->at(j)->id();
			OUT<<","<<vpid;
		}
		OUT<<endl;
	}
	OUT<<"END_NODE_LIST"<<endl;
	OUT<<"DEAD_VERT_LIST"<<endl;
	for (int i=0;i<vn;i++)
	{
		if (vertexList[i]->get_in_cell()) continue;
		int fpid = vertexList[i]->getfp()->id();
		int rpid = vertexList[i]->getrp()->id();
		int abid = vertexList[i]->getab()->id();
		int baid = vertexList[i]->getba()->id();
		int dvid = vertexList[i]->id();
		OUT<<dvid<<","<<fpid<<","<<rpid<<","<<abid<<","<<baid;
		OUT<<endl;
	}
	OUT<<"END_DEAD_VERT_LIST"<<endl;
	OUT<<"EDGE_PAIR_LIST"<<endl;
	for (int i=0;i<en;i++)
	{
		if (edgeList[i]->get_attach())
		{
			int eid = edgeList[i]->id();
			int enid = edgeList[i]->get_Neighbor()->id();
			int pid = edgeList[i]->get_Pair()->id();
			OUT<<eid<<","<<enid<<","<<pid<<endl;
		}
		else
		{
			int eid = edgeList[i]->id();
			OUT<<eid<<",-1,-1"<<endl;
		}
	}
	OUT<<"END_EDGE_PAIR_LIST"<<endl;
	OUT<<"MIGR_CELL_LIST"<<endl;
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Migrate())
		{
			int cid = cellList[i]->id();
			int pid = cellList[i]->get_migrate_filling_cells()->size();
			OUT<<cid<<","<<pid<<";";
			if (pid>0)
			{
				for (int j=0;j<pid;j++)
				{
					rear_fc *FC = cellList[cid]->get_migrate_filling_cells()->at(j);
					int cfid = FC->getC()->id();
					int initial_c = FC->get_initial();
					double fccx = FC->get_center(0);
					double fccy = FC->get_center(1);
					double angle = FC->get_angle();
					int fcn = (int)FC->get_NCs()->size();
					OUT<<cfid<<","<<initial_c<<","<<fccx<<","<<fccy<<","<<angle<<","<<fcn;
					for (int k=0;k<fcn;k++)
					{
						int cncid = FC->get_NCs()->at(k)->id();
						OUT<<","<<cncid;
					}
					OUT<<";";
				}
				OUT<<endl;
			}
		}
	}
	OUT<<"END_MIGR_CELL_LIST"<<endl;
	} /* end if (full) */
	OUT.close();
	t_e_c = clock();
	cout<<"Print out data file done: "<<(t_e_c-t_b_c)*0.001<<" seconds!"<<endl;
}

#endif
