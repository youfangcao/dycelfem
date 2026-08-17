// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     para.cpp            ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef PARA_CPP
#define PARA_CPP

#include <cmath>
#include "para.h"

using namespace std;

map<int,double> pressure_edge_class;
map<int,double> tension_edge_class;
map<int,double> mass_edge_class;
map<int,double> adhesion_pair_class;
double adhesion_pair_array[20][20];
map<int,double> friction_pair_class;
map<int,double> growth_rate_class;
map<int,double> growth_scale_factor_class;
map<int,double> color_class_R;
map<int,double> color_class_G;
map<int,double> color_class_B;
map<int,double> elastic_MU;
map<int,double> elastic_LAMBDA;
map<int,int>    fixed_cells;

void construct_map_of_para()  // initial value
{
	pressure_edge_class[0] = 4.0*tan(PI/SN);
	pressure_edge_class[1] = 2.0*tan(PI/SN);
	pressure_edge_class[2] = 2.0*tan(PI/SN);
	pressure_edge_class[3] = 2.0*tan(PI/SN);
	pressure_edge_class[4] = 2.0*tan(PI/SN);
	pressure_edge_class[5] = 2.0*tan(PI/SN);
	pressure_edge_class[6] = 2.0*tan(PI/SN);
	///////////////////////////
	tension_edge_class[0] = 2.0;
	tension_edge_class[1] = 1.0;
	tension_edge_class[2] = 1.0;
	tension_edge_class[3] = 1.0;
	tension_edge_class[4] = 1.0;
	tension_edge_class[5] = 1.0;
	tension_edge_class[6] = 1.0;
	///////////////////////////
	for (int i=0;i<20;i++)
	{
		for (int j=0;j<20;j++)
		{
			adhesion_pair_array[i][j] = 0;
		}
	}
	adhesion_pair_array[0][0] = 20; // ECM-ECM    0-0
	adhesion_pair_array[0][1] = 20; // ECM-BM     0-1
	adhesion_pair_array[0][2] = 2.3; // ECM-f      0-2
	adhesion_pair_array[0][3] = 9.8; // ECM-k      0-3
	adhesion_pair_array[0][4] = 1; // ECM-FC      0-4
	adhesion_pair_array[1][1] = 20; // BM-BM      1-1
	adhesion_pair_array[1][2] = 2; // BM-f      1-2
	adhesion_pair_array[1][3] = 6; // BM-k      1-3
	adhesion_pair_array[1][4] = 1; // BM-FC      1-4
	adhesion_pair_array[2][2] = 0; // f-f      2-2
	adhesion_pair_array[2][3] = 0; // f-k      2-3
	adhesion_pair_array[2][4] = 2; // f-FC      2-4
	adhesion_pair_array[3][3] = 3.2; // k-k      3-3
	adhesion_pair_array[3][4] = 0; // k-FC      3-4
	adhesion_pair_array[4][4] = 2; // FC-FC      4-4
	for (int i=1;i<20;i++)
	{
		for (int j=0;j<i;j++)
		{
			adhesion_pair_array[i][j] = adhesion_pair_array[j][i];
		}
	}
	///////////////////////////
	mass_edge_class[0] = 0.01;
	mass_edge_class[1] = 0.01;
	mass_edge_class[2] = 0.01;
	mass_edge_class[3] = 0.01;
	mass_edge_class[4] = 0.01;
	mass_edge_class[5] = 0.01;
	mass_edge_class[6] = 0.01;
	///////////////////////////
	growth_rate_class[0] = 0;
	growth_rate_class[1] = 0;
	growth_rate_class[2] = 0;
	growth_rate_class[3] = 0;
	growth_rate_class[4] = 0;
	growth_rate_class[5] = 0;
	growth_rate_class[6] = 0;
	//////////////////////////////////
	growth_scale_factor_class[0] = 0.004; // KGF
	growth_scale_factor_class[1] = 0.004; // EGF
	//////////////////////////////////
}

#endif

