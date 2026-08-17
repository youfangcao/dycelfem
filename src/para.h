// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     para.h              ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef PARA_H
#define PARA_H

#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

const double PI = 3.14159265;    // PI value
const int SN = 20;               // initial size number of polygon
const double ANGLE = 360;        // 360 degrees
const double B_ANGLE = 60;       // for pick correction
const double CR = 10;            // initial radius: 10 micrometers(10E-6 meter)
const double TM = 0.1;           // time step
const double BD = 2*PI*CR/SN;    // edge size of initialized polygon
const double AR = 2*CR*CR*2.598; // double size of the hexagon inside the shpere
const double AR0 = CR*CR*2.598;	 // resting cell area
const double GAR = PI*CR*CR/10;  // incremental area
const double RDT = 24.0;         // cell cycle time: 24 hours
//const double MU = 1.91;
//const double LAMBDA = 2.87; // 5 and 0.3
const double MU = 6711.4;
const double LAMBDA = 328859.1; // E: 2E4 Pa; v: 0.49;
const double MU2_LAMBDA = 2*MU + LAMBDA;
const double MU1_LAMBDA = MU + LAMBDA;
const double MU_virtual = 0.033;
const double LAMBDA_virtual = 1.63; // 0.1 and 0.49
const double MU_soft = 0.033; 
const double LAMBDA_soft = 1.63; // 0.1 and 0.49
const double MU_virtual2_LAMBDA_virtual = 2*MU_virtual + LAMBDA_virtual;
const double MU_virtual1_LAMBDA_virtual = MU_virtual + LAMBDA_virtual;
const double MU_soft2_LAMBDA_soft = 2*MU_soft + LAMBDA_soft;
const double MU_soft1_LAMBDA_soft = MU_soft + LAMBDA_soft;
const double pressure_edge_soften = 0.0*tan(PI/SN);
const double tension_edge_soften = 0.0;
const double Boltzmann_Kelvin = 4.28; // product of Boltzmann constant and Kelvin temparature: 10E-21 Joule(kg meter^2/second^2)

extern map<int,double> pressure_edge_class;
extern map<int,double> tension_edge_class;
extern map<int,double> mass_edge_class;
extern double adhesion_pair_array[20][20];
extern map<int,double> adhesion_pair_class;
extern map<int,double> friction_pair_class;
extern map<int,double> growth_rate_class;
extern map<int,double> growth_scale_factor_class;
extern map<int,double> color_class_R;
extern map<int,double> color_class_G;
extern map<int,double> color_class_B;
extern map<int,double> elastic_MU;
extern map<int,double> elastic_LAMBDA;
extern map<int,int>    fixed_cells;

void construct_map_of_para();
void UserInitialParas();

#endif
