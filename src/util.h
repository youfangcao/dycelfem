// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     util.h              ***
 ***    Author:   Jieling Zhao,       ***
 ***              Youfang Cao         ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef UTIL_H
#define UTIL_H

#include <sstream>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <ctime>
#include <ctype.h>
#include "para.h"

using namespace std;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
void str_split(string input_str, vector <string> & results, string delim);
void matrix_multiplication(double D1[3][3], double D2[3][3], double D3[3][3]);
void matrix_vector_multiplication(double D1[3][3], double D2[3], double D3[3]);
double determinant_matrix_3(double D[3][3]);
double determinant_matrix_2(double D[2][2]);
double clock_angle(double n1x, double n1y, double n2x, double n2y);
double vector2angle(double nx, double ny);
double get_lame_constant(double E, double v, int i); // E: Young's modulus; v: poisson's ratio; i: 0:mu, 1:lambda
void matrix_transpose(double D1[3][3], double D2[3][3]);
bool CLIPt(double denom, double num, double &tE, double &tL);
void clipping(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &zmin, bool &accept);
double max_val2(double a, double b);
double max_val3(double a, double b, double c);
int line_segment_intersection(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
int powerint(int base, int index);
double U_to_N(); // generate normal distribution from uniform distribution
double SOR_solver(int sm_idt, double converge, double w, int steps, int &flagtime, int *SOR_N, int *JA, double *VA, double *F, double *Xn, double *X, double *SOR_II_V);
double Euler (vector<int>* ia, vector<int>* ja, vector<double>* a, vector<double>* v, double t);
double RK4 (vector<int>* ia, vector<int>* ja, vector<double>* a, vector<double>* v, double t);

#endif