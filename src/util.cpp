// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     util.cpp            ***
 ***    Author:   Jieling Zhao,       ***
 ***              Youfang Cao         ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef UTIL_CPP
#define UTIL_CPP

#include "util.h"

using namespace std;

void str_split(string input_str, vector <string> & results, string delim)
{
	size_t cutAt;
	while( (cutAt = input_str.find_first_of(delim)) != input_str.npos )
	{
		if(cutAt > 0)
		{
			results.push_back(input_str.substr(0,cutAt));
		}
		input_str = input_str.substr(cutAt+1);
	}
	if(input_str.length() > 0)
	{
		results.push_back(input_str);
	}
}

void matrix_multiplication(double D1[3][3], double D2[3][3], double D3[3][3])
{
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			double element = 0;
			for (int k=0;k<3;k++)
			{
				element = element + D1[i][k]*D2[k][j];
			}
			D3[i][j] = element;
		}
	}
}

void matrix_vector_multiplication(double D1[3][3], double D2[3], double D3[3])
{
	for (int i=0;i<3;i++)
	{
		double element = 0;
		for (int j=0;j<3;j++)
		{
			element = element + D1[i][j]*D2[j];
		}
		D3[i] = element;
	}
}

double determinant_matrix_3(double D[3][3])
{
	/******************
	     00 01 02
	 D:  10 11 12
	     20 21 22
	******************/
	double value = D[0][0]*(D[1][1]*D[2][2] - D[2][1]*D[1][2]) -
				   D[1][0]*(D[0][1]*D[2][2] - D[2][1]*D[0][2]) +
				   D[2][0]*(D[0][1]*D[1][2] - D[1][1]*D[0][2]);
	return value;
}

double determinant_matrix_2(double D[2][2])
{
	/******************
	 D:  00 01
	     10 11
	******************/
	double value = D[0][0]*D[1][1] - D[1][0]*D[0][1];
	return value;
}

double clock_angle(double n1x, double n1y, double n2x, double n2y)
{
	double angle = 0;
	double angle1 = 0,angle2 = 0;
	
	if      (n1x>=0 && n1y>=0)
	{angle1 = asin(n1y)*180.0/PI;}
	else if (n1x>=0 && n1y<0)
	{angle1 = asin(n1y)*180.0/PI + 360.0;}
	else if (n1x<0 && n1y>=0)
	{angle1 = 180.0 - asin(n1y)*180.0/PI;}
	else if (n1x<0 && n1y<0)
	{angle1 = 180.0 - asin(n1y)*180.0/PI;}

	if      (n2x>=0 && n2y>=0)
	{angle2 = asin(n2y)*180.0/PI;}
	else if (n2x>=0 && n2y<0)
	{angle2 = asin(n2y)*180.0/PI + 360.0;}
	else if (n2x<0 && n2y>=0)
	{angle2 = 180.0 - asin(n2y)*180.0/PI;}
	else if (n2x<0 && n2y<0)
	{angle2 = 180.0 - asin(n2y)*180.0/PI;}

	angle = angle2 - angle1;
	if (angle<0) {angle = 360.0 + angle;}
	
	return angle;
}

double vector2angle(double nx, double ny)
{
	double angle = atan2(ny, nx) * 180.0 / PI;
	if (angle < 0)
	{
		angle += 360;
	}
	return angle;
}

double get_lame_constant(double E, double v, int i)
{
	/**********************
	   E: Young's modulus 
       v: poisson's ratio
	   i: 0 -> mu
	      1 -> lambda
	**********************/
	double value = 0;
	if (i==0)
	{
		value = E/2/(1+v);
	}
	else if (i==1)
	{
		value = E*v/(1+v)/(1-2*v);
	}
	return value;
}

void matrix_transpose(double D1[3][3], double D2[3][3])
{
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			D2[i][j] = D1[j][i];
		}
	}
}

bool CLIPt(double denom, double num, double &tE, double &tL)
{
	bool valid = true;
	double t;
	if (denom>0)
	{
		t = num/denom;
		if (t>tL) {valid = false;}
		else if (t<tL && t>tE) {tE = t;}
	}
	else if (denom<0)
	{
		t = num/denom;
		if (t<tE) {valid = false;}
		else if (t>tE && t<tL){tL = t;}
	}
	else if (num>0) {valid = false;}
	return valid;
}

void clipping(double &x0, double &y0, double &z0, double &x1, double &y1, double &z1, double &zmin, bool &accept)
{
	double tmin = 0, tmax = 1;
	double dx = x1 - x0;
	double dz = z1 - z0;
	accept = false;
	if (CLIPt(-dx-dz,x0+z0,tmin,tmax))
	{
		if (CLIPt(dx-dz,-x0+z0,tmin,tmax))
		{
			double dy = y1 - y0;
			if (CLIPt(dy-dz,-y0+z0,tmin,tmax))
			{
				if (CLIPt(-dy-dz,y0+z0,tmin,tmax))
				{
					if (CLIPt(-dz,z0-zmin,tmin,tmax))
					{
						if (CLIPt(dz,-z0-1,tmin,tmax))
						{
							accept = true;
							if (tmax<1)
							{
								x1 = x0 + tmax*dx;
								y1 = y0 + tmax*dy;
								z1 = z0 + tmax*dz;
							}
							if (tmin>0)
							{
								x0 = x0 + tmin*dx;
								y0 = y0 + tmin*dy;
								z0 = z0 + tmin*dz;
							}
						}
					}
				}
			}
		}
	}
}

double max_val2(double a, double b)
{
	double value;
	if (a>=b) {value=a;} 
	else {value=b;}
	return value;
}

double max_val3(double a, double b, double c)
{
	double value;
	if      (a>=b && a>=c) {value=a;} 
	else if (b>=a && b>=c) {value=b;}
	else if (c>=a && c>=b) {value=c;}
	return value;
}

int line_segment_intersection(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
	int index = -1;                                  // -1: no intersection
	double a = x1 - x2;
	double b = x4 - x3;
	double c = y1 - y2;
	double d = y4 - y3;
	double e = x4 - x2;
	double f = y4 - y2;
	double B = b*c - a*d;
	double tB = e*c - a*f;
	double sB = b*f - e*d;
	double t = tB/B;
	double s = sB/B;
	
	if ((s>0 && s<1) && (t>0 && t<1))
	{
		if      (s<0.5 && t< 0.5) {index = 3;}
		else if (s>0.5 && t> 0.5) {index = 4;}
		else if (s>0.5 && t<=0.5) {index = 5;}
		else if (s<=0.5 && t>0.5) {index = 6;}
	}
	else if ((s==0 || s==1) && (t==0 || t==1))
	{
		if      (s ==1 && t == 1) {index = 1;}
		else if (s ==0 && t == 0) {index = 2;}
		else if (s ==0 && t == 1) {index = 7;}
		else if (s ==1 && t == 0) {index = 8;}
	}
	return index;
}

int powerint(int base, int index)
{
	int value = 1;
	if (index>0)
	{
		for (int i=0;i<index;i++)
		{
			value = value*base;
		}
	}
	return value;
}

double U_to_N()
{
	/******************************************************
	 Convert universal distribution to normal distribution
	               Box-Muller transform
	                 U(0,1) => N(0,1)
	 ======================================================
	 (1) U1 ~ U[0,1]; U2 ~ U[0,1];
	 (2) theta = 2*PI*U2; rho = sqrt(-2*log(U1));
	 (3) Z1 = rho*cos(theta);
	     Z2 = rho*sin(theta);
	 both Z1 and Z2 are normal variables
	*******************************************************/
	double value = 0;
	double x1 = ((double)rand()+1.0)/((unsigned long)RAND_MAX+1);
	double x2 = ((double)rand()+1.0)/((unsigned long)RAND_MAX+1);
	double theta = x2;
	double rho = sqrt(-2*log(x1));
	value = rho*cos(theta*PI*2);
	return value;
}

double SOR_solver(int sm_idt, double converge, double w, int steps, int &flagtime, int *SOR_N, int *JA, double *VA, double *F, double *Xn, double *X, double *SOR_II_V)
{
	/*************************************
	 SOR input:
	 ====================================
	 sm_idt: size of solution
	 converge: convergence threshold
	 steps: number of steps up to
	 flagtime: running steps
	 JA: index of non-zero element
	 VA: value of non-zero element
	 SOR_N: range of indices of each row
	 Xn, X: solution at each step
	 SOR_II_V: value of diagnal element
	*************************************/
	double diff = 100;
	double max_X = 0;
	while (diff>converge && flagtime<steps)
	{
		flagtime++;
		int n_col_1 = 0;
		int n_col_2 = SOR_N[0];
		double gamma = 0;
		for (int i=n_col_1;i<n_col_2;i++)
		{
			int ja = JA[i];
			double aij = VA[i];
			gamma += aij*Xn[ja];
		}
		Xn[0] = (1-w)*Xn[0] + w/SOR_II_V[0]*(F[0] - gamma);
		for (int i=1;i<sm_idt;i++)
		{
			n_col_2 = SOR_N[i];
			n_col_1 = SOR_N[i-1];
			gamma = 0;
			for (int j=n_col_1;j<n_col_2;j++)
			{
				int ja = JA[j];
				double aij = VA[j];
				gamma += aij*Xn[ja];
			}
			Xn[i] = (1-w)*Xn[i] + w/SOR_II_V[i]*(F[i] - gamma);
		}
		diff = 0;
		for (int i=0;i<sm_idt;i++)
		{
			double dif = abs(Xn[i] - X[i]);
			if (diff<dif) {diff = dif;}
			if (abs(X[i]) > max_X) { max_X = abs(X[i]); }
			X[i] = Xn[i];
		}
		diff /= max_X;
	}
	return diff;
}

double Euler (vector<int>* ia, vector<int>* ja, vector<double>* a, vector<double>* v, double t)
{
	/***********************************************************
	  Linear Euler solver for finite element diffusion equation
	  dx/dt = Ax
	  YOUFANG CAO
	***********************************************************/
	double dt = 0.1;
	double tau = 0;
	double err = 1;
	double testv = 0;
	int    nnegv = 0;
	int    n = v->size();
	int    nnz = ia->size ();
	int    ns = 0;
	double maxv = 0;
	double maxa = 0;
	for (int i = 0; i < n; i++)
	{
		if (v->at(i) > maxv) { maxv = v->at(i); }
	}
	if (maxv < 1e-4) { return 0.0; }
	for (int i = 0; i < nnz; i++)
	{
		if (abs(a->at(i)) > maxa) { maxa = abs(a->at(i)); }
	}
	if (maxa < 1e-8) { return 0.0; }

	double* w = new double[n];
	while (tau < t)
	{
		for (int i = 0; i < n; i++) { w[i] = 0.0; }
		for (int i = 0; i < nnz; i++)
		{
			w[ia->at(i)] += a->at(i) * v->at(ja->at(i)) * dt;
		}
		/**/nnegv = 0;
		for (int i = 0; i < n; i++)
		{
			testv = v->at(i) + w[i];
			if (testv < 0)
			{
				nnegv ++;
				break;
			}
		}
		if (nnegv == 0)
		{
			err = 0;
			for (int i = 0; i < n; i++)
			{
				v->at(i) += w[i];
				err = max (err, w[i]);
			}
		}
		/**/else if (dt > 1e-4)
		{
			dt = dt/10.0;
			cout << "          WARNING: Concentration < 0, big time step suspected, change to smaller dt=" << dt << endl;
			continue;
		}
		else
		{
			cout << "          WARNING: dt too big, Euler not done." << endl;
			break;
		}
		tau += dt;
		ns ++;
	}
	delete[] w;
	return err;
}

double RK4 (vector<int>* ia, vector<int>* ja, vector<double>* a, vector<double>* v, double t)
{
	/***********************************************************
	  Linear Runge-Kutta solver for finite element diffusion equation
	  dx/dt = Ax
	  YOUFANG CAO
	***********************************************************/
	double dt = 0.1;
	double dthalf = 0.5 * dt;
	double h6 = dt / 6.0;
	double err = 0.0;
	double tau = 0;
	double testv = 0;
	int    nnegv = 0;
	int    n = v->size();
	int    nnz = ia->size ();
	int    ns = 0;
	double maxv = 0;
	double maxa = 0;
	for (int i = 0; i < n; i++)
	{
		if (v->at(i) > maxv) { maxv = v->at(i); }
	}
	if (maxv < 1e-4) { return 0.0; }
	for (int i = 0; i < nnz; i++)
	{
		if (abs(a->at(i)) > maxa) { maxa = abs(a->at(i)); }
	}
	if (maxa < 1e-8) { return 0.0; }


	double *w1;
	double *w2;
	double *w3;
	double *w4;
	double *w0;
	double *delta;

	w1 =  new double[n];
    w2 =  new double[n];
    w3 =  new double[n];
    w4 =  new double[n];
    w0 =  new double[n];
	delta =  new double[n];
	
	while(tau < t)
	{
		for (int i=0; i<n; i++) 
		{
			w0[i] = w1[i] = w2[i] = w3[i] = w4[i] = 0.0;
			delta[i] = 0.0;
		}

		for (int i=0; i<nnz; i++)
		{
			w1[ia->at(i)] += a->at(i) * v->at(ja->at(i)); 
		}

		for (int i=0; i<n; i++)
		{
			w0[i] = v->at(i) + dthalf * w1[i];
		}

		for (int i=0; i<nnz; i++)
		{
			w2[ia->at(i)] += a->at(i) * w0[ja->at(i)];
		}

		for (int i=0; i<n; i++)
		{
			w0[i] = v->at(i) + dthalf * w2[i];
		}

		for (int i=0; i<nnz; i++)
		{
			w3[ia->at(i)] += a->at(i) * w0[ja->at(i)];
		}

		for (int i=0; i<n; i++)
		{
			w0[i] = v->at(i) + dt * w3[i];
		}

		for (int i=0; i<nnz; i++)
		{
			w4[ia->at(i)] += a->at(i) * w0[ja->at(i)];
		}

		nnegv = 0;
		err = 0;
		for (int i = 0; i < n; i++)
		{
			delta[i] = h6 * (w1[i] + 2*w2[i] + 2*w3[i] + w4[i]);
			testv = v->at(i) + delta[i];
			if (testv < 0)
			{
				nnegv ++;
				break;
			}
			if (abs(delta[i]) > err)
			{
				err = abs(delta[i]);
			}
		}

		if (nnegv == 0)
		{
			for (int i = 0; i < n; i++)
			{
				v->at(i) += delta[i];
			}
		}
		else if (nnegv > 0)
		{
			if (dt > 1e-4)
			{
				dt = dt/10.0;
				cout << "          WARNING: Concentration < 0, big time step suspected, change to smaller dt=" << dt << endl;
				continue;
			}
			else
			{
				cout << "          WARNING: dt too small, RK4 not done." << endl;
				break;
			}
		}

		tau += dt;
	}

	delete[] w0;
	delete[] w1;
	delete[] w2;
	delete[] w3;
	delete[] w4;

	return err;
}




#endif

