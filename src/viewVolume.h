// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     viewVolume.h        ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef VIEWVOLUME_H
#define VIEWVOLUME_H

#include <vector>

using namespace std;

class viewVolume;
class point2D;

class viewVolume {
	private:
		// Members
		double uMin;
		double uMax;
		double vMin;
		double vMax;
		double uMin_default;
		double uMax_default;
		double vMin_default;
		double vMax_default;
		
		point2D *vpn;
		point2D *vup;
		point2D *vrp;
		point2D *prp;
	
	public:
		// Constructor/Destructor
		viewVolume();
		~viewVolume();
		
		// Methods
		void setUMin(double u)  {uMin = u;}
		void setUMax(double u)  {uMax = u;}
		void setVMin(double v)  {vMin = v;}
		void setVMax(double v)  {vMax = v;}
		void setUMin_default(double u)  {uMin_default = u;}
		void setUMax_default(double u)  {uMax_default = u;}
		void setVMin_default(double v)  {vMin_default = v;}
		void setVMax_default(double v)  {vMax_default = v;}
		void setVPN(point2D *p) {vpn = p;}
		void setVUP(point2D *p) {vup = p;}
		void setVRP(point2D *p) {vrp = p;}
		void setPRP(point2D *p) {prp = p;}
	
		double getUMin()  {return uMin;}
		double getUMax()  {return uMax;}
		double getVMin()  {return vMin;}
		double getVMax()  {return vMax;}
		double getUMin_default()  {return uMin_default;}
		double getUMax_default()  {return uMax_default;}
		double getVMin_default()  {return vMin_default;}
		double getVMax_default()  {return vMax_default;}
		point2D* getVPN() {return vpn;}
		point2D* getVUP() {return vup;}
		point2D* getVRP() {return vrp;}
		point2D* getPRP() {return prp;}
};

class point2D {
	private:
		double u;
		double v;
	public:
		// Constructor/Destructor
		point2D(double x, double y) {u=x;v=y;}
		~point2D() {};

		// Methods
		double x() {return u;}
		double y() {return v;}
		void setx(double x) {u = x;}
		void sety(double y) {v = y;}
};

#endif




