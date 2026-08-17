// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     delaunay.h          ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <string>
#include <vector>
#include "point3D.h"

using namespace std;

class CXYZW;
class CEdge;
class CCell;
class CSimplex;
class C0Simplex;
class C1Simplex;
class C2Simplex;
class CComplex;

class CXYZW
{
	protected:
		double X;
		double Y;
		double W;
		int N;
		void *Label;
		C0Simplex *S0;
		point3D *p;
		bool B;

	public:
		CXYZW() {N=0;B=false;};
		~CXYZW();
		int n() {return N;}
		double x() {return X;}
		double y() {return Y;}
		double w() {return W;}
		point3D* P() {return p;} 
		bool BN() {return B;}
		void set_BN(bool valid) {B=valid;}
		void set_p(point3D* P) {p = P;}
		void* label() const {return Label;}
		void set(double x,double y,double w) {X=x;Y=y;W=w;}
		void set(double x,double y,double w,int n,void *l) {set(x,y,w);N=n;Label=l;}
		bool operator==(CXYZW &A) const {return (double)X==(double)A.X && (double)Y==(double)A.Y;}
	friend class CD3DW;
	friend class CEdge;
	friend class CCell;
	friend class CSimplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
	friend class CComplex;
};

class CEdge
{
	protected:
		CXYZW *A,*B;
		CCell *Cell;
		CEdge *Rev,*Enx,*Fnx;
		char Flag; // 1:in stack;2:deleted
		int N;
		double center[2];
		double alphavalue;
		C1Simplex *S1;
		C2Simplex *S2;
		bool Bn;

	public:
		CEdge(CXYZW *a,CXYZW *b) {A=a;B=b;Flag=0;Bn=false;}
		~CEdge();
		void setA(CXYZW* a) {A = a;}
		void setB(CXYZW* b) {B = b;}
		const CXYZW* a() const {return A;}
		const CXYZW* b() const {return B;}
		void set_center();
		bool BN() {return Bn;}
		void set_BN(bool valid) {Bn=valid;}
	friend class CD3DW;
	friend class CCell;
	friend class CSimplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
	friend class CComplex;
};

class CCell
{
	protected:
		CXYZW *VV[3];
		CEdge *ab,*bc,*ca;
		CCell *CC[3];
		double center[2];
		double alphavalue;
		bool B;
		int type; // 0: in, 1: out

	public:
		CCell() {CC[0]=0;CC[1]=0;CC[2]=0;B=false;type = 1;}
		~CCell();
		void set_center();
		CXYZW* get_VV(int i) {return VV[i];}
		CEdge* getab() {return ab;}
		CEdge* getbc() {return bc;}
		CEdge* getca() {return ca;}
		double get_center(int i) {return center[i];}
		double get_alphav() {return alphavalue;}
		bool BN() {return B;}
		void set_BN(bool valid) {B=valid;}
		int get_type() {return type;}
		void set_type(int i) {type = i;}
		bool inside_complex(); // for the first step of inner mesh generation
	friend class CD3DW;
	friend class CEdge;
	friend class CSimplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
	friend class CComplex;
};

class CD3DW
{
	protected:
		std::vector<CEdge*> Es; // edge list
		std::vector<CXYZW*> Vs; // vertex list
		std::vector<CCell*> Cs; // triangle list
		std::vector<CEdge*> St; // edge list for flip process
		CXYZW *V0,*V1,*V2;
		
	public:
		CD3DW(double range[4]); // the boundary vertices, removed after triangulation is done
		~CD3DW();
		CXYZW* operator[](int n) {return Vs[n];}
		vector<CCell*> *get_Cs() {return &Cs;}
		vector<CXYZW*> *get_Vs() {return &Vs;}
		void add(double x, double y, double w, int n, void *l);	// each time, randomly toss one point into the box
		void add_resample(double x, double y, double w, int n, point3D *l);	// each time, randomly toss one point into the box
		void add_interior(double l);
		void add_interior_debug(double l);
		inline CEdge* locate(CXYZW *V);                         // locate the triangle which contains the inserted point
		inline bool judge_on_edge(CCell *C,CXYZW *V);                  // judge if point is on the edge of one triangle
		inline CEdge* get_on_edge(CCell *C,CXYZW *V);          // if point is on the edge, return the edge
		inline double judge_inside(CCell *C,CXYZW *V);
		inline bool inside(CCell *C,CXYZW *V);                  // test if the point is inside one triangle
		inline bool insphere(CCell *C,CXYZW *V);                // test if the point is inside the circumcircle of the triangle
		inline void flip(CEdge *ab);                            // flip the edge
		inline void flip_resample(CEdge *ab);                   // flip the edge
	friend class CComplex;
};

class CSimplex
{
	protected:
		double alphavalue;
		double center[2];
		CEdge *E;
	public:
		~CSimplex() {};
		virtual CXYZW* operator[](int)=0;
		virtual int size() = 0;
		double x(int i) {return operator[](i)->X;}
		double y(int i) {return operator[](i)->Y;}
		double r(int i) {return operator[](i)->W;}
		double alpha() {return alphavalue;}
		double Center(int i) {return center[i];}
		int n(int i) {return operator[](i)->N;}
		void* label(int i) {return operator[](i)->Label;}

	friend class CComplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
};

class C0Simplex: public CSimplex
{
	public:
		virtual ~C0Simplex() {};
		CXYZW *A;
		int a;
		int size() {return 1;}
		CXYZW* operator[](int);
	friend class CComplex;
};

class C1Simplex: public CSimplex
{
	public:
		virtual ~C1Simplex() {};
		CXYZW *A,*B;
		int a,b;
		int size() {return 2;}
		CXYZW* operator[](int);
		void geometry();       // alphavalue calculation
	friend class CComplex;
};

class C2Simplex: public CSimplex
{
	public:
		virtual ~C2Simplex() {};
		CXYZW *A,*B,*C;
		int a,b,c;
		int size() {return 3;}
		CXYZW* operator[](int);
		void geometry();       // alphavalue calculation
	friend class CComplex;
};

class CComplex
{
	public:
		CComplex(CD3DW*D);
		~CComplex();
		std::vector<CSimplex*> P;
		void sort(CSimplex *A);
};

#endif

