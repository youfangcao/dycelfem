// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     dbReader.h          ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef DBREADER_H
#define DBREADER_H

#include <string>
#include <vector>
#include <ctime>
#include "util.h"
#include "viewVolume.h"
#include "biology.h"
#include "point3D.h"
#include "delaunay.h"
#include "para.h"

using namespace std;

class Biology;

class dbReader {
	private:
		// Members
		int initial;
		viewVolume *vv;
		vector<point3D*> vertexList;
		vector<edge*> edgeList;
		vector<cell*> cellList;
		CComplex *Complex;
		vector<interpair*> collisionpairList;
		vector<node*> nodeList;
		Biology *biology;
		
	public:
		// Constructor/Destructor
		dbReader();
		~dbReader() {};
	
		// Methods
		viewVolume* getViewVolume();
		vector<point3D*> *V() {return &vertexList;}
		vector<edge*> *E() {return &edgeList;}
		vector<cell*> *C() {return &cellList;}
		CComplex *get_Complex() {return Complex;}
		vector<interpair*> *CL() {return &collisionpairList;}
		vector<node*> *NL() {return &nodeList;}
		Biology *get_biology() {return biology;}

		void read(char* fN, int number);						// read input file
		void delaunay_triangulation(double range[4]);			// delaunay triangulation for initialization
		int detect_collision(int number);                      // cell collision detection
		int intersection_type(point3D *P1, point3D *P2);
		double cell_distance(cell *C1, cell *C2);
		void collision_refresh(interpair *P);
		void collision_refresh_initial();													// only for initial collision detection step
		void cell_division(cell *C, int life, int type_Co, int type_Cn, double div_angle);  // INTERFACE: cell division: Co: mother cell; Cn: daughter cell; div_angle: division angle
		void cell_death_response(cell *C);													// only for single cell death
		void cell_migrate_response_0(cell *C, double angle, vector<int>* d, double F);		// angle: the migration direction; F: migration force;
		void cell_migrate_response_1(cell *C, double angle, vector<int>* d, double F);		// angle: the migration direction; F: migration force;
		cell* cell_migrate_angle(cell *C, double angle);                                    // return the closest cell on the angle direction 
		void cell_death_correction();														// cell death manipulation
		void collision_candidate(interpair *P,vector<point3D*> & P1,vector<point3D*> & P2);
		int collision_candidate2(interpair *P,vector<point3D*> & P1,vector<point3D*> & P2);
		void boundary_merge_1(point3D *S1, point3D *T1, point3D *SS1, point3D *TT1, int Pair_id);								// already adjacent
		void boundary_merge_0(point3D *S1, point3D *T1, point3D *S2, point3D *T2, int Pair_id);                                 // first time
		void boundary_merge_4(point3D *S1, point3D *T1, point3D *SS1, point3D *TT1, int Pair_id);								// already adjacent
		void boundary_merge_3(point3D *S1, point3D *T1, point3D *S2, point3D *T2, int Pair_id);	
		void boundary_merge_2(point3D *S1, point3D *T1, point3D *S2, point3D *T2);									// only for line segments between the two endpoints
		void boundary_merge_2_refresh(point3D **S1, point3D **T1, point3D **S2, point3D **T2, cell *CR1, cell *CR2);// check collapsed case when one cell is penetrating another one
		int boundary_merge_2_check(interpair *P, cell *C1, cell *C2, point3D *s1, point3D *s2, point3D **ps1, point3D **ps2); // switch: 0,1,2                                // check collapsed case when one cell is penetrating another one
		void boundary_merge_refresh(point3D *S1,point3D *T1, point3D *S2, point3D *T2);								// endpoints reassignment
		void collision_response(interpair *P);                   // collision response
		void collision_response_m(interpair *P);                 // collision response merge
		int collision_response_m2(interpair *P, int step, int round);                // collision response merge
		int collision_response_m3(interpair *P, interpair *Pr, int step, int round); // collision response merge: Pr: reference pair
		int edge_response_m(edge *E);                // edge length response:merge the edge into single vertex
		void edge_response_d(edge *E);                // edge length response: right divide from the middle
		void edge_response_b0(edge *E);               // break the edge pair; (2 points)
		void edge_response_b1(edge *E);               // break the edge pair; (3 points)
		void point_response_b0(point3D *V);           // break the point pair; (2 points)
		void point_response_b1(point3D *V);           // break the point pair; (3 points)
		void point_response_b2(point3D *V);           // break the point pair; (3 points)
		void cell_dynamic(int number);                // cell growth dynamic process, SIMPLE
		int cell_dynamic_FEM(int number);            // cell growth dynamic process, FEM
		int cell_edge_response(int time);            // cell edge resample 
		int cell_pick_response(int time);            // cell pick correction
		void cell_angle_response();                   // cell angle resample
		void cell_mesh_generation(int time);          // to generate triangular mesh for each cell
		void node_merge(point3D* P1, point3D* P2);    // to merge the nodes (for stiffness matrix)
		void node_remove_0(point3D* P);               // to remove one vertex from the node and the vertex is still in the cell
		void node_remove_1(point3D* P);               // to remove one vertex from the node and the vertex is not in the cell
		int pick_correction(point3D* P);             // to correct geometric shape
		void cell_interior_resample(cell* C, int time);						// mesh generation of cell
		void cell_simple_relaxation(cell *C, int time);						// simple relaxation
		bool cell_self_collision(cell *C, int time);						// cell self-collision detection
		int migrating_cell_rear(cell *C, int cell_type, int time);         // the rear side of migrating cell starts to grow to fill up the hole
		void cell_friction_response(cell *C, int time);						// friction force detect for migrating cell
		int cell_setupbiobehaviors(int Time);
		void node_test();                             // only to test if node related strategies work well
		void print_matlab_file();                     // for matlab manipulation
		void print_cell_file(int time);               // for non-sphere like shape cells
	
	friend class interpair;
	friend class mergepair;
};

#endif
