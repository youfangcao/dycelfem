// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     point3D.h           ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef POINT3D_H
#define POINT3D_H

#include <vector>
#include "para.h"
#include "util.h"

using namespace std;

class point3D;
class edge;
class triangle;
class cell;
class mergepair;
class interpair;
class node;
class le_sm;
class rear_fc;

class point3D {
	private:
		// Members
		int index;
		double u;
		double v;
		double u_3[4];
		double v_3[4]; // coordinates of last 3 steps stored 
		bool merge;
		bool in_cell;
		vector<point3D*> mp; // merged points
		point3D *fp;         // adjacent point ccw
		point3D *rp;         // adjacent point cw
		edge *ab;            // a->b ccw
		edge *ba;            // b->a cw
		cell *C;
		double Tension_ab[2];
		double Pressure_ab[2];
		double Tension_ba[2];
		double Pressure_ba[2];
		double Pressure_increment[2];
		double single_Force[2];
		double Force[2];
		double P_in_ab[2]; // pressure respect to edge ab
		double P_in_ba[2]; // pressure respect to edge ba
		double P_in[2];    // pressure due to incremental area: elastic model
		double Dp[2];      // displacement for simple relaxation
		double Mass;
		double A_ij[4]; // stiffness matrix: A_11, A_12, A_21, A_22
		double A_ij_single[4]; // stiffness matrix for single cell: A_11, A_12, A_21, A_22
		bool node_check; // this is to build up the global stiffness matrix
		bool cell_boundary;
		bool burry; // this is for interior points, true if it is wrapped by other interior points
		bool end_point;
		bool migrate_mark;
		int pop_out; // point having tendency to pop out: 0=>static; 1=>grow; -1=>shrink;
		double adhesion_break[2]; // when the edge is broken apart, the breaking force is stored for next step
		double stored_apart_force; // the force stored for adhesion break
		double stored_force_vis[2]; // the stored force for display
		double stored_friction_force; // the force stored for friction sliding
		double elastic_force[2]; // the internal force recovery for each step
		double elastic_force_migration[2]; // the internal force recovery in response to migration each step
		double migrate_response_force[2]; // when cell C migrates against cell Cn, the pressure force added on the neighbor cells on counter direction of migration
		double migration_rate; // for migration cells
		double single_growth_rate; // only for one single point
		node* node_id; // the id of temporal node vector
		vector<triangle*> inner_at;
		vector<le_sm*> node_sm; // the element sm for one edge
		int node_sm_n;
		int SM_ID;
		int SM_ID_single;
		vector<le_sm*> node_sm_single; // the element sm for one edge (only in one single cell)
		double pre_dp_mag; // the expected displacement magnitude

	public:
		// Constructor/Destructor
		point3D(int id, double x, double y);
		~point3D();

		// Methods
		int id() {return index;}
		double x() {return u;}
		double y() {return v;}
		double x_old(int i) {return u_3[i];}
		double y_old(int i) {return v_3[i];}
		vector<point3D*> *getmp() {return &mp;}
		bool get_in_cell() {return in_cell;}
		bool get_node_check() {return node_check;}
		bool get_cell_boundary() {return cell_boundary;}
		bool get_burry() {return burry;}
		bool get_end_point() {return end_point;}
		bool get_migrate_mark() {return migrate_mark;}
		int get_pop_out() {return pop_out;}
		double get_adhesion_break(int i) {return adhesion_break[i];}
		double get_stored_apart_force() {return stored_apart_force;}
		double get_stored_force_vis(int i) {return stored_force_vis[i];}
		double get_stored_friction_force() {return stored_friction_force;}
		double get_elastic_force(int i) {return elastic_force[i];}
		double get_elastic_force_migration(int i) {return elastic_force_migration[i];}
		double get_migrate_response_force(int i) {return migrate_response_force[i];}
		double get_migration_rate() {return migration_rate;}
		double get_single_growth_rate() {return single_growth_rate;}
		void set_node_check(bool valid) {node_check = valid;}
		void set_cell_boundary(bool valid) {cell_boundary = valid;}
		void set_burry(bool valid) {burry = valid;}
		void set_end_point(bool valid) {end_point = valid;}
		void set_migrate_mark(bool valid) {migrate_mark = valid;}
		void set_pop_out(int valid) {pop_out = valid;}
		void set_adhesion_break(double x, double y) {adhesion_break[0] = x;adhesion_break[1] = y;}
		void set_stored_apart_force(double v) {stored_apart_force = v;}
		void set_stored_force_vis(double x, double y) {stored_force_vis[0] = x;stored_force_vis[1] = y;}
		void set_stored_friction_force(double v) {stored_friction_force = v;}
		void set_elastic_force(double x, double y) {elastic_force[0] = x;elastic_force[1] = y;}
		void set_elastic_force_migration(double x, double y) {elastic_force_migration[0] = x;elastic_force_migration[1] = y;}
		void set_migration_rate(double MR) {migration_rate = MR;}
		void set_single_growth_rate(double SR) {single_growth_rate = SR;}
		void set_migrate_response_force(double x, double y) {migrate_response_force[0] = x;migrate_response_force[1] = y;}
		node* get_node_id() {return node_id;}
		void set_node_id(node* N) {node_id = N;}
		void set_in_cell(bool valid) {in_cell = valid;}
		point3D* getfp() {return fp;}
		point3D* getrp() {return rp;}
		edge* getab() {return ab;}
		edge* getba() {return ba;}
		cell* Cell() {return C;}
		double getTension_ab(int i) {return Tension_ab[i];}
		double getTension_ba(int i) {return Tension_ba[i];}
		double getPressure_ab(int i) {return Pressure_ab[i];}
		double getPressure_ba(int i) {return Pressure_ba[i];}
		double getPressure_increment(int i) {return Pressure_increment[i];}
		double getP_in_ab(int i) {return P_in_ab[i];}
		double getP_in_ba(int i) {return P_in_ba[i];}
		double getP_in(int i) {return P_in[i];}
		double getForce(int i) {return Force[i];}
		double get_single_Force(int i) {return single_Force[i];}
		double getDp(int i) {return Dp[i];}
		double getMass() {return Mass;}
		double getA_ij(int i) {return A_ij[i];}
		double getA_ij_single(int i) {return A_ij_single[i];}
		void setid(int id) {index = id;}
		void setx(double x) {u=x;}
		void sety(double y) {v=y;}
		void setx_old(double x, int i) {u_3[i] = x;}
		void sety_old(double y, int i) {v_3[i] = y;}
		void pushmp(point3D* P);
		bool checkmp(point3D *P);
		void clearmp();
		void removemp(point3D* P);
		void setfp(point3D* P) {fp = P;}
		void setrp(point3D* P) {rp = P;}
		void setab(edge* E) {ab = E;}
		void setba(edge* E) {ba = E;}
		void setCell(cell* Cell) {C = Cell;}
		void setForce(double Fx, double Fy) {Force[0] = Fx;Force[1] = Fy;}
		void set_single_Force(double Fx, double Fy) {single_Force[0] = Fx;single_Force[1] = Fy;}
		void setDp(double vx, double vy) {Dp[0] = vx;Dp[1] = vy;}
		void setP_in_ab(double Fx, double Fy) {P_in_ab[0] = Fx;P_in_ab[1] = Fy;}
		void setP_in_ba(double Fx, double Fy) {P_in_ba[0] = Fx;P_in_ba[1] = Fy;}
		void setP_in(double Fx, double Fy) {P_in[0] = Fx;P_in[1] = Fy;}
		void setTension_ab(double Tx, double Ty) {Tension_ab[0] = Tx;Tension_ab[1] = Ty;}
		void setTension_ba(double Tx, double Ty) {Tension_ba[0] = Tx;Tension_ba[1] = Ty;}
		void setPressure_ab(double Px, double Py) {Pressure_ab[0] = Px;Pressure_ab[1] = Py;}
		void setPressure_ba(double Px, double Py) {Pressure_ba[0] = Px;Pressure_ba[1] = Py;}
		void setPressure_increment();
		void setMass(double M) {Mass = M;}
		void setA_ij(double A11, double A12, double A21, double A22) {A_ij[0] = A11;A_ij[1] = A12;A_ij[2] = A21;A_ij[3] = A22;}
		void setA_ij_single(double A11, double A12, double A21, double A22) {A_ij_single[0] = A11;A_ij_single[1] = A12;A_ij_single[2] = A21;A_ij_single[3] = A22;}
		vector<triangle*> *get_inner_at() {return &inner_at;}
		void push_inner_at(triangle* T) {inner_at.push_back(T);}
		void clear_inner_at();
		double point_angle();
		vector<le_sm*> *get_node_sm() {return &node_sm;}
		inline int get_node_sm_n() {return node_sm_n;}
		inline void set_node_sm_n() {node_sm_n = (int)node_sm.size();}
		void push_node_sm(le_sm* L) {node_sm.push_back(L);}
		void clear_node_sm();
		inline int get_SM_ID() {return SM_ID;}
		inline void set_SM_ID(int n) {SM_ID = n;}
		inline int get_SM_ID_single() {return SM_ID_single;}
		inline void set_SM_ID_single(int n) {SM_ID_single = n;}
		vector<le_sm*> *get_node_sm_single() {return &node_sm_single;}
		void push_node_sm_single(le_sm* L) {node_sm_single.push_back(L);}
		void clear_node_sm_single();
		double get_pre_dp_mag() {return pre_dp_mag;}
		void set_pre_dp_mag(double va) {pre_dp_mag = va;}

	friend class edge;
	friend class cell;
	friend class triangle;
	friend class interpair;
	friend class node;
	friend class le_sm;
};

class edge {
	private:
		// Members
		int index;
		point3D* P1;
		point3D* P2;
		cell* C;
		edge* Neighbor;
		interpair* Pair;
		bool attach;
		bool in_cell;
		double normal[2];
		double virtual_stiffness_matrix[4]; // for p1-p2
		double virtual_triangle_stiffness_matrix[6][6]; // for p1:p2:virtual_point
		double soften_triangle_stiffness_matrix[6][6];

	public:
		// Constructor/Destructor
		edge(int id, point3D* A, point3D* B, cell* Cell);
		~edge();

		// Methods
		int id() {return index;}
		point3D* p1() {return P1;}
		point3D* p2() {return P2;}
		void reset_p1(point3D *p1) {P1 = p1;}
		void reset_p2(point3D *p2) {P2 = p2;}

		cell* Cell() {return C;}
		void reset_Cell(cell *NC) {C = NC;} 
		interpair* get_Pair() {return Pair;}
		void set_Pair(interpair *PA) {Pair = PA;}
		edge* get_Neighbor() {return Neighbor;}
		bool get_attach() {return attach;}
		bool get_in_cell() {return in_cell;}
		double getnormal(int i);
		double getlength();
		void set_Neighbor(edge* Edge) {Neighbor = Edge;}
		void remove_Neighbor() {Neighbor = NULL;}
		void set_attach(bool valid) {attach = valid;}
		void set_in_cell(bool valid) {in_cell = valid;}
		void setnormal() {normal[0] = P2->y() - P1->y();normal[1] = P1->x() - P2->x();}
		bool edge_twist();
		bool judge_vert_out_side(double x, double y); // return true if the test point is outside: one the same direction of the normal
		void set_virtual_stiffness_matrix();          // add one virtual point A to make the matrix regular, add the stiffness to the global matrix
		void set_virtual_stiffness_matrix_value_only(); // add one virtual point A to make the matrix regular, not to add the stiffness to the global matrix
		void set_soften_stiffness_matrix();
		double get_virtual_stiffness_matrix(int i) {return virtual_stiffness_matrix[i];}
		double get_virtual_triangle_stiffness_matrix(int i, int j) {return virtual_triangle_stiffness_matrix[i][j];}
		double get_soften_triangle_stiffness_matrix(int i, int j) {return soften_triangle_stiffness_matrix[i][j];}
	
	friend class point3D;
	friend class cell;
	friend class interpair;
};

class triangle {
	private:
		// Members
		int index;
		point3D* A;
		point3D* B;
		point3D* C;
		double alphavalue;
		double center[2];
		double lame[2]; // 0: MU, 1: LAMBDA
		bool pick;
		double Stiffness_matrix[6][6];

	public:
		// Constructor/Destructor
		triangle(int id, point3D *a, point3D *b, point3D *c);
		~triangle();

		// Methods
		inline point3D* getA() {return A;}
		inline point3D* getB() {return B;}
		inline point3D* getC() {return C;}
		double alpha() {return alphavalue;}
		double get_center(int i) {return center[i];}
		double get_lame(int i) {return lame[i];}
		bool get_pick() {return pick;}
		void setA(point3D* a) {A = a;}
		void setB(point3D* b) {B = b;}
		void setC(point3D* c) {C = c;}
		void set_alpha(double v) {alphavalue = v;}
		void set_center(double x, double y) {center[0] = x;center[1] = y;}
		void set_lame(int i, double v) {lame[i] = v;}
		void set_pick(bool valid) {pick = valid;}
		void set_Stiffness_matrix();
		inline double get_Stiffness_matrix(int i, int j) {return Stiffness_matrix[i][j];}

	friend class point3D;
	friend class cell;
};

class cell {
	private:
		// Members
		int index;
		double r; // radius
		double center[2];
		point3D* AABB[4]; // 0:xmin, 1:ymin, 2:xmax, 3:ymax
		bool burried;
		bool circle;
		bool dead;
		bool set_dead;
		bool migrate;    // 0: static; 1: migrate;
		bool migrate_relax;
		bool sense_mech; // 0: no mechanical signal transitted; 1: mechanical signal transitted;
		bool sense_migr; // 0: no migration signal transitted; 1: migration signal transitted;
		double sense_mech_vec[2]; // mechanical signal strength and orientation: determine the polarity orientation
		double sense_migr_vec[2]; // mechanical signal strength and orientation in response with migration
		double sense_migr_angle; // mechanical signal direction: keep for 3 steps
		int mark_number; // 0: normal; 1: migrate; 2:soften
		int mAngle;  // migration angle
		bool soften;     // the cell is soften
		bool migrate_response; // the cells on the counter direction of migration
		int cell_status; // -1: shrink; 0: static; 1: grow
		vector<int> state; // vector for cell state
		int behavior; // 0: waiting for behavior determination; 1: static; 2: grow & division; 3: apoptosis; 4: migration
		int bhsteps[5]; // Track number of steps of above 5 different behaviors.
		vector<edge*> sides;
		vector<triangle*> inner_t;
		vector<point3D*> inner_p;
		vector<point3D*> inner_p_slip; // when there is a narrow corner in the cell
		int inner_p_n;
		double area;
		double area_p[10]; // areas of last 2 steps; 0: last step; 1: last step of last step
		double initial_area;
		double tension_coef;
		double pressure_coef;
		double mass_coef;
		int cell_type;   // 0: at static status, 2: stem cell, 3: 
		int life;        // how many steps it has already endured
		double color[3]; // cell color: R, G, B (between 0 and 1)
		double diffusion_color[3];
		double lame[2];  // 0: MU, 1: LAMBDA
		double growth_rate;
		vector<interpair*> pairs;
		bool interior_refresh;
		int interior_time;
		double migrate_angle;
		point3D* migrate_p[2]; // 0: sp 1:tp
		vector<rear_fc*> migrate_dead_filling_cells;

	public:
		// Constructor/Destructor
		cell(int id);
		~cell();

		//Methods
		int id() {return index;}
		void set_r(double R) {r = R;}
		void set_center(double x, double y) {center[0]=x;center[1]=y;}
		void center_refresh();
		void setup_tcoef(double C) {tension_coef = C;}
		void setup_pcoef(double C) {pressure_coef = C;}
		void setup_mcoef(double M) {mass_coef = M;}
		void set_area();
		void set_original_area(double A) {area = A;}
		void set_area_p(int i, double v) {area_p[i] = v;}
		void set_initial_area(double A) {initial_area = A;}
		double shape_based_division_angle();
		double get_area() {return area;}
		double get_area_p(int i) {return area_p[i];}
		double get_initial_area() {return initial_area;}
		double get_tcoef() {return tension_coef;}
		double get_pcoef() {return pressure_coef;}
		double get_mass() {return mass_coef;}
		double get_center(int i) {return center[i];} 
		int get_cell_state(int i) {return state[i];}
		void push_cell_state(int v) {state.push_back(v);}
		int get_cell_statelen() {return state.size();}
		int get_cell_status() {return cell_status;}
		void set_cell_state(unsigned int i, int v) { if (i<state.size()) {state[i]=v;} else {state.push_back(v);}}
		void set_cell_status(int i) {cell_status = i;}
		int get_cell_type() {return cell_type;}    // INTERFACE: cell type
		int get_life() {return life;}			   // INTERFACE: cell age
		double get_color(int i) {return color[i];} // i=0,1,2
		double get_diffusion_color(int i) {return diffusion_color[i];}
		double get_lame(int i) {return lame[i];}
		void set_lame(int i, double v) {lame[i] = v;}
		void set_cell_color(double r, double g, double b) {color[0] = r;color[1] = g;color[2] = b;}
		void set_cell_diffusion_color(double r, double g, double b) {diffusion_color[0] = r;diffusion_color[1] = g;diffusion_color[2] = b;}
		double get_GR() {return growth_rate;}    // INTERFACE: growth rate
		void set_GR(double v) {growth_rate = v;} // INTERFACE: growth rate
		void set_cell_type(int i);				 // INTERFACE: cell type
		void set_life(int i) {life = i;}		 // INTERFACE: cell age
		void set_AABB(point3D*xmin, point3D*ymin, point3D*xmax, point3D*ymax) {AABB[0]=xmin;AABB[1]=ymin;AABB[2]=xmax;AABB[3]=ymax;}
		void set_AABB_s(point3D*p,int i) {AABB[i]=p;}
		void refresh_AABB();
		point3D* get_AABB(int i) {return AABB[i];}
		point3D* get_migrate_p(int i) {return migrate_p[i];}
		double R() {return r;}
		bool Burry() {return burried;}
		bool Circle() {return circle;}
		bool Dead() {return dead;}
		bool Set_Dead() {return set_dead;}
		bool Migrate() {return migrate;}
		bool get_migrate_relax() {return migrate_relax;}
		bool get_sense_mech() {return sense_mech;}
		bool get_sense_migr() {return sense_migr;}
		double get_sense_mech_vec(int i) {return sense_mech_vec[i];}
		double get_sense_migr_vec(int i) {return sense_migr_vec[i];}
		double get_sense_migr_angle() {return sense_migr_angle;}
		int Mark_number() {return mark_number;}
		bool Soften() {return soften;}
		bool Migrate_response() {return migrate_response;}
		void setup_Burry();
		void setup_Dead(bool D) {dead = D;}
		void set_set_dead(bool S) {set_dead = S;} // INTERFACE: cell death
		void set_migrate(bool M) {migrate = M;}   // INTERFACE: set migratable
		void set_migrate_relax(bool M) {migrate_relax = M;}
		void set_sense_mech(bool S) {sense_mech = S;}
		void set_sense_migr(bool S) {sense_migr = S;}
		void set_sense_mech_vec(double x, double y) {sense_mech_vec[0] = x;sense_mech_vec[1] = y;}
		void set_sense_migr_vec(double x, double y) {sense_migr_vec[0] = x;sense_migr_vec[1] = y;}
		void set_sense_migr_angle(double a) {sense_migr_angle = a;}
		void setup_Circle();
		void setup_Mark_number(int S) {mark_number = S;}
		void set_soften(bool S) {soften = S;}
		void set_migrate_response(bool M) {migrate_response = M;}
		double get_migrate_angle() {return migrate_angle;}
		void set_migrate_angle(double v) {migrate_angle = v;}
		void set_migrate_p(point3D *sp, point3D *tp) {migrate_p[0] = sp;migrate_p[1] = tp;}
		vector<edge*> *get_sides() {return &sides;}
		vector<triangle*> *get_inner_t() {return &inner_t;}
		vector<point3D*> *get_inner_p() {return &inner_p;}
		vector<point3D*> *get_inner_p_slip() {return &inner_p_slip;}
		int get_inner_p_n() {return inner_p_n;}
		void set_inner_p_n() {inner_p_n = (int)inner_p.size();}
		void clear_inner_t();
		void clear_inner_p();
		void clear_inner_p_slip();
		void pushtriangle(triangle *Tria) {inner_t.push_back(Tria);}
		void pushinnerpoi(point3D *Poin) {inner_p.push_back(Poin);}
		void pushinnerpoislip(point3D *Poin) {inner_p_slip.push_back(Poin);}
		void pushside(edge *Edge) {sides.push_back(Edge);}
		void removeside(edge *Edge);
		vector<interpair*> *get_pairs() {return &pairs;}
		void pushpair(interpair *Pair) {pairs.push_back(Pair);}
		void removepair(interpair *Pair); 
		bool get_interior_refresh() {return interior_refresh;}
		void set_interior_refresh(bool valid) {interior_refresh = valid;}
		int get_interior_time() {return interior_time;}
		void set_interior_time(int time) {interior_time = time;}
		bool test_point_inside(double x, double y);
		vector<rear_fc*> *get_migrate_filling_cells() {return &migrate_dead_filling_cells;}
		void push_filling_cell(rear_fc *FC) {migrate_dead_filling_cells.push_back(FC);}
		void remove_filling_cell(cell *DC);
		void clear_filling_cells();
		int get_behavior() {return behavior;}
		void set_behavior(int i) {behavior = i;}
		int get_bhsteps(int i) {return bhsteps[i];}
		void set_bhsteps(int i, int steps) {bhsteps[i] = steps;}

	friend class edge;
	friend class point3D;
	friend class triangle;
	friend class interpair;
	friend class rear_fc;
};

class mergepair
{
	private:
		int index;
		/************************
		         t1  s2
		     o--<--oo--<--o
				   ||
			  C1   ||   C2
				   ||
				   ||
			 o-->--oo-->--o
			     s1  t2
		************************/
		point3D* s1;
		point3D* t1; // from s1 to t1: ccw
		point3D* s2;
		point3D* t2; // from s2 to t2: ccw
		double F_s1;
		double F_t1;
		double F_s2;
		double F_t2; // stored force against adhesion to break edge apart
		double Fr_s1;
		double Fr_t1; // stored friction force to slide edge
		int Fr_dir;   // 1:fp:s1->t1, 0:rp:s1<-t1

	public:
		mergepair(int n, point3D* A, point3D* B, point3D* C, point3D* D);
		~mergepair();
		point3D* S1() {return s1;}
		point3D* T1() {return t1;}
		point3D* S2() {return s2;}
		point3D* T2() {return t2;}
		double get_F_s1() {return F_s1;}
		double get_F_t1() {return F_t1;}
		double get_F_s2() {return F_s2;}
		double get_F_t2() {return F_t2;}
		double get_Fr_s1() {return Fr_s1;}
		double get_Fr_t1() {return Fr_t1;}
		int get_Fr_dir() {return Fr_dir;}
		void set_F_s1(double fs1) {F_s1 = fs1;}
		void set_F_t1(double ft1) {F_t1 = ft1;}
		void set_F_s2(double fs2) {F_s2 = fs2;}
		void set_F_t2(double ft2) {F_t2 = ft2;}
		void set_Fr_s1(double frs1) {Fr_s1 = frs1;}
		void set_Fr_t1(double frt1) {Fr_t1 = frt1;}
		void set_Fr_dir(int d) {Fr_dir = d;}
		void setS1(point3D* A) {s1 = A;}
		void setT1(point3D* B) {t1 = B;}
		void setS2(point3D* C) {s2 = C;}
		void setT2(point3D* D) {t2 = D;}
		
	friend class point3D;
	friend class interpair;
};

class interpair
{
	private:
		int index;
		cell* I1;
		cell* I2;
		bool attach;
		bool redundant;      // the pair is redundant
		double adhesion;     // adhesion between cells
		double friction;     // friction between cells
		vector<mergepair*> P;

	public:
		// Constructor/Destructor
		interpair(int n, cell* A, cell* B);
		~interpair();

		//Methods
		int id() {return index;}
		cell* get_I1() {return I1;}
		cell* get_I2() {return I2;}
		bool get_attach() {return attach;}
		bool get_redundant() {return redundant;}
		bool check_detection();                  // if the boundary is burried, there is no need to do detection
		double get_adhesion() {return adhesion;} // INTERFACE: adhesion
		double get_friction() {return friction;} // INTERFACE: friction
		void set_adhesion(double A) {adhesion = A;} // INTERFACE: adhesion [0, 20]
		void set_friction(double F) {friction = F;} // INTERFACE: friction [0, 20]
		vector<mergepair*> *MP() {return &P;}
		void set_attach(bool valid) {attach = valid;}
		void set_redundant(bool valid) {redundant = valid;}
		void push_MP(mergepair* MP) {P.push_back(MP);}
		void remove_MP(mergepair* MP);
		void clear_MP() {P.clear();}

	friend class edge;
	friend class cell;
	friend class point3D;
	friend class mergepair;
};

class node // for global stiffness matrix
{
	private:
		int index;
		int sm_id; // the index in stiffmatrix
		bool in_global;
		bool joint;
		vector<point3D*> P;

	public:
		// Constructor/Destructor
		node(int n);
		~node();

		// Method
		int id() {return index;}
		int SM_ID() {return sm_id;}
		inline void set_SM_ID(int N) {sm_id = N;}
		inline bool get_in_global() {return in_global;}
		inline bool get_joint() {return joint;}
		inline void set_in_global(bool valid) {in_global = valid;}
		inline void set_joint(bool valid) {joint = valid;}
		void pushvertex(point3D* A);
		void removevertex(point3D* A);
		inline vector<point3D*> *getvertex() {return &P;}

	friend class point3D;
	friend class le_sm;
};

class le_sm
{
	private:
		int index;
		point3D* N1;
		point3D* N2;
		double stiffness_matrix_12[4];
		double stiffness_matrix_21[4];

	public:
		le_sm(int n);
		~le_sm();

		int id() {return index;}
		point3D* getN1() {return N1;}
		point3D* getN2() {return N2;}
		inline void setN1(point3D* n1) {N1 = n1;}
		inline void setN2(point3D* n2) {N2 = n2;}
		inline double get_stiffness_matrix_12(int i) {return stiffness_matrix_12[i];}
		inline double get_stiffness_matrix_21(int i) {return stiffness_matrix_21[i];}
		inline void set_sm_12(double v1, double v2, double v3, double v4) {stiffness_matrix_12[0] = v1;stiffness_matrix_12[1] = v2;stiffness_matrix_12[2] = v3;stiffness_matrix_12[3] = v4;}
		inline void set_sm_21(double v1, double v2, double v3, double v4) {stiffness_matrix_21[0] = v1;stiffness_matrix_21[1] = v2;stiffness_matrix_21[2] = v3;stiffness_matrix_21[3] = v4;}
	friend class point3D;
	friend class node;
	friend class triangle;
};

class rear_fc
{
	private:
		int index;
		double angle;
		double center[2];
		bool initial;
		cell* C;
		vector<cell*> NCs;

	public:
		rear_fc(int n);
		~rear_fc();

		int id() {return index;}
		double get_angle() {return angle;}
		double get_center(int i) {return center[i];}
		bool get_initial() {return initial;}
		void set_initial(bool V) {initial = V;}
		cell* getC() {return C;}
		vector<cell*> *get_NCs() {return &NCs;}
		void set_center(double x, double y) {center[0] = x;center[1] = y;}
		void set_angle(double A) {angle = A;}
		void set_C(cell *C1) {C = C1;}
		void push_NCs(cell *C1);
		void remove_NCs(cell *C1);
		void clear_NCs();

	friend class cell;
};

#endif
