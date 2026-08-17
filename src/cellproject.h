// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/***************************************
 ***    Project:  Cell Growth         ***
 ***    File:     cellproject.h       ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef CELLPROJECT_H
#define CELLPROJECT_H


#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "dbReader.h"
#include "util.h"

//#ifdef WINDOWS
#ifndef _WIN32
    #include <unistd.h>
    #define GetCurrentDir getcwd
	#define fopen_s(A,B,C) *A=fopen(B,C)
#endif
//#ifndef WINDOWS
#ifdef _WIN32
    #include <direct.h>
    #define GetCurrentDir _getcwd
#endif

#ifdef __APPLE__
	#include <GLUT/glut.h>
#endif
#ifndef __APPLE__
	#include <GL/glut.h>
#endif

using namespace std;

class cellproject
{
	private:
		// Members
		int width;
		int height;
		bool pause;
		int view_type;
		bool show_index;     // show cell index
		bool show_text;      // show molecule name
		bool show_direction; // show migration direction
		bool show_force;     // show adhesion
		int diffusion_number;
		int number;        // the time step
		int mouse;         // 0:null, 1:left, 2:right
		int mouses;        // 0:release, 1:press
		int mousec[2];     // window coordinates
		int mousec_tmp[2]; // window coordinates at motion
		dbReader *reader;
		
	public:
		// Constructor/Destructor
		cellproject(int w, int h);
		~cellproject() {};

		// Configuration helpers
		static string outputDir();                                   // config output.dir, default "out"
		static bool   writeDefaultConfig(const char *path);           // annotated template
		static bool   prepareOutputDir();                             // create output.dir up front
		static void   buildBMPHeader(unsigned char *h, int w, int hgt);
		void          autofitView();                                  // fit view to the tissue

		// OpenGL Methods
		void init();
		void display();
		void keyboard(unsigned char key, int x, int y);
		void specialKey(int key, int x, int y);
		void mouseButton(int button, int state, int x, int y);
		void reshape(int w, int h);
		void motion(int x, int y);
		//void idleFunc();
		void timerFunc(int value);

		// Project Methods
		dbReader* get_reader();
		/*** parameter input ***/
		void readDatabase(char* fileName);
		/*** cell growth ***/
		int collision_detection();
		int cell_biology();
		int cell_physics();
		int cell_edge_length();
		int cell_pick_check();
		void cell_inner_mesh();
		void cell_inner_mesh_debug();
		void cell_node_test();
		void cell_angle();
		int get_number() {return number;}
		void matlab_format();
		void cell_file_format();
		/*** visualization ***/
		void visualmatrix(double M[3][3]);
		void drawinnermesh(double M[3][3]);
		void drawtriangles(double M[3][3]);
		void drawindex(double M[3][3]);
		void drawforcedegree(double M[3][3]);
		void drawNumber(point2D *position, int number, double color[3]);
		void drawText(point2D *position, string *Text_str, double color[3]);
		//void printoutscreen(int number, int id);
		void printoutscreen(int number, int id, int viewtype, int diffnum);
		void drawBresenhamCircle(point2D *center,int radius, double color[3]);
		void Dataprojection();
};

#endif

