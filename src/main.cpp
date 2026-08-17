// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     main.cpp            ***
 ***    Author:   Jieling Zhao        ***
 ***                                  ***
 ***    Created on March 19, 2012     ***
 ***************************************/

#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <string>
#include "cellproject.h"
#include <common/dycelfem_config.h>

#define WIDTH  800
#define HEIGHT 600

using namespace std;

cellproject *Project = NULL;

void init()												{Project->init();}
void display()											{Project->display();}
void keyboard(unsigned char key, int x, int y)			{Project->keyboard(key, x, y);}
void specialKey(int key, int x, int y)					{Project->specialKey(key, x, y);}
void mouseButton(int button, int state, int x, int y)	{Project->mouseButton(button, state, x, y);}
void reshape(int w, int h)                              {Project->reshape(w,h);}
void motion(int x, int y)                               {Project->motion(x,y);}
//void idleFunc()                                         {Project->idleFunc();}
void timerFunction(int value)							{Project->timerFunc(value);}

int main(int argc, char** argv) {
	/* check command line arguments for database file name */
	if(argc<2) {
		cerr << " ****************************************************\n"
			 << " **                                                **\n"
			 << " **         CELL GROWTH DYNAMIC SIMULATION         **\n"
			 << " **                 (DyCelFEM)                     **\n"
			 << " ****************************************************\n\n"
			 << " USAGE: dycelfem <run.cfg>          (recommended)\n"
			 << "        dycelfem <celldatabase>     (legacy, uses built-in defaults)\n"
			 << "-----------------------------------------------------\n"
			 << " A run.cfg describes the whole simulation: input file, SBML models,\n"
			 << " number of steps, view and output settings. See examples/ and\n"
			 << " run 'dycelfem --write-default-config <file>' for an annotated template.\n\n"
			 << " Input file should include at least one cell position coordinate.\n"
			 << " The content should follow the format of the sample file.\n\n"
			 << " Keyboard and mouse input control:\n"
			 << "  -> Press arrow keys or use mouse to shift images around the screen\n"
			 << "  -> Press 'i' to zoom in\n"
			 << "  -> Press 'o' to zoom out\n"
			 << "  -> Press 's' to start, pause and resume the simulation\n"
			 << "  -> Press 'q' to quit the simulation\n"
			 << "  -> Press 'p' to print out the screen\n"
			 << "  -> Press 'c' to print out the point-based inputfile\n"
			 << "  -> Press '0' to swtich to solid mode\n"
			 << "  -> Press '1' to swtich to inner-mesh mode\n"
			 << "  -> Press '2' to swtich between molecule diffusion mode\n"
			 << "  -> Press '3' to active/deactive cell index\n"
			 << "  -> Press '4' to active/deactive migration mark\n"
			 << "  -> Press '5' to active/deactive adhesion degree\n"
			 <<endl;
		return 0; 
	}

	/* --write-default-config <file>: emit an annotated template and exit. */
	if (string(argv[1]) == "--write-default-config")
	{
		const char *out = (argc > 2) ? argv[2] : "run.cfg";
		if (cellproject::writeDefaultConfig(out)) {
			cout << "Wrote annotated configuration template to " << out << endl;
			return 0;
		}
		cerr << "ERROR: could not write " << out << endl;
		return 1;
	}

	/* argv[1] is either a run configuration or, for backward compatibility,
	   a bare cell database. */
	string database;
	if (dycfg::looksLikeConfig(argv[1]))
	{
		dycfg::load(argv[1]);
		dycfg::dumpEffective();
		database = dycfg::resolve(dycfg::getStr("input", ""));
		if (database.empty()) {
			cerr << "ERROR: " << argv[1] << " does not define 'input'." << endl;
			return 1;
		}
	}
	else
	{
		database = argv[1];
		cout << "No configuration file given; using built-in defaults." << endl;
	}

	if (!cellproject::prepareOutputDir()) return 1;

	int win_w = dycfg::getInt("view.width",  WIDTH);
	int win_h = dycfg::getInt("view.height", HEIGHT);

	/* create project and read input file */
	Project = new cellproject(win_w, win_h);
	Project->readDatabase(const_cast<char*>(database.c_str()));
	/* deal with any GLUT command Line options */
	glutInit(&argc, argv);
	/* create an output window */
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	// set up the physical window size
	glutInitWindowSize(win_w, win_h);
	// set the name of the window and try to create it
	glutCreateWindow("Dynamic Finite Element Modeling of Cell Growth");
	/* set up the logical graphics space */
	init();
	/* Receive keyboard inputs */
	glutKeyboardFunc(keyboard);
	/* Receive keyboard inputs */
	glutSpecialFunc(specialKey);
	/* Receive mouse button input*/
	glutMouseFunc(mouseButton);
	/* assign the display function */
	glutDisplayFunc(display);
	/* assign the idle function */
	glutReshapeFunc(reshape);
	/* reshape the window size */
	glutMotionFunc(motion);
	/* set the motion callback */
	//glutIdleFunc(idleFunc);
	glutTimerFunc(0,timerFunction,0);

	glutMainLoop();
	return (0);
}

#endif


