// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/***************************************
***    Project:  Cell Growth         ***
***    File:     cellproject.cpp     ***
***    Author:   Jieling Zhao        ***
***                                  ***
***    Created on March 19, 2012     ***
***************************************/

#ifndef CELLPROJECT_CPP
#define CELLPROJECT_CPP

#include "cellproject.h"
#include "biology.h"
#include <common/dycelfem_config.h>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>

using namespace std;

/* Where the .bmp frames go, relative to the working directory. On the original
   Windows setup this directory was created by hand before a run; when missing,
   fopen() returned NULL and the old code walked straight into fwrite() and
   segfaulted. It is now created on demand. */
string cellproject::outputDir()
{
	return dycfg::getStr("output.dir", "out");
}

static bool dycelfem_ensure_dir(const string &dir)
{
	struct stat st;
	if (stat(dir.c_str(), &st) == 0) return true;
#ifdef _WIN32
	int rc = _mkdir(dir.c_str());
#else
	int rc = mkdir(dir.c_str(), 0755);
#endif
	return (rc == 0 || stat(dir.c_str(), &st) == 0);
}

/* Called once at startup, before anything writes: the per-step printout files
   are written by dbReader and the frames by printoutscreen(), and whichever
   ran first used to fail on a missing directory. */
bool cellproject::prepareOutputDir()
{
	string dir = outputDir();
	if (dir.empty() || dir == ".") return true;
	if (dycelfem_ensure_dir(dir)) return true;
	cerr << "ERROR: cannot create output directory '" << dir << "'." << endl;
	cerr << "       Check output.dir in the configuration and the permissions "
	     << "of the working directory." << endl;
	return false;
}

/* Number of steps to run. DYCELFEM_MAXSTEPS is still honoured for the older
   scripts that set it. */
static int dycelfem_maxsteps()
{
	const char *legacy = getenv("DYCELFEM_MAXSTEPS");
	if (legacy != NULL && legacy[0] != '\0')
	{
		int n = atoi(legacy);
		if (n > 0) return n;
	}
	int n = dycfg::getInt("run.steps", 400);
	return (n > 0) ? n : 400;
}

/* Annotated configuration template, written by --write-default-config. */
bool cellproject::writeDefaultConfig(const char *path)
{
	FILE *f = NULL;
	fopen_s(&f, path, "w");
	if (f == NULL) return false;
	fprintf(f,
"# ===========================================================================\n"
"#  DyCelFEM run configuration\n"
"#\n"
"#  Every setting can also be given as an environment variable, which wins:\n"
"#  key 'run.steps' -> DYCELFEM_RUN_STEPS. Relative paths are resolved against\n"
"#  the directory holding THIS file, so a run directory can live anywhere.\n"
"#\n"
"#     dycelfem run.cfg\n"
"# ===========================================================================\n"
"\n"
"# Cell database describing the initial tissue.\n"
"input = init7140.txt.trimmed.skinf10.wound200\n"
"\n"
"[sbml]\n"
"# Reaction networks, one per cell type, plus the shared 'common' model that\n"
"# carries the diffusion constants (D_*), D0 and DeltaT.\n"
"common       = mcommon.xml\n"
"macrophage   = mmacrophage.xml\n"
"fibroblast   = mfibroblast.xml\n"
"keratinocyte = mkeratinocyte.xml\n"
"ecm          = mecm.xml\n"
"bm           = mbm.xml\n"
"clot         = mclot.xml\n"
"\n"
"[run]\n"
"# Number of simulation steps. One step is DeltaT minutes (60 in mcommon.xml),\n"
"# so 300 steps ~ 12.5 days of wound healing.\n"
"steps = 300\n"
"\n"
"# Start simulating immediately instead of waiting for the 's' key.\n"
"# 's' still pauses and resumes an autostarted run.\n"
"autostart = 1\n"
"\n"
"# Quit when the step count is reached, so a batch pipeline can continue.\n"
"# Set to 0 to leave the window open for inspection.\n"
"exit_when_done = 1\n"
"\n"
"# 0 = INITIALIZATION (prepare/relax a tissue; no chemistry at all)\n"
"# 1 = MODELING       (the actual simulation: reactions + diffusion)\n"
"# 2 = DELETING AND RELAXING CELLS\n"
"mode = 1\n"
"\n"
"# Inter-cell diffusion coefficient.\n"
"#   edgelen        D*D0*commedgelen        reproduces the published results\n"
"#   edgelen_distsq D*D0*commedgelen/distsq\n"
"#   distsq         D*D0/distsq             needs D0 raised by ~3e3 to match\n"
"diffusion = edgelen\n"
"\n"
"# 0 = seed the RNG from the clock (original behaviour, not reproducible)\n"
"# any positive integer = deterministic run\n"
"seed = 0\n"
"\n"
"[output]\n"
"# Directory for printout*.txt and the .bmp frames. Created if missing.\n"
"dir = out\n"
"\n"
"# Write a printout / an image every N steps (0 disables).\n"
"data_every  = 1\n"
"image_every = 1\n"
"\n"
"# Also render one image per molecular species each time an image is written.\n"
"species_images = 1\n"
"\n"
"# full   - include VERT/EDGE/NODE lists; the file can be used to restart a run\n"
"# simple - CELL_LIST and PAIR_LIST only; ~7x smaller, all the analysis needs\n"
"sections = full\n"
"\n"
"# Optional 54-byte BMP header template. Leave empty to generate a correct\n"
"# header (recommended); the old code required a hand-made reference.bmp.\n"
"bmp_reference =\n"
"\n"
"[view]\n"
"width  = 800\n"
"height = 600\n"
"\n"
"# Fit the view to the tissue instead of using the u/v window baked into the\n"
"# input file (which frames a small patch far from the wound).\n"
"autofit = 1\n"
"margin  = 0.03\n"
"\n"
"# Zoom about the view centre once it has been fitted. 1 = the whole tissue is\n"
"# guaranteed visible (aspect padding can leave empty space above and below);\n"
"# 2 = twice as close, with the edges cropped. Use view.center_x/_y to move the\n"
"# centre, e.g. onto the wound at 0,0.\n"
"zoom = 1.0\n"
"#center_x = 0\n"
"#center_y = 0\n"
"\n"
"# Draw the numeric cell index on every cell. Useful for debugging, unreadable\n"
"# on a 7000-cell tissue.\n"
"show_index = 0\n");
	fclose(f);
	return true;
}

/********************************
***     OPENGL FUNCTIONS      ***
********************************/

void timerFunction(int value);

// constructor
cellproject::cellproject(int w, int h)
{
	width = w;
	height = h;
	/* Start running immediately unless the config says otherwise; pressing 's'
	   still pauses and resumes. */
	pause = !dycfg::getBool("run.autostart", true);
	view_type = 0;
	/* Cell index labels are useful for debugging but clutter every frame of a
	   7000-cell tissue, so they are off unless asked for. */
	show_index = dycfg::getBool("view.show_index", false);
	show_text = false;
	show_direction = true;
	show_force = false;
	diffusion_number = 0;
	number = 0;
	mouse = 0;
	mouses = 0;
	mousec[0] = 0;
	mousec[1] = 0;
	mousec_tmp[0] = 0;
	mousec_tmp[1] = 0;
	reader = NULL;
}

// init function - intitializes OpenGL window
void cellproject::init() {
	glMatrixMode(GL_PROJECTION); // Matrix for projection transformation;
	glLoadIdentity();
	gluOrtho2D(0, width, 0, height); // Define a 2d orthographic projection matrix

	glClearColor(1.0, 1.0, 1.0, 1.0); // set background color to white
	glShadeModel(GL_SMOOTH); // specify color model
}

// display function - draw all pixels to screen
void cellproject::display(void) {	
	/* clear the screen to the clear color */
	glClearColor(1.0, 1.0, 1.0, 1.0); // set background color to white
	glClear(GL_COLOR_BUFFER_BIT);

	/* draw everything */
	Dataprojection();

	/* swap buffers */	
	glutSwapBuffers();
}

// keyboard function - called when user has keyboard input
void cellproject::keyboard(unsigned char key, int x, int y)
{
	switch (key) 
	{
	case 'q':
	case 'Q':
		cout << "Quit the program!" <<endl;
		exit(0);
		break;
	case 'o':
	case 'O':
		{
			//cout<< "Zoom out!" <<endl;
			double x_min = reader->getViewVolume()->getUMin()*1.05;
			double y_min = reader->getViewVolume()->getVMin()*1.05;
			double x_max = reader->getViewVolume()->getUMax()*1.05;
			double y_max = reader->getViewVolume()->getVMax()*1.05;
			reader->getViewVolume()->setUMin(x_min);
			reader->getViewVolume()->setVMin(y_min);
			reader->getViewVolume()->setUMax(x_max);
			reader->getViewVolume()->setVMax(y_max);
		}
		break;
	case 'i':
	case 'I':
		{
			//cout<< "Zoom in!" <<endl;
			double x_min = reader->getViewVolume()->getUMin()/1.05;
			double y_min = reader->getViewVolume()->getVMin()/1.05;
			double x_max = reader->getViewVolume()->getUMax()/1.05;
			double y_max = reader->getViewVolume()->getVMax()/1.05;
			reader->getViewVolume()->setUMin(x_min);
			reader->getViewVolume()->setVMin(y_min);
			reader->getViewVolume()->setUMax(x_max);
			reader->getViewVolume()->setVMax(y_max);
		}
		break;
	case 's':
	case 'S':
		{
			pause = !pause;
			if (pause)
			{
				cout <<"Growing stopped!"<<endl;
			}
			else 
			{
				cout <<"Resume growing!"<<endl;
			}
		}
		break;
	case 'p':
	case 'P':
		{
			//printoutscreen(number,-1);
			printoutscreen(number, -1, view_type, diffusion_number);
		}
		break;
	case 'W':
	case 'w':
		{
			matlab_format();
		}
		break;
	case 'C':
	case 'c':
		{
			cell_file_format();
		}
		break;
	case '0':
		{
			view_type = 0;
			show_text = false;
		}
		break;
	case '1':
		{
			view_type = 1;
			show_text = false;
		}
		break;
	case '2':
		{
			view_type = 2;
			show_text = true;
			diffusion_number++;
			int dn = reader->get_biology()->getNumRealSpecies();
			diffusion_number = diffusion_number%dn;
		}
		break;
	case '3':
		{
			show_index = !show_index;
			show_text = false;
		}
		break;
	case '4':
		{
			show_direction = !show_direction;
		}
		break;
	case '5':
		{
			show_force = !show_force;
		}
		break;
	}
}

// special keyboard function - called when user has special keys input
void cellproject::specialKey(int key, int x, int y) {
	switch (key) 
	{
	case GLUT_KEY_LEFT:
		{
			//cout << "Move left" <<endl;
			double x_vrp = reader->getViewVolume()->getVRP()->x() + 1;
			double y_vrp = reader->getViewVolume()->getVRP()->y();
			reader->getViewVolume()->getVRP()->setx(x_vrp);
			reader->getViewVolume()->getVRP()->sety(y_vrp);
		}
		break;
	case GLUT_KEY_RIGHT:
		{
			//cout << "Move right" <<endl;
			double x_vrp = reader->getViewVolume()->getVRP()->x() - 1;
			double y_vrp = reader->getViewVolume()->getVRP()->y();
			reader->getViewVolume()->getVRP()->setx(x_vrp);
			reader->getViewVolume()->getVRP()->sety(y_vrp);
		}
		break;
	case GLUT_KEY_UP:
		{
			//cout << "Move upward" <<endl;
			double x_vrp = reader->getViewVolume()->getVRP()->x();
			double y_vrp = reader->getViewVolume()->getVRP()->y() - 1;
			reader->getViewVolume()->getVRP()->setx(x_vrp);
			reader->getViewVolume()->getVRP()->sety(y_vrp);
		}
		break;
	case GLUT_KEY_DOWN:
		{
			//cout << "Move downward" <<endl;
			double x_vrp = reader->getViewVolume()->getVRP()->x();
			double y_vrp = reader->getViewVolume()->getVRP()->y() + 1;
			reader->getViewVolume()->getVRP()->setx(x_vrp);
			reader->getViewVolume()->getVRP()->sety(y_vrp);
		}
		break;
	}
}

// mouse button function - called when user clicks any mouse button
void cellproject::mouseButton(int button, int state, int x, int y) {
	switch(button){
		case GLUT_LEFT_BUTTON:
		{
			if (state==GLUT_DOWN)
			{
				mouse = 1;
				mouses = 1;
				//cout<<"mouse left button pressed at: "<<x<<" "<<y<<endl;
				mousec[0] = x;
				mousec[1] = y;
				mousec_tmp[0] = mousec[0];
				mousec_tmp[1] = mousec[1];
			}
			if (state==GLUT_UP)
			{
				mouse = 0;
				mouses = 0;
				//cout<<"mouse left button released at: "<<x<<" "<<y<<endl;
			}
		}
		break;
		case GLUT_RIGHT_BUTTON:
		{
			if (state==GLUT_DOWN)
			{
				mouse = 2;
				mouses = 1;
				//cout<<"mouse right button pressed at: "<<x<<" "<<y<<endl;
				mousec[0] = x;
				mousec[1] = y;
				mousec_tmp[0] = mousec[0];
				mousec_tmp[1] = mousec[1];
			}
			if (state==GLUT_UP)
			{
				mouse = 0;
				mouses = 0;
				//cout<<"mouse right button released at: "<<x<<" "<<y<<endl;
			}
		}
		break;
	}
}

// idle function - called when no other function is running
/*void cellproject::idleFunc() {
	glutPostRedisplay(); // redraw display function
}*/

// reshape function - reshape the viewport screen
void cellproject::reshape(int w, int h)
{
	glViewport(0,0,(GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, (GLdouble) w, 0.0, (GLdouble) h);
	width = w;
	height = h;
}

// mouse button function - mouse track
void cellproject::motion(int x, int y)
{
	if      (mouse==1 && mouses==1)
	{
		//cout<<"mouse left button placed at: "<<x<<" "<<y<<endl;
		mousec_tmp[0] = x;
		mousec_tmp[1] = y;
		double x_scale = (reader->getViewVolume()->getUMax()-reader->getViewVolume()->getUMin())/width;
		double y_scale = (reader->getViewVolume()->getVMax()-reader->getViewVolume()->getVMin())/height;
		double x_vrp = reader->getViewVolume()->getVRP()->x() - x_scale*(x - mousec[0]);
		double y_vrp = reader->getViewVolume()->getVRP()->y() + y_scale*(y - mousec[1]);
		reader->getViewVolume()->getVRP()->setx(x_vrp);
		reader->getViewVolume()->getVRP()->sety(y_vrp);
		mousec[0] = x;
		mousec[1] = y;
	}
	else if (mouse == 2 && mouses==1)
	{
		//cout<<"mouse right button placed at: "<<x<<" "<<y<<endl;
		mousec_tmp[0] = x;
		mousec_tmp[1] = y;
	}
}

// timer function
void cellproject::timerFunc(int value)
{
	if (!pause)
	{
		if (number<dycelfem_maxsteps())
		{
			number++;
			cout<<"Step "<<number<<endl;
			cell_inner_mesh();
			int error_id = -1;
			int cb = cell_biology();
			if (cb==-1)
			{
				int cp = cell_physics();
				if (cp==-1)
				{
					int cd = collision_detection();
					if (cd==-1)
					{
						int cpc = cell_pick_check();
						if (cpc==-1)
						{
							int cel = cell_edge_length();
							if (cel>=0) {error_id = cel;}
						}
						else {error_id = cpc;}
					}
					else {error_id = cd;}
				}
				else {error_id = cp;}
				//cell_node_test(); // not needed, only for test
				//if (number==29) {cell_file_format();} // not needed, only for test
			}
			else {error_id = cb;}
			cout<<endl;
			/** print out pics: one reference bmp file is needed in advance **/
			if (error_id==-1)
			{
				int data_every  = dycfg::getInt("output.data_every", 1);
				int image_every = dycfg::getInt("output.image_every", 1);
				if (data_every > 0 && number % data_every == 0)
				{
					// Output cell file
					cell_file_format();
				}
				if (image_every > 0 && number % image_every == 0)
				{
					// Output images
					int vtype = view_type;
					int dnumber = diffusion_number;
					bool stext = show_text;
					view_type = 0;
					glClearColor(1.0, 1.0, 1.0, 1.0); // set background color to white
					glClear(GL_COLOR_BUFFER_BIT);
					Dataprojection();
					glutSwapBuffers();
					//printoutscreen(number,-1);
					printoutscreen(number, -1, 0, 0);
					if (dycfg::getBool("output.species_images", true))
					{
						view_type = 2;
						show_text = true;
						int dn = reader->get_biology()->getNumRealSpecies();
						//cout << "    ##### SPECIES: " << dn << endl;
						for (int i=0;i<dn;i++)
						{
							diffusion_number = i;
							glClearColor(1.0, 1.0, 1.0, 1.0); // set background color to white
							glClear(GL_COLOR_BUFFER_BIT);
							Dataprojection();
							glutSwapBuffers();
							//printoutscreen(number,-1);
							printoutscreen(number, -1, 2, i);
						}
					}
					view_type = vtype;
					diffusion_number = dnumber;
					show_text = stext;
				}
			}
			else
			{
				//printoutscreen(number,error_id);
				printoutscreen(number, error_id, view_type, diffusion_number);
				cout<<"  -> Bug detected! Please check out cell "<<error_id<<" at output"<<number<<".error.bmp!"<<endl;
				exit(0);
			}
			if (number>=dycelfem_maxsteps())
			{
				pause = true;
				cout<<"Simulation finished: "<<number<<" steps."<<endl;
				/* Batch runs should terminate so a pipeline can continue; set
				   run.exit_when_done = 0 to keep the window open for inspection. */
				if (dycfg::getBool("run.exit_when_done", true))
				{
					exit(0);
				}
				cout<<"Window left open (run.exit_when_done = 0). Press 'q' to quit."<<endl;
			}
			//pause = true;
		}
	}
	glutPostRedisplay(); // redraw display function
	glutTimerFunc(16,timerFunction,value);
}

/********************************
***     PROJECT FUNCTIONS     ***
********************************/

// read database text file
dbReader* cellproject::get_reader()
{
	return reader;
}

void cellproject::readDatabase(char* fileName) {
	// free memory if previously allocated
	if(reader != NULL) delete reader;

	/* run.seed = 0 keeps the original clock-seeded, non-reproducible behaviour;
	   any positive value gives a deterministic run. */
	int seed = dycfg::getInt("run.seed", 0);
	if (seed > 0) { srand((unsigned int)seed);   cout<<"RNG seed: "<<seed<<" (deterministic)"<<endl; }
	else          { srand((unsigned int)clock()); cout<<"RNG seed: clock (not reproducible)"<<endl; }

	// create new reader
	dbReader *R = new dbReader();
	reader = R;
	reader->read(fileName,get_number());

	if (dycfg::getBool("view.autofit", true)) { autofitView(); }
}

/* Fit the viewing window to the tissue.
 *
 * visualmatrix() maps world point p to screen as
 *     screen_x = width /(uMax-uMin) * ( Rx . (p - VRP) )
 *     screen_y = height/(vMax-vMin) * ( Ry . (p - VRP) )
 * with Rx = normalised VUP and Ry = normalised VPN. So the visible region is
 * the box whose Rx-extent is [0, uMax-uMin] and Ry-extent is [0, vMax-vMin],
 * anchored at VRP. Fitting therefore means: measure the tissue extent along
 * Rx and Ry, pad it, widen the short axis so the aspect ratio matches the
 * window (otherwise the picture is stretched), then set the spans and place
 * VRP at the low corner.
 *
 * The archived input files carry u -100,100 / v -75,75 with VRP 200,100 - a
 * 200x150 world window on a tissue ~1700x950 across, which is why the default
 * view showed a tiny patch of normal skin far from the wound.
 */
void cellproject::autofitView()
{
	viewVolume *vv = reader->getViewVolume();
	if (vv == NULL || vv->getVPN() == NULL || vv->getVUP() == NULL) return;

	double dy = sqrt(vv->getVPN()->x()*vv->getVPN()->x() + vv->getVPN()->y()*vv->getVPN()->y());
	double dx = sqrt(vv->getVUP()->x()*vv->getVUP()->x() + vv->getVUP()->y()*vv->getVUP()->y());
	if (dx == 0 || dy == 0) return;
	double Rx[2] = { vv->getVUP()->x()/dx, vv->getVUP()->y()/dx };
	double Ry[2] = { vv->getVPN()->x()/dy, vv->getVPN()->y()/dy };

	bool   first = true;
	double uLo = 0, uHi = 0, vLo = 0, vHi = 0;
	for (int i = 0; i < (int)reader->V()->size(); i++)
	{
		point3D *p = reader->V()->at(i);
		if (p == NULL) continue;
		double u = Rx[0]*p->x() + Rx[1]*p->y();
		double v = Ry[0]*p->x() + Ry[1]*p->y();
		if (first) { uLo = uHi = u; vLo = vHi = v; first = false; }
		else {
			if (u < uLo) uLo = u;  if (u > uHi) uHi = u;
			if (v < vLo) vLo = v;  if (v > vHi) vHi = v;
		}
	}
	if (first) { cout<<"View autofit skipped: no vertices."<<endl; return; }

	double uSpan = uHi - uLo, vSpan = vHi - vLo;
	if (uSpan <= 0 || vSpan <= 0) return;

	double margin = dycfg::getNum("view.margin", 0.03);
	uLo -= uSpan*margin; uHi += uSpan*margin;
	vLo -= vSpan*margin; vHi += vSpan*margin;
	uSpan = uHi - uLo; vSpan = vHi - vLo;

	/* widen the short axis so pixels stay square */
	double want = (double)width/(double)height;
	if (uSpan/vSpan < want) {
		double grow = vSpan*want - uSpan;
		uLo -= grow/2; uHi += grow/2; uSpan = uHi - uLo;
	} else {
		double grow = uSpan/want - vSpan;
		vLo -= grow/2; vHi += grow/2; vSpan = vHi - vLo;
	}

	/* Aspect correction pads the short axis, which for a wide flat tissue
	   leaves a lot of empty space above and below. view.zoom scales the window
	   about its centre: 1 = whole tissue guaranteed visible, 2 = twice as
	   close (and the edges fall outside the window). view.center_x/_y move the
	   centre in world coordinates, e.g. onto the wound. */
	double cu = 0.5*(uLo + uHi);
	double cv = 0.5*(vLo + vHi);
	if (dycfg::has("view.center_x") || dycfg::has("view.center_y"))
	{
		double cx = dycfg::getNum("view.center_x", cu*Rx[0] + cv*Ry[0]);
		double cy = dycfg::getNum("view.center_y", cu*Rx[1] + cv*Ry[1]);
		cu = Rx[0]*cx + Rx[1]*cy;
		cv = Ry[0]*cx + Ry[1]*cy;
	}
	double zoom = dycfg::getNum("view.zoom", 1.0);
	if (zoom <= 0) zoom = 1.0;
	uSpan /= zoom; vSpan /= zoom;
	uLo = cu - uSpan/2;
	vLo = cv - vSpan/2;

	vv->setUMin(0); vv->setUMax(uSpan);
	vv->setVMin(0); vv->setVMax(vSpan);
	vv->setUMin_default(0); vv->setUMax_default(uSpan);
	vv->setVMin_default(0); vv->setVMax_default(vSpan);
	vv->getVRP()->setx(uLo*Rx[0] + vLo*Ry[0]);
	vv->getVRP()->sety(uLo*Rx[1] + vLo*Ry[1]);

	cout<<"View autofit: "<<uSpan<<" x "<<vSpan<<" world units (zoom "<<zoom<<"), VRP ("
	    <<vv->getVRP()->x()<<", "<<vv->getVRP()->y()<<")"<<endl;
}

int cellproject::collision_detection()
{
	int flag_cd = reader->detect_collision(get_number());
	return flag_cd;
}

int cellproject::cell_biology()
{
	int flag_cb = reader->cell_setupbiobehaviors(get_number());
	return flag_cb;
}

int cellproject::cell_physics()
{
	int value = -1;
	switch(pause) 
	{
		// pause == 0 --> resume the growing process
		case 0:
			{
				value = reader->cell_dynamic_FEM(get_number());
			}
			break;
		// pause == 1 --> freeze the growing process 
		case 1:
			{
				cout <<"Stop growing!"<<endl;
			}
			break;
	}
	return value;
}

int cellproject::cell_edge_length()
{
	int flag_cel = reader->cell_edge_response(get_number());
	return flag_cel;
}

int cellproject::cell_pick_check()
{
	int flag_cpr = reader->cell_pick_response(get_number());
	return flag_cpr;
}

void cellproject::cell_inner_mesh()
{
	reader->cell_mesh_generation(get_number());
}

void cellproject::cell_inner_mesh_debug()
{
	//reader->cell_interior_resample_debug(reader->C().at(7),number);
}

void cellproject::cell_node_test()
{
	reader->node_test();
}

void cellproject::cell_angle()
{
	reader->cell_angle_response();
}

void cellproject::matlab_format()
{
	reader->print_matlab_file();
}

void cellproject::cell_file_format()
{
	reader->print_cell_file(get_number());
}

void cellproject::drawText(point2D *position, string *Text_str, double color[3])
{
	int size = Text_str->size();
	glRasterPos2i(position->x(),position->y());
	for (int i=0; i<size; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, Text_str->at(i));
	}
}

void cellproject::printoutscreen(int number, int id, int viewtype, int diffnum)
{
	//cout << "     #### VT -- DN -- NUM: " << viewtype << " -- " << diffnum << " -- " << number;
	ostringstream temp;
	temp<<number;
	string tmps = temp.str();
	string cwd("");
	string outname = "";
	string Molecule = "";
	// Get the current working directory:
	cwd.assign(GetCurrentDir( NULL, 0 ));
	string odir = outputDir();
	if (odir.empty() || odir[0] != '/') odir = cwd + "/" + odir;
	dycelfem_ensure_dir(odir);
	if (viewtype == 0)
	{
		outname = odir+"/output_"+tmps+".bmp";
	}
	else if (viewtype==2 && show_text)
	{
		Molecule = reader->get_biology()->getSpeciesName(diffnum);
		outname = odir+"/output"+Molecule+"_"+tmps+".bmp";
	}
	//cout << "; " << outname << endl;
	glClearColor(1.0, 1.0, 1.0, 1.0); // set background color to white
	glClear(GL_COLOR_BUFFER_BIT);
	Dataprojection();
	glutSwapBuffers();
	if (id>=0)
	{
		outname = odir+"/output"+tmps+".error.bmp";
		double x_min = reader->getViewVolume()->getUMin_default();
		double y_min = reader->getViewVolume()->getVMin_default();
		double x_max = reader->getViewVolume()->getUMax_default();
		double y_max = reader->getViewVolume()->getVMax_default();
		reader->getViewVolume()->setUMin(x_min);
		reader->getViewVolume()->setVMin(y_min);
		reader->getViewVolume()->setUMax(x_max);
		reader->getViewVolume()->setVMax(y_max);
		double cx = reader->C()->at(id)->get_center(0);
		double cy = reader->C()->at(id)->get_center(1);
		//////////
		double Rx[2];
		double Ry[2];
		double dy = sqrt(reader->getViewVolume()->getVPN()->x()*reader->getViewVolume()->getVPN()->x()+reader->getViewVolume()->getVPN()->y()*reader->getViewVolume()->getVPN()->y());
		Ry[0] = reader->getViewVolume()->getVPN()->x()/dy;
		Ry[1] = reader->getViewVolume()->getVPN()->y()/dy;
		double dx = sqrt(reader->getViewVolume()->getVUP()->x()*reader->getViewVolume()->getVUP()->x()+reader->getViewVolume()->getVUP()->y()*reader->getViewVolume()->getVUP()->y());
		Rx[0] = reader->getViewVolume()->getVUP()->x()/dx;
		Rx[1] = reader->getViewVolume()->getVUP()->y()/dx;
		double scale_x = width/(x_max - x_min);
		double scale_y = height/(y_max - y_min);
		double a1 = scale_x*Rx[0];
		double a2 = scale_y*Ry[0];
		double b1 = scale_x*Rx[1];
		double b2 = scale_y*Ry[1];
		double c1 = scale_x*(Rx[0]*cx + Rx[1]*cy) - width/2;
		double c2 = scale_y*(Ry[0]*cx + Ry[1]*cy) - height/2;
		double vrp_x = (c1*b2 - c2*b1)/(a1*b2 - a2*b1);
		double vrp_y = (c1*a2 - c2*a1)/(b1*a2 - b2*a1);
		//////////
		reader->getViewVolume()->getVRP()->setx(vrp_x);
		reader->getViewVolume()->getVRP()->sety(vrp_y);
		view_type = 1;
		show_index = true;
		//////////
		glClearColor(1.0, 1.0, 1.0, 1.0); // set background color to white
		glClear(GL_COLOR_BUFFER_BIT);
		Dataprojection();
		glutSwapBuffers();
	}
	GLubyte BMP_Header[54];
	GLint viewPort[4] = {0};
	glGetIntegerv(GL_VIEWPORT,viewPort);
	GLint PixelDataLength = viewPort[2]*viewPort[3]*3;
	GLubyte * pPixelData = (GLubyte*)malloc(PixelDataLength);
	if (PixelDataLength==0) exit(0);
	FILE * pWriting = NULL;
	/* An unrecognised viewtype leaves outname empty; fopen("") then returns
	   NULL just like a missing directory would. */
	if (!outname.empty()) fopen_s(&pWriting,outname.c_str(),"wb");
	if (pWriting==NULL)
	{
		cerr<<"ERROR: cannot create image file "
		    <<(outname.empty() ? string("(no name: unhandled viewtype)") : outname)<<endl;
		free(pPixelData);
		return;
	}

	/* The original code copied the first 54 bytes of a hand-made reference.bmp
	   and patched only width/height at offset 0x12, which left bfSize and
	   biSizeImage stale in every image it ever wrote. The header is built from
	   scratch here, so no reference.bmp is needed and the size fields are
	   correct. Set output.bmp_reference to a file to restore the old
	   copy-a-template behaviour. */
	string ref = dycfg::getStr("output.bmp_reference", "");
	FILE * pRefer = NULL;
	if (!ref.empty()) fopen_s(&pRefer,dycfg::resolve(ref).c_str(),"rb");

	if (pRefer!=NULL)
	{
		if (fread(BMP_Header,sizeof(BMP_Header),1,pRefer) != 1)
		{
			cerr<<"WARNING: "<<ref<<" is shorter than a 54-byte BMP header; "
			    <<"generating one instead."<<endl;
			fclose(pRefer);
			pRefer = NULL;
		}
	}
	if (pRefer==NULL)
	{
		buildBMPHeader(BMP_Header, viewPort[2], viewPort[3]);
	}
	else
	{
		fclose(pRefer);
	}

	glReadPixels(viewPort[0],viewPort[1],viewPort[2],viewPort[3],GL_BGR_EXT,GL_UNSIGNED_BYTE,pPixelData);
	fwrite(BMP_Header,sizeof(BMP_Header),1,pWriting);
	fseek(pWriting,0x0012,SEEK_SET);
	fwrite(&viewPort[2],sizeof(viewPort[2]),1,pWriting);
	fwrite(&viewPort[3],sizeof(viewPort[3]),1,pWriting);
	fseek(pWriting,0,SEEK_END);
	fwrite(pPixelData,PixelDataLength,1,pWriting);
	fclose(pWriting);
	free(pPixelData);
}

/* 54-byte BITMAPFILEHEADER + BITMAPINFOHEADER for a bottom-up 24-bit BGR
   image, which is exactly what glReadPixels(GL_BGR_EXT) hands us. */
void cellproject::buildBMPHeader(unsigned char *h, int w, int hgt)
{
	unsigned int rowbytes  = ((unsigned int)w * 3 + 3) & ~3u;
	unsigned int imagesize = rowbytes * (unsigned int)hgt;
	unsigned int filesize  = 54 + imagesize;
	for (int i = 0; i < 54; i++) h[i] = 0;
	h[0] = 'B'; h[1] = 'M';
	h[2] = (unsigned char)( filesize        & 0xFF);
	h[3] = (unsigned char)((filesize >>  8) & 0xFF);
	h[4] = (unsigned char)((filesize >> 16) & 0xFF);
	h[5] = (unsigned char)((filesize >> 24) & 0xFF);
	h[10] = 54;                    /* pixel data offset */
	h[14] = 40;                    /* BITMAPINFOHEADER size */
	h[18] = (unsigned char)( (unsigned int)w        & 0xFF);
	h[19] = (unsigned char)(((unsigned int)w >>  8) & 0xFF);
	h[20] = (unsigned char)(((unsigned int)w >> 16) & 0xFF);
	h[21] = (unsigned char)(((unsigned int)w >> 24) & 0xFF);
	h[22] = (unsigned char)( (unsigned int)hgt        & 0xFF);
	h[23] = (unsigned char)(((unsigned int)hgt >>  8) & 0xFF);
	h[24] = (unsigned char)(((unsigned int)hgt >> 16) & 0xFF);
	h[25] = (unsigned char)(((unsigned int)hgt >> 24) & 0xFF);
	h[26] = 1;                     /* planes */
	h[28] = 24;                    /* bits per pixel */
	h[34] = (unsigned char)( imagesize        & 0xFF);
	h[35] = (unsigned char)((imagesize >>  8) & 0xFF);
	h[36] = (unsigned char)((imagesize >> 16) & 0xFF);
	h[37] = (unsigned char)((imagesize >> 24) & 0xFF);
	h[38] = 0x13; h[39] = 0x0B;    /* 2835 px/m ~ 72 dpi */
	h[42] = 0x13; h[43] = 0x0B;
}

void cellproject::drawNumber(point2D *position, int number, double color[3])
{
	if (number<10)
	{
		glColor3f(color[0],color[1],color[2]);
		glRasterPos2i(position->x(),position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+number);
	}
	else if (number>=10 && number<100)
	{
		int dig1 = number/10;
		int dig2 = number%10;
		glColor3f(color[0],color[1],color[2]);
		glRasterPos2i(position->x()-4,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig1);
		glRasterPos2i(position->x()+4,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig2);
	}
	else if (number>=100 && number<1000)
	{
		int dig1 = number/100;
		int dig2 = (number - dig1*100)/10;
		int dig3 = number%10;
		glColor3f(color[0],color[1],color[2]);
		glRasterPos2i(position->x()-8,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig1);
		glRasterPos2i(position->x(),position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig2);
		glRasterPos2i(position->x()+8,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig3);
	}
	else if (number>=1000 && number<10000)
	{
		int dig1 = number/1000;
		int dig2 = (number - dig1*1000)/100;
		int dig3 = (number - dig1*1000 - dig2*100)/10;
		int dig4 = number%10;
		glColor3f(color[0],color[1],color[2]);
		glRasterPos2i(position->x()-12,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig1);
		glRasterPos2i(position->x()-4,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig2);
		glRasterPos2i(position->x()+4,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig3);
		glRasterPos2i(position->x()+12,position->y());
		glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig4);
	}
	else
	{
		int numbert = number;
		double logupper = log(number*1.00);
		double logbase = log(10.00);
		int t_log = logupper/logbase;
		for (int i=t_log;i>=0;i--)
		{
			int base = powerint(10,i);
			int dig = numbert/base;
			numbert = numbert - base;
			int pos = 4*t_log - i*8;
			glRasterPos2i(position->x()+pos,position->y());
			glutBitmapCharacter(GLUT_BITMAP_8_BY_13,48+dig);
		}
	}
}

void cellproject::drawBresenhamCircle(point2D *center, int radius, double color[3])
{
	glBegin(GL_POINTS);
	glColor3f(color[0],color[1],color[2]);
	int x1 = center->x();
	int y1 = center->y();
	int x=0;
	int y=radius;
	int d=1-radius;
	int dE=3;
	int dSE=-2*radius+5;
	glVertex2i(x1,y+y1);
	glVertex2i(x1,-y+y1);
	glVertex2i(y+x1,y1);
	glVertex2i(-y+x1,y1);
	while (y>x)
	{
		if (d<0)
		{
			d=d+dE;
			dE=dE+2;
			dSE=dSE+2;
		}
		else 
		{
			d=d+dSE;
			dE=dE+2;
			dSE=dSE+4;
			y--;
		}
		x++;
		glVertex2i(x+x1,y+y1);
		glVertex2i(x+x1,-y+y1);
		glVertex2i(-x+x1,y+y1);
		glVertex2i(-x+x1,-y+y1);
		glVertex2i(y+x1,x+y1);
		glVertex2i(y+x1,-x+y1);
		glVertex2i(-y+x1,x+y1);
		glVertex2i(-y+x1,-x+y1);
	}
	if (y==x)
	{
		glVertex2i(x+x1,x+y1);
		glVertex2i(x+x1,-x+y1);
		glVertex2i(-x+x1,x+y1);
		glVertex2i(-x+x1,-x+y1);
	}
	glEnd();
}

void cellproject::visualmatrix(double M[3][3])
{
	// Step 1: translate VRP to the origin
	double T[3][3];
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			if (i==j) {T[i][j]=1;}
			else if (j==2 && i==0) {T[i][j]=reader->getViewVolume()->getVRP()->x()*(-1);}
			else if (j==2 && i==1) {T[i][j]=reader->getViewVolume()->getVRP()->y()*(-1);}
			else {T[i][j]=0;}
		}
	}
	// Step 2: rotate VRC
	double R[3][3];
	double Rx[2];
	double Ry[2];
	
	double dy = sqrt(reader->getViewVolume()->getVPN()->x()*reader->getViewVolume()->getVPN()->x()+reader->getViewVolume()->getVPN()->y()*reader->getViewVolume()->getVPN()->y());
	Ry[0] = reader->getViewVolume()->getVPN()->x()/dy;
	Ry[1] = reader->getViewVolume()->getVPN()->y()/dy;
	
	double dx = sqrt(reader->getViewVolume()->getVUP()->x()*reader->getViewVolume()->getVUP()->x()+reader->getViewVolume()->getVUP()->y()*reader->getViewVolume()->getVUP()->y());
	Rx[0] = reader->getViewVolume()->getVUP()->x()/dx;
	Rx[1] = reader->getViewVolume()->getVUP()->y()/dx;

	R[0][0] = Rx[0];R[0][1] = Rx[1];R[0][2] = 0;
	R[1][0] = Ry[0];R[1][1] = Ry[1];R[1][2] = 0;
	R[2][0] = 0;    R[2][1] = 0;    R[2][2] = 1;

	// Step 3: Scale intot the canonical viewport
	double S[3][3];
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			if      (i==j && i==0) {S[i][j] = width/(reader->getViewVolume()->getUMax()-reader->getViewVolume()->getUMin());}
			else if (i==j && i==1) {S[i][j] = height/(reader->getViewVolume()->getVMax()-reader->getViewVolume()->getVMin());}
			else if (i==j && i==2) {S[i][j] = 1;}
			else {S[i][j]=0;}
		}
	}
	double Temp[3][3];
	matrix_multiplication(R,T,Temp);
	matrix_multiplication(S,Temp,M);
}

void cellproject::drawinnermesh(double M[3][3])
{
	for (int i=0;i<(int)reader->C()->size();i++)
	{
		if (reader->C()->at(i)->Dead()) continue;
		for (int j=0;j<(int)reader->C()->at(i)->get_inner_t()->size();j++)
		{
			double x1 = reader->C()->at(i)->get_inner_t()->at(j)->getA()->x();
			double y1 = reader->C()->at(i)->get_inner_t()->at(j)->getA()->y();
			double x2 = reader->C()->at(i)->get_inner_t()->at(j)->getB()->x();
			double y2 = reader->C()->at(i)->get_inner_t()->at(j)->getB()->y();
			double x3 = reader->C()->at(i)->get_inner_t()->at(j)->getC()->x();
			double y3 = reader->C()->at(i)->get_inner_t()->at(j)->getC()->y();
			double v1[3];
			double v2[3];
			double v3[3];
			double v1_temp[3];
			double v2_temp[3];
			double v3_temp[3];
			for (int k=0;k<3;k++)
			{
				if      (k==0) {v1[k]=x1;v2[k]=x2;v3[k]=x3;}
				else if (k==1) {v1[k]=y1;v2[k]=y2;v3[k]=y3;}
				else if (k==2) {v1[k]=1; v2[k]=1; v3[k]=1;}
			}
			matrix_vector_multiplication(M,v1,v1_temp);
			matrix_vector_multiplication(M,v2,v2_temp);
			matrix_vector_multiplication(M,v3,v3_temp);
			glLineWidth(1.0);
			glBegin(GL_LINES);
			glColor3f(0.8,0.8,0.8);
			glVertex2f(v1_temp[0],v1_temp[1]);
			glVertex2f(v2_temp[0],v2_temp[1]);
			glVertex2f(v1_temp[0],v1_temp[1]);
			glVertex2f(v3_temp[0],v3_temp[1]);
			glVertex2f(v2_temp[0],v2_temp[1]);
			glVertex2f(v3_temp[0],v3_temp[1]);
			glEnd();
		}
		for (int j=0;j<(int)reader->C()->at(i)->get_inner_p()->size();j++)
		{
			double x1 = reader->C()->at(i)->get_inner_p()->at(j)->x();
			double y1 = reader->C()->at(i)->get_inner_p()->at(j)->y();
			double v1[3];
			double v1_temp[3];
			for (int k=0;k<3;k++)
			{
				if      (k==0) {v1[k]=x1;}
				else if (k==1) {v1[k]=y1;}
				else if (k==2) {v1[k]=1;}
			}
			matrix_vector_multiplication(M,v1,v1_temp);
			double color[3];
			color[0]=0.3;color[1]=0.9;color[2]=0.1;
			point2D *P = new point2D(v1_temp[0],v1_temp[1]);
			drawBresenhamCircle(P,4,color);
			delete P;
		}
		for (int j=0;j<(int)reader->C()->at(i)->get_sides()->size();j++)
		{
			double x1 = reader->C()->at(i)->get_sides()->at(j)->p1()->x();
			double y1 = reader->C()->at(i)->get_sides()->at(j)->p1()->y();
			double x2 = reader->C()->at(i)->get_sides()->at(j)->p2()->x();
			double y2 = reader->C()->at(i)->get_sides()->at(j)->p2()->y();
			double v1[3];
			double v2[3];
			double v1_temp[3];
			double v2_temp[3];
			for (int k=0;k<3;k++)
			{
				if      (k==0) {v1[k]=x1;v2[k]=x2;}
				else if (k==1) {v1[k]=y1;v2[k]=y2;}
				else if (k==2) {v1[k]=1;v2[k]=1;}
			}
			matrix_vector_multiplication(M,v1,v1_temp);
			matrix_vector_multiplication(M,v2,v2_temp);
			glLineWidth(1.0);
			glBegin(GL_LINES);
			glColor3f(0.1,0.1,0.1);
			glVertex2f(v1_temp[0],v1_temp[1]);
			glVertex2f(v2_temp[0],v2_temp[1]);
			glEnd();
			//// interested cell to track ////
			if (i==6763 || i==1236 || i==3610)
			{
				double color[3];
				color[0]=0.9;color[1]=0.1;color[2]=0.7;
				point2D *P = new point2D(v1_temp[0],v1_temp[1]);
				drawBresenhamCircle(P,3,color);
				delete P;
			}
			/*if (reader->C()->at(i)->get_sides()->at(j)->p1()->id()==129773 ||
				reader->C()->at(i)->get_sides()->at(j)->p1()->id()==131082)
			{
				double color[3];
				color[0]=0.9;color[1]=0.1;color[2]=0.4;
				point2D *P = new point2D(v1_temp[0],v1_temp[1]);
				drawBresenhamCircle(P,3,color);
				delete P;
			}*/
		}
	}
}

void cellproject::drawtriangles(double M[3][3])
{
	for (int i=0;i<(int)reader->C()->size();i++)
	{
		if (reader->C()->at(i)->Dead()) continue;
		double colort[3];
		double colorc[3];
		colort[0] = reader->C()->at(i)->get_color(0);
		colort[1] = reader->C()->at(i)->get_color(1);
		colort[2] = reader->C()->at(i)->get_color(2);
		//colorc[0] = 1.0;
		//colorc[1] = 1.0;
		//colorc[2] = 1.0;
		colorc[0] = colort[0];
		colorc[1] = colort[1];
		colorc[2] = colort[2];
		if (view_type==2)
		{
			double pigment = 1.0 - min(100.0,(double)reader->C()->at(i)->get_cell_state(diffusion_number))/100.0;
			if (pigment<=0.01) {pigment = 0.1;}
			colort[0] = pigment;
			colort[1] = pigment;
			colort[2] = pigment;
			colorc[0] = pigment;
			colorc[1] = pigment;
			colorc[2] = pigment;
		}
		double xc = reader->C()->at(i)->get_center(0);
		double yc = reader->C()->at(i)->get_center(1);
		double vc[3];
		double vc_temp[3];
		for (int k=0;k<3;k++)
		{
			if      (k==0) {vc[k]=xc;}
			else if (k==1) {vc[k]=yc;}
			else if (k==2) {vc[k]=1;}
		}
		matrix_vector_multiplication(M,vc,vc_temp);
		for (int j=0;j<(int)reader->C()->at(i)->get_sides()->size();j++)
		{
			double colorl[3];
			if (!reader->C()->at(i)->get_sides()->at(j)->get_attach() &&
				 reader->C()->at(i)->get_cell_type()==0) // ECM
			{
				colorl[0] = 0.8;
				colorl[1] = 0.8;
				colorl[2] = 0.8;
				/*if (number==21)
				{
					if (VA[ind + j*2]>1000000 || VA[ind + j*2]<-1000000 ||
						VA[ind + j*2 + 1]>1000000 || VA[ind + j*2 + 1]<-1000000 || 
						VA[ind + nnzr*2 + 1 + j*2]>1000000 || VA[ind + nnzr*2 + 1 + j*2]<-1000000 ||
						VA[ind + nnzr*2 + 1 + j*2 + 1]>1000000 || VA[ind + nnzr*2 + 1 + j*2 + 1]<-1000000)
					{
						if (SMl->getN1()->get_cell_boundary())
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN1()->getab()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						else
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_SM_ID()<<" of cell "<<SMl->getN1()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						if (SMl->getN2()->get_cell_boundary())
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN2()->getab()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<": ";
							if (SMl->getN2()->getmp()->size()>0) {cout<<SMl->getN2()->getmp()->at(0)->getab()->Cell()->id()<<endl;}
							else {cout<<endl;}
						}
						else
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_SM_ID()<<" of cell "<<SMl->getN2()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<endl;
						}
						cout<<" triangle NULL N1-N2: 00: "<<SMl->get_stiffness_matrix_12(0)<<endl;
						cout<<" triangle NULL N1-N2: 01: "<<SMl->get_stiffness_matrix_12(1)<<endl;
						cout<<" triangle NULL N1-N2: 10: "<<SMl->get_stiffness_matrix_12(2)<<endl;
						cout<<" triangle NULL N1-N2: 11: "<<SMl->get_stiffness_matrix_12(3)<<endl;
					}
				}*/
			}
			else if (reader->C()->at(i)->get_sides()->at(j)->get_attach() && 
					 reader->C()->at(i)->get_sides()->at(j)->get_Pair()->get_I1()->get_cell_type()==0 &&
					 reader->C()->at(i)->get_sides()->at(j)->get_Pair()->get_I2()->get_cell_type()==0) // ECM
			{
				colorl[0] = 0.8;
				colorl[1] = 0.8;
				colorl[2] = 0.8;
				/*if (number==21)
				{
					if (VA[ind + j*2]>1000000 || VA[ind + j*2]<-1000000 ||
						VA[ind + j*2 + 1]>1000000 || VA[ind + j*2 + 1]<-1000000 || 
						VA[ind + nnzr*2 + 1 + j*2]>1000000 || VA[ind + nnzr*2 + 1 + j*2]<-1000000 ||
						VA[ind + nnzr*2 + 1 + j*2 + 1]>1000000 || VA[ind + nnzr*2 + 1 + j*2 + 1]<-1000000)
					{
						if (SMl->getN1()->get_cell_boundary())
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN1()->getab()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						else
						{
							cout<<" triangle NULL point N1: "<<SMl->getN1()->get_SM_ID()<<" of cell "<<SMl->getN1()->Cell()->id()<<" mpn: "<<SMl->getN1()->getmp()->size()<<endl;
						}
						if (SMl->getN2()->get_cell_boundary())
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_node_id()->SM_ID()<<" of cell "<<SMl->getN2()->getab()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<": ";
							if (SMl->getN2()->getmp()->size()>0) {cout<<SMl->getN2()->getmp()->at(0)->getab()->Cell()->id()<<endl;}
							else {cout<<endl;}
						}
						else
						{
							cout<<" triangle NULL point N2: "<<SMl->getN2()->get_SM_ID()<<" of cell "<<SMl->getN2()->Cell()->id()<<" mpn: "<<SMl->getN2()->getmp()->size()<<endl;							
						}
						cout<<" triangle NULL N1-N2: 00: "<<SMl->get_stiffness_matrix_12(0)<<endl;
						cout<<" triangle NULL N1-N2: 01: "<<SMl->get_stiffness_matrix_12(1)<<endl;
						cout<<" triangle NULL N1-N2: 10: "<<SMl->get_stiffness_matrix_12(2)<<endl;
						cout<<" triangle NULL N1-N2: 11: "<<SMl->get_stiffness_matrix_12(3)<<endl;
					}
				}*/
			}
			else
			{
				colorl[0] = 0.0;
				colorl[1] = 0.0;
				colorl[2] = 0.0;
			}
			double x1 = reader->C()->at(i)->get_sides()->at(j)->p1()->x();
			double y1 = reader->C()->at(i)->get_sides()->at(j)->p1()->y();
			double x2 = reader->C()->at(i)->get_sides()->at(j)->p2()->x();
			double y2 = reader->C()->at(i)->get_sides()->at(j)->p2()->y();
			double v1[3];
			double v2[3];
			double v1_temp[3];
			double v2_temp[3];
			for (int k=0;k<3;k++)
			{
				if      (k==0) {v1[k]=x1;v2[k]=x2;}
				else if (k==1) {v1[k]=y1;v2[k]=y2;}
				else if (k==2) {v1[k]=1;v2[k]=1;}
			}
			matrix_vector_multiplication(M,v1,v1_temp);
			matrix_vector_multiplication(M,v2,v2_temp);

			glBegin(GL_TRIANGLES);
			glColor3f(colort[0],colort[1],colort[2]);
			glVertex2f(v1_temp[0],v1_temp[1]);
			glColor3f(colort[0],colort[1],colort[2]);
			glVertex2f(v2_temp[0],v2_temp[1]);
			glColor3f(colorc[0],colorc[1],colorc[2]);
			glVertex2f(vc_temp[0],vc_temp[1]);
			glEnd();

			glLineWidth(2.0);
			glBegin(GL_LINES);
			glColor3f(colorl[0],colorl[1],colorl[2]);
			glVertex2f(v1_temp[0],v1_temp[1]);
			glVertex2f(v2_temp[0],v2_temp[1]);
			glEnd();
		}
		if (reader->C()->at(i)->Mark_number()==1 && show_direction)
		{
			double mangle = reader->C()->at(i)->get_migrate_angle();
			double dx = cos(mangle*PI/180);
			double dy = sin(mangle*PI/180);
			double p0x = dx*24 + vc_temp[0];
			double p0y = dy*24 + vc_temp[1];
			double p1x = dx*8 + vc_temp[0] + (-dy*4);
			double p1y = dy*8 + vc_temp[1] + (dx*4);
			double p2x = dx*8 + vc_temp[0] + (dy*4);
			double p2y = dy*8 + vc_temp[1] + (-dx*4);
			double p3x = dx*10 + vc_temp[0];
			double p3y = dy*10 + vc_temp[1];
			double colort1 = 0.1;
			double colort2 = 0.9;
			glBegin(GL_TRIANGLES);
			glColor3f(colort1,colort1,colort1);
			glVertex2f(p0x,p0y);
			glColor3f(colort1,colort1,colort1);
			glVertex2f(p1x,p1y);
			glColor3f(colort1,colort1,colort1);
			glVertex2f(p3x,p3y);
			glColor3f(colort1,colort1,colort1);
			glVertex2f(p0x,p0y);
			glColor3f(colort1,colort1,colort1);
			glVertex2f(p2x,p2y);
			glColor3f(colort1,colort1,colort1);
			glVertex2f(p3x,p3y);
			glEnd();
		}
	}
}

void cellproject::drawindex(double M[3][3])
{
	for (int i=0;i<(int)reader->C()->size();i++)
	{
		if (reader->C()->at(i)->Dead()) continue;
		double v1[3];
		double v1_temp[3];
		for (int k=0;k<3;k++)
		{
			if      (k==0) {v1[k]=reader->C()->at(i)->get_center(0);}
			else if (k==1) {v1[k]=reader->C()->at(i)->get_center(1);}
			else if (k==2) {v1[k]=1;}
		}
		if (i>=0 && reader->C()->at(i)->Mark_number()==1) // migrate
		{
			matrix_vector_multiplication(M,v1,v1_temp);
			double color[3];
			color[0] = 0.9;color[1] = 0.1;color[2] = 0.1;
			point2D *P = new point2D(v1_temp[0],v1_temp[1]);
			drawNumber(P,i,color);
			delete P;
			P = NULL;
		}
		else if (i>=0 && reader->C()->at(i)->Mark_number()==0)
		{
			matrix_vector_multiplication(M,v1,v1_temp);
			double color[3];
			color[0] = 0.5;color[1] = 0.5;color[2] = 0.5;
			point2D *P = new point2D(v1_temp[0],v1_temp[1]);
			drawNumber(P,i,color);
			delete P;
			P = NULL;
		}
		else if (i>=0 && reader->C()->at(i)->Mark_number()==2) // soften
		{
			matrix_vector_multiplication(M,v1,v1_temp);
			double color[3];
			color[0] = 0.1;color[1] = 0.1;color[2] = 0.9;
			point2D *P = new point2D(v1_temp[0],v1_temp[1]);
			drawNumber(P,i,color);
			delete P;
			P = NULL;
		}
		else if (i>=0 && reader->C()->at(i)->Mark_number()==3) // for user track
		{
			matrix_vector_multiplication(M,v1,v1_temp);
			double color[3];
			color[0] = 0.1;color[1] = 0.8;color[2] = 0.3;
			point2D *P = new point2D(v1_temp[0],v1_temp[1]);
			drawNumber(P,i,color);
			delete P;
			P = NULL;
		}
	}
}

void cellproject::drawforcedegree(double M[3][3])
{
	for (int i=0;i<(int)reader->C()->size();i++)
	{
		if (reader->C()->at(i)->Dead()) continue;
		for (int j=0;j<(int)reader->C()->at(i)->get_sides()->size();j++)
		{
			if (reader->C()->at(i)->get_sides()->at(j)->p1()->get_pop_out()<=0)
			{
				double x1 = reader->C()->at(i)->get_sides()->at(j)->p1()->x();
				double y1 = reader->C()->at(i)->get_sides()->at(j)->p1()->y();
				double sf_x = reader->C()->at(i)->get_sides()->at(j)->p1()->get_stored_force_vis(0);
				double sf_y = reader->C()->at(i)->get_sides()->at(j)->p1()->get_stored_force_vis(1);
				double force_v = sqrt(sf_x*sf_x + sf_y*sf_y);
				if (force_v>=1)
				{
					double dx = sf_x/force_v;
					double dy = sf_y/force_v;
					double v1[3];
					double v1_temp[3];
					for (int k=0;k<3;k++)
					{
						if      (k==0) {v1[k]=x1;}
						else if (k==1) {v1[k]=y1;}
						else if (k==2) {v1[k]=1;}
					}
					matrix_vector_multiplication(M,v1,v1_temp);
					double p0x = dx*16 + v1_temp[0];
					double p0y = dy*16 + v1_temp[1];
					double p1x = dx*8 + v1_temp[0] + (-dy)*2;
					double p1y = dy*8 + v1_temp[1] + (dx)*2;
					double p2x = dx*8 + v1_temp[0] + (dy)*2;
					double p2y = dy*8 + v1_temp[1] + (-dx)*2;
					double p3x = dx*8 + v1_temp[0];
					double p3y = dy*8 + v1_temp[1];
					double p4x = v1_temp[0] - 3;
					double p4y = v1_temp[1] + 3;
					double p5x = v1_temp[0] + 3;
					double p5y = v1_temp[1] + 3;
					double p6x = v1_temp[0] - 3;
					double p6y = v1_temp[1] - 3;
					double p7x = v1_temp[0] + 3;
					double p7y = v1_temp[1] - 3;
					double colora[3];
					colora[0] = 1.0;
					if (force_v>=10) {colora[1] = 0;colora[2] = 0;}
					else
					{
						colora[1] = 1 - force_v/10.0;
						colora[2] = 1 - force_v/10.0;
					}
					glBegin(GL_TRIANGLES);
					/*glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p0x,p0y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p1x,p1y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p2x,p2y);*/
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p4x,p4y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p5x,p5y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p6x,p6y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p5x,p5y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p6x,p6y);
					glColor3f(colora[0],colora[1],colora[2]);
					glVertex2f(p7x,p7y);
					glEnd();
				}
			}
		}
	}
	/*double colorb[3];
	colorb[0] = 1;colorb[1] = 1;colorb[2] = 1;
	glBegin(GL_TRIANGLES);
	//
	glColor3f(colorb[0],colorb[1],colorb[2]);
	glVertex2f(width-80,0);
	glColor3f(colorb[0],colorb[1],colorb[2]);
	glVertex2f(width-80,20);
	glColor3f(colorb[0],colorb[1],colorb[2]);
	glVertex2f(width,20);
	glColor3f(colorb[0],colorb[1],colorb[2]);
	glVertex2f(width,20);
	glColor3f(colorb[0],colorb[1],colorb[2]);
	glVertex2f(width,0);
	glColor3f(colorb[0],colorb[1],colorb[2]);
	glVertex2f(width-80,0);
	//
	glColor3f(1,0.9,0.9);
	glVertex2f(width-70,5);
	glColor3f(1,0.9,0.9);
	glVertex2f(width-70,15);
	glColor3f(1,0,0);
	glVertex2f(width-10,15);
	glColor3f(1,0,0);
	glVertex2f(width-10,15);
	glColor3f(1,0,0);
	glVertex2f(width-10,5);
	glColor3f(1,0.9,0.9);
	glVertex2f(width-70,5);
	//
	glEnd;
	double colort[3];
	colort[0] = 0.1;colort[1] = 0.1;colort[2] = 0.1;
	point2D *P1 = new point2D(width-80,10);
	drawNumber(P1,1,colort);
	delete P1;
	P1 = NULL;
	point2D *P2 = new point2D(width-10,10);
	drawNumber(P2,10,colort);
	delete P2;
	P2 = NULL;*/
}

void cellproject::Dataprojection()
{
	double M[3][3];
	visualmatrix(M);

	/////////////////////////////////////////////
	// project objects in WC to viewing window //
	/////////////////////////////////////////////

	if (reader->C()->size()>0)
	{
		/********************************
		  view_type:
		  0: cell type
		  1: inner mesh
  		  2: concentration
		  3: index
		  4: migrate point
		********************************/
		if (view_type==0 || view_type==2)
		{
			if (show_text)
			{
				string Molecule = reader->get_biology()->getSpeciesName(diffusion_number);
				int sizeM = (int)Molecule.size();
				double colort[3];colort[0] = 0.1;colort[1] = 0.2;colort[2] = 0.1;
				int poss = 9*sizeM;
				if (sizeM<=5) {poss = 18*sizeM;}
				else if (sizeM>5 && sizeM<=10) {poss = 15*sizeM;}
				else if (sizeM>10 && sizeM<=20) {poss = 12*sizeM;}
				point2D *TP = new point2D(width-poss,height-20);
				drawText(TP,&Molecule,colort);
				delete TP;
			}
			drawtriangles(M);
		}
		else if (view_type==1)
		{
			drawinnermesh(M);
		}
		if (show_index)
		{
			drawindex(M);
		}
		if (show_force)
		{
			drawforcedegree(M);
		}
	}
	/* selection box */
	if (mouse==2 && mouses==1) // mouse right button
	{
		glLineWidth(1.0);
		glBegin(GL_LINES);
		glColor3f(0.4,0.2,0.1);
		glVertex2f(mousec[0],height - mousec[1]);
		glVertex2f(mousec[0],height - mousec_tmp[1]);
		glVertex2f(mousec[0],height - mousec_tmp[1]);
		glVertex2f(mousec_tmp[0],height - mousec_tmp[1]);
		glVertex2f(mousec_tmp[0],height - mousec_tmp[1]);
		glVertex2f(mousec_tmp[0],height - mousec[1]);
		glVertex2f(mousec_tmp[0],height - mousec[1]);
		glVertex2f(mousec[0],height - mousec[1]);
		glEnd();
	}
}

#endif

