// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     biology.cpp         ***
 ***    Author:   Youfang Cao,        ***
 ***              Jieling Zhao        ***
 ***                                  ***
 ***    Created on May 20, 2013       ***
 ***************************************/

#ifndef BIOLOGY_CPP
#define BIOLOGY_CPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include "dbReader.h"
#include "biology.h"
#include "para.h"
#include <libBioModel/Network.h>
#include <libBioModel/Stoichiometry.h>
#include <libBioModel/State.h>
#include <libBioModel/GaussianRNG.h>
#include <sbml/Model.h>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <common/dycelfem_config.h>

#define MIN_FBEHAVIOR_PERSISTENCE 5

/* cell_setupbiobehaviors() below is switched by a MODE variable:
     0 = INITIALIZATION       prepare/relax a tissue, no chemistry at all
     1 = MODELING             the actual simulation: reactions + diffusion
     2 = DELETING AND RELAXING
   Committed trunk r517 had MODE = 1.  Both archived working copies were left on
   MODE = 0 by an uncommitted edit, which is why they moved cells around but
   never ran any reaction or diffusion step.  Default restored to 1; override
   with DYCELFEM_MODE=0 or =2 to reach the tissue-preparation workflows. */
/* Load one SBML model into a Network, with a clear error rather than a crash
   when the file is missing: readSBML() happily returns a document with no
   model, and the very next getNumSpecies() call would dereference NULL. */
static void dycelfem_load_network(Network *nw, const char *cfgkey, const char *fallback)
{
	std::string file = dycfg::resolve(dycfg::getStr(cfgkey, fallback));
	std::ifstream probe(file.c_str());
	if (!probe.good())
	{
		std::cerr << "ERROR: SBML model file not found: " << file << std::endl;
		std::cerr << "       (config key '" << cfgkey << "'). Every run needs all seven"
		          << " model files." << std::endl;
		exit(1);
	}
	probe.close();
	std::cout << "  SBML " << cfgkey << " <- " << file << std::endl;
	if (nw->initializeNetwork(const_cast<char*>(file.c_str())) != 0)
	{
		std::cerr << "ERROR: failed to build the reaction network from " << file << std::endl;
		exit(1);
	}
}

/* Inter-cell diffusion kernel selector. See the switch in
   cell_setupbiobehaviors(). DYCELFEM_DIFFUSION = edgelen | edgelen_distsq | distsq */
enum { DIFF_DISTSQ = 0, DIFF_EDGELEN = 1, DIFF_EDGELEN_DISTSQ = 2 };

static int dycelfem_diffusion_kernel()
{
	const char *legacy = getenv("DYCELFEM_DIFFUSION");   /* pre-config env name */
	std::string v = (legacy && legacy[0]) ? std::string(legacy)
	                                      : dycfg::getStr("run.diffusion", "edgelen");
	if (v == "distsq")         return DIFF_DISTSQ;
	if (v == "edgelen_distsq") return DIFF_EDGELEN_DISTSQ;
	return DIFF_EDGELEN;                                  // reproduces the published runs
}

static int dycelfem_mode()
{
	const char *legacy = getenv("DYCELFEM_MODE");        /* pre-config env name */
	int n = (legacy && legacy[0]) ? atoi(legacy) : dycfg::getInt("run.mode", 1);
	return (n == 0 || n == 1 || n == 2) ? n : 1;
}

using namespace std;

Biology::Biology()
{
	nw_macrophage = new Network();
	nw_fibroblast = new Network();
	nw_keratinocyte = new Network();
	nw_ecm = new Network();
	nw_bm = new Network();
	nw_clot = new Network();
	nw_common = new Network();

	/* Model files come from the run configuration ([sbml] section); the
	   defaults are the filenames the code used to hardcode. Paths are resolved
	   relative to the config file, so a run directory can live anywhere. */
	dycelfem_load_network(nw_macrophage,   "sbml.macrophage",   "mmacrophage.xml");
	dycelfem_load_network(nw_fibroblast,   "sbml.fibroblast",   "mfibroblast.xml");
	dycelfem_load_network(nw_keratinocyte, "sbml.keratinocyte", "mkeratinocyte.xml");
	dycelfem_load_network(nw_ecm,          "sbml.ecm",          "mecm.xml");
	dycelfem_load_network(nw_bm,           "sbml.bm",           "mbm.xml");
	dycelfem_load_network(nw_clot,         "sbml.clot",         "mclot.xml");
	dycelfem_load_network(nw_common,       "sbml.common",       "mcommon.xml");

	nklaws_macrophage = nw_macrophage->getNumReactions();
	nklaws_fibroblast = nw_fibroblast->getNumReactions();
	nklaws_keratinocyte = nw_keratinocyte->getNumReactions();
	nklaws_ecm = nw_ecm->getNumReactions();
	nklaws_bm = nw_bm->getNumReactions();
	nklaws_clot = nw_clot->getNumReactions();
	nklaws_common = nw_common->getNumReactions();

	kls_macrophage = new KineticLaw*[nklaws_macrophage];
	kls_fibroblast = new KineticLaw*[nklaws_fibroblast];
	kls_keratinocyte = new KineticLaw*[nklaws_keratinocyte];
	kls_ecm = new KineticLaw*[nklaws_ecm];
	kls_bm = new KineticLaw*[nklaws_bm];
	kls_clot = new KineticLaw*[nklaws_clot];
	kls_common = new KineticLaw*[nklaws_common];

	for (int k=0; k<nklaws_macrophage; k++)
	{
		kls_macrophage[k] = nw_macrophage->getModel()->getReaction(abs(nw_macrophage->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}
	for (int k=0; k<nklaws_fibroblast; k++)
	{
		kls_fibroblast[k] = nw_fibroblast->getModel()->getReaction(abs(nw_fibroblast->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}
	for (int k=0; k<nklaws_keratinocyte; k++)
	{
		kls_keratinocyte[k] = nw_keratinocyte->getModel()->getReaction(abs(nw_keratinocyte->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}
	for (int k=0; k<nklaws_ecm; k++)
	{
		kls_ecm[k] = nw_ecm->getModel()->getReaction(abs(nw_ecm->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}
	for (int k=0; k<nklaws_bm; k++)
	{
		kls_bm[k] = nw_bm->getModel()->getReaction(abs(nw_bm->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}
	for (int k=0; k<nklaws_clot; k++)
	{
		kls_clot[k] = nw_clot->getModel()->getReaction(abs(nw_clot->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}
	for (int k=0; k<nklaws_common; k++)
	{
		kls_common[k] = nw_common->getModel()->getReaction(abs(nw_common->getStoichiometricMatrix()->getReactionAttribute(k)) - 1)->getKineticLaw();
	}

	// Actual area of a cell, 30 micrometers diameter assumed.
	AREA0 = 7.07E-6; // cm^2/cell. pi*15*15
	// Basal cell growth rates
	fg0 = 0.1 * AR; // fibroblast
	kg0 = 0.1 * AR; // keratinocyte
	// Basal dimension for diffusion rate
	D0 = nw_common->getModel()->getParameter("D0")->getValue();
	// Time Step for simulation
	DeltaT = nw_common->getModel()->getParameter("DeltaT")->getValue();

	// Mutually Exclusive Cell behaviors
	// 0. Waiting for Behavior Determination; 1. Do nothing; 2. Growth/Division; 3. Apoptosis; 4.Migration.
	Prob_DoNothing = 0.0;
	Prob_GrowthDiv = 0.0;
	Prob_Apoptosis = 0.0;
	Prob_Migration = 0.0;

	// 0. fibrin clot; 1. collagen; 2. PDGF; 3. IL-1; 4. TGF-alpha; 5. EGF; 
	// 6. TGF-beta1; 7. TGF-beta3; 8. FGF2; 9. MMP; 10. TNF-alpha; 11. KGF(FGF7); 12. procollagen
	N_SPECIES = 20;
	nrealSPECIES = nw_common->getModel()->getNumSpecies();

	/* 	cell type:
		0: ECM (virtual cell type)
		1: BM (basement membrane, virtual cell type)
		2: Fibroblast
		3: Keratinocyte
		4: Clot (virtual cell type)
		5: Hypodermis
		6: Macrophage
	*/
	N_CELLTYPES = 6;

		/***************************************
		cell type:
		0: ECM (virtual cell type, not real cell)
		1: BM (basement membrane)
		2: fibroblast
		3: keratinocyte
		4: clot
		5: Hypodermis
		6: Macrophage
		---------------------------------------
		adhesion type:
		0: ECM-ECM      0-0
		1: ECM-BM       0-1
		2: ECM-f        0-2
		3: ECM-k        0-3
		4: ECM-FC       0-4
		5: BM-BM        1-1
		6: BM-f         1-2
		7: BM-k         1-3
		8: BM-FC        1-4
		9: f-f          2-2
		10: f-k          2-3
		11: f-FC         2-4
		12: k-k          3-3
		13: k-FC         3-4
		14: FC-FC        4-4
		16: MP-ECM      0-6
		17: MP-BM       1-6
		18: MP-f        2-6
		19: MP-k        3-6
		20: MP-FC       4-6
		21: MP-MP       6-6
		***************************************/


	// 0. fibrin clot; 1. collagen; 2. PDGF; 3. IL-1; 4. TGF-alpha; 5. EGF; 
	// 6. TGF-beta1; 7. TGF-beta3; 8. FGF2; 9. MMP; 10. TNF-alpha; 11. KGF(FGF7); 
	// Diffusion rates STORED IN nw_fibroblast NETWORK
	D = new double[nrealSPECIES];
	D[0]  = nw_common->getModel()->getParameter("D_WoundSignal")->getValue();   // 1.0E-9; cm^2/s, 0. wound signal
	D[1]  = nw_common->getModel()->getParameter("D_clot")->getValue();		 // 0; cm^2/s, 1. fibrin clot
	D[2]  = nw_common->getModel()->getParameter("D_collagen")->getValue();   // 1.0E-9; cm^2/s, 2. collagen
	D[3]  = nw_common->getModel()->getParameter("D_PDGF")->getValue();    // 2.78E-8; cm^2/s, 3. PDGF
	D[4]  = nw_common->getModel()->getParameter("D_KGF")->getValue();  // 9.5E-6; cm^2/s, 4.KGF(FGF7)
	D[5]  = nw_common->getModel()->getParameter("D_IL1")->getValue(); // 5.18E-7; cm^2/s, 5. EGF (Thorne et al. 2004 J Neurophysiol 92:3471-3481.)
	D[6]  = nw_common->getModel()->getParameter("D_TGFb1")->getValue(); // 2.94E-7; cm^2/s, 6. TGF-beta1
	D[7]  = nw_common->getModel()->getParameter("D_Procollagen")->getValue(); // cm^2/s, 7. procollagen
	D[8]  = nw_common->getModel()->getParameter("D_MMP")->getValue(); // 8.0E-9; cm^2/s, 9. MMP-1 (Saffarian et al. 2004 Science 306:108-111.)
	
	SpeciesNames = new string[nrealSPECIES];
	for (int i = 0; i<nrealSPECIES; i++)
	{
		SpeciesNames[i] = nw_common->getSpeciesID(i);
	}
}

int Biology::getNumSpecies ()
{
	return N_SPECIES;
}

int Biology::getNumRealSpecies ()
{
	return nrealSPECIES;
}

int Biology::getNumCellTypes ()
{
	return N_CELLTYPES;
}

double Biology::getDeltaT()
{
	return DeltaT;
}

double Biology::getD0()
{
	return D0;
}

double Biology::getD(int k)
{
	return D[k];
}

string Biology::getSpeciesName(int i)
{
	return SpeciesNames[i];
}

int Biology::UpdateCellState (cell* pcell)
{
	Network* nw = 0;
	KineticLaw** klaws = 0;
	int nklaws = 0;
	int celltype = pcell->get_cell_type();
	int nsteps = 0;
	switch (celltype)
	{
	case 0:
		nw = nw_ecm;
		klaws = kls_ecm;
		nklaws = nklaws_ecm;
		break;
	case 1:
		nw = nw_bm;
		klaws = kls_bm;
		nklaws = nklaws_bm;
		break;
	case 2:
		nw = nw_fibroblast;
		klaws = kls_fibroblast;
		nklaws = nklaws_fibroblast;
		break;
	case 3:
		nw = nw_keratinocyte;
		klaws = kls_keratinocyte;
		nklaws = nklaws_keratinocyte;
		break;
	case 4:
		nw = nw_clot;
		klaws = kls_clot;
		nklaws = nklaws_clot;
		break;
	case 5:
		return 0;
		break;
	case 6:
		nw = nw_macrophage;
		klaws = kls_macrophage;
		nklaws = nklaws_macrophage;
		break;
	}
	double* rates = new double[nklaws];
	double* rates_sqrt = new double[nklaws];
	double* gaussrns = new double[nklaws];
	double* state = new double[nrealSPECIES];
	double* d1_state = new double[nrealSPECIES];
	double* d2_state = new double[nrealSPECIES];
	State* currstate = new State(nrealSPECIES);
	double t = 0;
	//double dt = 0.1;
	double dt = 6.0;
	double dt_sqrt = sqrt(dt);
	double tau1 = 0;
	double tau2 = 0;
	double rates0 = 0;
	double cum = 0.0;
    double cri = 0.0;
    int f = -1;

	for (int i=0; i<nrealSPECIES; i++)
	{
		currstate->setStateValue(i, (unsigned short)pcell->get_cell_state (i));
		state[i] = (double)pcell->get_cell_state (i);
		d1_state[i] = 0.0;
		d2_state[i] = 0.0;
	}

	for (int k=0; k<nklaws; k++)
	{
		rates[k] = 0.0;
		gaussrns[k] = 0.0;
	}

	// SSA
	//while (t < DeltaT)
	//{
	//	double r1 = randn_trig();
	//	rates0 = 0.0;
	//	for (int k=0; k<nklaws; k++)
	//	{
	//		rates[k] = nw->getReactionRate(klaws[k], k, currstate);
	//		rates0 += rates[k];
	//	}
	//	if (rates0 == 0)
	//	{
	//		cout << "Reaction rates sum to zero." << endl;
	//		return 0;
	//	}
	//	tau1 = ((double)rand()+1.0)/((unsigned long)RAND_MAX+1);
	//	tau2 = ((double)rand()+1.0)/((unsigned long)RAND_MAX+1);
	//	dt = -log(tau1)/rates0;
	//	cum = 0.0;
	//	cri = tau2 * rates0;
	//	f = -1;
 //       while ((cum < cri) && (f < nklaws))
 //       {
 //               ++f;
 //               cum += rates[f];
 //       }
 //       /* Now, I got which reaction to go */
 //       if (f == -1)
 //       {
 //               cout << "What happened?" << endl;
 //               return -1;
 //       }
 //       /* Update the system state X and time t */
	//	nw->getStoichiometricMatrix()->updateStateWithReaction(f, currstate);
	//	t += dt;
	//	nsteps ++;
	//}
	// END SSA

	// SDE
	while (t < DeltaT) // DeltaT=60mins defined in mcommon.xml.
	{
		//cout << "time=" << t << ", step=" << nsteps << ", dT=" << DeltaT << endl;
		for (int k=0; k<nklaws; k++)
		{
			rates[k] = nw->getReactionRate(klaws[k], k, currstate);
			rates_sqrt[k] = sqrt(rates[k]);
			gaussrns[k] = randn_trig();
		}
		
		//if (celltype == 6)
		//{
		//	cout << "    >>>> Mphi PDGF SynRate=" << rates[0] << endl;
		//}


		Stoichiometry* st = nw->getStoichiometricMatrix();
		for (int i = 0; i<nrealSPECIES; i++)
		{
			d1_state[i] = 0.0;
			d2_state[i] = 0.0;
			for (int j = 0; j<nklaws; j++)
			{
				d1_state[i] += st->getValue(i,j) * rates[j];
				d2_state[i] += st->getValue(i,j) * rates_sqrt[j] * gaussrns[j];
			}
			state[i] += d1_state[i] * dt + d2_state[i] * dt_sqrt;
			if (state[i] < 0.0) { state[i] = 0.0; }
	        /* Update the system state X and time t */
			currstate->setStateValue(i, (unsigned short)state[i]);
		}
		t += dt;
		nsteps ++;
	}
	// END SDE

	
	for (int i=0; i<nrealSPECIES; i++)
	{
		pcell->set_cell_state(i, currstate->getStateValue(i));
	}


	delete rates;
	delete currstate;
	delete rates_sqrt;
	delete gaussrns;
	delete state;
	delete d1_state;
	delete d2_state;

	return nsteps;
}

int Biology::setValueProbGrowthDiv(double p)
{
	Prob_GrowthDiv = p;
	return 0;
}

int Biology::setValueProbMigration(double p)
{
	Prob_Migration = p;
	return 0;
}

int Biology::setValueProbApoptosis(double p)
{
	Prob_Apoptosis = p;
	return 0;
}

int Biology::setValueProbDoNothing(double p)
{
	Prob_DoNothing = p;
	return 0;
}

double Biology::getValueProbGrowthDiv()
{
	return Prob_GrowthDiv;
}

double Biology::getValueProbMigration()
{
	return Prob_Migration;
}

double Biology::getValueProbApoptosis()
{
	return Prob_Apoptosis;
}

double Biology::getValueProbDoNothing()
{
	return Prob_DoNothing;
}

double Biology::CalcValueProbGrowthDiv(cell* pcell)
{
	Network* nw = 0;
	KineticLaw** klaws = 0;
	int nklaws = 0;
	int celltype = pcell->get_cell_type();
	double prob_growthdiv = 0;
	switch (celltype)
	{
	case 0:
		nw = nw_ecm;
		klaws = kls_ecm;
		nklaws = nklaws_ecm;
		break;
	case 2:
		nw = nw_fibroblast;
		klaws = kls_fibroblast;
		nklaws = nklaws_fibroblast;
		break;
	case 3:
		nw = nw_keratinocyte;
		klaws = kls_keratinocyte;
		nklaws = nklaws_keratinocyte;
		break;
	case 6:
		nw = nw_macrophage;
		klaws = kls_macrophage;
		nklaws = nklaws_macrophage;
		break;
	default:
		return prob_growthdiv;
	}

	State* currstate = new State(pcell->get_cell_statelen());

	for (int i=0; i<pcell->get_cell_statelen(); i++)
	{
		currstate->setStateValue(i, pcell->get_cell_state (i));
	}
	
	prob_growthdiv = nw->getValueOfParameter("Prob_GrowthDiv", 0, currstate);
	//TEST
	//if (celltype == 2)
	//{
	//	cout << "   >>>>> >>>>> >>>>> " << prob_growthdiv << " ---- " << currstate->getStateValue(3) << endl;
	//}
	//
	delete currstate;
	return prob_growthdiv;
}

double Biology::CalcValueProbMigration(cell* pcell)
{
	Network* nw = 0;
	KineticLaw** klaws = 0;
	int nklaws = 0;
	int celltype = pcell->get_cell_type();
	double prob_migration = 0;
	switch (celltype)
	{
	case 0:
		nw = nw_ecm;
		klaws = kls_ecm;
		nklaws = nklaws_ecm;
		break;
	case 2:
		nw = nw_fibroblast;
		klaws = kls_fibroblast;
		nklaws = nklaws_fibroblast;
		break;
	case 3:
		nw = nw_keratinocyte;
		klaws = kls_keratinocyte;
		nklaws = nklaws_keratinocyte;
		break;
	case 6:
		nw = nw_macrophage;
		klaws = kls_macrophage;
		nklaws = nklaws_macrophage;
		break;
	default:
		return prob_migration;
	}

	State* currstate = new State(pcell->get_cell_statelen());

	for (int i=0; i<pcell->get_cell_statelen(); i++)
	{
		currstate->setStateValue(i, pcell->get_cell_state (i));
	}
	
	prob_migration = nw->getValueOfParameter("Prob_Migration", 0, currstate);
	//TEST
	//if (celltype == 2)
	//{
	//	cout << "   >>>>> >>>>> >>>>> " << prob_migration << " ---- " << currstate->getStateValue(3) << endl;
	//}
	//
	delete currstate;
	return prob_migration;
}

double Biology::CalcValueProbApoptosis(cell* pcell)
{
	Network* nw = 0;
	KineticLaw** klaws = 0;
	int nklaws = 0;
	int celltype = pcell->get_cell_type();
	double prob_apoptosis = 0;
	switch (celltype)
	{
	case 0:
		nw = nw_ecm;
		klaws = kls_ecm;
		nklaws = nklaws_ecm;
		break;
	case 2:
		nw = nw_fibroblast;
		klaws = kls_fibroblast;
		nklaws = nklaws_fibroblast;
		break;
	case 3:
		nw = nw_keratinocyte;
		klaws = kls_keratinocyte;
		nklaws = nklaws_keratinocyte;
		break;
	case 6:
		nw = nw_macrophage;
		klaws = kls_macrophage;
		nklaws = nklaws_macrophage;
		break;
	default:
		return prob_apoptosis;
	}

	State* currstate = new State(pcell->get_cell_statelen());

	for (int i=0; i<pcell->get_cell_statelen(); i++)
	{
		currstate->setStateValue(i, pcell->get_cell_state (i));
	}
	
	prob_apoptosis = nw->getValueOfParameter("Prob_Apoptosis", 0, currstate);
	//TEST
	//if (celltype == 2)
	//{
	//	cout << "   >>>>> >>>>> >>>>> " << prob_apoptosis << " ---- " << currstate->getStateValue(3) << endl;
	//}
	//

	delete currstate;
	return prob_apoptosis;
}

double Biology::CalcValueProbDoNothing(cell* pcell)
{
	Network* nw = 0;
	KineticLaw** klaws = 0;
	int nklaws = 0;
	int celltype = pcell->get_cell_type();
	double prob_donothing = 0;
	switch (celltype)
	{
	case 0:
		nw = nw_ecm;
		klaws = kls_ecm;
		nklaws = nklaws_ecm;
		break;
	case 2:
		nw = nw_fibroblast;
		klaws = kls_fibroblast;
		nklaws = nklaws_fibroblast;
		break;
	case 3:
		nw = nw_keratinocyte;
		klaws = kls_keratinocyte;
		nklaws = nklaws_keratinocyte;
		break;
	case 6:
		nw = nw_macrophage;
		klaws = kls_macrophage;
		nklaws = nklaws_macrophage;
		break;
	default:
		return prob_donothing;
	}

	State* currstate = new State(pcell->get_cell_statelen());

	for (int i=0; i<pcell->get_cell_statelen(); i++)
	{
		currstate->setStateValue(i, pcell->get_cell_state (i));
	}
	
	prob_donothing = nw->getValueOfParameter("Prob_DoNothing", 0, currstate);
	//TEST
	//if (celltype == 2)
	//{
	//	cout << "   >>>>> >>>>> >>>>> " << prob_donothing << " ---- " << currstate->getStateValue(3) << endl;
	//}
	//

	delete currstate;
	return prob_donothing;
}

double Biology::CalcValueMigForce(cell* pcell)
{
	Network* nw = 0;
	KineticLaw** klaws = 0;
	int nklaws = 0;
	int celltype = pcell->get_cell_type();
	double migforce = 0;
	switch (celltype)
	{
	case 0:
	case 1:
	case 4:
	case 5:
	case 6:
		return 0.0;
	case 2:
		nw = nw_fibroblast;
		klaws = kls_fibroblast;
		nklaws = nklaws_fibroblast;
		break;
	case 3:
		nw = nw_keratinocyte;
		klaws = kls_keratinocyte;
		nklaws = nklaws_keratinocyte;
		break;
	default:
		return migforce;
	}

	State* currstate = new State(pcell->get_cell_statelen());

	for (int i=0; i<pcell->get_cell_statelen(); i++)
	{
		currstate->setStateValue(i, pcell->get_cell_state (i));
	}
	
	migforce = nw->getValueOfParameter("MigForce", 0, currstate);
	//TEST
	//if (celltype == 2)
	//{
	//	cout << "   >>>>> >>>>> >>>>> " << prob_growthdiv << " ---- " << currstate->getStateValue(3) << endl;
	//}
	//
	delete currstate;
	return migforce;
}


void UserInitialParas()
{
	// Adhesion as pressure
	// 1 unit = 4nN/micrometer^2
	adhesion_pair_array[0][0] = 20;  // ECM-ECM   0-0
	adhesion_pair_array[0][1] = 20;  // ECM-BM    0-1
	adhesion_pair_array[0][2] = 2.3; // ECM-f     0-2 (3nN/micrometer^2 * Pi) / 4nN/micrometer^2 [R. Ananthakrishnan, et al. Int. J. Biol. Sci. 2007, 3:303.]
	adhesion_pair_array[0][3] = 9.8; // ECM-k     0-3 (12.7nN/micrometer^2 * Pi) / 4nN/micrometer^2 [O. Roure et al., PNAS., 2005 102(7):2390.]
	adhesion_pair_array[0][4] = 1;   // ECM-FC    0-4
	adhesion_pair_array[1][1] = 20;  // BM-BM     1-1
	adhesion_pair_array[1][2] = 2;   // BM-f      1-2
	adhesion_pair_array[1][3] = 6;   // BM-k      1-3
	adhesion_pair_array[1][4] = 1;   // BM-FC     1-4
	adhesion_pair_array[2][2] = 0;   // f-f       2-2 (3.9E-3nN/micrometer^2 * Pi) / 4nN/micrometer^2 [B. Hinz, et al., Mol. Biol. Cell, 2004 15:4310.]
	adhesion_pair_array[2][3] = 0;   // f-k       2-3
	adhesion_pair_array[2][4] = 2;   // f-FC      2-4
	adhesion_pair_array[3][3] = 3.2; // k-k       3-3 (4.1nN/micrometer^2 * Pi) / 4nN/micrometer^2 [Y. Chu, et al., J. Cell Biol., 2004 167(6):1183.]
	adhesion_pair_array[3][4] = 0;   // k-FC      3-4
	adhesion_pair_array[4][4] = 2;   // FC-FC     4-4
	adhesion_pair_array[0][6] = 2.3; // ECM-MP    0-6
	adhesion_pair_array[1][6] = 0;   // BM-MP     1-6
	adhesion_pair_array[2][6] = 0;   // f-MP      2-6
	adhesion_pair_array[3][6] = 0;   // k-MP      3-6
	adhesion_pair_array[4][6] = 0;   // FC-MP     4-6
	for (int i=1;i<20;i++)
	{
		for (int j=0;j<i;j++)
		{
			adhesion_pair_array[i][j] = adhesion_pair_array[j][i];
		}
	}
	
	////// COLORING //////
	// ECM
	color_class_R[0] = 0.9;
	color_class_G[0] = 0.9;
	color_class_B[0] = 0.9;
	// BM
	color_class_R[1] = 0.9;
	color_class_G[1] = 0.66;
	color_class_B[1] = 0.15;
	// Fibroblast
	color_class_R[2] = 1.0;
	color_class_G[2] = 0.0;
	color_class_B[2] = 0.4;
	// Keratinocyte
	color_class_R[3] = 0.0;
	color_class_G[3] = 0.4;
	color_class_B[3] = 1.0;
	// Fibrin Clot
	color_class_R[4] = 0.2;
	color_class_G[4] = 1.0;
	color_class_B[4] = 0.2;
	// OTHERS
	color_class_R[5] = 0.1;
	color_class_G[5] = 0.1;
	color_class_B[5] = 0.1;
	// Macrophage
	color_class_R[6] = 1.0;
	color_class_G[6] = 0.8;
	color_class_B[6] = 0.0;

	////// LAME CONSTANT ////////
	elastic_MU[0] = 1.91; // ECM; Use: get_lame_constant(E, v, 0)
	elastic_MU[1] = 1.91; // BM
	elastic_MU[2] = 1.91; // Fibroblast
	elastic_MU[3] = 1.91; // Keratinocyte
	elastic_MU[4] = 1.91; // Fibrin Clot
	elastic_MU[5] = 1.91; // OTHERS
	elastic_MU[6] = 1.91; // Macrophage
	//
	elastic_LAMBDA[0] = 2.87; // ECM; Use: get_lame_constant(E, v, 1)
	elastic_LAMBDA[1] = 2.87; // BM
	elastic_LAMBDA[2] = 2.87; // Fibroblast
	elastic_LAMBDA[3] = 2.87; // Keratinocyte
	elastic_LAMBDA[4] = 2.87; // Fibrin Clot
	elastic_LAMBDA[5] = 2.87; // OTHERS
	elastic_LAMBDA[6] = 2.87; // Macrophage
	//
	fixed_cells[0] = 0;
	fixed_cells[1] = 0;
	fixed_cells[2] = 0;
	fixed_cells[3] = 0;
	fixed_cells[4] = 0;
	fixed_cells[5] = 1;  // FIX SPATIAL POSTION OF TYPE 5 CELLS
	fixed_cells[6] = 0;
}

bool SortNeighborCellsByAngles (vector<int> i, vector<int> j) { return (i[1]<j[1]); }

int dbReader::cell_setupbiobehaviors(int Time)
{
	// INTERFACE: MAIN ENTRY
	int MODE = dycelfem_mode(); // SETUP A MODELING MODE, 0: INITIALIZATION; and 1: MODELING; and 2: DELETING AND RELAXING CELLS.
	cell *pcell0 = 0, *pcell1 = 0;
	point3D *s1 = 0;
	point3D *t1 = 0;
	point3D *s2 = 0;
	point3D *t2 = 0;
	point3D *pp = 0;
	// Definition of the diffusion matrix.
	vector<int> ia;
	vector<int> ja;
	vector< vector<double> > a;
	vector< vector<double> > v;
	for (int j = 0; j < biology->getNumRealSpecies(); j++)
	{
		vector<double> a0;
		a.push_back(a0);
		vector<double> v0;
		v.push_back(v0);
	}
	vector<int> idx2cellid;
	vector<int> cellid2idx;
	vector<double> vii;
	for (int j = 0; j < biology->getNumRealSpecies(); j++)
	{
		vii.push_back (0.0);
	}
	int iind = 0;
	// 
	double commedgelen = 0;
	double *dx = new double[biology->getNumRealSpecies()];
	int ra = 0;
	int num_neighbors = 0;
	vector<int> neighbors;
	vector<double> cytokines;
	int cn = (int)cellList.size();
	
	// 0: INITIALIZATION MODE
	if (MODE == 0)
	{
		cn = (int)cellList.size();
		for (int i=0;i<cn;i++)
		{
			if (cellList[i]->Dead()) continue; // skip dead cells.

			if (cellList[i]->get_cell_type () == 0) // ECM
			{
				//if (abs(cellList[i]->get_center(0)) < 200) { cellList[i]->set_GR (10); }
				if (cellList[i]->get_area()>=AR)
				{
					ra = (int)cellList[i]->shape_based_division_angle();
					cout<<"  -> ";
					cell_division(cellList[i], Time, 0, 0, ra);
				}
			}
		}
		cn = (int)cellList.size();
		for (int i=0;i<cn;i++)
		{
			if (cellList[i]->Dead()) continue; // skip dead cells.

			if (cellList[i]->get_cell_type () == 0) // ECM
			{
				//if (abs(cellList[i]->get_center(0)) < 200) { cellList[i]->set_GR (10); }
				cellList[i]->set_GR (GAR);
			}
		}
		return -1;
	}

	// 2: RELAXATION MODE
	else if (MODE == 2)
	{
		cout << "RELAXATION MODE." << endl;
		double max_area = 0;
		double min_area = 10000;
		double area = 0;
		cn = (int)cellList.size();
		for (int i=0;i<cn;i++)
		{
			if (cellList[i]->Dead()) continue; // skip dead cells.

			cellList[i]->set_GR (0);

			/*if (cellList[i]->get_cell_type () == 2)
			{
				cellList[i]->set_GR (GAR*5);
			}
			else
			{
				cellList[i]->set_GR (0);
			}*/
			if (cellList[i]->get_area()>=AR)
			{
				int ra = (int)cellList[i]->shape_based_division_angle();
				// All divisions are self-renewal type
				cell_division(cellList[i], Time, cellList[i]->get_cell_type(), cellList[i]->get_cell_type(), ra);
			}
		}
		cn = (int)cellList.size();
		for (int i=0;i<cn;i++)
		{
			if (cellList[i]->Dead()) continue; // skip dead cells.

			area = cellList[i]->get_area();
			if (area > max_area) max_area = area;
			if (area < min_area) min_area = area;
		}
		cout << "Area Diff: " << max_area - min_area << ", MAX_AREA: " << max_area << ", MIN_AREA: " << min_area << endl;
		return -1;
	}


	// 1: MODELING MODE
	cout << "Biology -> Start" << endl;
	cout << "  0. Building Cell-ID maps, ";
	for (int i=0;i<cn;i++) // build the maps first.
	{
		if (cellList[i]->Dead()) { cellid2idx.push_back(-1); continue; } // skip dead cells.
		else
		{
			idx2cellid.push_back(i);
			cellid2idx.push_back(iind);
			iind ++;
			for (int j = 0; j < biology->getNumRealSpecies(); j++)
			{
				v[j].push_back(0.0);
			}
		}
	}
	cout << "Done." << endl;
	
	// SIMULATION STARTS HERE
	// STEP 1: REACTIONS
	cout << "  1. Updating Cell States by Reactions, ";
	time_t time_0 = clock();
	int nsteps = 0;
	for (int i=0; i<cn; i++)
	{
		if (cellList[i]->Dead()) { continue; } // skip dead cells.
		//if (cellList[i]->get_cell_type() <= 4)
		//{
		//cout << ">>>cell reactions: " << i << endl;
			nsteps += biology->UpdateCellState(cellList[i]);
		//}
	}
	time_t time_1 = clock();
	cout << "Done. " << nsteps << " total steps in " << (time_1-time_0)*0.001 << " seconds." << endl;


	// STEP 2: Build the diffusion matrix
	for (int i=0; i<cn; i++)
	{
		if (cellList[i]->Dead()) { continue; } // skip dead cells.

		// Construct the diffusion matrix
		for (int j = 0; j < biology->getNumRealSpecies(); j++)
		{
			v[j][cellid2idx[cellList[i]->id()]] = (double)cellList[i]->get_cell_state(j);
		}
		// Each cells is a container of multiple growth factors and cytokines, which can diffuse through extracellular spaces.
		num_neighbors = cellList[i]->get_pairs()->size();
		for (int j = 0; j < biology->getNumRealSpecies(); j++)	{ vii[j] = 0.0; }
		for (int j = 0; j < num_neighbors; j++)
		{
			if (cellList[i]->get_pairs()->at(j)->get_redundant()) continue;
			if (cellList[i]->get_pairs()->at(j)->MP()->size() > 0)
			{
				pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
				pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
				if (pcell0 != cellList[i])
				{
					pcell1 = pcell0;
					pcell0 = cellList[i];
				}

				///// NEIGHBOR? /////
				commedgelen = 0.0;
				s1 = cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1();
				t1 = cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1();
				s2 = cellList[i]->get_pairs()->at(j)->MP()->at(0)->S2();
				t2 = cellList[i]->get_pairs()->at(j)->MP()->at(0)->T2();
				if (!s1->get_in_cell() ||
					!t1->get_in_cell() ||
					!s2->get_in_cell() ||
					!t2->get_in_cell())
				{
					cout<<"  -> error: wrong pair endpoints: "<<cellList[i]->get_pairs()->at(j)->id()<<": "<<cellList[i]->get_pairs()->at(j)->get_I1()->id()<<", "<<cellList[i]->get_pairs()->at(j)->get_I2()->id()<<endl;
					return i;
				}
				pp=s1;
				edge *edge1 = NULL;
				int flag = 0;
				while (pp!=t1)
				{
					pp=pp->getfp(); // getfp() is function in class point3D to get the point next to p in ccw direction.
					edge1 = pp->getab(); // getab() is funciton in class point3D to get the edge starting from p in ccw direction.
					commedgelen += edge1->getlength(); // calculate the length of edge1.
					flag++;
					if (flag==100)
					{
						cout<<"  -> error: wrong pair endpoints: "<<cellList[i]->get_pairs()->at(j)->id()<<": "<<cellList[i]->get_pairs()->at(j)->get_I1()->id()<<", "<<cellList[i]->get_pairs()->at(j)->get_I2()->id()<<endl;
						return i;
					}
				}
				if (s1==t1) {commedgelen = 0.1;}
				int idx_i = 0;
				int idx_j = 0;
				if (commedgelen > 0) ///// NEIGHBOR! /////
				{
					double distsq = 0;
					distsq =  (pcell0->get_center(0) - pcell1->get_center(0)) * (pcell0->get_center(0) - pcell1->get_center(0));
					distsq += (pcell0->get_center(1) - pcell1->get_center(1)) * (pcell0->get_center(1) - pcell1->get_center(1));
					idx_i = cellid2idx[pcell1->id()];
					idx_j = cellid2idx[i];
					ia.push_back (idx_i);
					ja.push_back (idx_j);
					for(int k=0; k<biology->getNumRealSpecies(); k++)
					{
						if (pcell1->get_cell_type() == 3 && (k == 0 || k == 1 || k == 2 || k == 7 || k == 8))
						{
							a[k].push_back(0.0);
						}
						else
						{
							/* Diffusion kernel, selectable via DYCELFEM_DIFFUSION.
							   The archived tree used "distsq"; committed trunk r517
							   used "edgelen".  D0 in mcommon.xml was calibrated for
							   one of them, so the choice changes the effective
							   diffusivity by ~3 orders of magnitude. */
							double value = 0;
							switch (dycelfem_diffusion_kernel())
							{
								case DIFF_EDGELEN:        // trunk r517
									value = biology->getD(k) * biology->getD0() * commedgelen;
									break;
								case DIFF_EDGELEN_DISTSQ:
									value = biology->getD(k) * biology->getD0() * commedgelen / distsq;
									break;
								default:                  // DIFF_DISTSQ, archived working copy
									value = biology->getD(k) * biology->getD0() / distsq;
									break;
							}
							a[k].push_back (value);
							vii[k] += value;
						}
					}
				}
			}
		}
		ia.push_back (cellid2idx[i]);
		ja.push_back (cellid2idx[i]);
		for(int k=0; k<biology->getNumRealSpecies(); k++)
		{
			a[k].push_back (-1 * vii[k]);
		}
	}
	cout << "  2. Building Diffusion Matrix, Done." << endl;
	// Diffusion matrix, DONE

	// STEP 3: SOLVING THE DIFFUSION EQUATION
	cout << "  3. Solving Diffusion Equation:" << endl;
	double err = 0;
	for (int i=0; i<biology->getNumRealSpecies(); i++)
	{
		//err = Euler (&ia, &ja, &a[i], &v[i], biology->getDeltaT());
		err = RK4 (&ia, &ja, &a[i], &v[i], biology->getDeltaT());
		cout << "       Diffusion Solved for " << biology->getSpeciesName(i) << ", err=" << err << endl;
	}
	cout << "     All Diffusions Done." << endl;

	// STEP 4: Update concentrations after diffusion
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) { continue; } // skip dead cells.
		cellList[i]->setup_Mark_number(0);

		for (int j=0; j<biology->getNumRealSpecies(); j++)
		{
			cellList[i]->set_cell_state (j, (int)v[j][cellid2idx[i]]);
		}
	}
	cout << "  4. Updating Cell States by Diffusion, Done." << endl;

	// IF FIBROBLASTS SYNTHESIZED ENOUGH PROCOLLAGEN, THEN THEY GENERATE NEW ECM ELEMENTS.
	// WILL FINISH IT LATER!!! CORRELATE ECM ELEMENT WITH COLLAGEN PRODUCTION
	/*for (int i=0; i<cn; i++)
	{
		if (cellList[i]->get_cell_type() == 2) // IF FIBROBLAST, 
		{
			int procol = cellList[i]->get_cell_state(7);
			if (procol > 200) // IF FIBROBLAST AND IT HAS ENOUGH PROCOLLAGEN
			{
				ra = (int)cellList[i]->shape_based_division_angle();
				cell_division(cellList[i], Time, cellList[i]->get_cell_type(), 0, ra);
				cellList[i]->center_refresh();
				cellList[cellList.size()-1]->center_refresh();
				cellList[cellList.size()-1]->set_cell_state(7, 100);
				cellList[i]->set_cell_state(7, procol - 100);
			}
		}
	}*/

	// STEP 5: Setup Cell Behaviors - A 4-way stochastic switch
	// (0) TBD; (1) Do nothing; (2) Growth/Division; (3) Apoptosis; (4) Migration.
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue; // skip dead cells.
		if (cellList[i]->get_behavior() != 0) continue; // Only make decisions for those waiting cells. 

		if (cellList[i]->get_cell_type () == 0 || 
			cellList[i]->get_cell_type () == 1 || 
			cellList[i]->get_cell_type () == 4 || 
			cellList[i]->get_cell_type () == 5) // BM, Clot, ECM, and hypodermis
		{
			cellList[i]->set_GR (0);  // They do not grow, do not migrate, do not die, and do nothing.
			// LOOK: NEED TO ADD MORE FANCY STUFF? TO CONTROL ECM AND CLOT BREAKDOWN BY MMPs?
		}
		else if (cellList[i]->get_cell_type () == 2 || 
			     cellList[i]->get_cell_type () == 3 ||
				 cellList[i]->get_cell_type () == 6) // Macrophage, Fibroblasts and Keratinocytes will do something
		{
			double sum = 0;
			biology->setValueProbDoNothing (biology->CalcValueProbDoNothing(cellList[i]));
			biology->setValueProbGrowthDiv (biology->CalcValueProbGrowthDiv(cellList[i]));
			biology->setValueProbMigration (biology->CalcValueProbMigration(cellList[i]));
			biology->setValueProbApoptosis (biology->CalcValueProbApoptosis(cellList[i]));
			sum = biology->getValueProbDoNothing() + biology->getValueProbGrowthDiv() + biology->getValueProbMigration() + biology->getValueProbApoptosis();
			biology->setValueProbDoNothing (biology->getValueProbDoNothing()/sum);
			biology->setValueProbGrowthDiv (biology->getValueProbGrowthDiv()/sum);
			biology->setValueProbMigration (biology->getValueProbMigration()/sum);
			biology->setValueProbApoptosis (biology->getValueProbApoptosis()/sum);
			double tau1 = ((double)rand()+1.0)/((unsigned long)RAND_MAX+1);
			// TEST
			/*if (cellList[i]->get_cell_type () == 2){
			// || cellList[i]->get_cell_type () == 3) {
			//if (i == 7105) {
			cout << "  >> Cell behavior prob:";
			printf(" Type=%d", cellList[i]->get_cell_type ());
			printf(" DoNothing=%0.2f", biology->getValueProbDoNothing()); 
			printf(" GrowthDiv=%0.2f", biology->getValueProbGrowthDiv()); 
			printf(" Migration=%0.2f", biology->getValueProbMigration()); 
			printf(" Apoptosis=%0.2f", biology->getValueProbApoptosis()); 
			cout << endl;}*/
			//
			if (tau1 >= 0 && tau1 < biology->getValueProbGrowthDiv())
			{
				// Growth and Division
				cellList[i]->set_behavior(2);
				cellList[i]->set_bhsteps(2,0);
			}
			else if (tau1 >= biology->getValueProbGrowthDiv() && tau1 < biology->getValueProbGrowthDiv() + biology->getValueProbMigration())
			{
				// Migration
				cellList[i]->set_behavior(4);
				cellList[i]->set_bhsteps(4,0);
			}
			else if (tau1 >= biology->getValueProbGrowthDiv() + biology->getValueProbMigration() && tau1 < 1 - biology->getValueProbDoNothing())
			{
				// Apoptosis
				cellList[i]->set_behavior(3);
				cellList[i]->set_bhsteps(3,0);
			}
			else if (tau1 >= 1 - biology->getValueProbDoNothing() && tau1 <= 1.0)
			{
				// Do Nothing
				cellList[i]->set_behavior(1);
				cellList[i]->set_bhsteps(1,0);
			}
		}
	}
	cout << "  5. Setting Up Cell Behaviors, Done." << endl;

	//cout << "    >>>>>> CELL 7105: " << cellList[7105]->get_behavior() << endl;
			
	// Step 6: Cell Apoptosis
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue; // skip dead cells.

		if (cellList[i]->get_cell_type () == 0) // A ECM transforms into a macrophage. Not now.
		{
			if (cellList[i]->get_behavior() == 3)
			{
				// Cell Apoptosis: ECM -> Macrophage previously. Now gone.
				//cellList[i]->set_cell_type(6);
				//cellList[i]->set_behavior(1);
				//cellList[i]->set_bhsteps(1,0);
			}
		}
		else if (cellList[i]->get_cell_type () == 3) // Keratinocytes 
		{
			if (cellList[i]->get_behavior() == 3)
			{
				// Cell Apoptosis
				cellList[i]->set_set_dead(1); // DO NOT USE setup_Dead() TO KILL A CELL.
			}
		}
		else if (cellList[i]->get_cell_type () == 2) // Fibroblasts
		{
			if (cellList[i]->get_behavior() == 3)
			{
				cellList[i]->set_cell_type(0);
			}
		}
		else if (cellList[i]->get_cell_type () == 6) // Macrophage
		{
			// ASSUME MACROPHAGE NOT TO DIE, SO COMMENT OUT FOLLOWING LINES
			/*if (cellList[i]->get_behavior() == 3)
			{
				// At least live for 5 steps.
				if (cellList[i]->get_bhsteps(1) > 5) 
				{
					cellList[i]->set_cell_type(0);
					cellList[i]->set_behavior(0);
					cellList[i]->set_bhsteps(0,0);
					cellList[i]->set_bhsteps(1,0);
				}
				else
				{
					cellList[i]->set_bhsteps(1, cellList[i]->get_bhsteps(1)+1);
				}
			}
			else if (cellList[i]->get_behavior() == 1)
			{
				cellList[i]->set_bhsteps(1, cellList[i]->get_bhsteps(1)+1);
			}*/
		}
	}
	cout << "  6. Doing Cell Apoptosis, Done." << endl;

	// Step 7: Cell Division
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue; // skip dead cells.

		if (cellList[i]->get_cell_type () == 2 || 
			cellList[i]->get_cell_type () == 3 ||
			cellList[i]->get_cell_type () == 6) // Macrophage, Fibroblasts and Keratinocytes
		{
			if (cellList[i]->get_behavior() == 2)
			{
				// CELL DIVISION AT DOUBLE SIZE
				if (cellList[i]->get_area()>=AR)
				{
					ra = (int)cellList[i]->shape_based_division_angle();
					cellList[i]->set_behavior(0);
					cellList[i]->set_bhsteps(0,0);
					cout<<"  -> ";
					// All divisions are self-renewal type
					cell_division(cellList[i], Time, cellList[i]->get_cell_type(), cellList[i]->get_cell_type(), ra);
					cellList[i]->center_refresh();
					cellList[cellList.size()-1]->center_refresh();
					for (int j=0; j<biology->getNumRealSpecies(); j++)
					{
						cellList[cellList.size()-1]->set_cell_state(j, (int)(cellList[i]->get_cell_state(j)/2.0));
						cellList[i]->set_cell_state(j, (int)(cellList[i]->get_cell_state(j)/2.0));
					}
				}
			}
		}

		//
		if (cellList[i]->get_cell_type() == 0)
		{
			if (cellList[i]->get_area() >= AR)
			{
				//cout << ">>>>TEST<<<< COL_DIV " << cellList[i]->get_area() << endl;
				ra = (int)cellList[i]->shape_based_division_angle();
				cell_division(cellList[i], Time, cellList[i]->get_cell_type(), cellList[i]->get_cell_type(), ra);
				cellList[i]->center_refresh();
				cellList[cellList.size()-1]->center_refresh();
				for (int j=0; j<biology->getNumRealSpecies(); j++)
				{
					cellList[cellList.size()-1]->set_cell_state(j, (int)(cellList[i]->get_cell_state(j)/2.0));
					cellList[i]->set_cell_state(j, (int)(cellList[i]->get_cell_state(j)/2.0));
				}
			}
		}
	}
	cout << "  7. Doing Cell Divisions, Done." << endl;

	// Step 8: Cell Growth
	cn = (int)cellList.size();
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue; // skip dead cells.

		if (cellList[i]->get_cell_type () == 2 || 
			cellList[i]->get_cell_type () == 3 ||
			cellList[i]->get_cell_type () == 6) // Macrophage, Fibroblasts and Keratinocytes
		{
			if (cellList[i]->get_behavior() == 2)
			{
				// Cell Growth and Division
				cellList[i]->set_GR(GAR*4.5);
				//if (cellList[i]->get_cell_type()==3) {cout<<"  -> Kera cell "<<i<<" is growing!"<<endl;}
				if (cellList[i]->get_bhsteps(2) > 5) // Grow for 10 steps.
				{
					cellList[i]->set_behavior(0);
					cellList[i]->set_bhsteps(0,0);
				}
				else
				{
					cellList[i]->set_bhsteps(2, cellList[i]->get_bhsteps(2)+1);
				}
			}
		}

		// 
		if (cellList[i]->get_cell_type() == 0)
		{
			if (cellList[i]->get_cell_state(2)/1000 > cellList[i]->get_area()/AR0) // ECM VOLUME SHOULD REFLECT ITS MASS OF COLLAGEN.
			{
				//cout << ">>>>TEST<<<< COL " << cellList[i]->get_cell_state(2) << " AREA " << cellList[i]->get_area() << endl;
				cellList[i]->set_GR(AR0 * (cellList[i]->get_cell_state(2)/1000 - cellList[i]->get_area()/AR0));
			}
		}/**/
	}
	cout << "  8. Doing Cell Growths, Done." << endl;

	// STEP 9: Cell Migration
	for (int i=0; i<cn; i++)
	{
		if (cellList[i]->Dead()) continue; // skip dead cells.
		if (cellList[i]->get_cell_type () == 2) // FIBROBLASTS MIGRATION
		{
			if (cellList[i]->get_behavior() == 4)
			{
				cellList[i]->set_GR(0.0);
				// Cell MIGRATION
				int MIGCYTOKINEID = 3;  // 3: PDGF
				int maxcytokine_id = -1;
				interpair* maxcytokine_pair = 0;
				double maxcytokine = 0;
				int num_neighbors = cellList[i]->get_pairs()->size();
				double cdist = 0.0;
				double cdist_maxcytokine = 0.0;
				vector<int> kct;
				kct.push_back(0); kct.push_back(1); kct.push_back(4);
				for (int j = 0; j < num_neighbors; j++)
				{
					if (cellList[i]->get_pairs()->at(j)->get_redundant()) continue;
					if (cellList[i]->get_pairs()->at(j)->MP()->size() > 0)  // Direct contact between cell i and other cell in neighbor j.
					{
						pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
						pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
						if (pcell0 != cellList[i])
						{
							pcell1 = pcell0;
							pcell0 = cellList[i];
						}
						// 3: PDGF
						if ((pcell1->get_cell_type() <= 1 || pcell1->get_cell_type() == 4) && pcell1->get_cell_state(MIGCYTOKINEID) > maxcytokine)
						{
							if (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()==cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1())
							{
								double b_x = pcell1->get_center(0);
								double b_y = pcell1->get_center(1);
								double a_x = pcell0->get_center(0);
								double a_y = pcell0->get_center(1);
								double angle_tmp = vector2angle((b_x - a_x), (b_y - a_y));
								cell *Cod_tmp = cell_migrate_angle(pcell0, angle_tmp);
								if (Cod_tmp==NULL ||
								   (Cod_tmp!=NULL && find(kct.begin(),kct.end(),Cod_tmp->get_cell_type())!=kct.end()))
								{
									maxcytokine_id = pcell1->id();
									maxcytokine_pair = cellList[i]->get_pairs()->at(j);
									maxcytokine = pcell1->get_cell_state(MIGCYTOKINEID);
									cdist_maxcytokine = 0.0;
								}
							}
							else
							{
								maxcytokine_id = pcell1->id();
								maxcytokine_pair = cellList[i]->get_pairs()->at(j);
								maxcytokine = pcell1->get_cell_state(MIGCYTOKINEID);
								cdist_maxcytokine = 0.0;
							}
						}
					}
					else if (cellList[i]->get_pairs()->at(j)->MP()->size() == 0)  // No direct contact between cell i and other cell in neighbor j.
					{
						pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
						pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
						if (pcell0 != cellList[i])
						{
							pcell1 = pcell0;
							pcell0 = cellList[i];
						}
						cdist = cell_distance (pcell0, pcell1);
						if (cdist < CR*3)
						{
							//if ((pcell1->get_cell_type() <= 1 || pcell1->get_cell_type() == 4) && pcell1->get_cell_state(MIGCYTOKINEID)/(1 + cdist) > maxcytokine)
							if ((pcell1->get_cell_type() <= 1 || pcell1->get_cell_type() == 4) && pcell1->get_cell_state(MIGCYTOKINEID) > maxcytokine)
							{
								double b_x = pcell1->get_center(0);
								double b_y = pcell1->get_center(1);
								double a_x = pcell0->get_center(0);
								double a_y = pcell0->get_center(1);
								double angle_tmp = vector2angle((b_x - a_x), (b_y - a_y));
								cell *Cod_tmp = cell_migrate_angle(pcell0, angle_tmp);
								if (Cod_tmp==NULL ||
								   (Cod_tmp!=NULL && find(kct.begin(),kct.end(),Cod_tmp->get_cell_type())!=kct.end()))
								{
									maxcytokine_id = pcell1->id();
									maxcytokine_pair = cellList[i]->get_pairs()->at(j);
									maxcytokine = pcell1->get_cell_state(MIGCYTOKINEID);
									cdist_maxcytokine = cdist;
								}
							}
						}
					}
				}
				if (maxcytokine_id >= 0)
				{
					double b_x = 0.0, b_y = 0.0, a_x = 0.0, a_y = 0.0;
					double angle = 0.0;
					double angle_ow = 0.0;
					if (cdist_maxcytokine == 0.0)
					{
						if (maxcytokine_pair->MP()->at(0)->S1()!=maxcytokine_pair->MP()->at(0)->T1())
						{
							b_x = (maxcytokine_pair->MP()->at(0)->S1()->x() + maxcytokine_pair->MP()->at(0)->T1()->x())/2.0;
							b_y = (maxcytokine_pair->MP()->at(0)->S1()->y() + maxcytokine_pair->MP()->at(0)->T1()->y())/2.0;
						}
						else
						{
							b_x = cellList[maxcytokine_id]->get_center(0);
							b_y = cellList[maxcytokine_id]->get_center(1);
						}
						a_x = cellList[i]->get_center(0);
						a_y = cellList[i]->get_center(1);
					}
					else if (cdist_maxcytokine > 0.0)
					{
						b_x = cellList[maxcytokine_id]->get_center(0);
						b_y = cellList[maxcytokine_id]->get_center(1);
						a_x = cellList[i]->get_center(0);
						a_y = cellList[i]->get_center(1);
					}
					angle = vector2angle ((b_x - a_x), (b_y - a_y));
					angle_ow = angle;
					cellList[i]->set_migrate(1);
					cellList[i]->set_migrate_relax(1);
					//double force = max(0.0, (maxcytokine - (double)cellList[i]->get_cell_state(MIGCYTOKINEID))); // Variable force depending on signal strength
					//double force = GAR*5; // Constant initial force
					double force = GAR*2.5; // Constant initial force
					cell *Cod = cell_migrate_angle(cellList[i], angle);
					if (cdist_maxcytokine == 0.0 && Cod!=NULL && Cod==cellList[maxcytokine_id])
					{
						cellList[maxcytokine_id]->set_set_dead(1);
						//cout<<"   * -> cell "<<i<<" migrates towards dying cell "<<maxcytokine_id<<endl;
					}
					else if (cdist_maxcytokine == 0.0 && Cod!=NULL && Cod!=cellList[maxcytokine_id])
					{
						Cod->set_set_dead(1);
						int cipn = (int)cellList[i]->get_pairs()->size();
						for (int kk=0;kk<cipn;kk++)
						{
							if ((cellList[i]->get_pairs()->at(kk)->get_I1()==cellList[i] &&
								 cellList[i]->get_pairs()->at(kk)->get_I2()==Cod) ||
								(cellList[i]->get_pairs()->at(kk)->get_I2()==cellList[i] &&
								 cellList[i]->get_pairs()->at(kk)->get_I1()==Cod))
							{
								b_x = (cellList[i]->get_pairs()->at(kk)->MP()->at(0)->S1()->x() + cellList[i]->get_pairs()->at(kk)->MP()->at(0)->T1()->x())/2.0;
								b_y = (cellList[i]->get_pairs()->at(kk)->MP()->at(0)->S1()->y() + cellList[i]->get_pairs()->at(kk)->MP()->at(0)->T1()->y())/2.0;
								angle_ow = vector2angle ((b_x - a_x), (b_y - a_y));
							}
						}
					}
					else if (cdist_maxcytokine>0 && Cod!=NULL && Cod!=cellList[maxcytokine_id])
					{
						Cod->set_set_dead(1);
						int cipn = (int)cellList[i]->get_pairs()->size();
						for (int kk=0;kk<cipn;kk++)
						{
							if ((cellList[i]->get_pairs()->at(kk)->get_I1()==cellList[i] &&
								 cellList[i]->get_pairs()->at(kk)->get_I2()==Cod) ||
								(cellList[i]->get_pairs()->at(kk)->get_I2()==cellList[i] &&
								 cellList[i]->get_pairs()->at(kk)->get_I1()==Cod))
							{
								b_x = (cellList[i]->get_pairs()->at(kk)->MP()->at(0)->S1()->x() + cellList[i]->get_pairs()->at(kk)->MP()->at(0)->T1()->x())/2.0;
								b_y = (cellList[i]->get_pairs()->at(kk)->MP()->at(0)->S1()->y() + cellList[i]->get_pairs()->at(kk)->MP()->at(0)->T1()->y())/2.0;
								angle_ow = vector2angle ((b_x - a_x), (b_y - a_y));
							}
						}
						//cout<<"   * -> cell "<<Cod->id()<<" on the way of cell "<<i<<" towards to cell "<<maxcytokine_id<<endl;
					}
					else if (cdist_maxcytokine>0 && Cod==NULL)
					{
						//cout<<"   * -> cell "<<i<<" migrates towards solid cell "<<maxcytokine_id<<" through NULL"<<endl;
					}
					cellList[i]->set_migrate_angle(angle_ow);
					cell_migrate_response_1(cellList[i], angle_ow, &kct,  2*force); // MIGRATION SPEED
					//cellList[i]->set_GR(1.236 * force);
					cellList[i]->set_GR(1.8*force);

					// Fibroblasts deposite new ECMs
					ra = (int)angle_ow + 90;
					if (ra < 0) { ra = ra + 360; }
					else if (ra > 360) { ra = ra - 360; }
					if (cellList[i]->get_area() > AR)// && cellList[i]->get_cell_state(7) > 200)
					{
						cell_division(cellList[i], Time, cellList[i]->get_cell_type(), cellList[i]->get_cell_type(), ra);
						cellList[i]->center_refresh();
						cellList[cellList.size()-1]->center_refresh();
						double dir_x = cellList[i]->get_center(0) - cellList[cellList.size()-1]->get_center(0);
						double dir_y = cellList[i]->get_center(1) - cellList[cellList.size()-1]->get_center(1);
						double an = vector2angle(dir_x, dir_y);
						//cout<<"angle: "<<angle_ow<<", an: "<<an<<", cell "<<i<<" migrating: "<<cellList[i]->Migrate()<<endl;
						/*double cos_angle_ow = cos(angle_ow*PI/180.0);
						double sin_angle_ow = sin(angle_ow*PI/180.0);
						double cos_an = cos(an*PI/180.0);
						double sin_an = sin(an*PI/180.0);
						double plot_an = cos_angle_ow*cos_an + sin_angle_ow*sin_an;*/
						double procoldensity = 0.0;
						if      (abs(an-angle_ow) < 90 || abs(an-angle_ow) > 270)
						{
							procoldensity = cellList[i]->get_cell_state (7) / 2.0;
							cellList[cellList.size()-1]->set_cell_type(0);
							cellList[cellList.size()-1]->set_cell_state(7, (int)procoldensity); // Deposit new collagen
							cellList[i]->set_cell_state(7, (int)procoldensity);
						}
						else if (abs(an-angle_ow) >= 90 && abs(an-angle_ow)<=270)
						{
							procoldensity = cellList[cellList.size()-1]->get_cell_state (7) / 2.0;
							cellList[i]->set_cell_type(0);
							cellList[i]->set_cell_state(7, (int)procoldensity); // Deposit new collagen
							cellList[cellList.size()-1]->set_cell_state(7, (int)procoldensity);
							cellList[i]->setup_Mark_number(0);
							cellList[cellList.size()-1]->setup_Mark_number(1);
							cellList[cellList.size()-1]->set_migrate_angle(angle_ow);
							if (cellList[i]->Migrate())
							{
								cellList[i]->set_migrate(0);
								if (cellList[cellList.size()-1]->get_migrate_p(0)!=NULL &&
									cellList[cellList.size()-1]->get_migrate_p(1)!=NULL)
								{
									cellList[cellList.size()-1]->set_migrate(1);
								}
							}
						}
					}
					//// END deposition of new ECMs

					if (cellList[i]->get_bhsteps(4) > MIN_FBEHAVIOR_PERSISTENCE) // At least migrate for 30 steps.
					{
						cellList[i]->set_behavior(0);
						cellList[i]->set_bhsteps(0,0);
					}
					else
					{
						cellList[i]->set_bhsteps(4, cellList[i]->get_bhsteps(4)+1);
					}
				}
			}
		}
		else if (cellList[i]->get_cell_type () == 3) // KERATINOCYTES MIGRATION
		{
			/* NEW K MIGRATION MODEL
			// NEW K MIGRATION USING MECHANICAL SIGNALS
			// MIGRATION SPEED COMES FROM CYTOKINE, BUT DIRECTION COMES FROM MECHANICAL SIGNALS
			if (cellList[i]->get_behavior() == 4)
			{

				// Cell MIGRATION
				cellList[i]->set_GR(0.0);
				int MIGCYTOKINEID = 4;  // 4: KGF
				// Step 1: Get the migration direction
				int flag_mig = -1;
				int killcell_id = -1;
				double cont_cytokine = -1;
				int num_neighbors = cellList[i]->get_pairs()->size();
				vector<vector<int> > cell_angles;
				vector<int>          id_angle;
				id_angle.push_back(0);
				id_angle.push_back(0);
				double migangle = 0;
				double cdist = 0.0;
				vector<int> kct;
				kct.push_back(4); // Clot
				//kct.push_back(0);
				// For keratinocytes, migration direction is determined by the boundary between collagen and clot.
				for (int j = 0; j < num_neighbors; j++)  // SCAN NEIGHBORING CELLS, BOTH CONTACTING AND NON-CONTACTING, RECORD ANGLES AND IDS
				{
					if (cellList[i]->get_pairs()->at(j)->get_redundant()) continue;
					if (cellList[i]->get_pairs()->at(j)->MP()->size() > 0)  // CONTACTING NEIGHBORING CELLS
					{
						pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
						pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
						if (pcell0 != cellList[i])
						{
							pcell1 = pcell0;
							pcell0 = cellList[i];
						}
						double a_x = pcell0->get_center(0);
						double a_y = pcell0->get_center(1);
						if (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()==cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1())
						{
							double b_x = (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()->x() + cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1()->x())/2.0;
							double b_y = (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()->y() + cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1()->y())/2.0;
							double angle = vector2angle ((b_x - a_x), (b_y - a_y));
							cell_angles.push_back(vector<int>(pcell1->id(), (int)angle));
						}
						else
						{
							double b_x = pcell1->get_center(0);
							double b_y = pcell1->get_center(1);
							double angle = vector2angle ((b_x - a_x), (b_y - a_y));
							id_angle[0] = pcell1->id();
							id_angle[1] = (int)angle;
							cell_angles.push_back(id_angle);
						}
					}
					else if (cellList[i]->get_pairs()->at(j)->MP()->size() == 0)  // NON-CONTACTING NEIGHBORING CELLS
					{
						pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
						pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
						if (pcell0 != cellList[i])
						{
							pcell1 = pcell0;
							pcell0 = cellList[i];
						}
						cdist = cell_distance (pcell0, pcell1);
						if (cdist < CR*3)  // NON-BLOCKED CELLS WITH 3 RADIA ARE CONSIDERED NON-CONTACTING NEIGHBORING CELLS
						{
							double a_x = pcell0->get_center(0);
							double a_y = pcell0->get_center(1);
							double b_x = pcell1->get_center(0);
							double b_y = pcell1->get_center(1);
							double angle = vector2angle ((b_x - a_x), (b_y - a_y));
							id_angle[0] = pcell1->id();
							id_angle[1] = (int)angle;
							cell_angles.push_back(id_angle);
						}
					}
				}
				sort(cell_angles.begin(), cell_angles.end(), SortNeighborCellsByAngles);  // ASCENDING SORT ACCORDING TO ANGLES
				cell_angles.push_back(cell_angles[0]);
				// FIND OUT THE BOUNDARY BETWEEN (ECMS,FIBROBLASTS,MACROPHAGES) AND (CLOTS,KERATINOCYTES)
				int num_clot = 0;
				for (unsigned int j=0; j<cell_angles.size()-1; j++)
				{
					// Find neighboring ECM and Clot cells.
					int celltype0 = cellList[cell_angles[j][0]]->get_cell_type();
					int celltype1 = cellList[cell_angles[j+1][0]]->get_cell_type();
					if (celltype0 == 4)
					{
						num_clot ++;
					}
					if ((celltype0 <= 2 || celltype0 == 6 ) && celltype1 == 4)  // IF NEXT NEIGHBOR IS CLOT, MIGRATE TO THAT DIRECTION
					{
						migangle = cell_angles[j+1][1];
						flag_mig = 1;
						killcell_id = cell_angles[j+1][0];
					}
					else if ((celltype1 <= 2 || celltype1 == 6 ) && (celltype0 == 4))  // IF NEXT NEIGHBOR IS CLOT, MIGRATE TO THAT DIRECTION
					{
						migangle = cell_angles[j][1];
						flag_mig = 1;
						killcell_id = cell_angles[j][0];
					}
				}

				if (num_clot == 0) // OTHERWISE, LET'S LOOK AT THE MECHANICAL SIGNAL
				{
					if (cellList[i]->get_sense_migr())
					{
						migangle = vector2angle(cellList[i]->get_sense_migr_vec (0), cellList[i]->get_sense_migr_vec (1));
						flag_mig = 1;
						killcell_id = -1;
					}
				}
				cout << "    >>>>>>>>>>>>>>> MIGANGLE <<<<<<: " << migangle << endl;
				// Start to migrate
				if (flag_mig == 1)
				{
					cellList[i]->set_migrate(1);
					cellList[i]->set_migrate_relax(1);
					cellList[i]->set_migrate_angle(migangle);
					cell *Cod = cell_migrate_angle(cellList[i],migangle);
					bool degtype = false;
					if (Cod!=NULL)
					{
						if (find(kct.begin(),kct.end(),Cod->get_cell_type())!=kct.end()) 
						{degtype = true;}
					}

					if (killcell_id >= 0) // IF THERE IS A CLOT CELL NEEDS TO BE DESTROYED, THEN DESTROY; OTHERWISE JUST MIGRATE.
					{
						for (int k=0;k<(int)cellList[i]->get_pairs()->size();k++)
						{
							if ((cellList[i]->get_pairs()->at(k)->get_I1()==cellList[i] &&
								 cellList[i]->get_pairs()->at(k)->get_I2()==cellList[killcell_id]) ||
								(cellList[i]->get_pairs()->at(k)->get_I2()==cellList[i] &&
								 cellList[i]->get_pairs()->at(k)->get_I1()==cellList[killcell_id]))
							{
								if (cellList[i]->get_pairs()->at(k)->MP()->size()==0)
								{
									cdist = cell_distance(cellList[i],cellList[killcell_id]);
								}
								else
								{
									cdist = 0.0;
								}
							}
						}
						if (cdist == 0.0 && Cod!=NULL && Cod==cellList[killcell_id])
						{ 
							if (degtype) {cellList[killcell_id]->set_set_dead(1);}
							//cout<<"    * -> cell "<<i<<" migrates towards dying cell "<<maxcytokine_id<<endl;
						}
						else if (cdist == 0.0 && Cod!=NULL && Cod!=cellList[killcell_id])
						{
							if (degtype) {Cod->set_set_dead(1);}
							//cout<<"    * -> cell "<<Cod->id()<<" on the way of cell "<<i<<" towards to cell "<<maxcytokine_id<<endl;
						}
						else if (cdist>0 && Cod!=NULL && Cod!=cellList[killcell_id])
						{
							if (degtype) {Cod->set_set_dead(1);}
							//cout<<"    * -> cell "<<Cod->id()<<" on the way of cell "<<i<<" towards to cell "<<maxcytokine_id<<endl;
						}
						else if (cdist>0 && Cod==NULL)
						{
							//cout<<"    * -> cell "<<i<<" migrates towards solid cell "<<maxcytokine_id<<" through NULL"<<endl;
						}
					}

					// Step 2: Define the migration force
					// DO NOT PUSH FRONT KERATINOCYTE TOO HARD
					double ffactor = 1.0;
					if (Cod!=NULL && Cod->get_cell_type() == 3) // CELLTYPE 3: KERATINOCYTES
					{
						ffactor = 0.2;
					}
					//double force = ffactor * max(0.0, min(GAR, (double)cellList[i]->get_cell_state(MIGCYTOKINEID)));
					double force = ffactor * biology->CalcValueMigForce(cellList[i]) * GAR;
					cell_migrate_response_1(cellList[i], migangle, &kct, force);
					//cellList[i]->set_GR(0.6 * force);
					cellList[i]->set_GR(force);
					//cout << "   *-> K MIGRATION FORCE: " << force << endl;
					//// KERATINOCYTE DIVIDE AT DOUBLE SIZE
					ra = (int)migangle + 90;
					if (ra < 0) { ra = ra + 360; }
					else if (ra > 360) { ra = ra - 360; }
					if (cellList[i]->get_area() > AR)
					{
						cell_division(cellList[i], Time, cellList[i]->get_cell_type(), cellList[i]->get_cell_type(), ra);
						cellList[i]->center_refresh();
						cellList[cellList.size()-1]->center_refresh();
					}
					//// END ////
					if (cellList[i]->get_bhsteps(4) > 5) // At least migrate for 5 steps.
					{
						cellList[i]->set_behavior(0);
						cellList[i]->set_bhsteps(0,0);
					}
					else
					{
						cellList[i]->set_bhsteps(4, cellList[i]->get_bhsteps(4)+1);
					}
				}
			}*/

			/* OLD K MIGRATION WITHOUT MECHANICAL SIGNALS */
			if (cellList[i]->get_behavior() == 4)
			{
				cellList[i]->set_GR(0.0);
				// Cell MIGRATION
				int MIGCYTOKINEID = 4;  // 4: KGF
				int maxcytokine_id = -1;
				double maxcytokine = -1;
				int num_neighbors = cellList[i]->get_pairs()->size();
				vector<vector<int> > cell_angles;
				vector<int>          id_angle;
				id_angle.push_back(0);
				id_angle.push_back(0);
				double migangle = 0;
				double cdist = 0.0;
				vector<int> kct;
				kct.push_back(4);
				//kct.push_back(0);
				// For keratinocytes, migration direction is determined by the boundary between collagen and clot.
				for (int j = 0; j < num_neighbors; j++)  // SCAN NEIGHBORING CELLS, BOTH CONTACTING AND NON-CONTACTING, RECORD ANGLES AND IDS
				{
					if (cellList[i]->get_pairs()->at(j)->get_redundant()) continue;
					if (cellList[i]->get_pairs()->at(j)->MP()->size() > 0)  // CONTACTING NEIGHBORING CELLS
					{
						pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
						pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
						if (pcell0 != cellList[i])
						{
							pcell1 = pcell0;
							pcell0 = cellList[i];
						}
						double a_x = pcell0->get_center(0);
						double a_y = pcell0->get_center(1);
						if (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()==cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1())
						{
							double b_x = (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()->x() + cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1()->x())/2.0;
							double b_y = (cellList[i]->get_pairs()->at(j)->MP()->at(0)->S1()->y() + cellList[i]->get_pairs()->at(j)->MP()->at(0)->T1()->y())/2.0;
							double angle = vector2angle ((b_x - a_x), (b_y - a_y));
							cell_angles.push_back(vector<int>(pcell1->id(), (int)angle));
						}
						else
						{
							double b_x = pcell1->get_center(0);
							double b_y = pcell1->get_center(1);
							double angle = vector2angle ((b_x - a_x), (b_y - a_y));
							id_angle[0] = pcell1->id();
							id_angle[1] = (int)angle;
							cell_angles.push_back(id_angle);
						}
					}
					else if (cellList[i]->get_pairs()->at(j)->MP()->size() == 0)  // NON-CONTACTING NEIGHBORING CELLS
					{
						pcell0 = cellList[i]->get_pairs()->at(j)->get_I1();
						pcell1 = cellList[i]->get_pairs()->at(j)->get_I2();
						if (pcell0 != cellList[i])
						{
							pcell1 = pcell0;
							pcell0 = cellList[i];
						}
						cdist = cell_distance (pcell0, pcell1);
						if (cdist < CR*3)  // NON-BLOCKED CELLS WITH 3 RADIA ARE CONSIDERED NON-CONTACTING NEIGHBORING CELLS
						{
							double a_x = pcell0->get_center(0);
							double a_y = pcell0->get_center(1);
							double b_x = pcell1->get_center(0);
							double b_y = pcell1->get_center(1);
							double angle = vector2angle ((b_x - a_x), (b_y - a_y));
							id_angle[0] = pcell1->id();
							id_angle[1] = (int)angle;
							cell_angles.push_back(id_angle);
						}
					}
				}
				sort(cell_angles.begin(), cell_angles.end(), SortNeighborCellsByAngles);  // ASCENDING SORT ACCORDING TO ANGLES
				cell_angles.push_back(cell_angles[0]);
				for (unsigned int j=0; j<cell_angles.size()-1; j++)  // FIND OUT THE BOUNDARY BETWEEN (ECMS,FIBROBLASTS,MACROPHAGES) AND (CLOTS,KERATINOCYTES)
				{
					// Find neighboring ECM and Clot cells.
					int celltype0 = cellList[cell_angles[j][0]]->get_cell_type();
					int celltype1 = cellList[cell_angles[j+1][0]]->get_cell_type();
					if ((celltype0 <= 2) && (celltype1 == 3 || celltype1 == 4 || celltype1 == 6))
					{
						if (celltype1 == 4 || celltype1 == 6)  // IF NEXT NEIGHBOR IS CLOT, MIGRATE TO THAT DIRECTION
						{
							migangle = cell_angles[j+1][1];
							maxcytokine = cellList[cell_angles[j+1][0]]->get_cell_state(MIGCYTOKINEID);
							maxcytokine_id = cell_angles[j+1][0];
							if (celltype1 == 6)
							{
								cellList[cell_angles[j+1][0]]->set_cell_type(4);
							}
						}
						else  // ELSE, LET'S COMPARE KGF GRADIENT
						{
							if (maxcytokine_id < 0)
							{
								if (cellList[cell_angles[j+1][0]]->get_cell_state(MIGCYTOKINEID) > maxcytokine)
								{
									migangle = cell_angles[j+1][1];
									maxcytokine = cellList[cell_angles[j+1][0]]->get_cell_state(MIGCYTOKINEID);
									maxcytokine_id = cell_angles[j+1][0];
								}
							}
							else
							{
								if (cellList[cell_angles[j+1][0]]->get_cell_state(MIGCYTOKINEID) > maxcytokine && cellList[maxcytokine_id]->get_cell_type() != 4)
								{
									migangle = cell_angles[j+1][1];
									maxcytokine = cellList[cell_angles[j+1][0]]->get_cell_state(MIGCYTOKINEID);
									maxcytokine_id = cell_angles[j+1][0];
								}
							}
						}
					}
					else if ((celltype1 <= 2) && (celltype0 == 3 || celltype0 == 4 || celltype0 == 6))
					{
						if (celltype0 == 4)  // IF NEXT NEIGHBOR IS CLOT, MIGRATE TO THAT DIRECTION
						{
							migangle = cell_angles[j][1];
							maxcytokine = cellList[cell_angles[j][0]]->get_cell_state(MIGCYTOKINEID);
							maxcytokine_id = cell_angles[j][0];
							if (celltype0 == 6)
							{
								cellList[cell_angles[j][0]]->set_cell_type(4);
							}
						}
						else
						{
							if (maxcytokine_id < 0)
							{
								if (cellList[cell_angles[j][0]]->get_cell_state(MIGCYTOKINEID) > maxcytokine)
								{
									migangle = cell_angles[j][1];
									maxcytokine = cellList[cell_angles[j][0]]->get_cell_state(MIGCYTOKINEID);
									maxcytokine_id = cell_angles[j][0];
								}
							}
							else
							{
								if (cellList[cell_angles[j][0]]->get_cell_state(MIGCYTOKINEID) > maxcytokine && cellList[maxcytokine_id]->get_cell_type() != 4)
								{
									migangle = cell_angles[j][1];
									maxcytokine = cellList[cell_angles[j][0]]->get_cell_state(MIGCYTOKINEID);
									maxcytokine_id = cell_angles[j][0];
								}
							}
						}
					}
				}
				if (maxcytokine_id >= 0)
				{
					cellList[i]->set_migrate(1);
					cellList[i]->set_migrate_relax(1);
					cellList[i]->set_migrate_angle(migangle);
					cell *Cod = cell_migrate_angle(cellList[i],migangle);
					bool degtype = false;
					if (Cod!=NULL)
					{
						if (find(kct.begin(),kct.end(),Cod->get_cell_type())!=kct.end()) 
						{degtype = true;}
					}
					for (int k=0;k<(int)cellList[i]->get_pairs()->size();k++)
					{
						if ((cellList[i]->get_pairs()->at(k)->get_I1()==cellList[i] &&
							 cellList[i]->get_pairs()->at(k)->get_I2()==cellList[maxcytokine_id]) ||
							(cellList[i]->get_pairs()->at(k)->get_I2()==cellList[i] &&
							 cellList[i]->get_pairs()->at(k)->get_I1()==cellList[maxcytokine_id]))
						{
							if (cellList[i]->get_pairs()->at(k)->MP()->size()==0)
							{
								cdist = cell_distance(cellList[i],cellList[maxcytokine_id]);
							}
							else
							{
								cdist = 0.0;
							}
						}
					}
					if (cdist == 0.0 && Cod!=NULL && Cod==cellList[maxcytokine_id])
					{ 
						if (degtype) {cellList[maxcytokine_id]->set_set_dead(1);}
						//cout<<"    * -> cell "<<i<<" migrates towards dying cell "<<maxcytokine_id<<endl;
					}
					else if (cdist == 0.0 && Cod!=NULL && Cod!=cellList[maxcytokine_id])
					{
						if (degtype) {Cod->set_set_dead(1);}
						//cout<<"    * -> cell "<<Cod->id()<<" on the way of cell "<<i<<" towards to cell "<<maxcytokine_id<<endl;
					}
					else if (cdist>0 && Cod!=NULL && Cod!=cellList[maxcytokine_id])
					{
						if (degtype) {Cod->set_set_dead(1);}
						//cout<<"    * -> cell "<<Cod->id()<<" on the way of cell "<<i<<" towards to cell "<<maxcytokine_id<<endl;
					}
					else if (cdist>0 && Cod==NULL)
					{
						//cout<<"    * -> cell "<<i<<" migrates towards solid cell "<<maxcytokine_id<<" through NULL"<<endl;
					}

					// DO NOT PUSH FRONT KERATINOCYTE TOO HARD
					double ffactor = 1.0;
					if (Cod!=NULL && Cod->get_cell_type() == 3)
					{
						ffactor = 0.2;
					}
					//double force = ffactor * max(0.0, min(GAR, (double)cellList[i]->get_cell_state(MIGCYTOKINEID)));
					//double force = ffactor * min(2.0*GAR, max(0.0, (maxcytokine - (double)cellList[i]->get_cell_state(MIGCYTOKINEID))));
					//cout << "    >>>> KMigForceFactor: " << biology->CalcValueMigForce(cellList[i]) << endl;
					double force = ffactor * biology->CalcValueMigForce(cellList[i]) * GAR * 5;
					cell_migrate_response_1(cellList[i], migangle, &kct, force);
					//cellList[i]->set_GR(0.6 * force);
					cellList[i]->set_GR(force);
					//cout << "   *-> K MIGRATION FORCE: " << force << endl;
					//// KERATINOCYTE DIVIDE AT DOUBLE SIZE
					ra = (int)migangle + 90;
					if (ra < 0) { ra = ra + 360; }
					else if (ra > 360) { ra = ra - 360; }
					if (cellList[i]->get_area() > AR)
					{
						cell_division(cellList[i], Time, cellList[i]->get_cell_type(), cellList[i]->get_cell_type(), ra);
						cellList[i]->center_refresh();
						cellList[cellList.size()-1]->center_refresh();
						for (int j=0; j<biology->getNumRealSpecies(); j++)
						{
							cellList[cellList.size()-1]->set_cell_state(j, (int)(cellList[i]->get_cell_state(j)/2.0));
							cellList[i]->set_cell_state(j, (int)(cellList[i]->get_cell_state(j)/2.0));
						}

					}
					//// END ////
					if (cellList[i]->get_bhsteps(4) > 5) // At least migrate for 5 steps.
					{
						cellList[i]->set_behavior(0);
						cellList[i]->set_bhsteps(0,0);
					}
					else
					{
						cellList[i]->set_bhsteps(4, cellList[i]->get_bhsteps(4)+1);
					}
				}
			}
			// OLD MIGRATION HAVE NO MECHANICAL SIGNALS
		}
	}
	cout << "  9. Doing Cell Migrations, Done." << endl;

	// Step 10: Do Nothing
	for (int i=0;i<cn;i++)
	{
		if (cellList[i]->Dead()) continue; // skip dead cells.

		if (cellList[i]->get_cell_type() == 4)
		{
			//if (cellList[i]->get_cell_state(2) > 100 && cellList[i]->get_cell_state(1) < 10)
			if (cellList[i]->get_cell_state(2) > 2 * cellList[i]->get_cell_state(1))
			{
				// If clot has more collagen than clot, then change cell type to ECM
				cellList[i]->set_cell_type(0);
			}
		}

		if (cellList[i]->get_cell_type () == 2 || 
			cellList[i]->get_cell_type () == 3 ||
			cellList[i]->get_cell_type () == 6) // Macrophage, Fibroblasts and Keratinocytes
		{
			if (cellList[i]->get_behavior() == 1)
			{
				// Cell Growth Rate zero
				cellList[i]->set_GR(0.0);
				if (cellList[i]->get_bhsteps(1) >= MIN_FBEHAVIOR_PERSISTENCE) // Stay dormant for 5 steps.
				{
					cellList[i]->set_behavior(0);
					cellList[i]->set_bhsteps(0,0);
				}
				else
				{
					cellList[i]->set_bhsteps(1, cellList[i]->get_bhsteps(1)+1);
				}
			}
		}
	}
	cout << "  10. Doing Nothing, Done." << endl;

	cout << "Biology -> Done." << endl;
	
	return -1;
}

#endif
