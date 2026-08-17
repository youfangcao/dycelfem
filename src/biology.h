// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************
 ***    Project:  Cell Growth         ***
 ***    File:     biology.h           ***
 ***    Author:   Youfang Cao,        ***
 ***              Jieling Zhao        ***
 ***                                  ***
 ***    Created on May 20, 2013       ***
 ***************************************/

#ifndef BIOLOGY_H
#define BIOLOGY_H

using namespace std;

class Network;
class Stoichiometry;
class KineticLaw;
class cell;

class Biology
{
	// Networks first
	Network* nw_macrophage;
	Network* nw_fibroblast;
	Network* nw_keratinocyte;
	Network* nw_ecm;
	Network* nw_bm;
	Network* nw_clot;
	Network* nw_common;

	int nklaws_macrophage;
	int nklaws_fibroblast;
	int nklaws_keratinocyte;
	int nklaws_ecm;
	int nklaws_bm;
	int nklaws_clot;
	int nklaws_common;

	KineticLaw** kls_macrophage;
	KineticLaw** kls_fibroblast;
	KineticLaw** kls_keratinocyte;
	KineticLaw** kls_ecm;
	KineticLaw** kls_bm;
	KineticLaw** kls_clot;
	KineticLaw** kls_common;

	// Actual area of a cell, 30 micrometers diameter assumed.
	double AREA0; // cm^2/cell.
	// Basal cell growth rates
	double fg0; // fibroblast
	double kg0; // keratinocyte
	// Basal dimension for diffusion rate
	double D0;
	// Time Step for simulation
	double DeltaT;

	// Mutually Exclusive Cell behaviors
	// 0. Do nothing; 1. Growth/Division; 2. Migration; 3. Apoptosis.
	double Prob_DoNothing;
	double Prob_GrowthDiv;
	double Prob_Migration;
	double Prob_Apoptosis;

	// 0. fibrin clot; 1. collagen; 2. PDGF; 3. IL-1; 4. TGF-alpha; 5. EGF; 
	// 6. TGF-beta1; 7. TGF-beta3; 8. FGF2; 9. MMP; 10. TNF-alpha; 11. KGF(FGF7); 
	int N_SPECIES;
	
	int nrealSPECIES;

	/* 	cell type:
		0: ECM (virtual cell type)
		1: BM (basement membrane, virtual cell type)
		2: Fibroblast
		3: Keratinocyte
		4: Fibrin clot (virtual cell type)
		5: Hypodermis
		6: Macrophage
	*/
	int N_CELLTYPES;

		/***************************************
		cell type:
		0: ECM (virtual cell type, not real cell)
		1: BM (basement membrane)
		2: fibroblast
		3: keratinocyte
		4: fibrin clot
		5: Hypodermis
		6: Macrophage
		***************************************/


	// 0. fibrin clot; 1. collagen; 2. PDGF; 3. IL-1; 4. TGF-alpha; 5. EGF; 
	// 6. TGF-beta1; 7. TGF-beta3; 8. FGF2; 9. MMP; 10. TNF-alpha; 11. KGF(FGF7); 
	// Diffusion rates
	double* D;

	string* SpeciesNames;


public:
	Biology();
	~Biology();

	int getNumSpecies ();
	int getNumRealSpecies ();
	int getNumCellTypes ();
	double getDeltaT();
	double getD0();
	double getD(int k);
	string getSpeciesName(int i);

	// int UpdateCellState (Network* nw, KineticLaw** klaws, int nklaws, cell* pcell, double dT);
	// int UpdateCellState (int celltype, cell* pcell, double dT);
	int UpdateCellState (cell* pcell); // celltype: 2: nw_fibroblast; 3: nw_keratinocyte; 0: nw_ecm; 1: nw_bm; 4: nw_clot.
	int setValueProbGrowthDiv(double p);
	int setValueProbMigration(double p);
	int setValueProbApoptosis(double p);
	int setValueProbDoNothing(double p);
	double getValueProbGrowthDiv();
	double getValueProbMigration();
	double getValueProbApoptosis();
	double getValueProbDoNothing();
	double CalcValueProbGrowthDiv(cell* pcell);
	double CalcValueProbMigration(cell* pcell);
	double CalcValueProbApoptosis(cell* pcell);
	double CalcValueProbDoNothing(cell* pcell);
	double CalcValueMigForce(cell* pcell);
	

};

#endif

