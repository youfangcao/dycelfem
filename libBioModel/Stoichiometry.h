// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#ifndef STOICHIOMETRY_H_
#define STOICHIOMETRY_H_

#include <vector>

using namespace std;

//class Model;
class State;
class SamplingState;
class ListOf;
class KineticLaw;

class Stoichiometry
{
	unsigned int  n_reactions; // Real number of reactions, reversible reactions counted as two. 
	unsigned int  n_species;
	string*       SpeciesID;
	double*       StoichiometricMatrix;
	vector<int*> ReversiblePairs;
	KineticLaw**  klaws;
	//ListOf*       ReactionModifier;
	unsigned int* ReactionModifier;
	int*          ReactionAttribute;	// if reaction is irreversible then its value is positive, else negative.
						// the absolute value is the the serial number of the reaction, 1 based.
	
public:
	Stoichiometry();
	~Stoichiometry();

	double getValue(unsigned int species_i, unsigned int reactions_j);
	unsigned int getModifierValue(unsigned int species_i, unsigned int reactions_j);
	unsigned int getNumSpecies ();
	int setNumSpecies(unsigned int n);
	int getSpeciesIndex(string Id);
	unsigned int getNumReactions ();
	int setNumReactions(unsigned int n);
	double* getStoichiometricMatrix();
	int setStoichiometricMatrix(double* matrix);
	int genStoichiometricMatrix(class Model* m);
	KineticLaw* getKineticLaw (unsigned int i);
	int getReactionAttribute(unsigned int i);

	int MapAllReversiblePairs();
	unsigned int getNumReversiblePairs();
	int getReversiblePairID(unsigned int i_pair, unsigned int idx);
	
	State* reactWithState(unsigned int i_reaction, State* astate);
	int    reactWithState(unsigned int i_reaction, State* astate, State* destState);
	//int    reactWithState(unsigned int i_reaction, State& astate, State& destState);
	int    updateStateWithReaction (unsigned int i_reaction, State* curstate);
	int    updateNewStateWithReaction (unsigned int i_reaction, State* before, State* after);
	int    canReact (unsigned int i_reaction, State* curstate);
	//int    canReact (unsigned int i_reaction, State& curstate);

	SamplingState* reactWithState(unsigned int i_reaction, SamplingState* astate);
	int    reactWithState(unsigned int i_reaction, SamplingState* astate, SamplingState* destState);
	int    updateStateWithReaction (unsigned int i_reaction, SamplingState* curstate);
	int    updateNewStateWithReaction (unsigned int i_reaction, SamplingState* before, SamplingState* after);
	int    canReact (unsigned int i_reaction, SamplingState* curstate);
	
	double* getStoichiometricVector (unsigned int i_reaction, int withModifiers);
	
	int printStoichiometricMatrix();

};





#endif /*STOICHIOMETRY_H_*/

