// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#ifndef NETWORK_H_
#define NETWORK_H_

/* CHOOSE ONE LENGTH SETUP */
/* global.h */
#include "StateType.h"

#include <string>

using namespace std;



class SparseMatrix;
class MatrixElement;
class CLParameters;
class KineticLaw;
class ASTNode;
class SBMLDocument;
class Model;
class Stoichiometry;
class State;
class SamplingState;
class DestinationSet;


class Network
{
	SBMLDocument*        d;
	Model*               m;
	int                  n_errors;
	STATETYPE*	     States;
	Stoichiometry*       stoichiometry;
	SparseMatrix*        TransMatrix;
	MatrixElement*       matrixEle;
	string*		     CompartmentsID;
	string*              SpeciesID;
	string*              RulesID;
	string*		     ParametersID;
	string*              ReactionsID;
	string* 	     LocalParamIDs;
	unsigned long        n_states;
	unsigned int         n_compartments;
	unsigned int         n_species;
	unsigned int         n_reactions;
	unsigned int         n_stoireactions;
	unsigned int         n_rules;
	unsigned int	     n_parameters;
	unsigned int         n_realStates;
	int*                 ReactionAttribute;
	unsigned int         n_assignmentRules;

	double* rates_a;
	double* rates_b;

public:
	// Constructors here.
	Network();

	int initializeNetwork(char* sbmlfile);

	int getNumErrors();

	int printErrors();

	Model* getModel();

	/*
	 * Define a useful function here
	 */
	int getSpeciesIndex(string Id);

	int getCompartmentIndex(string Id);

	int getRateRuleIndex(string Id);

	string getSpeciesID(unsigned int sindex);

	string getReactionID(unsigned int sindex);

	string getRuleID(unsigned int rindex);

	string getParameterID(unsigned int pindex);

	//int validateTransMatrix(MatrixElement* TMatrix, unsigned int n_states);

	//int writeStateSpace();

	//int writeTransMatrix();

	/*
	 * reaction index: sindex - stoichiometry matrix index
	 */
	double getReactionRate(const KineticLaw* klaw, unsigned int sindex, State* curstate);
	double getReactionRate(const KineticLaw* klaw, unsigned int sindex, SamplingState* curstate);

	double computeASTNodeValue(const ASTNode* node, unsigned int sindex, State* curstate);

	double getValueOfParameter(const char* parameter_id, unsigned int rindex, State* curstate);

	int    isParameter(const char* parameter_id, unsigned int rindex);

	double getValueOfRateConstant(unsigned int rindex);

	void initializeRateArrays();

	SparseMatrix* getTransitionMatrix();

	void setTransitionMatrix(SparseMatrix* tmatrix);

	Stoichiometry* getStoichiometricMatrix();

	STATETYPE* getStateSpace();

	unsigned int getNumCompartments();

	unsigned int getNumReactions();

	unsigned int getNumStoiReactions();

	unsigned int getNumSpecies();

	unsigned int getNumRules();

	unsigned int getNumParameters();

	unsigned int getNumRealStates();

	unsigned int getNumStates();

	void setNumRealStates(unsigned long num);

	void setStateSpace(STATETYPE* states);

	void setNumStates(unsigned long n_states);
	
};

#endif /*NETWORK_H_*/


