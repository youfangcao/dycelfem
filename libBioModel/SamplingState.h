// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#ifndef SAMPLINGSTATE_H_
#define SAMPLINGSTATE_H_

/* CHOOSE ONE LENGTH SETUP */
/* State.h, StateSpace.h */
#include "StateType.h"
#include "fstream"

using namespace std;

class State;

class SamplingState
{
private:
	unsigned long sn; // serial number in state space, 0-based.
	unsigned int n_species;
	unsigned int n_stoireactions;
	STATETYPE* state;
	double* rates_a;
	double* rates_b;
	double* gamma;
	double  a0;
	double  b0;
	double* pfn;
	double* pbn;
	unsigned long nVisits;
	SamplingState* prev;
	SamplingState* next;

public:
	SamplingState();
	SamplingState(SamplingState* astate);
	SamplingState(SamplingState* astate, unsigned int n_stoireactions);
	SamplingState(State* astate, unsigned int n_stoireactions);
	//SamplingState(unsigned int n_species);
	SamplingState(unsigned int n_species, unsigned int n_stoireactions);
	SamplingState(int* values, unsigned int n_species, unsigned int n_stoireactions);
	~SamplingState();
	
	int initialize (unsigned int n_stoireactions);
	int CopyfromState(SamplingState* astate);
	//int CopyfromState(State* astate);
	int CopyfromState(State* astate, unsigned int n_stoireactions);
	int PrintState();
	int PrintState_Full();
	int PrintState_Full(ofstream& fout);
	int PrintState_nonl();
	int PrintState(ofstream &fout);
	STATETYPE getStateValue(unsigned int i); // zero-based
	unsigned int  getNumSpecies();
	int setStateValue(unsigned int i, STATETYPE value);
	int setState(STATETYPE* astate, unsigned int n_species);
	SamplingState* getPrev();
	SamplingState* getNext();
	int setPrev(SamplingState* pt);
	int setNext(SamplingState* pt);
	int setSN (unsigned long index);
	unsigned long getSN ();
	int setRatesA (int j, double v);
	int setRatesB (int j, double v);
	int setGamma  (int j, double v);
	double getRatesA (int j);
	double getRatesB (int j);
	double getGamma  (int j);
	double geta0 ();
	int seta0 (double v);
	int calca0 ();
	double getb0 ();
	int setb0 (double v);
	int setpFn(unsigned int i, double v);
	int setpBn(unsigned int i, double v);
	double getpFn(unsigned int i);
	double getpBn(unsigned int i);
	unsigned long getnVisits ();
	int setnVisits (unsigned long v);
	int increase1Visit ();
	double calcB0 ();
	unsigned int getNumStoiReactions();
	int setNumStoiReactions(unsigned int n_stoireactions);
	
	SamplingState& operator= (const SamplingState &pt);
	bool operator== (const SamplingState &other);

	//void delete();
	
	// Special members for 2Dwalk boundary problem.
	//int moveStateOnB (SamplingState* p, SamplingState* p1, int distance);
};



#endif /*SAMPLINGSTATE_H_*/


