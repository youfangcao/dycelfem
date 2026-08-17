// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#ifndef STATE_H_
#define STATE_H_

/* CHOOSE ONE LENGTH SETUP */
/* State.h, StateSpace.h */
#include "StateType.h"
#include "fstream"

using namespace std;

class SamplingState;

class State
{
private:
	unsigned long sn; // serial number in state space, 0-based.
	unsigned int n_species;
	STATETYPE* state;
	State*        prev;
	State*        next;
	
public:
	State();
	State(State* astate);
	State(SamplingState* astate);
	State(unsigned int n_species);
	State(int* values, unsigned int n_species);
	~State();
	
	int CopyfromState(State* astate);
	int PrintState();
	int PrintState_nonl();
	int PrintState(ofstream &f);
	STATETYPE getStateValue(unsigned int i); // zero-based
	unsigned int  getNumSpecies();
	int setStateValue(unsigned int i, STATETYPE value);
	int setState(STATETYPE* astate, unsigned int n_species);
	State* getPrev();
	State* getNext();
	int setPrev(State* pt);
	int setNext(State* pt);
	int setSN (unsigned long index);
	unsigned long getSN ();
	int initialize ();
	
	State& operator= (const State &pt);
	bool operator== (const State &other);

	//void delete();
};



#endif /*STATE_H_*/
