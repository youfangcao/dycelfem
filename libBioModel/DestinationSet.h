// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#ifndef DESTINATIONSET_H_
#define DESTINATIONSET_H_

#include "StateType.h"
#include <vector>

using namespace std;

class State;
class SamplingState;
class SLinkedListS;

class DestinationSet
{
private:
	unsigned int  n_species;
	unsigned int  n_effspecies;  /* those species that are not -1 in destfile */
	long**        Destination;    /* 2D array of destination sets */
	int*          effspecies;
	
public:
	DestinationSet();
	~DestinationSet();
	
	int readDestinationSet (char* destfile);
	double getDistanceFromState (State* astate);
	double getDistanceBtwTwoStates (State* astate, State* bstate);

	double getDistanceFromState (SamplingState* astate);
	double getDistanceBtwTwoStates (SamplingState* astate, SamplingState* bstate);

	int getDistance2Edge (SamplingState* astate);

	double getDistanceFromPath (SLinkedListS* path);
	double getPathDistanceToDest(SamplingState* astate, SamplingState** path, unsigned int length);
	double getPathDistanceToBegin(SamplingState* X0, SamplingState* astate, SamplingState** path, unsigned int length);

	vector<int>* geteffspecies();
};

#endif /*DESTINATIONSET_H_*/
