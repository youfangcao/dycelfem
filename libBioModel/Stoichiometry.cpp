// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#include <iostream>
#include <string>
#include "Stoichiometry.h"
#include "State.h"
#include "SamplingState.h"
#include "sbml/SBMLTypes.h"
#include "debug.h"


using namespace std;

Stoichiometry::Stoichiometry()
{
	this->n_reactions = 0;
	this->n_species = 0;
	this->SpeciesID = 0;
	this->StoichiometricMatrix = 0;
	ReactionModifier = 0;
	this->ReactionAttribute = 0;
	this->klaws = 0;
	this->ReversiblePairs.clear();
}

Stoichiometry::~Stoichiometry()
{}

double Stoichiometry::getValue(unsigned int species_i, unsigned int reactions_j)
{
	return this->StoichiometricMatrix[species_i * n_reactions + reactions_j];
}

unsigned int Stoichiometry::getModifierValue(unsigned int species_i, unsigned int reactions_j)
{
	return ReactionModifier[(species_i + 1) * n_reactions + reactions_j];
}

unsigned int Stoichiometry::getNumSpecies ()
{
	return this->n_species;
}

int Stoichiometry::setNumSpecies (unsigned int n)
{
	this->n_species = n;
	return 1;
}

int Stoichiometry::getSpeciesIndex(string Id)
{
	int index = -1;
	unsigned int i = 0;
	
	for (i = 0; i < n_species; i++)
	{
		if (SpeciesID[i] == Id)
		{
			index = i;
			break;
		}
	}
	
	return index;
}

unsigned int Stoichiometry::getNumReactions ()
{
	return this->n_reactions;
}

int Stoichiometry::setNumReactions (unsigned int n)
{
	this->n_reactions = n;
	return 1;
}


double* Stoichiometry::getStoichiometricMatrix()
{
	return this->StoichiometricMatrix;
}


int Stoichiometry::setStoichiometricMatrix(double* matrix)
{
	// set a stoichiometric matrix given species number and reaction number.
	unsigned int scale = this->n_species * this->n_reactions;
	this->StoichiometricMatrix = new double[scale];
	for (unsigned int i = 0; i < scale; i++)
	{
		this->StoichiometricMatrix[i] = matrix[i];
	}
	return 1;
}


int Stoichiometry::getReactionAttribute(unsigned int i)
{
	return this->ReactionAttribute[i];
}


int Stoichiometry::MapAllReversiblePairs()
{
	double* reactionvec0 = 0;
	double* reactionvec1 = 0;
	int* pair = 0;
	int     ispair = 1;
	int     paircount = 0;
	unsigned int* visited = new unsigned int[this->n_reactions];

	if (!this->StoichiometricMatrix)
	{
		cout << "WRONG call sequence." << endl;
		cout << "NULL stoichiometry matrix." << endl;
		return -1;
	}

	for (unsigned int i = 0; i < this->n_reactions; i++)
	{
		visited[i] = 0;
	}

	for (unsigned int i = 0; i < this->n_reactions; i++)
	{
		if (visited[i] == 1) { continue; }

		paircount = 0;
		reactionvec0 = getStoichiometricVector(i, 0);
		for (unsigned int j = 0; j < this->n_species; j ++)
		{
			reactionvec0[j] = -1 * reactionvec0[j];
		}
		for (unsigned int k = i+1; k < this->n_reactions; k ++)
		{
			reactionvec1 = getStoichiometricVector(k, 0);
			ispair = 1;
			for (unsigned int j = 0; j < this->n_species; j ++)
			{
				if (reactionvec0[j] != reactionvec1[j]) ispair = 0;
			}
			if (ispair == 1 && visited[i] == 0 && visited[k] == 0)
			{
				pair = new int[2];
				pair[0] = (int)i; pair[1] = (int)k;
				ReversiblePairs.push_back(pair);
				visited[i] = 1;
				visited[k] = 1;
				paircount ++;
			}
			//delete reactionvec1;
		}
		if (paircount == 0)
		{
			pair = new int[2];
			pair[0] = (int)i; pair[1] = -1;
			ReversiblePairs.push_back(pair);
			visited[i] = 1;
		}
		if (paircount > 1) { cout << "I don't feel good: triple reversible for reaction " << i << "?" << endl; }
		//delete reactionvec0;
	}

	for (unsigned int i = 0; i < ReversiblePairs.size(); i ++)
	{
		cout << ReversiblePairs[i][0] << ":" << ReversiblePairs[i][1] << endl;
	}

	//delete visited;
	return 0;
}

unsigned int Stoichiometry::getNumReversiblePairs()
{
	return this->ReversiblePairs.size();
}

int Stoichiometry::getReversiblePairID(unsigned int i_pair, unsigned int idx)
{
	return this->ReversiblePairs[i_pair][idx];
}


/*
 * Generate Stoichiometrix Matrix
 */
int Stoichiometry::genStoichiometricMatrix(Model* m)
{
	Reaction*          cur_reaction = 0;
	SpeciesReference*  cur_refspecies = 0;
	unsigned int       i=0, j=0, k=0, n=0;
	unsigned int       n_reactants = 0;
	unsigned int       n_products = 0;
	
	/*
	 * build stoichiometic matrix according to reactions in the model.
	 */
	/* Count in reversible reactions. */
	this->n_species = m->getNumSpecies();
	this->n_reactions = m->getNumReactions();
	
	int num_reversibleReaction = 0;
	
	for (i = 0; i < this->n_reactions; i++)
	{
		if (m->getReaction(i)->getReversible())
			num_reversibleReaction ++;
	}
	this->n_reactions += num_reversibleReaction;

	//printf ("%d, %d\n", this->n_species, this->n_reactions);

	// for debug
	//printf ("     >>>>>>>>>sss444.\n");
	
	/*
	 * get a modifier list vector for each reaction.
	 * store all modifier information and modifier number in this 2D array.
	 * The first line is the modifier number of each reaction.
	 */
	ReactionModifier = new unsigned int[(n_species + 1) * n_reactions];
	for (i = 0; i < n_species + 1; i++)
	{
		for (j = 0; j < n_reactions; j++)
		{
			ReactionModifier[i * n_reactions + j] = 0;
		}
	}
	
	/*
	 * get SpeciesID list from the model.
	 */
	this->SpeciesID = new string[this->n_species];
	for (i = 0; i < this->n_species; i++)
	{
		this->SpeciesID[i] = m->getSpecies(i)->getId();
		//printf ("%s\n", SpeciesID[i].data());
	}
	
	/*
	 * Initialize klaws
	 */
	this->klaws = new KineticLaw*[this->n_reactions];

	// for debug
	//printf ("     >>>>>>>>>sss555.\n");
	
	/* 
	 * allocate memory space for stoichiometric matrix. 
	 */
	this->StoichiometricMatrix = new double[this->n_species * this->n_reactions];
	this->ReactionAttribute = new int[this->n_reactions];

	unsigned int cur_StoichiometryIndex = 0;
	
	unsigned int totalReactionsInModel = m->getNumReactions(); /* the number of reactions without reversible reactions. */
	for (j = 0; j < totalReactionsInModel; j++)
	{
		/*
		 * initialize the Stoichiometry
		 */
		for(i = 0; i < this->n_species; i++)
		{
			this->StoichiometricMatrix[i * this->n_reactions + cur_StoichiometryIndex] = 0;
		}
		
		/*
		 * acquire current reaction from the model
		 */
		cur_reaction = m->getReaction(j);

		/*
		 * Deal with irreversible reaction.
		 */

		/*
		 * get modifier speciesID from current reaction.
		 */
		//ListOf RM = cur_reaction->getListOfModifiers();
		ReactionModifier[cur_StoichiometryIndex] = cur_reaction->getNumModifiers();
		if (ReactionModifier[cur_StoichiometryIndex] > 0)
		{
			for (k = 0; k < ReactionModifier[cur_StoichiometryIndex]; k++)
			{
				ReactionModifier[ (getSpeciesIndex(cur_reaction->getModifier(k)->getSpecies()) + 1) * n_reactions + cur_StoichiometryIndex] += 1;
			}
		}
		// for test
		//printf ("Reaction %d has %d modifiers.\n", j+1, cur_reaction->getNumModifiers());
		//unsigned int nmodi = cur_reaction->getNumModifiers();
		//for (k = 0; k < nmodi; k++)
		//{
		//	printf ("   Modifier %d is %s.\n", k, cur_reaction->getModifier(k)->getSpecies().c_str());
		//}
		//
		
		/* 
		 * dealing with reactants first
		 */
		n_reactants = cur_reaction->getNumReactants();
		//printf ("%d\n", n_reactants);
		for(k = 0; k < n_reactants; k++)
		{
			cur_refspecies = cur_reaction->getReactant(k);
			string id = cur_refspecies->getSpecies();
			
			/* find out the index number of this ID in the SpeciesID array. */
			int index = -1;
			for (n = 0; n < this->n_species; n++)
			{
				if (this->SpeciesID[n] == id)
				{
					index = n;
					break;
				}
			}
			if (index == -1)
			{
				printf ("Error: Unknown %dth reactant of reaction: %s\n", k+1, cur_reaction->getId().data());
				return -1;
			}
			/* OK, here we have got the species index, 
			 * we can use it in building the SMatrix.
			 */
			this->StoichiometricMatrix[index * this->n_reactions + cur_StoichiometryIndex] += -1 * cur_refspecies->getStoichiometry(); /* times by -1, because they are reactants. */
			//printf ("%f\t%f\n", cur_refspecies->getStoichiometry(), StoichiometricMatrix[index * n_reactions + cur_StoichiometryIndex]);
		}
		
		/* 
		 * then, dealing with products now.
		 */
		n_products = cur_reaction->getNumProducts();
		for(k = 0; k < n_products; k++)
		{
			cur_refspecies = cur_reaction->getProduct(k);
			string id = cur_refspecies->getSpecies();
			
			/* find out the index number of this ID in the SpeciesID array. */
			int index = -1;
			for (n = 0; n < this->n_species; n++)
			{
				if(this->SpeciesID[n] == id)
				{
					index = n;
					break;
				}
			}
			if (index == -1)
			{
				printf ("Error: Unknown %dth product of reaction: %s\n", k+1, cur_reaction->getId().data());
				return -1;
			}

			/* OK, here we have got the index, 
			 * we can use this index to set the value of StoichiometricMatrix.
			 */
			this->StoichiometricMatrix[index * this->n_reactions + cur_StoichiometryIndex] += cur_refspecies->getStoichiometry(); /* they are products. */
		}
		
		/*
		 * Record the attribute of current reaction to discriminate reversible and irreversible reactions.
		 */
		this->ReactionAttribute[cur_StoichiometryIndex] = j + 1;
		this->klaws[cur_StoichiometryIndex] = cur_reaction->getKineticLaw();

		cur_StoichiometryIndex ++;


		//printf ("   <><><>here<><><>%s\n", cur_reaction->getId().data());
		
		/*********************************************\
		 * If this reaction is reversible, then deal with the reverse reaction -- add the reaction to Stoichiometry below.
		\*********************************************/
		if (cur_reaction->getReversible())
		{
			/*
			 * initialize the Stoichiometry
			 */
			for(i = 0; i < this->n_species; i++)
			{
				this->StoichiometricMatrix[i * this->n_reactions + cur_StoichiometryIndex] = 0;
			}
			
			/*
			 * get modifier speciesID from current reaction.
			 */
			//ReactionModifier[cur_StoichiometryIndex] = cur_reaction->getListOfModifiers();
			ReactionModifier[cur_StoichiometryIndex] = cur_reaction->getNumModifiers();
			if (ReactionModifier[cur_StoichiometryIndex] > 0)
			{
				for (k = 0; k < ReactionModifier[cur_StoichiometryIndex]; k++)
				{
					ReactionModifier[ (getSpeciesIndex(cur_reaction->getModifier(k)->getSpecies()) + 1) * n_reactions + cur_StoichiometryIndex] += 1;
				}
			}
		
			/* 
			 * dealing with reactants first
			 * ACTUALLY, THEY ARE PRODUCTS HERE.
			 */
			n_reactants = cur_reaction->getNumReactants();
			for(k = 0; k < n_reactants; k++)
			{
				cur_refspecies = cur_reaction->getReactant(k);
				string id = cur_refspecies->getSpecies();
				
				/* find out the index number of this ID in the SpeciesID array. */
				int index = -1;
				for (n = 0; n < this->n_species; n++)
				{
					if(this->SpeciesID[n] == id)
					{
						index = n;
						break;
					}
				}
				if (index == -1)
				{
					printf ("Error: Unknown %dth reactant of reaction: %s\n", k+1, cur_reaction->getId().data());
					return -1;
				}
				/* OK, here we have got the species index, 
				 * we can use it in building the SMatrix.
				 */
				this->StoichiometricMatrix[index * this->n_reactions + cur_StoichiometryIndex] += /*-1 * */ cur_refspecies->getStoichiometry(); /* they are products now. */
			}
			
			/* 
			 * then, dealing with products now.
			 * ACTUALLY, THEY ARE REACTANTS HERE.
			 */
			n_products  = cur_reaction->getNumProducts();
			for(k = 0; k < n_products; k++)
			{
				cur_refspecies = cur_reaction->getProduct(k);
				string id = cur_refspecies->getSpecies();
				
				/* find out the index number of this ID in the SpeciesID array. */
				int index = -1;
				for (n = 0; n < this->n_species; n++)
				{
					if(this->SpeciesID[n] == id)
					{
						index = n;
						break;
					}
				}
				if (index == -1)
				{
					printf ("Error: Unknown %dth product of reaction: %s\n", k+1, cur_reaction->getId().data());
					return -1;
				}
				/* OK, here we have got the index, 
				 * we can use it in building the SMatrix.
				 */
				this->StoichiometricMatrix[index * this->n_reactions + cur_StoichiometryIndex] += -1 * cur_refspecies->getStoichiometry(); /* times by -1, they are reactants now. */
			}
			
			/*
		 	* Remember the attribute of current reaction to discriminate reversible and irreversible reactions.
		 	*/
			this->ReactionAttribute[cur_StoichiometryIndex] = -1 * (j + 1);
			this->klaws[cur_StoichiometryIndex] = cur_reaction->getKineticLaw();

			cur_StoichiometryIndex ++;
		}
		

	}


	/*
	 * print out SMatrix for debug
	 * 
	 * comment out if below code block is not necessary any more.
	 * 
	 */
	/*
	if (clpara->verbose)
	{
		printf ("Stoichiometric matrix:\n");
		for (i = 0; i < n_species; i++)
		{
			for (j = 0; j < n_reactions; j++)
			{
				printf ("%d  ", (int)StoichiometricMatrix[i * n_reactions + j]);
			}
			printf ("\n");
		}
	}*/
	
	return 1;
}
////

//int Stoichiometry::reactWithState(unsigned int i_reaction, State& astate, State& destState)
//{
//	int value = 0;
//	unsigned int i = 0, j = 0;
//	cout << "INSIDE REACT" << endl;
//	astate.PrintState();
//	destState.PrintState();
//	
//	// initialize destState
//	for (i = 0; i < this->n_species; i++)
//	{
//			 astate.PrintState();
//		destState.setStateValue(i, 0);
//	}
//			 astate.PrintState();
//	cout << "INSIDE REACT 1" << endl;
//	
//	for (i = 0; i < this->n_species; i++)
//	{
//			 astate.PrintState();
//		if (this->StoichiometricMatrix[i * this->n_reactions + i_reaction] < 0)
//		{
//			 astate.PrintState();
//			 if (astate.getStateValue(i) == 0 || (int)(astate.getStateValue(i)) < abs((int)(this->StoichiometricMatrix[i * this->n_reactions + i_reaction])))
//			 {
//			 astate.PrintState();
//	cout << "WHAT HAPPENED HERE." << endl;
//				return 0;
//			 }
//		}		
//	}
//	cout << "INSIDE REACT E" << endl;
//	
//	// check Modifier list vector.
//	if (ReactionModifier[i_reaction] != 0)
//	{
//		for (j = 0; j < n_species; j++)
//		{
//			if (astate.getStateValue(j) < ReactionModifier[(j + 1) * n_reactions + i_reaction])
//			{
//				return 0;
//			}
//		}
//	}
//	cout << "INSIDE REACT EN" << endl;
//	
//	for (i = 0; i < this->n_species; i++)
//	{
//		value = (int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction] + astate.getStateValue(i);
//		destState.setStateValue(i, value);
//	}
//	destState.PrintState();
//	cout << "INSIDE REACT END" << endl;
//	
//	return 1;
//}


int Stoichiometry::reactWithState(unsigned int i_reaction, State* astate, State* destState)
{
	int value = 0;
	unsigned int i = 0, j = 0;
	
	// initialize destState
	for (i = 0; i < this->n_species; i++)
	{
		destState->setStateValue(i, 0);
	}
	
	for (i = 0; i < this->n_species; i++)
	{
		if (this->StoichiometricMatrix[i * this->n_reactions + i_reaction] < 0)
		{
			 if (astate->getStateValue(i) == 0 || (int)(astate->getStateValue(i)) < abs((int)(this->StoichiometricMatrix[i * this->n_reactions + i_reaction])))
			 {
				return 0;
			 }
		}		
	}
	
	// check Modifier list vector.
	if (ReactionModifier[i_reaction] != 0)
	{
		for (j = 0; j < n_species; j++)
		{
			if (astate->getStateValue(j) < ReactionModifier[(j + 1) * n_reactions + i_reaction])
			{
				return 0;
			}
		}
	}
	
	for (i = 0; i < this->n_species; i++)
	{
		value = (int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction] + astate->getStateValue(i);
		if (destState != NULL)
		{
			destState->setStateValue(i, value);
		}
		else
		{
			printf ("Variable(State) not initialized.\n");
			return -1;
		}
	}
	
	return 1;
}

State* Stoichiometry::reactWithState(unsigned int i_reaction, State* astate)
{
	int value = 0;
	unsigned int i = 0, j = 0;
	
	for (i = 0; i < this->n_species; i++)
	{
		if (this->StoichiometricMatrix[i * this->n_reactions + i_reaction] < 0)
		{
			 if (astate->getStateValue(i) == 0 || (int)(astate->getStateValue(i)) < abs((int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction]))
			 {
				return NULL;
			 }
		}		
	}
	
	// check Modifier list vector.
	if (ReactionModifier[i_reaction] != 0)
	{
		for (j = 0; j < n_species; j++)
		{
			if (astate->getStateValue(j) < ReactionModifier[(j + 1) * n_reactions + i_reaction])
			{
				return 0;
			}
		}
	}
	
	
	State* s = new State(this->n_species);
	for (i = 0; i < this->n_species; i++)
	{
		value = (int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction] + astate->getStateValue(i);
		
		s->setStateValue(i, value);
	}
	
	return s;
}

int Stoichiometry::reactWithState(unsigned int i_reaction, SamplingState* astate, SamplingState* destState)
{
	int value = 0;
	unsigned int i = 0, j = 0;
	
	// initialize destState
	for (i = 0; i < this->n_species; i++)
	{
		destState->setStateValue(i, 0);
	}
	
	for (i = 0; i < this->n_species; i++)
	{
		if (this->StoichiometricMatrix[i * this->n_reactions + i_reaction] < 0)
		{
			 if (astate->getStateValue(i) == 0 || (int)(astate->getStateValue(i)) < abs((int)(this->StoichiometricMatrix[i * this->n_reactions + i_reaction])))
			 {
				return 0;
			 }
		}		
	}
	
	// check Modifier list vector.
	if (ReactionModifier[i_reaction] != 0)
	{
		for (j = 0; j < n_species; j++)
		{
			if (astate->getStateValue(j) < ReactionModifier[(j + 1) * n_reactions + i_reaction])
			{
				return 0;
			}
		}
	}
	
	for (i = 0; i < this->n_species; i++)
	{
		value = (int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction] + astate->getStateValue(i);
		if (destState != NULL)
		{
			destState->setStateValue(i, value);
		}
		else
		{
			printf ("Variable(State) not initialized.\n");
			return -1;
		}
	}
	
	return 1;
}

SamplingState* Stoichiometry::reactWithState(unsigned int i_reaction, SamplingState* astate)
{
	int value = 0;
	unsigned int i = 0, j = 0;
	
	for (i = 0; i < this->n_species; i++)
	{
		if (this->StoichiometricMatrix[i * this->n_reactions + i_reaction] < 0)
		{
			 if (astate->getStateValue(i) == 0 || (int)(astate->getStateValue(i)) < abs((int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction]))
			 {
				return NULL;
			 }
		}		
	}
	
	// check Modifier list vector.
	if (ReactionModifier[i_reaction] != 0)
	{
		for (j = 0; j < n_species; j++)
		{
			if (astate->getStateValue(j) < ReactionModifier[(j + 1) * n_reactions + i_reaction])
			{
				return 0;
			}
		}
	}
	
	
	SamplingState* s = new SamplingState(this->n_species, astate->getNumStoiReactions());
	for (i = 0; i < this->n_species; i++)
	{
		value = (int)this->StoichiometricMatrix[i * this->n_reactions + i_reaction] + astate->getStateValue(i);
		
		s->setStateValue(i, value);
	}
	
	return s;
}

int Stoichiometry::updateStateWithReaction (unsigned int i_reaction, State* curstate)
{
	unsigned int i = 0;
	int before = 0;

	for (i=0; i<this->n_species; i++)
	{
		before = curstate->getStateValue(i);
		curstate->setStateValue(i, before+(int)(this->StoichiometricMatrix[i*this->n_reactions+i_reaction]));
	}

	return 1;
}

int Stoichiometry::updateNewStateWithReaction (unsigned int i_reaction, State* before, State* after)
{
	unsigned int i = 0;

	for (i=0; i<this->n_species; i++)
	{
		after->setStateValue(i, before->getStateValue(i)+(int)(this->StoichiometricMatrix[i*this->n_reactions+i_reaction]));
	}

	return 1;
}

int Stoichiometry::canReact (unsigned int i_reaction, State* curstate)
{
	unsigned int i = 0;
	vector<double> s;

	for (i=0; i<this->n_species; i++)
	{
		s.push_back (this->StoichiometricMatrix[i*this->n_reactions+i_reaction]);
		s[i] *= -1;
		if (s[i] < 0) s[i] = 0;
		s[i] += (double)ReactionModifier[(i + 1) * n_reactions + i_reaction];
	}

	for (i=0; i<this->n_species; i++)
	{
		if (curstate->getStateValue(i) < s[i])
		{
			return 0;
		}
	}


	return 1;
}

int Stoichiometry::updateStateWithReaction (unsigned int i_reaction, SamplingState* curstate)
{
	unsigned int i = 0;
	int before = 0;

	for (i=0; i<this->n_species; i++)
	{
		before = curstate->getStateValue(i);
		curstate->setStateValue(i, before+(int)(this->StoichiometricMatrix[i*this->n_reactions+i_reaction]));
	}

	return 1;
}

int Stoichiometry::updateNewStateWithReaction (unsigned int i_reaction, SamplingState* before, SamplingState* after)
{
	unsigned int i = 0;

	for (i=0; i<this->n_species; i++)
	{
		after->setStateValue(i, before->getStateValue(i)+(int)(this->StoichiometricMatrix[i*this->n_reactions+i_reaction]));
	}

	return 1;
}

int Stoichiometry::canReact (unsigned int i_reaction, SamplingState* curstate)
{
	unsigned int i = 0;
	vector<double> s;

	for (i=0; i<this->n_species; i++)
	{
		s.push_back (this->StoichiometricMatrix[i*this->n_reactions+i_reaction]);
		s[i] *= -1;
		if (s[i] < 0) s[i] = 0;
		s[i] += (double)ReactionModifier[(i + 1) * n_reactions + i_reaction];
	}

	for (i=0; i<this->n_species; i++)
	{
		if (curstate->getStateValue(i) < s[i])
		{
			return 0;
		}
	}


	return 1;
}

//int Stoichiometry::canReact (unsigned int i_reaction, State& curstate)
//{
//	unsigned int i = 0;
//	vector<double> s;
//	cout << "canReact" << endl;
//	curstate.PrintState();
//
//	for (i=0; i<this->n_species; i++)
//	{
//		s.push_back (this->StoichiometricMatrix[i*this->n_reactions+i_reaction]);
//		s[i] *= -1;
//		if (s[i] < 0) s[i] = 0;
//		s[i] += (double)ReactionModifier[(i + 1) * n_reactions + i_reaction];
//	}
//
//	for (i=0; i<this->n_species; i++)
//	{
//		if (curstate.getStateValue(i) < s[i])
//		{
//			return 0;
//		}
//	}
//
//
//	return 1;
//}

int Stoichiometry::printStoichiometricMatrix()
{
	printf ("Stoichiometric matrix:\n");
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		for (unsigned int j = 0; j < this->n_reactions; j++)
		{
			printf ("% d  ", (int)this->StoichiometricMatrix[i * this->n_reactions + j]);
		}
		printf ("\n");
	}
	printf ("Reaction Modifiers:\n");
	for (unsigned int i = 0; i < n_species; i++)
	{
		if (i == 0)
		{
			printf ("Number of modifiers in each reaction:\n");
			for (unsigned int j = 0; j < this->n_reactions; j++)
			{
				printf ("%d ", ReactionModifier[i * this->n_reactions + j]);
			}
			printf ("\n\n");
		}
		
		for (unsigned int j = 0; j < this->n_reactions; j++)
		{
			printf ("%d ", this->ReactionModifier[(i+1) * this->n_reactions + j]);
		}
		printf ("\n");
	}
	//printf ("\n");
	
	return 0;
}

double* Stoichiometry::getStoichiometricVector (unsigned int i_reaction, int withModifiers)
{
	double* stoichiometry = new double[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		stoichiometry[i] = 0;
	}

	for (unsigned int i = 0; i < this->n_species; i++)
	{
		stoichiometry[i] = this->StoichiometricMatrix[i * this->n_reactions + i_reaction];
		if (withModifiers > 0)
		{
			 stoichiometry[i] -= ReactionModifier[(i + 1) * this->n_reactions + i_reaction];
		}
	}
	

	return stoichiometry;
}

KineticLaw* Stoichiometry::getKineticLaw (unsigned int i)
{
	return klaws[i];
}

/* EOF */


