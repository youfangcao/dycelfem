// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois

#include <iostream>
#include "SamplingState.h"
#include "State.h"
#include "debug.h"

using namespace std;


SamplingState::SamplingState()
{
	this->sn = 0;
	this->n_species = 0;
	this->n_stoireactions = 0;
	this->state = 0;
	this->prev = 0;
	this->next = 0;
	this->rates_a = 0;
	this->rates_b = 0;
	this->pfn = 0;
	this->pbn = 0;
	this->gamma = 0;
	this->nVisits = 0;
	a0 = 0;
	b0 = 0;
}

SamplingState::SamplingState(SamplingState* astate)
{
	this->sn = astate->sn;
	this->n_species = astate->n_species;
	this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate->state[i];
	}
	this->prev = astate->prev;
	this->next = astate->next;
	this->n_stoireactions = astate->n_stoireactions;
	this->rates_a = new double[this->n_stoireactions];
	this->rates_b = new double[this->n_stoireactions];
	this->gamma   = new double[this->n_stoireactions];
	this->nVisits = astate->nVisits;
	if (astate->pfn)
	{
		this->pfn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pfn[i] = astate->pfn[i];
		}
	}
	else
	{
		this->pfn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pfn[i] = 0;
		}
	}
	if (astate->pbn)
	{
		this->pbn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pbn[i] = astate->pbn[i];
		}
	}
	else
	{
		this->pbn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pbn[i] = 0;
		}
	}

	if ((astate->rates_a) && (astate->rates_b) && (astate->gamma))
	{
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			rates_a[i] = astate->rates_a[i];
			rates_b[i] = astate->rates_b[i];
			gamma[i] = astate->gamma[i];
		}
		a0 = astate->a0;
		b0 = astate->b0;
	}
	else if (!(astate->rates_a) && !(astate->rates_b) && !(astate->gamma))
	{
		this->rates_a = new double[this->n_stoireactions];
		this->rates_b = new double[this->n_stoireactions];
		this->gamma   = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			rates_a[i] = 0;
			rates_b[i] = 0;
			gamma[i] = 0;
		}
		a0 = 0;
		b0 = 0;
	}
}

SamplingState::SamplingState(SamplingState* astate, unsigned int n_stoireactions)
{
	this->sn = astate->sn;
	this->n_species = astate->n_species;
	this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate->state[i];
	}
	this->prev = astate->prev;
	this->next = astate->next;
	this->n_stoireactions = n_stoireactions;
	this->nVisits = astate->nVisits;
	if (astate->pfn)
	{
		this->pfn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pfn[i] = astate->pfn[i];
		}
	}
	else
	{
		this->pfn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pfn[i] = 0;
		}
	}
	if (astate->pbn)
	{
		this->pbn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pbn[i] = astate->pbn[i];
		}
	}
	else
	{
		this->pbn = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			pbn[i] = 0;
		}
	}

	if (!rates_a && !rates_b && !gamma && (astate->rates_a) && (astate->rates_b) && (astate->gamma))
	{
		this->rates_a = new double[this->n_stoireactions];
		this->rates_b = new double[this->n_stoireactions];
		this->gamma   = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			rates_a[i] = astate->rates_a[i];
			rates_b[i] = astate->rates_b[i];
			gamma[i] = astate->gamma[i];
		}
		a0 = astate->a0;
		b0 = astate->b0;
	}
	else if (rates_a && rates_b && gamma && (astate->rates_a) && (astate->rates_b) && (astate->gamma))
	{
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			rates_a[i] = astate->rates_a[i];
			rates_b[i] = astate->rates_b[i];
			gamma[i] = astate->gamma[i];
		}
		a0 = astate->a0;
		b0 = astate->b0;
	}
	else if (!(astate->rates_a) && !(astate->rates_b) && !(astate->gamma))
	{
		this->rates_a = new double[this->n_stoireactions];
		this->rates_b = new double[this->n_stoireactions];
		this->gamma   = new double[this->n_stoireactions];
		for (unsigned int i = 0; i < this->n_stoireactions; i++)
		{
			rates_a[i] = 0;
			rates_b[i] = 0;
			gamma[i] = 0;
		}
		a0 = 0;
		b0 = 0;
	}
}

SamplingState::SamplingState(State* astate, unsigned int n_stoireactions)
{
	this->sn = astate->getSN();
	this->n_species = astate->getNumSpecies();
	this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate->getStateValue(i);
	}
	this->prev = 0;
	this->next = 0;
	this->n_stoireactions = n_stoireactions;
	this->rates_a = new double[this->n_stoireactions];
	this->rates_b = new double[this->n_stoireactions];
	this->gamma   = new double[this->n_stoireactions];
	this->pfn = new double[this->n_stoireactions];
	this->pbn = new double[this->n_stoireactions];
	this->nVisits = 0;
	for (unsigned int i = 0; i < this->n_stoireactions; i++)
	{
		rates_a[i] = 0;
		rates_b[i] = 0;
		gamma[i] = 0;
		pfn[i] = 0;
		pbn[i] = 0;
	}
	a0 = b0 = 0;
}

/*
SamplingState::SamplingState(unsigned int n_species)
{
	this->sn = 0;
	this->n_species = n_species;
	state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = 0;
	}
	this->prev = 0;
	this->next = 0;
	this->n_stoireactions = 0;
	this->rates_a = 0;
	this->rates_b = 0;
	this->gamma   = 0;
	a0 = b0 = 0;
}
*/

SamplingState::SamplingState(unsigned int n_species, unsigned int n_stoireactions)
{
	this->sn = 0;
	this->n_species = n_species;
	this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = 0;
	}
	this->prev = 0;
	this->next = 0;
	this->n_stoireactions = n_stoireactions;
	this->rates_a = new double[this->n_stoireactions];
	this->rates_b = new double[this->n_stoireactions];
	this->gamma   = new double[this->n_stoireactions];
	this->pfn = new double[this->n_stoireactions];
	this->pbn = new double[this->n_stoireactions];
	this->nVisits = 0;
	for (unsigned int i = 0; i < this->n_stoireactions; i++)
	{
		rates_a[i] = 0;
		rates_b[i] = 0;
		gamma[i] = 0;
		pfn[i] = 0;
		pbn[i] = 0;
	}
	a0 = b0 = 0;
}

SamplingState::SamplingState(int* values, unsigned int n_species, unsigned int n_stoireactions)
{
	this->sn = 0;
	this->n_species = n_species;
	this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = values[i];
	}
	this->prev = 0;
	this->next = 0;
	this->n_stoireactions = n_stoireactions;
	this->rates_a = new double[n_stoireactions];
	this->rates_b = new double[n_stoireactions];
	this->gamma   = new double[n_stoireactions];
	this->pfn = new double[this->n_stoireactions];
	this->pbn = new double[this->n_stoireactions];
	this->nVisits = 0;
	for (unsigned int i = 0; i < n_stoireactions; i++)
	{
		rates_a[i] = 0;
		rates_b[i] = 0;
		gamma[i] = 0;
		pfn[i] = 0;
		pbn[i] = 0;
	}
	a0 = b0 = 0;
}


SamplingState::~SamplingState()
{
	if (state != 0)   delete[] state;
	if (rates_a != 0) delete[] rates_a;
	if (rates_b != 0) delete[] rates_b;
	if (gamma != 0)   delete[] gamma;
	if (pfn!= 0)   delete[] pfn;
	if (pbn!= 0)   delete[] pbn;
}

int SamplingState::CopyfromState(SamplingState* astate)
{
	this->sn = astate->sn;
	this->n_species = astate->n_species;
	if (!this->state) this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate->state[i];
	}
	this->prev = astate->prev;
	this->next = astate->next;

	this->n_stoireactions = astate->n_stoireactions;
	if (!this->rates_a) this->rates_a = new double[this->n_stoireactions];
	if (!this->rates_b) this->rates_b = new double[this->n_stoireactions];
	if (!this->gamma)   this->gamma   = new double[this->n_stoireactions];
	if (!this->pfn)     this->pfn     = new double[this->n_stoireactions];
	if (!this->pbn)     this->pbn     = new double[this->n_stoireactions];
	for (unsigned int i = 0; i < this->n_stoireactions; i++)
	{
		this->rates_a[i] = astate->rates_a[i];
		this->rates_b[i] = astate->rates_b[i];
		this->gamma[i] = astate->gamma[i];
		this->pfn[i] = astate->pfn[i];
		this->pbn[i] = astate->pbn[i];
	}
	this->a0 = astate->a0;
	this->b0 = astate->b0;
	this->nVisits = astate->nVisits;

	return 1;
}

/*
int SamplingState::CopyfromState(State* astate)
{
	this->sn = astate->getSN();
	this->n_species = astate->getNumSpecies();
	this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate->getStateValue(i);
	}
	this->prev = 0;
	this->next = 0;
	this->n_stoireactions = 0;
	this->rates_a = 0;
	this->rates_b = 0;
	this->gamma   = 0;
	a0 = 0;
	b0 = 0;

	return 1;
}
*/

int SamplingState::CopyfromState(State* astate, unsigned int n_stoireactions)
{
	this->sn = astate->getSN();
	this->n_species = astate->getNumSpecies();
	if(!this->state) this->state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate->getStateValue(i);
	}
	this->prev = 0;
	this->next = 0;
	this->n_stoireactions = n_stoireactions;
	if (!this->rates_a) this->rates_a = new double[this->n_stoireactions];
	if (!this->rates_b) this->rates_b = new double[this->n_stoireactions];
	if (!this->gamma)   this->gamma   = new double[this->n_stoireactions];
	if (!this->pfn)     this->pfn     = new double[this->n_stoireactions];
	if (!this->pbn)     this->pbn     = new double[this->n_stoireactions];
	for (unsigned int i = 0; i < this->n_stoireactions; i++)
	{
		this->rates_a[i] = 0;
		this->rates_b[i] = 0;
		this->gamma[i] = 0;
		this->pfn[i] = 0;
		this->pbn[i] = 0;
	}
	this->a0 = 0;
	this->b0 = 0;
	this->nVisits = 0;

	return 1;
}

int SamplingState::PrintState()
{
	for (unsigned int j=0; j<this->n_species; j++)
        {
		cout << this->getStateValue(j);
		if (j < this->n_species-1) cout << " ";
	}
	cout << endl;
	/*
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		if (rates_a && rates_b && gamma)
			printf ("a[%d]=%e, b[%d]=%e, gamma[%d]=%f\n", j, rates_a[j], j, rates_b[j], j, gamma[j]);
	}
	cout << "a0=" << a0 << ", b0=" << b0 << endl;
	*/
	//cout << "A add: " << rates_a << ", B add: " << rates_b << ", G add: " << gamma << endl;
	/*
	*/

	return 1;
}

int SamplingState::PrintState_Full(ofstream& fout)
{
	for (unsigned int j=0; j<this->n_species; j++)
        {
		fout << this->getStateValue(j);
		if (j < this->n_species-1) fout << " ";
	}
	fout << ", ";
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		//if (rates_a && rates_b && gamma)
		//	printf ("a[%d]=%.3g  b[%d]=%.3g  g[%d]=%.3f; ", j, this->rates_a[j], j, this->rates_b[j], j, this->gamma[j]);
		if (rates_a && rates_b && gamma)
			fout << "a" << dec << j+1 << "=" << scientific << this->rates_a[j] << "  b" << dec << j+1 << "=" << scientific << this->rates_b[j] << "  g" << dec << j+1 << "=" << scientific << this->gamma[j] << ";  ";
	}
	fout << scientific << "a0=" << this->a0 << ", b0=" << this->b0 << ", nV=" << this->nVisits;
	
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		fout << ", pfn" << dec << j+1 << "=" << scientific << this->pfn[j];
	}
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		fout << ", pbn" << dec << j+1 << "=" << scientific << this->pbn[j];
	}
	fout << endl;
	/*
	*/
	//fout << "A add: " << rates_a << ", B add: " << rates_b << ", G add: " << gamma << endl;
	/*
	*/

	return 1;
}

int SamplingState::PrintState_Full()
{
	for (unsigned int j=0; j<this->n_species; j++)
        {
		cout << this->getStateValue(j);
		if (j < this->n_species-1) cout << " ";
	}
	cout << ", ";
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		if (rates_a && rates_b && gamma)
			cout << "a" << dec << j+1 << "=" << scientific << this->rates_a[j] << "  b" << dec << j+1 << "=" << scientific << this->rates_b[j] << "  g" << dec << j+1 << "=" << scientific << this->gamma[j] << ";  ";
	}
	cout << scientific << "a0=" << this->a0 << ", b0=" << this->b0 << ", nV=" << this->nVisits;
	
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		cout << ", pfn" << dec << j+1 << "=" << scientific << this->pfn[j];
	}
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		cout << ", pbn" << dec << j+1 << "=" << scientific << this->pbn[j];
	}
	cout << endl;

	return 1;
}

int SamplingState::setpFn(unsigned int i, double v)
{
	this->pfn[i] = v;
	return 1;
}

int SamplingState::setpBn(unsigned int i, double v)
{
	this->pbn[i] = v;
	return 1;
}

double SamplingState::getpFn(unsigned int i)
{
	return this->pfn[i];
}

double SamplingState::getpBn(unsigned int i)
{
	return this->pbn[i];
}

int SamplingState::PrintState(ofstream& fout)
{
	for (unsigned int j=0; j<this->n_species; j++)
        {
		fout << this->getStateValue(j);
		if (j < this->n_species-1) fout << " ";
	}
	fout << endl;
	/*
	for (unsigned int j=0; j<this->n_stoireactions; j++)
        {
		if (rates_a && rates_b && gamma)
			printf ("a[%d]=%e, b[%d]=%e, gamma[%d]=%f\n", j, rates_a[j], j, rates_b[j], j, gamma[j]);
	}
	cout << "a0=" << a0 << ", b0=" << b0 << endl;
	*/
	//cout << "A add: " << rates_a << ", B add: " << rates_b << ", G add: " << gamma << endl;
	/*
	*/

	return 1;
}

int SamplingState::PrintState_nonl()
{
	for (unsigned int j=0; j<this->n_species; j++)
        {
		cout << this->getStateValue(j);
		if (j < this->n_species-1) cout << " ";
	}

	return 1;
}


STATETYPE SamplingState::getStateValue(unsigned int i)
{
	return this->state[i];
}

unsigned int SamplingState::getNumSpecies()
{
	return this->n_species;
}

int SamplingState::setStateValue(unsigned int i, STATETYPE value)
{
	this->state[i] = value;
	return 0;
}


int SamplingState::setState(STATETYPE* astate, unsigned int n_species)
{
	this->sn = 0;
	this->n_species = n_species;
	state = new STATETYPE[this->n_species];
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = astate[i];
	}
	return 1;
}

double SamplingState::geta0 ()
{
	return a0;
}

int SamplingState::seta0 (double v)
{
	a0 = v;
	return 1;
}

int SamplingState::calca0 ()
{
	a0 = 0;
	for (unsigned int i = 0; i < this->n_stoireactions; i++)
	{
		a0 += rates_a[i];
	}
	return 1;
}

double SamplingState::getb0 ()
{
	return b0;
}

int SamplingState::setb0 (double v)
{
	b0 = v;
	return 1;
}

unsigned long SamplingState::getnVisits ()
{
	return nVisits;
}

int SamplingState::setnVisits (unsigned long v)
{
	this->nVisits = v;
	return 1;
}

int SamplingState::increase1Visit ()
{
	this->nVisits ++;
	return 1;
}

double SamplingState::calcB0 ()
{
	this->b0 = 0.0;
	for (unsigned int i=0; i<n_stoireactions; i++)
	{
		this->b0 += this->rates_b[i];
	}

	return this->b0;
}

SamplingState* SamplingState::getPrev()
{
	return this->prev;
}

SamplingState* SamplingState::getNext()
{
	return this->next;
}

int SamplingState::setPrev(SamplingState* pt)
{
	this->prev = pt;
	return 0;
}

int SamplingState::setNext(SamplingState* pt)
{
	this->next = pt;
	return 0;
}

int SamplingState::setSN (unsigned long index)
{
	this->sn = index;
	return 1;
}

unsigned long SamplingState::getSN ()
{
	return sn;
}

int SamplingState::initialize (unsigned int n_stoireactions)
{
	this->sn = 0;
	if (!state)
	{
		state = new STATETYPE[this->n_species];
	}
	for (unsigned int i = 0; i < this->n_species; i ++)
	{
		this->state[i] = 0;
	}
	this->prev = NULL;
	this->next = NULL;
	this->n_stoireactions = n_stoireactions;

	if (!rates_a)
	{
		rates_a = new double[this->n_stoireactions];
	}
	if (!rates_b)
	{
		rates_b = new double[this->n_stoireactions];
	}
	if (!gamma)
	{
		gamma = new double[this->n_stoireactions];
	}
	for (unsigned int i = 0; i < this->n_stoireactions; i ++)
	{
		this->rates_a[i] = 0;
		this->rates_b[i] = 0;
		this->gamma[i] = 0;
	}

	return 1;
}

int SamplingState::setRatesA (int j, double v)
{
	rates_a[j] = v;
	return 1;
}

int SamplingState::setRatesB (int j, double v)
{
	rates_b[j] = v;
	return 1;
}

int SamplingState::setGamma  (int j, double v)
{
	gamma[j] = v;
	return 1;
}

double SamplingState::getRatesA (int j)
{
	return rates_a[j];
}

double SamplingState::getRatesB (int j)
{
	return rates_b[j];
}

double SamplingState::getGamma  (int j)
{
	return gamma[j];
}

unsigned int SamplingState::getNumStoiReactions()
{
	return n_stoireactions;
}

int SamplingState::setNumStoiReactions(unsigned int n_stoireactions)
{
	this->n_stoireactions = n_stoireactions;
	return 1;
}

SamplingState& SamplingState::operator= (const SamplingState &pt)
{
	this->sn = pt.sn;
	this->n_species = pt.n_species;
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		this->state[i] = pt.state[i];
	}
	
	this->prev = pt.prev;
	this->next = pt.next;
	this->n_stoireactions = pt.n_stoireactions;
	this->n_species = pt.n_species;

	if (this->n_stoireactions > 0)
	{
		if (pt.rates_a && pt.rates_b && pt.gamma)
		{
			if (!(this->rates_a) && !(this->rates_b) && !(this->gamma))
			{
				this->rates_a = new double[this->n_stoireactions];
				this->rates_b = new double[this->n_stoireactions];
				this->gamma = new double[this->n_stoireactions];
			}
			for (unsigned int i = 0; i < this->n_stoireactions; i++)
			{
				this->rates_a[i] = pt.rates_a[i];
				this->rates_b[i] = pt.rates_b[i];
				this->gamma[i] = pt.gamma[i];
			}
			this->a0 = pt.a0;
			this->b0 = pt.b0;
		}
		if (pt.pfn)
		{
			if (!(this->pfn))
			{
				this->pfn = new double[this->n_stoireactions];
			}
			for (unsigned int i = 0; i < this->n_stoireactions; i++)
			{
				this->pfn[i] = pt.pfn[i];
			}
		}
		if (pt.pbn)
		{
			if (!(this->pbn))
			{
				this->pbn = new double[this->n_stoireactions];
			}
			for (unsigned int i = 0; i < this->n_stoireactions; i++)
			{
				this->pfn[i] = pt.pfn[i];
			}
		}
	}
	else if (this->n_stoireactions == 0)
	{
		this->rates_a = 0;
		this->rates_b = 0;
		this->gamma = 0;
		this->pfn = 0;
		this->pbn = 0;
		this->a0 = 0;
		this->b0 = 0;
		cout << "Zero number of reactions in SamplingState." << endl;
	}
	
	return *this;
}

bool SamplingState::operator== (const SamplingState &other)
{
	bool equal = true;
	if (this->n_species != other.n_species)
	{
		equal = false;
	}
	
	for (unsigned int i = 0; i < this->n_species; i++)
	{
		if (this->state[i] != other.state[i])
		{
			equal = false;
		}
	}
	
	return equal;
}

/*
void SamplingState::delete ()
{
	if (this == NULL) return;
	free (this->sn);
	free (this->n_species);
	free (this->next);
	free (this->state);
}
*/

// Special members for 2Dwalk boundary problem.
//int SamplingState::moveStateOnB (SamplingState* p, SamplingState* p1, int distance)
//{
//}




/* EOF */


