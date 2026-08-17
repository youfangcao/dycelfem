// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
#ifndef MATRIXELEMENT_H_
#define MATRIXELEMENT_H_

#define VALUETYPE double

#include <vector>

using namespace std;

/*
 * Define a class for building the linked list
 * used for storing of non-zero elements of sparse
 * transition matrix.
 */
class MatrixElement
{
private:
	unsigned long   i; // to state i
	unsigned long   j; // from state j
	VALUETYPE       value; // transition rate
	vector<unsigned int>*        r; // by which reaction
	MatrixElement* next;

	
public:
	MatrixElement();
	MatrixElement(unsigned int i, unsigned int j, VALUETYPE v);
	MatrixElement(unsigned int i, unsigned int j, VALUETYPE v, vector<unsigned int>& rv);

	virtual ~MatrixElement();

	unsigned int SetValues (unsigned int i, unsigned int j, VALUETYPE v, unsigned int k);

	unsigned long getI();
	unsigned long getJ();
	VALUETYPE getValue();
	vector<unsigned int>* getR();
	MatrixElement* getNext();
	void setNext(MatrixElement* next);
};


#endif /*MATRIXELEMENT_H_*/

