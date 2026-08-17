// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/*
 * Compute the steady state of the chemical master equation
 * through Runge-Kutta 4
 * like this: dx=A*x*dt
 *
 * RK4A.cpp
 * 9/5/2010
 * 
 * By Youfang Cao
 *
 */


#include <cmath>
#include <iostream>
#include <ctime>
#include "anyoption.h"

using namespace std;

double Xnext (double dt, int n, int nnz, int* irow, int* jcol, double* A, double* xin, double* xout);

int main(int argc, char* argv[])
{
  // Defining variables of A matrix file.

  char* Afile;
  char* Outfile;

  // Defining variables needed to store A in CSC format.

  int     n;     // Dimension of matrix.
  int     nnz;   // Number of nonzero elements in A.
  int*    irow;  // Row index of all nonzero elements of A.
  int*    jcol;  // Col index of all nonzero elements of A.
  double* A;     // Nonzero elements of A.

  double** X = 0;

  double* PSS = 0;

  long tic, tac;
  time_t start_t, end_t;
  double time_tot = 0;
  double time_used = 0;
  double dt = 0.01;
	double tol = 1E-15;

  tic = clock();
  start_t = time(NULL);

        /* CREATE AN OBJECT TO HANDLE COMMAND LINE */
        AnyOption *opt = new AnyOption();

        /* SET THE USAGE/HELP   */
        opt->addUsage( "" );
        opt->addUsage( "Usage: " );
        opt->addUsage( "" );
        opt->addUsage( " -h  --help  			Prints this help " );
        opt->addUsage( " -t  --matrix <tmfile> 		Transition matrix file" );
        opt->addUsage( " -T  --time 100 		Time interval" );
        opt->addUsage( " -d  --dt 0.01 			Time step" );
        opt->addUsage( " -e  --tol 1e-8			Error tolerance" );
        opt->addUsage( " -o  --output			Output file" );
        opt->addUsage( "" );

        /* SET THE OPTION STRINGS/CHARACTERS */
        opt->setFlag(  "help", 'h' );
        opt->setOption(  "matrix", 't' );
        opt->setOption(  "time", 'T' );
        opt->setOption(  "dt", 'd' );
        opt->setOption(  "tol", 'e' );
        opt->setOption(  "output", 'o' );

	/* go through the command line and get the options  */
        opt->processCommandArgs( argc, argv );

        /* GET THE VALUES */
        if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) 
	{
                opt->printUsage();
		delete opt;
		return 1;
	}
	if( opt->getValue( 'T' ) != NULL  || opt->getValue( "time" ) != NULL  )
		time_tot = atof(opt->getValue( 'T' ));
	if( opt->getValue( 'd' ) != NULL  || opt->getValue( "dt" ) != NULL  )
		dt = atof(opt->getValue( 'd' ));
	if( opt->getValue( 'e' ) != NULL  || opt->getValue( "tol" ) != NULL  )
		tol = atof(opt->getValue( 'e' ));
	if( opt->getValue( 't' ) != NULL  || opt->getValue( "matrix" ) != NULL  )
		Afile = opt->getValue( 't' );
	else
	{
		printf ("Please specify the matrix filename.\n\n");
		opt->printUsage();
		return 1;
	}
	if( opt->getValue( 'o' ) != NULL  || opt->getValue( "output" ) != NULL  )
		Outfile = opt->getValue( 'o' );
	else
	{
		printf ("Please specify the output filename.\n\n");
		opt->printUsage();
		return 1;
	}
	//cout << endl;

  FILE* fp = 0;
  fp = fopen (Afile, "r");
  fscanf (fp, "%d %d\n", &n, &nnz);
  //cout << "n = " << n << ", nnz =  " << nnz << endl;
  // Initialize the CSC matrix.
  irow = new int[nnz];
  jcol = new int[nnz];
  A = new double[nnz];
  // Allocate memory space for the initial vector.
  X = new double*[2];
  X[0] = new double[n];
  X[1] = new double[n];
  PSS = new double[n];  // solution.


  // Define three temporary arrays to read matrix of Expokit format first.
  // And then transform them to CSC format.
  char ch;
  for (long i = 0; i < nnz; i++)
  {
    fscanf (fp, "%d %d %lg", irow+i, jcol+i, A+i);
    //cout << irow[i] << " " << A[i] << endl;
    fscanf (fp, "%c", &ch);
  }

  // Initial vector.
  //InitVec = new (nothrow) double(n);
  X[0][0] = 1.0;
  X[1][0] = 0.0;
  double suminit = 0;
  for (int i = 0; i < n; i++)
  {
	X[0][i] = abs( double(rand())/double(RAND_MAX) );
	suminit += X[0][i];
  	//X[0][i] = 0.0;
  	X[1][i] = 0.0;
    //cout << InitVec[i] << endl;
  }
  for (int i = 0; i < n; i++)
  {
  	//X[0][i] /= suminit;
  	X[0][i] = 1.0/(double)n;
  }
  //cout << "Initial vector initialized..." << endl;

	double err = 1;
	int n0 = 0, n1 = 1;
	int tmp = 0;
	int counter = 0;

  if (time_tot == 0)
  {
	while (err > tol)
	{
		err = Xnext (dt, n, nnz, irow, jcol, A, X[n0], X[n1]);
		cout << counter << " " << err << endl;
		tmp = n0; n0 = n1; n1 = tmp;
		counter ++;
	}
  }
  else
  {
  	time_used = 0;
	while (time_used <= time_tot)
	{
		err = Xnext (dt, n, nnz, irow, jcol, A, X[n0], X[n1]);
		cout << "step=" << counter << ", t=" << time_used << ", err=" << err << endl;
		tmp = n0; n0 = n1; n1 = tmp;
		counter ++;
		time_used += dt;
	}
  }

	for (int i = 0; i < n; i++)
	{
		PSS[i] = X[n0][i];
	}
    

  // Printing computing time used.
  tac = clock();
  end_t = time(NULL);
  cout << "Total CPU time used: "  << (tac-tic)/CLOCKS_PER_SEC << " seconds" << endl;
  cout << "Total time used: " << end_t - start_t << " seconds" << endl;

  // Printing eigenvalues.

    //ev_one = 0;
    // end
    double ev_sum = 0;
    for (int i = 0; i < n; i++)
    {
      ////if (abs(EigVec[ev_one*n+i]) < ev_tol) EigVec[ev_one*n+i] = 0.0; // Set very tiny numbers to zero.
      //printf ("%0.10g\n", EigVec[ev_one*n+i]);
      ev_sum += PSS[i];
    }
    cout << endl << endl;

  if (time_tot == 0)
    cout << "Probability distribution for steady state" << endl;
  else
    cout << "Probability distribution for t=" << time_tot << endl;
    
  FILE* fpo = 0;
  fpo = fopen (Outfile, "w");
    fprintf (fpo, "1.0\n"); 
    for (int i = 0; i < n; i++)
    {
      fprintf (fpo, "%0.8E\n", PSS[i]/ev_sum);
    }

  // Release memory space.
  delete[] irow;
  delete[] jcol;
  delete[] A;
  //delete[][] X;
  delete[] PSS;

  fclose (fp);
  fclose (fpo);

} // main


double Xnext (double dt, int n, int nnz, int* irow, int* jcol, double* A, double* xin, double* xout)
{
	double *k1;
	double *k2;
	double *k3;
	double *k4;
	double *k0;

	k1 =  new double[n];
        k2 =  new double[n];
        k3 =  new double[n];
        k4 =  new double[n];
        k0 =  new double[n];
	

	for (int i=0; i<n; i++)
	{
		xout[i] = 0.0;
		k0[i] = k1[i] = k2[i] = k3[i] = k4[i] = 0.0;
	}

	for (int i=0; i<nnz; i++)
	{
		k1[irow[i]] += A[i] * xin[jcol[i]];
		//cout << xout[irow[i]] << endl;
	}

	for (int i=0; i<n; i++)
	{
		k0[i] = xin[i] + 0.5 * dt * k1[i];
	}

	for (int i=0; i<nnz; i++)
	{
		k2[irow[i]] += A[i] * k0[jcol[i]];
		//cout << xout[irow[i]] << endl;
	}

	for (int i=0; i<n; i++)
	{
		k0[i] = xin[i] + 0.5 * dt * k2[i];
	}

	for (int i=0; i<nnz; i++)
	{
		k3[irow[i]] += A[i] * k0[jcol[i]];
		//cout << xout[irow[i]] << endl;
	}

	for (int i=0; i<n; i++)
	{
		k0[i] = xin[i] + dt * k3[i];
	}

	for (int i=0; i<nnz; i++)
	{
		k4[irow[i]] += A[i] * k0[jcol[i]];
		//cout << xout[irow[i]] << endl;
	}

	double h6 = dt / 6.0;
	double err = 0.0;
	for (int i = 0; i < n; i++)
	{
		xout[i] = xin[i] + h6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		double delta = abs(xout[i] - xin[i]);
		if (delta > err)
		{
			err = delta;
		}
	}

	/*for (int i = 0; i < n; i++)
	{
		err += (xout[i] - xin[i]) * (xout[i] - xin[i]);
	}
	err = sqrt(err);
	*/

	delete[] k0;
	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;

	return err;
}



