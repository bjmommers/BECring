/*******************************************************************************
Sherman Morrison header
*******************************************************************************/

//#pragma once in place of header guard
#pragma once

#include<cmath>
#include<complex>


int cyctridag
//Sherman-Morrison method for tridiagonal systems with cyclic boundary conds
(vector< complex<double> > &zadiag, vector< complex<double> > &zalower, 
vector< complex<double> > &zaupper, vector< complex<double> > &zrrhs, 
vector< complex<double> > &zxsoln, int n, vector< complex<double> > &uvec,
vector< complex<double> > &vvec);


int tridag
//Tridiagonal system solver based on numerical recipes 3rd edition
(vector< complex<double> > &zadiag, vector< complex<double> > &zalower, 
vector< complex<double> > &zaupper, vector< complex<double> > &zrrhs, 
vector< complex<double> > &zxsoln, int n);


/********************************DEFINITIONS***********************************/

/*
int osoevo
// Optimised Split-Operator time evolution, per Clinton Jackson's work
(
vector< complex<double> > &zadiag,
vector< complex<double> > &zalower,
vector< complex<double> > &zaupper,
vector< complex<double> > &zrrhs, 
vector< complex<double> > &zxsoln,
int n, 
vector< complex<double> > &uvec,
vector< complex<double> > &vvec
)

{
    
}
*/

/*----------------------------------------------------------------------------*/

int cyctridag
//Sherman-Morrison method for tridiagonal systems with cyclic boundary conds
(
vector< complex<double> > &zadiag,
vector< complex<double> > &zalower,
vector< complex<double> > &zaupper,
vector< complex<double> > &zrrhs, 
vector< complex<double> > &zxsoln,
int n, 
vector< complex<double> > &uvec,
vector< complex<double> > &vvec
)

{

	// initialise variables, array pointers
	int ierry, ierrz, j; //for tridag calls
	vector< complex<double> > zysoln, zzsoln;
	complex<double> zvdoty, zvdotz;
	//working vectors
	zysoln.resize(n);
	zzsoln.resize(n);
	zvdoty = complex<double> {0.0, 0.0};
	zvdotz = complex<double> {0.0, 0.0};

	//for y, zrrhs is rhs
	ierry = tridag(zadiag,zalower,zaupper,zrrhs,zysoln,n);

	// for z, uvec is rhs
	ierrz = tridag(zadiag,zalower,zaupper,uvec,zzsoln,n);

	// check all is well before proceeding
	if((ierry != 0) || (ierrz != 0)){
		printf("tridag failed! exiting...\n");
		exit(EXIT_FAILURE);
	}

	// dot products for final solution
		zvdoty = zysoln[0] + (vvec[n-1]*zysoln[n-1]);
		zvdotz = zzsoln[0] + (vvec[n-1]*zzsoln[n-1]);

	// evaluate final solution
	for(j=0;j<n;j++){
		zxsoln[j] = zysoln[j] - (zzsoln[j]*(zvdoty/(1.0+zvdotz)));
	}

// solution is returned in zxsoln array
return 0;
}

/*----------------------------------------------------------------------------*/

int tridag
//Tridiagonal system solver based on numerical recipes 3rd edition
(
vector< complex<double> > &zadiag,
vector< complex<double> > &zalower, 
vector< complex<double> > &zaupper,
vector< complex<double> > &zrrhs, 
vector< complex<double> > &zxsoln,
int n
)

{ // modified to handle complex data
	int j;
	complex<double> zbet {0.0, 0.0};
	vector< complex<double> > zgam;

	/* if the following happens rethink your equations */
	if (abs(zadiag[0]) == 0.0) return -1;    // Error

	/* one vector of workspace is needed */
//	zgam = malloc(n*sizeof(complex<double>));
	zgam.resize(n);

	zbet=zadiag[0];
	zxsoln[0]=zrrhs[0]/zbet;

	for (j=1;j<n;j++) {
		zgam[j]=zaupper[j-1]/zbet;
		zbet=zadiag[j]-zalower[j]*zgam[j];
		if (abs(zbet) == 0.0)
        {
            zxsoln[j] = {0.0,0.0}; //error will kill wavefunction
        }//return -1;   // Pivoting Error
		zxsoln[j]=(zrrhs[j]-zalower[j]*zxsoln[j-1])/zbet;
	}
	
	for (j=n-2;j>=0;j--) {
		zxsoln[j]-= zgam[j+1]*zxsoln[j+1];  // Backsubstitution
	}

return 0;
}

