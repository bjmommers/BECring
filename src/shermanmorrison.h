/*******************************************************************************
Sherman Morrison header
*******************************************************************************/

//#pragma once in place of header guard
#pragma once

#include<cmath>
#include<complex>
#include "bec-headers.h"


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

int osoevo
/* Optimised Split-Operator time evolution, per Clinton Jackson's work
   Wrapper function that incorporates updatepotential, prepSMvecs, calcrrhs, cyctridag
   While this requires 3 solves per timestep (=6 tridag calls), it is fourth-order accurate
*/
(
//from updatepotential()
vector<double> &x,
vector<double> &potential,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2, //two-component
int impot,
int imm,
double fweight,
int rt,
double gii, 
double gij,
double omega1, 
double omega2, 
double t, 
double potamp,
double sigma,
double gdepth,
int husimitype,
double imprintamp,
double impoffset,
double impsigma,
int imprinttype,
int imptimeprofile,
double omegaou,
double printoffset,
//from prepSMvecs
vector< complex<double> > &zadiag, 
vector< complex<double> > &zaupper, 
vector< complex<double> > &zalower,
vector< complex<double> > &zuvec,
vector< complex<double> > &zvvec,
complex<double> alpha,
complex<double> beta,
complex<double> gamma,
//from calcrrhs()
vector< complex<double> > &zrrhs,
vector< complex<double> > &zoldpsi,
int Nx,
//from cyctridag()
vector< complex<double> > &zxsoln,
int n,
//additional required parameters
double tstep
)

{
    /* What this function needs to do:
     * 1: updatepotential() using the average potential
     * 2: prepSMvecs() for the first iteration of time evolution (forward by s\lambda)
     * 3: calcrrhs() 
     * 4: cyctridag() to evolve forward in time
     * 5: prepSMvecs() for the second iteration of time evolution (back by (1-2s)\lambda)
     * 6: calcrrhs()
     * 7: cyctridag() to evolve back in time
     * 8: prepSMvecs() for final iteration of time evolution (forward by s\lambda to t_f)
     * 9: calcrrhs()
     * 10: cyctridag()
     */
    
    int ret; //return value for checking
    double s = 1.0/(2.0-pow(2.0,1.0/3.0)); //s-value from Clinton's work
    double alt_tstep = tstep;
    
    if(t == 0.0) { alt_tstep = 0.0; } // for imaginary time
    
    
    //1: updatepotential using half-timestep (gives average potential over the duration of the timestep)
//    updatepotential(x, potential, psi1, psi2, impot, imm, fweight, rt, gii, gij, omega1, omega2, t+(alt_tstep/2.0), potamp, sigma, gdepth, husimitype, imprintamp, impoffset, impsigma, imprinttype, imptimeprofile, omegaou, printoffset);
//    updatepotential(x, potential, psi1, psi2, impot, imm, fweight, rt, gii, gij, omega1, omega2, t, potamp, sigma, gdepth, husimitype, imprintamp, impoffset, impsigma, imprinttype, imptimeprofile, omegaou, printoffset);
    
    //2: prepSMvecs for first iteration: multiply alpha/beta/gamma by s
    // this works because they're all linear in (\Delta t), so we can scale time easily
    prepSMvecs(zadiag,zaupper,zalower,zuvec,zvvec,potential,alpha*s,beta*s,gamma*s);
//    prepSMvecs(zadiag,zaupper,zalower,zuvec,zvvec,potential,alpha,beta,gamma);
    
    //3: calculate right hand side for first step
    calcrrhs(zrrhs, zoldpsi, potential, alpha*s, beta*s, gamma*s, Nx);
//    calcrrhs(zrrhs, zoldpsi, potential, alpha, beta, gamma, Nx);

    //4: perform first time evolution step
    ret = cyctridag (zadiag, zalower, zaupper, zrrhs, zxsoln, n, zuvec, zvvec);
    if(ret != 0) { return 1;} //check if cyctridag threw an error
    
    //5: prepSMvecs for backward time evolution
    prepSMvecs(zadiag,zaupper,zalower,zuvec,zvvec,potential,alpha*(1.0-2.0*s),beta*(1.0-2.0*s),gamma*(1.0-2.0*s));
    
    //6: calcrrhs
//    calcrrhs(zrrhs, zoldpsi, potential, alpha*(1.0-2.0*s), beta*(1.0-2.0*s), gamma*(1.0-2.0*s), Nx);
    calcrrhs(zrrhs, zxsoln, potential, alpha*(1.0-2.0*s), beta*(1.0-2.0*s), gamma*(1.0-2.0*s), Nx);
    
    //7: cyctridag for backward time evolution
    ret = cyctridag (zadiag, zalower, zaupper, zrrhs, zxsoln, n, zuvec, zvvec);
    if(ret != 0) { return 1;} //check if cyctridag threw an error
    
    //8: prepSMvecs for final iteration
    prepSMvecs(zadiag,zaupper,zalower,zuvec,zvvec,potential,alpha*s,beta*s,gamma*s);
    
    //9: calcrrhs
//    calcrrhs(zrrhs, zoldpsi, potential, alpha*s, beta*s, gamma*s, Nx);
    calcrrhs(zrrhs, zxsoln, potential, alpha*s, beta*s, gamma*s, Nx);
    
    //10: cyctridag
    ret = cyctridag (zadiag, zalower, zaupper, zrrhs, zxsoln, n, zuvec, zvvec);
return ret;
}
    

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

