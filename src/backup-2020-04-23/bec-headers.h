/*******************************************************************************
NLSE solver header file
*******************************************************************************/

//#pragma once in place of header guard
#pragma once

#include<complex>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<vector>
#include<functional>
#include<cstdarg>
#include<valarray>
#include<limits>
 
const double PI = 3.141592653589793238460;
const double hbar = 1.054571817e-34; //Joule seconds
const double rbmass = 1.4431609e-25;   //kg
const double kboltz = 1.381e-23;     //Joule/Kelvin

const double globtol = __DBL_EPSILON__;
 
//constant definitions
#define PI 3.14159265358979323846264338

//set default namespace
using namespace std;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

//const complex<double> I(0.0,1.0);

/************************FUNCTION DECLARATIONS*********************************/

int loadinput
//Take input from file, place into double array
(
string filename, vector<string> &inputlist
);

void loadmulticoldata
// load phase data into multi-column vector
(
string filename,
vector<double> &phaseinput
);

void extractrotation
// calculates rotation rate from phase data
(
vector<double> &phasevec,
vector<double> &gradvec,
vector<double> &fourierlist,
double tstep,
double threshold,
int Nt,
int nfreq,
int imm
);

void multicolfp
// print vector to file in columns
(
string filename,
vector<double> &data,
int ncolumns,
string colhead
);

double calcoverlap
// calculate overlap of 2 waavefunctions
(
vector< complex<double> > wf1,
vector< complex<double> > wf2,
int Nx,
double xstep
);


void setxgrid
//Create a grid from upper/lower bounds and spacing
(
vector<double> &x, double upper, double lower, int Nx
);


int setinitialGaussian
// calculate initiali Gaussian waveform
// without kick (added if no imag time set)
(
vector< complex<double> > &waveform, vector<double> &x, double x0, 
double sigma0, double xstep, double upper, double lower, int Nx
);

void addkick
// apply momentum kick to wavefunction
(
 vector< complex<double> > &waveform, vector<double> &x, 
 double k0, int Nx
);

void addmodekick
// apply quantisedmomentum kick to wavefunction
(
 vector< complex<double> > &waveform,
 vector<double> &x,
 double x0,
 double k0, 
 double tmphase,
 int Nx
);

void addhusimipotential
// calculates Husimi potential term and adds it to potential vector
(
vector<double> &potential, vector<double> &x, double t, double potamp,
double omega1, double omega2, double gdepth, double sigma, int husimitype
);

void addtdpotential
// calculates time-dependent imprinting potential and adds it to potentialvec
(
vector<double> &potentialvec,
vector<double> &x,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2,
double t,
double sigma,
int impot,
int imm,
double fweight,
double gii,
double gij,
double imprintamp,
int imptimeprofile,
double impoffset,
double impsigma,
double omegaou,
double printoffset
);

double calcenergy
// calculate energy
(
vector< complex<double> > &waveform, vector<double> &potential, 
vector<double> &x, double xstep, double g, double omega
);

double calcenergy
// calculate energy, needs rotation energy
(vector< complex<double> > &waveform1, vector< complex<double> > &waveform2,
vector<double> &potential, vector<double> &x, double xstep, double gself,
double gint
);

double calclabenergy
// calculate energy
(
vector< complex<double> > &waveform, vector<double> &potential, 
vector<double> &x, double xstep, double g
);

double calclabenergy
// calculate energy 
(vector< complex<double> > &waveform1, vector< complex<double> > &waveform2,
vector<double> &potential, vector<double> &x, double xstep, double gself,
double gint
);

double calclz
//calculate expectation value of rotation operator Lz
(
vector< complex<double> > &waveform,
vector<double> &x,
double xstep
);

void calcstats
// calculate generic wavefunction statistics and place in statlist vector
(
vector< complex<double> > &waveform, vector< complex<double> > &oldwaveform,
vector<double> &potential, vector<double> &statlist, vector<double> &x, 
double xstep, int t, double timestep, int renorm, double omegaou
);


void calcstats
// calculate generic wavefunction statistics and place in statlist vector
(
vector< complex<double> > &waveform, vector< complex<double> > &oldwaveform,
vector<double> &potential, vector<double> &statlist, vector<double> &x, 
double xstep, int t, double timestep, int renorm, int modulus, double omegaou
);

void aggregatestats
// combine stats for both components into a single vector for printing
(
vector<double> &statvec1,
vector<double> &statvec2,
vector<double> &aggstatvec
);

void statfileprint
// print data contained in statlist vector to file
(
ofstream &outfile,
vector<double> &statlist
);


void wffileprint
// print whole vectors to file as columns
(
ofstream outfile, vector<vector<double> * > veclist, int vecsize
);


void calcrrhs
// Calculate right hand side vector for Sherman-Morrison method w/ Husimi 
// potential
(
vector< complex<double> > &zrrhs, vector< complex<double> > &zoldpsi,
vector<double> &potential, complex<double> zalpha, complex<double> zbeta,
complex<double> zgamma, int Nx
);


void updatediagonal
// update the values of Sherman-Morrison matrix diagonal terms for
// time-dependent potential
(
vector< complex<double> > &zadiag, vector<double> &potential, 
complex<double> &zalpha, complex<double> &zbeta, int Nx
);


void twocompfileprint
// print two-component waveforms to file
(ofstream &outfile,string filename,vector<double> &xgrid,
vector< complex<double> > &waveform1,vector< complex<double> > &waveform2,
vector<double> &potential1,vector<double> &potential2
);

/*
void updatepotential
//re-calculates the potential function potentialfunc at time t
(vector<double> &x, vector<double> &potentialvec,
vector< complex<double> > &psi1, int impot, int imm, int rt, double gii = 0.0, 
double omega1 = 0.0, double omega2 = 0.0, double t = 0.0, double potamp = 1.0,
double sigma = 1.0, double gdepth = 1.0, int husimitype = 1
);
*/

void updatepotential
//re-calculates the potential function potentialfunc at time t; two component
(
vector<double> &x,
vector<double> &potentialvec,
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
double printoffset
);

void addnonlinearity 
//one component, calculates interaction potential, adds to vector
(
vector<double> &potential, vector< complex<double> > &psi1, double gii 
);

void addnonlinearity //two components, calculates interaction potential
(
vector<double> &potential,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2,
double gii,
double gij
);

void calcVimprint //calculate imprinting potential
(
vector<double> &x,
vector<double> &potentialvec,
vector< complex<double> > &psi1, //one component
int impot,
int imm,
double fweight,
double gii,
double potamp,
double printoffset
);


void calcVimprint //calculate imprinting potential
(
vector<double> &x,
vector<double> &potentialvec,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2, //two-component
int impot,
int imm,
double fweight,
double gii, 
double gij,
double potamp,
double printoffset
);

double heaviside //step function
(
double x
);

double sgn
// implementation of sgn() function for square wave
(
double x
);

void prepSMvecs
// prepare A matrix tridiagonals, u and v vectors for Sherman-Morrison
(
vector< complex<double> > &zadiag, vector< complex<double> > &zaupper, 
vector< complex<double> > &zalower, vector< complex<double> > &zuvec,
vector< complex<double> > &zvvec, vector<double> &potential,
complex<double> alpha, complex<double> beta, complex<double> gamma
);

void prepSMvecs
// prepare A matrix tridiagonals, u and v vectors for Sherman-Morrison
// overloaded version without upper, lower vecs (use this for BEC component 2)
// debug: may need to remove/comment out
(
vector< complex<double> > &zadiag, vector< complex<double> > &zuvec,
vector< complex<double> > &zvvec, vector<double> &potential,
complex<double> alpha, complex<double> beta, complex<double> gamma
);

void updatediag
// update potential-dependent matrix diagonal
(
vector< complex<double> > &zadiag,
vector<double> &potential,
complex<double> zalpha,
complex<double> zbeta
);

complex<double> fixedmodefourierint
// perform Fourier transform of density at fixed mode number by integration
(vector< complex<double> > &wavefunction, vector<double> &x, double m);

complex<double> fixedmodedensityfourierint
// perform Fourier transform of density at fixed wavenumber by integration
(vector< complex<double> > &wavefunction,vector<double> &x,double m);

complex<double> fixedkfourierint
// placeholder: Fourier transform at fixed wavenumber
(vector< complex<double> > &wavefunction, vector<double> &x, double k);

double complexphase
// calculate the phase angle of a complex quantity using datan
// mostly exists to avoid calling fourierint funcs twice on the same data
(complex<double> cnumber);

double datan
// atan implementation
(double x, double y);

double dvecavg
// returns average value of given vector<double>
(vector<double> &vec);

void col2vec
// copies the contents of nth matrix column into empty vector
(
vector<double> &source,
vector<double> &dest,
int ncolumns, int column, int header);

void vecmax
// finds and returns indices of local maxima 
(vector<double> &vec, vector<int> &maxima);

string makeheader
// creates file column headers from list of frequencies
(vector<double> freqvec);

void setomegaou
// set oscillator units scaling frequency for unit conversions
(double &omegaou,double ztrapfreq,double rtrapfreq,int ouselect,double omega1,
double husimisigma,double husimigdepth,int husimitype
);

void setoscunits
//convert quantities to osc units
( double omegaou, double &xupper, double &xlower, double &tmax, double &tstep);

vector<double> threemodeoptimise
// gradient descent function,  only optimises on a,delta,a0
(// can't pass by reference here
 double a0,     //constant background
 double a,      //average modulation amplitude
 double delta,  //mode splitting
 double l,      //angular momentum quantum number
 double omega,  //rotation rate
 double alpha,  //imprinting potential amplitude
 double res,    //search resolution
 double gii     //interaction strength
);

double threemodeHexpect
// hamiltonian expectation value function for three-mode model optimisation
// ignores potential constant V0 (since it's constant)
(
 double a0,     //constant background
 double a,      //average modulation amplitude
 double delta,  //mode splitting
 double l,      //angular momentum quantum number
 double omega,  //rotation rate
 double alpha,  //imprinting potential amplitude
 double gii     //interaction strength
);

double threemodeHexpect_manual
// expectation value of Hamiltonian for wavefunction constructed using
// given parameters
//
// 
// TESTING ONLY FOR NOW, SOME VALUES HARDCODED
// IF VALID, REWORK THIS FUNCTION
(
 double a0,     //constant background
 double a,      //average modulation amplitude
 double delta,  //mode splitting
 double l,      //angular momentum quantum number
 double omega,  //rotation rate
 double alpha   //imprinting potential amplitude
);

vector< complex<double> > calcthreemode
// returns a three-mode model state
(
 vector<double> &x,
 double a0,
 double a,
 double delta,
 double l
);

/************************FUNCTION DEFINITIONS**********************************/

int loadinput
	(
	string filename,
	vector<string> &inputlist
	)
// load input parameters from file into double array
	{
		ifstream infile;
		string line;

		infile.open(filename.c_str());
		while(getline(infile,line,'\n'))
		{
			if((line.size() == 0) || (line.at(0) == '#'))
				{/*do nothing*/}
			else
				{
					inputlist.push_back(line);
				}
		}

	return 0;
	}

/*----------------------------------------------------------------------------*/

void loadstatefromfile
// loads initial state from file containing only data
(
string statefilename,
vector< complex<double> > &statevector
)

{
    ifstream infile;
    string line;
    stringstream tempstream;
    vector< complex<double> > inputlist;
    complex<double> tempcomplex;
    unsigned int j;

    cerr << "Loading initial state from file" << endl;
    infile.open(statefilename.c_str());
    cerr << "Initial state file opened successfully" << endl;
    while(getline(infile,line,'\n'))
    {
        tempstream << line;
        tempstream >> tempcomplex;
        inputlist.push_back(tempcomplex); //does this correctly convert to dcomplex?
    }

    if (inputlist.size() != statevector.size()) 
    {
        throw "Input state file size does not match grid size!";
        // does this terminate?
    }
    
    
    for(j=0;j<statevector.size();j++)
    {
        statevector[j] = inputlist[j];
    }
    //may need to print both vectors for debugging

return;
}

/*----------------------------------------------------------------------------*/

void savestatetofile
// saves state to file containing only data
(
string statefilename,
vector< complex<double> > &statevector
)

{
    ofstream outfile;
    unsigned int j;

    outfile.open(statefilename.c_str());
    for(j=0;j<(statevector.size()-1);j++)
    {
        outfile << statevector[j] << "\n";
    }
    //write last line without trailing newline?
    outfile << statevector[statevector.size()-1] << endl;

return;
}

/*----------------------------------------------------------------------------*/

void loadmulticoldata
(
string filename,
vector<double> &phaseinput
)

{
    ifstream phasefile;
    string line;
    string data;

    phasefile.open(filename.c_str());
    while(getline(phasefile,line,'\n'))
    {
        istringstream buff; //needs to be declared here to reset flags each loop

        buff.str(line);  // load line into stringstream to use another getline

        if(line.at(0) == '#') 
        {/*do nothing*/}
        else
        {
            while(getline(buff,data,'\t'))
            {
//                cout << "data: " << data << endl;
                phaseinput.push_back(atof(data.c_str()));
            } 

            buff.str("");

        }

    }

}

/*----------------------------------------------------------------------------*/

void extractrotation
(
vector<double> &phasevec,
vector<double> &gradvec,
vector<double> &fourierlist,
double tstep,
double threshold,
int Nt,
int nfreq,
int imm
)

{
    int j,k;
    double gradient;
    
    for(j=0;j<Nt;j++) //k=0 (i.e. time) column
    {
        gradvec[j*(nfreq+1)] = phasevec[j*(nfreq+1)];
    }

     
    for(k=0;k<nfreq;k++)
    {
        //gradient at edges using forward/backward space finite diff
        gradvec[k+1] = (phasevec[(nfreq+1)+(k+1)] - phasevec[k+1])/(tstep*fourierlist[k]);
        if(abs(gradvec[k+1]) > threshold) {gradvec[k+1] = 0.0;}

        gradvec[(Nt-1)*(nfreq+1)+(k+1)] = (phasevec[(Nt-1)*(nfreq+1)+(k+1)] - phasevec[(Nt-2)*(nfreq+1)+(k+1)])/(tstep*fourierlist[k]);
        if(abs(gradvec[(Nt-1)*(nfreq+1)+(k+1)]) > threshold) {gradvec[(Nt-1)*(nfreq+1)+(k+1)] = 0.0;}

         
        for(j=1;j<Nt-1;j++)
        {
            gradient = (phasevec[((j+1)*(nfreq+1))+(k+1)] - phasevec[((j-1)*(nfreq+1))+(k+1)])
                        /(2.0*tstep);
            if(abs(gradient) > threshold) // threshold is arbitrary, just need to avoid asymptotes
            {
                gradient = 0.0;
            }
            else
            {
                gradvec[(j*(nfreq+1))+(k+1)] = gradient/fourierlist[k];
            }
        }
    }

/*
 * Replace this function with the below?
 *
    vector<double> results;
    double oneosc;
    double radius;
    double numosc;
    int j;
    int index;
*/
    //determine stroboscopic times from analytics
    /*
     * Here we use:
     * t = (n \pi \hbar)/(E-E0)
     *   = (n \pi) / (l^2/2R^2)
     * with n=1 to give the oscilation period
     */
/*
    //calculate the ring radius, and therefore the oscillation time
    radius = (xmax-xmin)/(2.0*PI);
    oneosc = PI/((imm*imm)/(2.0*radius*radius));

    //determine the number of full oscillations in the sim, size results vector
    numosc = floor(tmax/oneosc);
    results.resize(numosc);
    timevec.resize(numosc);


    for(j=0;j<numosc;j++)
    {
        //todo: include both odd and even measurements (should have 2n+1 total measurements,
        //      where n is the number of oscillations
        results[j] = phasevec[(int) (((j+1)*oneosc)/tstep)] - phasevec[(int) ((j*oneosc)/tstep)];
        timevec[j] = (j+1)*oneosc;
    }
*/


return;
}

/*----------------------------------------------------------------------------*/

void multicolfp
(
string filename,
vector<double> &data,
int ncolumns,
string colhead
)

{
    int j,k;
    ofstream outfile;

    outfile.open(filename);

    outfile << colhead << endl;

    for(j=0;j< (int) data.size();j++)
    {
        outfile << setw(8) << data[j];

        if((((j+1)%(ncolumns)) == 0) && (j>0))
        {
            outfile << '\n';
        }
        else
        {
            outfile << '\t';
        }
    }

return;
}

/*----------------------------------------------------------------------------*/

double calcoverlap
// calculate overlap of 2 waavefunctions
(
vector< complex<double> > wf1,
vector< complex<double> > wf2,
int Nx,
double xstep
)

{
    int j;
    double sum = 0.0;

    for(j=0;j<Nx;j++)
    {
        sum += fmin(pow(abs(wf1[j]),2.0),pow(abs(wf2[j]),2.0))*xstep;
    }

return sum;
}

/*----------------------------------------------------------------------------*/

void setxgrid
//Create a grid from upper/lower bounds and spacing
(
vector<double> &x, double upper, double lower, int Nx
)

{
	int j;
	double spacing = (upper-lower)/((double) Nx);

	for(j=0;j<Nx;j++)
	{
		x[j] = (lower + (j * spacing));
	}
	
	cout << "x grid runs from " << x[0] << " to " << x[Nx-1] << endl;
}

/*----------------------------------------------------------------------------*/

int setinitialGaussian
// calculate initial Gaussian waveform
(
vector< complex<double> > &waveform,
vector<double> &x,
double x0,
double sigma0,
double xstep,
double upper,
double lower,
int Nx
)

{
	int j;
	double pnorm = 0.0;
	double dx;

    if(abs(sigma0) < globtol) 
    { // setting width to zero instead produces flat waveform
         
        //pragma omp for
	    for(j=0;j<Nx;j++)
        {
            waveform[j] = 1.0;
            pnorm += xstep*pow(abs(waveform[j]),2.0);
        }
    }
    else
    {

         
        //pragma omp for reduction(+:pnorm)
	    for(j=0;j<Nx;j++)
	    {
    		dx = fmin(x[j]-x0, (upper-lower)-(x[j]-x0));
    		waveform[j] = sqrt(1.0/(sqrt(PI)*sigma0))
    					*exp((-(pow(dx,2.0)))
    					/(2.0*pow(sigma0,2.0)));
//    					*exp(k0 * x[j] * complex<double> (0.0,1.0)); //momentum kick
    		pnorm += xstep*pow(abs(waveform[j]),2.0);
    	}
    }

// normalisation check
	if(abs(pnorm - 1.0) > globtol) // is this dangerous?!
	{
        //pragma omp parallel
        //pragma omp for
		for(j=0;j<Nx;j++)
		{
			waveform[j] /= sqrt(pnorm);
		}
	}


return 0;
}

/*----------------------------------------------------------------------------*/

void initialisethreemode
// calculate three-mode model initial state
(
 vector< complex<double> > &waveform,
 vector<double> &x,
 double tmphase,
 double bgmode,   //a0
 double modprop,  //a
 double modebias, //delta
 int immode, //l
 double omega,
 double alpha,
 double res,
 double gii
)

{  
    int j;
    int Nx = x.size();
    double xstep = x[1] - x[0];
    double pnorm = 0.0;
    double l = (double) immode;
    double initial_a0, initial_a, initial_delta;
    complex<double> i = {0.0,1.0};
    vector<double> params;
    
    cout << "Inititalising three-mode state" << endl;

    //todo: calculate initial guess values
    initial_a0 = bgmode;
    initial_a  = modprop;
    initial_delta = modebias;

    params.resize(6);
    // do a nonlinear search first
    params = threemodeoptimise(initial_a0,initial_a,initial_delta,immode,omega,alpha,res,0.0);
    // use the nonlinear result for the full linear result
    if(gii != 0.0)
    {
        params = threemodeoptimise(params[0],params[1],params[2],immode,omega,alpha,res,gii);
    }

    double new_a0 = params[0];
    double new_a = params[1];
    double new_delta = params[2];

    /*
    new_a0 = bgmode;
    new_a = (2.0/pow(l,2.0))*pow(l*omega,2.0) + 0.25*alpha*new_a0;
    new_delta = ((2.0*omega)/(double) l)*new_a;
    */


//    cout << "new vars: " << params[0] << "\t" << params[1] << "\t" << params[2] << endl;

//    waveform = calcthreemode(x,new_a0,new_a,new_delta,immode);

    for(j=0;j<Nx;j++)
    {
        waveform[j] = (1.0/sqrt(2.0*PI*(new_a0*new_a0 + 2.0*new_a*new_a 
                                                      + 2.0*new_delta*new_delta) ) ) 
                      *(new_a0 + (new_a + new_delta)*exp( i*l*2.0*PI*((double) j/(double) Nx)) 
                               + (new_a - new_delta)*exp(-i*l*2.0*PI*((double) j/(double) Nx)));
        pnorm += pow(abs(waveform[j]),2.0)*((2.0*PI)/(double) Nx);
                      
    }

    cout << "Three mode norm: " << pnorm << endl;

    if(abs(pnorm - 1.0) > globtol)
    {
        cout << "Normalising three-mode state with factor 1/" << sqrt(pnorm) << endl;
        for(j=0;j<Nx;j++)
        {
            waveform[j] /= sqrt(pnorm);
        }
    }


//    for(j=0;j<Nx;j++)
//    {  // e^ikx where k=number of modes
//        waveform[j] = modprop*((0.5+modebias)*exp(i*(mode*2.0*PI*((x[j]-x[0])/(x[Nx-1]-x[0]) + tmphase)))
//                    + (0.5-modebias)*exp(-i*mode*2.0*PI*((x[j]-x[0])/(x[Nx-1]-x[0])))) // no offset
//                    + (1.0-modprop);
//                    
//        waveform[j] = modprop*((0.4+modebias)*exp(i*(mode*2.0*PI*(((double) j / (double) Nx)  + tmphase)))
//                    + (0.4-modebias)*exp(-i*mode*2.0*PI*((double) j / (double) Nx))
//                    + (0.1-modebias)*exp(i*2.0*mode*2.0*PI*((double) j / (double) Nx))
//                    + (0.1+modebias)*exp(-i*2.0*mode*2.0*PI*((double) j / (double) Nx)) ) // no offset
//                    + (1.0-modprop);
//        pnorm += pow(abs(waveform[j]),2.0);
//    }
//
////    cout << pnorm << endl;
//
//	if(abs(pnorm-1.0) > globtol) 
//	{
//        cout << "Normalising three-mode state with factor 1/" << sqrt(pnorm) << endl;
//		for(j=0;j<Nx;j++)
//		{
//			waveform[j] /= sqrt(pnorm);
//		}
//	}

return;
}

/*----------------------------------------------------------------------------*/

void addkick
// apply momentum kick to wavefunction
(
 vector< complex<double> > &waveform,
 vector<double> &x,
 double k0, 
 int Nx
)

{
    int j;
    const complex<double> i = {0.0, 1.0}; 

    for(j=0;j<Nx;j++)
    {
        waveform[j] *= exp(i*k0*x[j]); //non-quantised k
    }

return;
}

/*----------------------------------------------------------------------------*/

void addmodekick
// apply quantisedmomentum kick to wavefunction
(
 vector< complex<double> > &waveform,
 vector<double> &x,
 double x0,
 double k0, 
 double tmphase,
 int Nx
)

{
    int j;
    const complex<double> i = {0.0, 1.0}; 
    double pnorm = 0.0;

    for(j=0;j<Nx;j++)
    {
        waveform[j] *= (exp(i*2.0*PI*k0*x[j]) + exp(-i*2.0*PI*k0*x[j] - i*tmphase)); //quantised k
        pnorm += pow(abs(waveform[j]),2.0);
    }

	if(abs(pnorm-1.0) > globtol) 
	{
		for(j=0;j<Nx;j++)
		{
			waveform[j] /= sqrt(pnorm);
		}
	}

return;
}

/*----------------------------------------------------------------------------*/

double calcenergy
// calculate energy
(
vector< complex<double> > &waveform,
vector<double> &potential,
vector<double> &x,
double xstep,
double g,
double omega //rotating frame 
)

{
	int j;
	int Nx = x.size();
	complex<double> energy = {0.0,0.0};
    complex<double> i = {0.0,1.0};
    
    xstep = (2.0*PI)/(double) Nx;

	// energy sensitive to cyclic boundary condition
    energy += conj(waveform[0])*xstep*
              (
              -0.5*((waveform[1] - 2.0*waveform[0] + waveform[Nx-1])/(pow(xstep,2.0)))
                + potential[0]*waveform[0] 
                - (0.5*g*pow(abs(waveform[0]),2.0)*waveform[0])
                - i*omega*((x[Nx-1]-x[0])/(2.0*PI))*((waveform[1] - waveform[Nx-1])/(2.0*xstep))
              );

    energy += conj(waveform[Nx-1])*xstep*
              (
              -0.5*((waveform[0]-2.0*waveform[Nx-1]+waveform[Nx-2])/(pow(xstep,2.0)))
              + potential[Nx-1]*waveform[Nx-1] 
              - (0.5*g*pow(abs(waveform[Nx-1]),2.0)*waveform[Nx-1])
              - i*omega*((x[Nx-1]-x[0])/(2.0*PI))*((waveform[0] - waveform[Nx-2])/(2.0*xstep))
              );
    	// rest of the array
	for(j=1;j<(Nx-1);j++)
	{
    energy += conj(waveform[j])*xstep*
              (
              -0.5*((waveform[j+1]-2.0*waveform[j]+waveform[j-1])/(pow(xstep,2.0)))
              + potential[j]*waveform[j] 
              - (0.5*g*pow(abs(waveform[j]),2.0)*waveform[j]) // -ve as potential[] includes interactions
              - i*omega*((x[Nx-1]-x[0])/(2.0*PI))*((waveform[j+1] - waveform[j-1])/(2.0*xstep))
              );
    }

return abs(energy);
}

/*----------------------------------------------------------------------------*/

double calcenergy
// calculate energy in the non-rotating lab frame
(
vector< complex<double> > &waveform1,
vector< complex<double> > &waveform2,
vector<double> &potential,
vector<double> &x,
double xstep,
double gself,
double gint
)

{
    cout << "wrong calcenergy() called" << endl;
	int j;
	int Nx = x.size();
	double energy = 0.0;

	// energy sensitive to cyclic boundary condition
	energy += 0.5*abs((waveform1[1] - waveform1[Nx-1])/(2.0*xstep))
				+ potential[0]*pow(abs(waveform1[0]),2.0)
				+ 0.5*gself*pow(abs(waveform1[0]),4.0)
				+ 0.5*gint*pow(abs(waveform2[0]),4.0);
	energy += 0.5*abs((waveform1[0] - waveform1[Nx-2])/(2.0*xstep))
				+ potential[Nx-1]*pow(abs(waveform1[Nx-1]),2.0)
				+ 0.5*gself*pow(abs(waveform1[Nx-1]),4.0)
				+ 0.5*gint*pow(abs(waveform2[Nx-1]),4.0);

	for(j=1;j<(Nx-1);j++)
	{
		energy += 0.5*abs((waveform1[j+1] - waveform1[j-1])/(2.0*xstep))
					+ potential[j]*pow(abs(waveform1[j]),2.0)
					+ 0.5*gself*pow(abs(waveform1[j]),4.0)
					+ 0.5*gint*pow(abs(waveform2[j]),4.0);
	}

return energy;
}

/*----------------------------------------------------------------------------*/

double calclabenergy
// calculate energy in the lab frame
(
vector< complex<double> > &waveform,
vector<double> &potential,
vector<double> &x,
double xstep,
double g
)

{
    cout << "calling calclabenergy()" << endl;
	int j;
	int Nx = x.size();
	double energy = 0.0;

	// energy sensitive to cyclic boundary condition
	energy += 0.5*abs((waveform[1] - 2.0*waveform[0] + waveform[Nx-1])/(pow(xstep,2.0)))
				+ potential[0]*pow(abs(waveform[0]),2.0)
				+ g*pow(abs(waveform[0]),4.0);
	energy += 0.5*abs((waveform[0] - 2.0*waveform[Nx-1] + waveform[Nx-2])/(pow(xstep,2.0)))
				+ potential[Nx-1]*pow(abs(waveform[Nx-1]),2.0)
				+ g*pow(abs(waveform[Nx-1]),4.0);

	for(j=1;j<(Nx-1);j++)
	{
		energy += 0.5*abs((waveform[j+1] - 2.0*waveform[j] + waveform[j-1])/(pow(xstep,2.0)))
					+ potential[j]*pow(abs(waveform[j]),2.0)
					+ g*pow(abs(waveform[j]),4.0);
	}

return energy;
}

/*----------------------------------------------------------------------------*/

double calclabenergy
// calculate energy in the Husimi potential
(
vector< complex<double> > &waveform1,
vector< complex<double> > &waveform2,
vector<double> &potential,
vector<double> &x,
double xstep,
double gself,
double gint
)

{
	int j;
	int Nx = x.size();
	double energy = 0.0;

	// energy sensitive to cyclic boundary condition
	energy += 0.5*abs((waveform1[1] - waveform1[Nx-1])/(2.0*xstep))
				+ potential[0]*pow(abs(waveform1[0]),2.0)
				+ 0.5*gself*pow(abs(waveform1[0]),4.0)
				+ 0.5*gint*pow(abs(waveform2[0]),4.0);
	energy += 0.5*abs((waveform1[0] - waveform1[Nx-2])/(2.0*xstep))
				+ potential[Nx-1]*pow(abs(waveform1[Nx-1]),2.0)
				+ 0.5*gself*pow(abs(waveform1[Nx-1]),4.0)
				+ 0.5*gint*pow(abs(waveform2[Nx-1]),4.0);

	for(j=1;j<(Nx-1);j++)
	{
		energy += 0.5*abs((waveform1[j+1] - waveform1[j-1])/(2.0*xstep))
					+ potential[j]*pow(abs(waveform1[j]),2.0)
					+ 0.5*gself*pow(abs(waveform1[j]),4.0)
					+ 0.5*gint*pow(abs(waveform2[j]),4.0);
	}

return energy;
}

/*----------------------------------------------------------------------------*/

double calclz
//calculate expectation value of rotation operator Lz
(
vector< complex<double> > &waveform,
vector<double> &x,
double xstep
)

{
    int j;
    complex<double> lzexpect = {0.0,0.0};
    complex<double> i = {0.0,1.0};
    double Nx = x.size();

    //terms at boundary
    lzexpect += conj(waveform[0]) * xstep * (-i*((x[Nx-1]-x[0])/(4.0*PI*xstep))*(waveform[1]-waveform[Nx-1]));
    lzexpect += conj(waveform[Nx-1]) * xstep * (-i*((x[Nx-1]-x[0])/(4.0*PI*xstep))*(waveform[0]-waveform[Nx-2]));
    
    for(j=1;j<(Nx-1);j++)
    {
        lzexpect += conj(waveform[j]) * xstep * (-i*((x[Nx-1]-x[0])/(4.0*PI*xstep))*(waveform[j+1]-waveform[j-1]));
    }

    if(abs(lzexpect) < 1e-12)
    {
        lzexpect = {0.0,0.0};
    }

return abs(lzexpect);
}

/*----------------------------------------------------------------------------*/

void calcstats
// calculate generic wavefunction statistics and place in statlist vector
(
vector< complex<double> > &waveform,
vector< complex<double> > &oldwaveform,
vector<double> &potential,
vector<double> &statlist,
vector<double> &x,
double xstep,
int t,
double timestep,
int renorm,
double omegaou
)

{

	int j;
	int Nx = (int) x.size();
	double norm = 0.0;
    double oldnorm = 0.0;
	double checknorm = 0.0;
	double avgx = 0.0;	  // <x>
	double avgxsq = 0.0;  // <x^2> for variance/stdev
	complex<double> overlap = {0.0,0.0}; // overlap of probability distributions
    double potexpect = 0.0; // expectation value of potential
    double lzexpect = 0.0; // expectation value of Lz
    complex<double> i = {0.0,1.0};
    double beta = sqrt(hbar/(rbmass*omegaou));

	for(j=0;j<Nx;j++)
	{
	norm += xstep * pow(abs(waveform[j]),2.0);
    oldnorm += xstep * pow(abs(oldwaveform[j]),2.0); //inefficient
	}
	
	if(renorm == 1) //renormalise if needed
	{
		for(j=0;j<Nx;j++)
		{
			waveform[j] *= (1.0/sqrt(norm));
			checknorm += pow(abs(waveform[j]),2.0)*xstep;
		}
		norm = checknorm;
	}

	for(j=0;j<Nx;j++)
	{
        avgx += x[j] * xstep * pow(abs(waveform[j]),2.0);
        avgxsq += pow(x[j],2.0) * xstep * pow(abs(waveform[j]),2.0);
        overlap += conj(oldwaveform[j])*waveform[j]*xstep;
        potexpect += xstep * potential[j];
    }

    lzexpect = calclz(waveform,x,xstep);

/*
    if(t==0)
    {
        cout << "Lz(t=0): " << lzexpect << endl;
    }
*/

    //output stats in arbitrary units
	statlist[0] = t;
	statlist[1] = timestep*t; 
	statlist[2] = norm;
	statlist[3] = avgx;
    statlist[4] = avgxsq;
	statlist[5] = sqrt(avgxsq - pow(avgx,2.0));
	statlist[6] = abs(overlap)/sqrt(oldnorm);
    statlist[7] = potexpect;
    statlist[8] = lzexpect; 

return;
}

/*----------------------------------------------------------------------------*/

void calcstats
// calculate generic wavefunction statistics and place in statlist vector
(
vector< complex<double> > &waveform,
vector< complex<double> > &oldwaveform,
vector<double> &potential,
vector<double> &statlist,
vector<double> &x,
double xstep,
int t,
double timestep,
int renorm,
int modulus,
double omegaou
)

{
/*
	int j;
	int Nx = (int) x.size();
	double norm = 0.0;
	double checknorm = 0.0;
	double avgx = 0.0;	  // <x>
	double avgxsq = 0.0;  // <x^2> for variance/stdev
	double overlap = 0.0; // overlap of probability distributions
    double potexpect = 0.0; // expectation value of potential
    double beta = sqrt(hbar/(rbmass*omegaou));
*/
    
    
	int j;
	int Nx = (int) x.size();
	double norm = 0.0;
    double oldnorm = 0.0;
	double checknorm = 0.0;
	double avgx = 0.0;	  // <x>
	double avgxsq = 0.0;  // <x^2> for variance/stdev
	complex<double> overlap = {0.0,0.0}; // overlap of probability distributions
    double potexpect = 0.0; // expectation value of potential
    double lzexpect = 0.0; // expectation value of Lz
    complex<double> i = {0.0,1.0};
    double beta = sqrt(hbar/(rbmass*omegaou));
    
    
	for(j=0;j<Nx;j++)
	{
	norm += xstep * pow(abs(waveform[j]),2.0);
    oldnorm += xstep * pow(abs(oldwaveform[j]),2.0); //inefficient
	}
	
	if(renorm == 1) //renormalise if needed
	{
		for(j=0;j<Nx;j++)
		{
			waveform[j] *= (1.0/sqrt(norm));
			checknorm += pow(abs(waveform[j]),2.0)*xstep;
		}
		norm = checknorm;
	
        /*
        if(norm != 1.0)
        {
            cout << "Renormalisation failed, re-renormalising" << endl;
            for(j=0;j<Nx;j++)
            {
                waveform[j] *= (1.0/sqrt(norm));
            }
        }
        */
	}

	for(j=0;j<Nx;j++)
	{
        avgx += x[j] * xstep * pow(abs(waveform[j]),2.0);
        avgxsq += pow(x[j],2.0) * xstep * pow(abs(waveform[j]),2.0);
        overlap += conj(oldwaveform[j])*waveform[j]*xstep;
        potexpect += xstep * potential[j];
    }

    lzexpect = calclz(waveform,x,xstep);
/*
	for(j=0;j<Nx;j++)
	{
	norm += xstep * pow(abs(waveform[j]),2.0);
	avgx += x[j] * pow(abs(waveform[j]),2.0) * xstep;
	avgxsq += pow(x[j],2.0) * pow(abs(waveform[j]),2.0) * xstep;
	overlap += fmin(pow(abs(waveform[j]),2.0),
					pow(abs(oldwaveform[j]),2.0))*xstep;
	}
	// account for the gap between x[Nx-1] and x[0]
	norm += xstep * pow(abs(waveform[0]),2.0);
	avgx += x[0] * pow(abs(waveform[0]),2.0) * xstep;
	avgxsq += pow(x[0],2.0) * pow(abs(waveform[0]),2.0) * xstep;
	overlap += fmin(pow(abs(waveform[0]),2.0),
					pow(abs(oldwaveform[0]),2.0))*xstep;
    potexpect += xstep * potential[j];

	if((norm != 1.0) && (renorm == 1)) //renormalise if needed
	{

		for(j=0;j<Nx;j++)
		{
			waveform[j] *= (1.0/sqrt(norm));
			checknorm += pow(abs(waveform[j]),2.0)*xstep;
		}
		norm = checknorm;
	}
*/

	statlist[0] = t-t/modulus;
	statlist[1] = timestep*(t-t/modulus);
	statlist[2] = norm;
	statlist[3] = avgx;
    statlist[4] = avgxsq;
	statlist[5] = sqrt(avgxsq - pow(avgx,2.0));
	statlist[6] = abs(overlap);
    statlist[7] = potexpect;
    

return;
}

/*----------------------------------------------------------------------------*/

void aggregatestats
// combine stats for both components into a single vector for printing
(
vector<double> &statvec1,
vector<double> &statvec2,
vector<double> &aggstatvec
)

{
    int j;
    
    for(j=0;j<(int) aggstatvec.size();j++)
    {
        if(j<(int) statvec1.size())
        {
            aggstatvec[j] = statvec1[j];
        }
        else if(j >= (int) statvec1.size())
        {
            aggstatvec[j] = statvec2[j-(statvec1.size()-2)];
        }
    }

return;
}

/*----------------------------------------------------------------------------*/

void statfileprint
// print data contained in statlist vector to file
/* If headers are required they must be printed in advance, this function is
designed to print one line  to file eg. to write summary stats at each time
step.
Stats are stored in the statlist array and are printed in index order
*/
(
ofstream &outfile,
vector<double> &statlist
)

{
	int j;
	int ncolumns = statlist.size();
	
	//write stats to file
	for(j=0;j<(ncolumns-1);j++)
	{
		outfile << showpoint << statlist[j]  << setw(6)  << "\t";
	}
	//print last value without trailing tab, end line & flush stream
		outfile << statlist[ncolumns-1] << endl;

return;
}

/*----------------------------------------------------------------------------*/

void wffileprint
// print whole vectors to file as columns
/* veclist is a vector of pointers to the vectors you want to print to file.
This function prints the vectors as columns in the order they appear in veclist.
You can have as many columns as needed, this will iterate over the entirety of
veclist printing values.
Note: all vectors should be the same length as specified by vecsize parameter.

Example for creating a veclist:
vector<vector<double> * > veclist;
veclist = {&vec1, &vec2, ...};
*/
(
ofstream outfile,
vector<vector<double> * > veclist,
int vecsize
)

{
	vector<double> vecselect;
	int j, k;	
	int nvecs = veclist.size();

	for(j=0;j<vecsize;j++)
	{
		for(k=0;k<(nvecs-1);k++)
		{
			vecselect = *veclist[k]; // select the desired vector
			outfile << vecselect[j] << "\t"; //print contents of vector
		}
		//print last column without trailing tab, end the line
		vecselect = *veclist[nvecs-1];
		outfile << vecselect[j] << endl;
	}
return;
}

/*----------------------------------------------------------------------------*/

void calcrrhs
// Calculate right hand side vector for Sherman-Morrison method
(
vector< complex<double> > &zrrhs,
vector< complex<double> > &zoldpsi,
vector<double> &potential,
complex<double> zalpha,
complex<double> zbeta,
complex<double> zgamma,
int Nx
)

{
	int j;

	if(zrrhs.size() != zoldpsi.size())
	{
		cout << "calcrrhshusimi() error: arrays must be same size!\n" << endl;
		exit(EXIT_FAILURE);
	}

	zrrhs[0] = (1.0-(zbeta*potential[0])-2.0*zalpha)*zoldpsi[0]
			   + ((zalpha+zgamma))*zoldpsi[1] + ((zalpha-zgamma))*zoldpsi[Nx-1];

    //pragma omp parallel for
	for(j=1;j<(Nx-1);j++)
	{
		zrrhs[j] = (1.0-(zbeta*potential[j])-2.0*zalpha)*zoldpsi[j]
			   + ((zalpha+zgamma))*zoldpsi[j+1] + ((zalpha-zgamma))*zoldpsi[j-1];
	}

	zrrhs[Nx-1] = (1.0-(zbeta*potential[Nx-1])-2.0*zalpha)*zoldpsi[Nx-1]
		  	  + ((zalpha+zgamma))*zoldpsi[0] + ((zalpha-zgamma))*zoldpsi[Nx-2];

return;
}

/*----------------------------------------------------------------------------*/

void updatediagonal
// update the values of Sherman-Morrison matrix diagonal terms for
// time-dependent potential
(
vector< complex<double> > &zadiag,
vector<double> &potential,
complex<double> &zalpha,
complex<double> &zbeta,
int Nx
)

{
	int j;

	if(zadiag.size() != potential.size())
	{
		cout << "husimiupdatediagonal() error: vectors must be same size!\n"
			 << endl;
		exit(EXIT_FAILURE);
	}


    //pragma omp parallel for
	for(j=0;j<Nx;j++){
			zadiag[j] = 1.0 + zbeta*potential[j] + 2.0*zalpha;
		}

		//re-apply cyclic corrections to matrix diagonal
		zadiag[Nx-1] -= ((zalpha*zalpha)/(-zadiag[0]));
		zadiag[0] += zadiag[0]; // b_0 - (-b_0)
}

/*----------------------------------------------------------------------------*/

void twocompfileprint
// print two-component waveforms to file
(
ofstream &outfile,
string filename,
vector<double> &xgrid,
vector< complex<double> > &waveform1,
vector< complex<double> > &waveform2,
vector<double> &potential1,
vector<double> &potential2
)

{
	unsigned int j;

	outfile.open(filename.c_str());
	outfile << "#x\t\t\t\t|psi1|^2\t\t\tRe(psi1)\t\t\tIm(psi1)\t\t\tV1(x,t)\t\t\t"
				"phase(psi)\t\t\t|psi2|^2\t\t\tRe(psi2)\t\t\tIm(psi2)\t\t\tV2(x,t)"
				"\t\t\tfreq\t\t\txrotating\n" 
			<< endl;

	for(j=0;j<xgrid.size();j++)
	{
		outfile << setw(8) 
                << xgrid.at(j) // 1
                << "\t\t" << pow(abs(waveform1[j]),2.0) // 2
			 	<< "\t\t" << real(waveform1[j])  // 3
                << "\t\t" << imag(waveform1[j])  // 4
                << "\t\t" << potential1[j] // 5
                << "\t\t" << atan2(imag(waveform1[j]),real(waveform1[j]))/PI // 6
				<< "\t\t" << pow(abs(waveform2[j]),2.0)  // 7
                << "\t\t" << real(waveform2[j])   // 8
                << "\t\t" << imag(waveform2[j]) // 9
				<< "\t\t" << potential2[j]   // 10
				<< "\n";
	}

    outfile.close();

return;
}

/*----------------------------------------------------------------------------*/

/*
void updatepotential
//re-calculates the potential function potentialfunc at time t
(
vector<double> &x,
vector<double> &potentialvec,
vector< complex<double> > &psi1, //one component
int impot,
int imm,
double fweight,
int rt,
double gii, 
double omega1, 
double omega2, 
double t, 
double potamp,
double sigma,
double gdepth,
int husimitype
)

{
	int j, Nx = x.size();

// This can probably be more efficient than it currently is
    
// first zero the potential to avoid issues
    for(j=0;j<Nx;j++)
    {
        potentialvec[j] = 0.0;
    }

// now construct the potential with each desired element

    for(j=0;j<Nx;j++)
    {
        // self-interaction (disabled by setting gii=0
        potentialvec[j] += calcpotential(psi1[j],gii);
        
        // Husimi potential (disabled by setting omega1=omega2=0
        switch(husimitype) 
        {
            case 2: //Gaussian trap
cout << "Gaussian" << endl;                
                potentialvec[j] += -gdepth*(exp(-pow(x[j]-(0.5*(x[Nx-1]-x[0])),2.0)/(2.0*sigma))) - potamp*(x[j]-(0.5*(x[Nx-1]-x[0])))*sin(omega2*t);
                // 
            break;

            default: //standard harmonic trap
cout << "Harmonic" << endl;                
                potentialvec[j] += 0.5*pow(omega1,2.0)*pow((x[j]-(0.5*(x[Nx-1]-x[0]))),2.0) - potamp*(x[j]-(0.5*(x[Nx-1]-x[0])))*sin(omega2*t);
                break;
        }
        
    }
     
    // imprinting done outside loop due to structure of calcVimprint()
    if((rt == 0) && (impot != 0) && (omega1 == 0.0)) //imprinting potential in imag time
    {
        calcVimprint(x,potentialvec,psi1,impot,imm,fweight,gii,potamp);
    }

}
*/

/*----------------------------------------------------------------------------*/

void updatepotential
//re-calculates the potential function potentialfunc at time t; two component
(
vector<double> &x,
vector<double> &potentialvec,
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
double printoffset
)

{
	int j, Nx = x.size();
    
    // first zero the potential
    for(j=0;j<Nx;j++)
    {
        potentialvec[j] = 0.0;
    }

    // now reconstruct the potential based on passed options

    // add nonlinearity
    addnonlinearity(potentialvec,psi1,psi2,gii,gij);

    // add Huismi potential
    if(husimitype != 0)
    {
        addhusimipotential(potentialvec,x,t,potamp,omega1,omega2,gdepth,sigma,husimitype);
    }

    //add imprinting potential in imag time
    if(((rt == 0) && (imprinttype != 1)) || ((rt == 1) && (imprinttype == 2))) 
    {
        calcVimprint(x,potentialvec,psi1,psi2,impot,imm,fweight,gii,gij,imprintamp,printoffset);
    }

    //add real time imprinting potential
    if((rt == 1) && (imprinttype == 1))  
    {
        addtdpotential(potentialvec,x,psi1,psi2,t,sigma,impot,imm,fweight,gii,gij,imprintamp,imptimeprofile,impoffset,impsigma,omegaou,printoffset);
    }

return;
}

/*----------------------------------------------------------------------------*/

void addtdpotential
// calculates time-dependent imprinting potential and adds it to potentialvec
(
vector<double> &potentialvec,
vector<double> &x,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2,
double t,
double sigma,
int impot,
int imm,
double fweight,
double gii,
double gij,
double imprintamp,
int imptimeprofile,
double impoffset,
double impsigma,
double omegaou,
double printoffset
)

{
    int j, Nx = x.size();
    double rampfactor; // has Gaussian profile in time, with constant total integral
    double newamp; //potential amplitude scaled by time-dependent Gaussian profile


    //put a switch here to change rampfactor for each profile type
    switch(imptimeprofile)
    {
        case 0:
            rampfactor = (1.0/impsigma)*exp(-(pow((t-impoffset),2.0))/(2.0*pow(impsigma,2.0)));
            break;
    }

    newamp = rampfactor*imprintamp;

    calcVimprint(x,potentialvec,psi1,psi2,impot,imm,fweight,gii,gij,newamp,printoffset);

return;
}

/*----------------------------------------------------------------------------*/

void addhusimipotential
// calculates Husimi potential term and adds it to potential vector
(
vector<double> &potential,
vector<double> &x,
double t,
double potamp,
double omega1,
double omega2,
double gdepth,
double sigma,
int husimitype
)

{
    int j, Nx = x.size();
    double dx = x[1] - x[0];

    switch(husimitype)
    {
        case 1: // harmonic trap (standard Husimi potential)
            for(j=0;j<Nx;j++)
            {
            potential[j] += 0.5*pow(omega1,2.0)*pow((x[j]-(x[0]+0.5*((x[Nx-1]+dx)-x[0]))),2.0) - potamp*(x[j]-(x[0]+0.5*((x[Nx-1]+dx)-x[0])))*sin(omega2*t);
            }
            break;
        case 2: // Gaussian trap
            for(j=0;j<Nx;j++)
            {
            potential[j] += -gdepth*(exp(-pow(x[j]-(x[0] + (0.5*(x[Nx-1]-x[0]))),2.0)/(2.0*sigma*sigma))) - potamp*(x[j]-(0.5*(x[Nx-1]-x[0])))*sin(omega2*t);
            }
            break;
    }


return;
}

/*----------------------------------------------------------------------------*/

void addnonlinearity //one component, calculates interaction potential
(
vector<double> &potential,
vector< complex<double> > &psi1,
double gii 
)

{
    unsigned int j;
    double result;

    for(j=0;j<(potential.size());j++)
    {
        potential[j] += gii*pow(abs(psi1[j]),2.0);
    }

}

/*----------------------------------------------------------------------------*/

void addnonlinearity //two components, calculates interaction potential
(
vector<double> &potential,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2,
double gii,
double gij
)

{
    unsigned int j;

    for(j=0;j<(potential.size());j++)
    {
        potential[j] += gii*pow(abs(psi1[j]),2.0) //self-interaction nonlinearity
					    + gij*pow(abs(psi2[j]),2.0); //two-component interactions
    }

return;
}

/*----------------------------------------------------------------------------*/

void calcVimprint //calculate imprinting potential and add to total potential
(
vector<double> &x,
vector<double> &potentialvec,
vector< complex<double> > &psi1, //one component
int impot,
int imm,
double fweight,
double gii,
double imprintamp,
double printoffset
)

{
    int j;
    int Nx = (int) x.size();
    double dx = x[1] - x[0];
    vector<double> normvec;
    double sum, polymod; //for numerically normalising the polynomial
    double xmod;
    double lwidth = x[0] + (0.5*(x[Nx-1]-x[0])/imm);
    double sinmod = (2.0 * lwidth)/(lwidth + (1.0/(imm*PI))*x[Nx-1]*pow(-cos(PI*lwidth*imm/x[Nx-1]),2.0));

    normvec.resize(Nx);

/*
    if(gii == 0.0)
    {
            printoffset = 0.0;
    }
    else
    {
            printoffset = 1.0;
    }
*/

    switch (impot) //which type of imprinting do you want?
    {
        case 1:
            //heaviside step-function propeller
            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += 0.5*imprintamp*(sgn(-cos((imm*2.0*PI)*((double) j/(double) Nx))) + printoffset);
            }
            break;
            
        case 2:
            //sinusoid
            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += 0.5*imprintamp*sinmod*(-cos(imm*2.0*PI*((double) j/(double) Nx)) + printoffset);
            }
            break;
            
        case 3:
            //Gaussian comb
            for (j=0;j<Nx;j++)
            {
                xmod = fmod(imm*(x[j]-x[0])/(x[Nx-1]-x[0]),1.0);
                potentialvec[j] += imprintamp*exp(-(pow(xmod - 0.5,2.0))/(2*pow(0.1,2.0)));
            }
            break;
            
        case 4:
            //Polynomial approximation to sine/square
            // need to normalise this numerically
            for(j=0;j<Nx;j++)
            {
                normvec[j] = 0.5*imprintamp*(sgn(-cos(2.0*PI*imm*((x[j]-x[0])/(x[Nx-1]-x[0]))))
                                    *(1.0-(pow(fmod((4.0*x[j]),2.0)-1.0,2.0*fweight)))+1.0);
                if(x[j] <= lwidth)
                {
                    sum += normvec[j]*dx;
                }
            }

            if(!((sum - (imprintamp*lwidth)) < 0.000001)) //should switch to global tolerance
            {
                polymod = (imprintamp*lwidth)/sum; //normalisation factor
            }

            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += 0.5*imprintamp*polymod*(sgn(-cos(2.0*PI*imm*((x[j]-x[0])/(x[Nx-1]-x[0]))))
                                    *(1.0-(pow(fmod((4.0*x[j]),2.0)-1.0,2.0*fweight)))+1.0);
            }
            break;
    }
    
}

/*----------------------------------------------------------------------------*/

void calcVimprint //calculate imprinting potential
(
vector<double> &x,
vector<double> &potentialvec,
vector< complex<double> > &psi1,
vector< complex<double> > &psi2, //two-component
int impot,
int imm,
double fweight,
double gii, 
double gij,
double potamp,
double printoffset
)

{
    int j;
    unsigned int Nx = x.size();
    double dx = x[1] - x[0];
    vector<double> normvec(Nx);
    double sum, polymod; //for numerically normalising the polynomial
    double xmod;
    double lwidth = x[0] + (0.5*(x[Nx-1]-x[0])/(double) imm);
    double sinmod = (2.0 * lwidth)/(lwidth + (1.0/(imm*PI))*x[Nx-1]*pow(-cos(PI*lwidth*imm/x[Nx-1]),2.0));
    
    double cx;
    double x0 = potamp; // temporary

//    normvec.resize(Nx);


    // printoffset centres the potential function around V=0 for
    // clearer plots
/*
    if(gii == 0.0)
    {
            printoffset = 0.0;
    }
    else
    {
            printoffset = 1.0;
    }
*/

    switch (impot) //which type of imprinting do you want?
    {
        case 1:
            //heaviside step-function propeller
            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += 0.5*potamp*(sgn(-cos((imm*2.0*PI)*((double) j/(double) Nx))) + printoffset);
            }
            break;
            
        case 2:
            //sinusoid: normalised for comparison with other imprinting functions
            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += 0.5*potamp*sinmod*(-cos((imm*2.0*PI)*((double) j/(double) Nx)) + printoffset);
            }
            break;
            
        case 3:
            //Gaussian comb
            for (j=0;j<Nx;j++)
            {
                xmod = fmod(imm*(x[j]-x[0])/(x[Nx-1]-x[0]),1.0);
                potentialvec[j] += potamp*exp(-(pow(xmod - 0.5,2.0))/(2*pow(0.1,2.0)));
            }
            break;
            
        case 4:
            //Polynomial approximation to sine/square

            for(j=0;j<Nx;j++)
            {
                normvec[j] = 0.5*potamp*(sgn(-cos(2.0*PI*imm*((double) j/(double) Nx)))
                                    *(1.0-(pow(fmod((4.0*x[j]),2.0)-1.0,2.0*fweight)))+1.0);
                if(x[j] <= lwidth)
                {
                    sum += normvec[j]*dx;
                }
            }

            // numerical normalisation
            if(!((sum - (potamp*lwidth)) < 0.000001))
            {
                polymod = (potamp*lwidth)/sum; //normalisation factor
            }

            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += 0.5*potamp*polymod*(sgn(-cos(2.0*PI*imm*((x[j]-x[0])/(x[Nx-1]-x[0]))))
                                    *(1.0-(pow(fmod((4.0*x[j]),2.0)-1.0,2.0*fweight)))+1.0);
            }
            break;

        case 5:
            // Harmonic trap

            for(j=0;j<Nx;j++)
            {
    		    cx = fmin(x[j]-x0, (x[Nx-1]-x[0])-(x[j]-x0)); // centred x coordinate
                potentialvec[j] += 0.5*pow((double) imm,2.0)*pow(cx,2.0);
            }
            break;
            
        case 6:
            //alternate sinusoid: 1-alpha cos(theta) , for analytics
            for (j=0;j<Nx;j++)
            {
                potentialvec[j] += (1.0 - potamp*cos(imm*((2.0*PI)*((double) j/(double) Nx))));
            }
//            cout << "imprintamp: " << potamp << endl;
//            cout << "imm       : " << imm << endl;
            break;
    }
    
}

/*----------------------------------------------------------------------------*/

double heaviside
(
double x
)

{
    if(x<0.0)
    {
        return 0.0;
    }
    else
    {
        return 1.0;
    }
}

/*----------------------------------------------------------------------------*/

double sgn
// implementation of signum function for square wave
(
double x
)

{
    return ((x > 0.0) - (x < 0.0));
}

/*----------------------------------------------------------------------------*/

void prepSMvecs
// prepare A matrix tridiagonals, u and v vectors for Sherman-Morrison
(
vector< complex<double> > &zadiag, 
vector< complex<double> > &zaupper, 
vector< complex<double> > &zalower,
vector< complex<double> > &zuvec,
vector< complex<double> > &zvvec,
vector<double> &potential,
complex<double> alpha,
complex<double> beta,
complex<double> gamma
)

{
	int j;
	int Nx = (int) zadiag.size();

    //pragma omp parallel for
	for(j=0;j<Nx;j++)
	{
		//diagonal
		zadiag[j]  = 1.0 + beta*potential[j] + 2.0*alpha;

		//upper & lower off-diagonals
	    zaupper[j]  = (-alpha - gamma);
		zalower[j]  = (-alpha + gamma);

		//fill u,v with zeroes, then replace first & last terms outside loop
		zuvec[j]  = {0.0,0.0};
		zvvec[j]  = {0.0,0.0};

	}


	zuvec[0] = -zadiag[0]; // u[0] = -b_0
	zuvec[Nx-1] = (-alpha - gamma); // bottom-left corner term

	zvvec[0] = {1.0, 0.0};
	zvvec[Nx-1] = ((-alpha + gamma))/(-zadiag[0]); 
		// NR(beta/gamma), where beta is upper-right correction term

	//corrections to matrix, NR eq 2.7.11
	zadiag[Nx-1] -= (-alpha + gamma)*((-alpha - gamma)/(-zadiag[0]));
	zadiag[0] += zadiag[0]; // b_0 - (-b_0)


    //debug
    //cout << "debug start" << endl;
    //cout << "===========" << endl;
    //cout << "alpha = " << alpha << endl;
    //cout << "beta = " << beta << endl;
    //cout << "gamma = " << gamma << endl;
    //cout << " " << endl;
    //cout << "zadiag[j] = " << zadiag[3] << endl;
    //cout << "zaupper[j] = " << zaupper[3] << endl;
    //cout << "zalower[j] = " << zalower[3] << endl;
    //cout << "zuvec[0] = " << zuvec[0] << endl;
    //cout << "zuvec[Nx-1] = " << zuvec[Nx-1] << endl;
    //cout << "zvvec[0] = " << zvvec[0] << endl;
    //cout << "zvvec[Nx-1] = " << zvvec[Nx-1] << endl;
    //cout << "zadiag[Nx-1] = " << zadiag[Nx-1] << endl;
    //cout << "zadiag[0] = " << zadiag[0] << endl;
    //cout << "debug end" << endl;
    //cout << "===========" << endl;

}

/*----------------------------------------------------------------------------*/

void prepSMvecs
// prepare A matrix tridiagonals, u and v vectors for Sherman-Morrison
// overloaded version without upper, lower vecs (use this for BEC component 2)
(
vector< complex<double> > &zadiag, 
vector< complex<double> > &zuvec,
vector< complex<double> > &zvvec,
vector<double> &potential,
complex<double> alpha,
complex<double> beta,
complex<double> gamma
)

{
	int j;
	int Nx = (int) zadiag.size();

    //pragma omp parallel for
	for(j=0;j<Nx;j++)
	{
		//diagonal
		zadiag[j]  = 1.0 + beta*potential[j] + 2.0*alpha;


		//fill u,v with zeroes, then replace first & last terms outside loop
		zuvec[j]  = {0.0,0.0};
		zvvec[j]  = {0.0,0.0};
	}

	zuvec[0] = -zadiag[0]; // u[0] = -b_0
	zuvec[Nx-1] = (-alpha - gamma); // bottom-left corner term

	zvvec[0] = {1.0, 0.0};
	zvvec[Nx-1] = ((-alpha + gamma))/(-zadiag[0]); 
		// NR(beta/gamma), where beta is upper-right correction term

	//corrections to matrix, NR eq 2.7.11
	zadiag[Nx-1] -= (-alpha + gamma)*((-alpha - gamma)/(-zadiag[0]));
	zadiag[0] += zadiag[0]; // b_0 - (-b_0)

}

/*----------------------------------------------------------------------------*/

void updatediag
// update potential-dependent matrix diagonal
// deprecated due to u,v vectors; which are dependent on diagonal and therefore
// the potential
(
vector< complex<double> > &zadiag,
vector<double> &potential,
complex<double> zalpha,
complex<double> zbeta
)

{
	int j;
	int Nx = (int) zadiag.size();

    //pragma omp parallel for
	for(j=0;j<Nx;j++)
	{
		zadiag[j] = 1.0 + zbeta*potential[j] + 2.0*zalpha;
	}
	zadiag[Nx-1] -= ((zalpha*zalpha)/(-zadiag[0]));
	zadiag[0] += zadiag[0]; // b_0 - (-b_0)

}

/*----------------------------------------------------------------------------*/

complex<double> fixedmodefourierint
// perform Fourier transform of wavefunction at fixed wavenumber by integration
(
vector< complex<double> > &wavefunction,
vector<double> &x,
double m
)

{
    const complex<double> i = {0.0,1.0};
    int j;
    int Nx = (int) x.size();
    complex<double> sum = {0.0,0.0};
    double xstep = x[1]-x[0];
    double L = (x[Nx-1] + xstep) - x[0];
    double angularstep = (2.0*PI)/(Nx);

    for(j=0;j<Nx;j++)
    {
        /*
            sum += ((1.0)/(2.0*PI))*sqrt(L)* // sqrt(L) factor normalises
                    wavefunction[j]*angularstep*
                    exp(-i*m*2.0*PI*(((double) j)/((double) Nx)));
        */
        
            sum += (1.0/(2.0*PI))*sqrt(L)* 
                wavefunction[j]*angularstep*
                exp(i*m*2.0*PI*((double) j/(double) Nx));
    }
return sum;
}

/*----------------------------------------------------------------------------*/

complex<double> fixedmodedensityfourierint
// perform Fourier transform of density at fixed wavenumber by integration
(
vector< complex<double> > &wavefunction,
vector<double> &x,
double m //mode number
)

{
    const complex<double> i = {0.0,1.0};
    int j;
    int Nx = (int) x.size();
    complex<double> sum = {0.0,0.0};
    double xstep = x[1]-x[0];
    double L = x[Nx-1] - x[0];
    double angularstep = (2.0*PI)/(Nx);

    for(j=0;j<Nx;j++)
    {
        /*
            sum += (1.0/(2.0*PI))*sqrt(L)* 
                conj(wavefunction[j])*wavefunction[j]*angularstep*
                exp(-i*m*2.0*PI*((x[j]-x[0])/(x[Nx-1]-x[0])));
        */
        
            sum += (1.0/(2.0*PI))*sqrt(L)* 
                conj(wavefunction[j])*wavefunction[j]*angularstep*
                exp(i*m*2.0*PI*((double) j/(double) Nx ));
//        cout << (double) j/(double) (Nx-1) << endl;
    }
//    cout << endl;
    
return sum;
}

/*----------------------------------------------------------------------------*/

complex<double> fixedkfourierint
// placeholder: Fourier transform at fixed wavenumber
(
vector< complex<double> > &wavefunction,
vector<double> &x,
double k
)

{
    const complex<double> i = {0.0,1.0};
    int j;
    int Nx = (int) x.size();
    complex<double> sum = {0.0,0.0};
    double xstep = x[1]-x[0];

    //pragma omp parallel
    for(j=0;j<Nx;j++)
    {
            // x runs from 0 to 2pi in the exponent
            sum += wavefunction[j]*xstep*exp((-i*k*((x[j]-x[0])/(x[Nx-1]))));
    }
return sum;
}


/*----------------------------------------------------------------------------*/

complex<double> fixedkdensityfourierint
// placeholder: Fourier transform density at fixed wavenumber
(
vector< complex<double> > &wavefunction,
vector<double> &x,
double k
)

{
    const complex<double> i = {0.0,1.0};
    int j;
    int Nx = (int) x.size();
    complex<double> sum = {0.0,0.0};
    double xstep = x[1]-x[0];

    //pragma omp parallel
    for(j=0;j<Nx;j++)
    {
            // x runs from 0 to 2pi in the exponent
            sum += pow(abs(wavefunction[j]),2.0)*xstep*exp((-i*k*((x[j]-x[0])/(x[Nx-1]))));
    }
return sum;
}
/*----------------------------------------------------------------------------*/

double complexphase
// calculate the phase angle of a complex quantity using datan
// mostly exists to avoid calling fourierint funcs twice on the same data
(
complex<double> cnumber
)

{
    return datan(real(cnumber),imag(cnumber));
}

/*----------------------------------------------------------------------------*/


double datan
// atan implementation based on Bromley code datannew()
(
double x,
double y
)

{
    double temp;

    if(abs(x) < globtol)
    {
        if(abs(y) < globtol)
        {
            return 0.0;
        }
        else
        {
            if(y > globtol) { return 0.5*PI; }
            if(y < globtol) { return 1.5*PI; }
        }
    }
    else
    {
        temp = atan(y/x);
        if((y > globtol) && (x < globtol)) { return temp + PI; }
        if((y < globtol) && (x < globtol)) { return temp + PI; }
        if((y < globtol) && (x > globtol)) { return temp + 2.0*PI; }
    }

return temp;
}

/*----------------------------------------------------------------------------*/

double dvecavg
// returns average value of given vector<double>
(
vector<double> &vec
)

{
    unsigned int j;
    double sum;

    for(j=0;j<vec.size();j++)
    {
        sum += vec[j];
    }

return sum/((double) vec.size());
}

/*----------------------------------------------------------------------------*/

void col2vec
// copies the contents of nth matrix column into empty vector
(
vector<double> &source,
vector<double> &dest,
int ncolumns, //number of columns in source
int column,   //column index (zero-based) to copy
int header    //specifies whether data has headings (ie from file)
)

{
    int j=0;
    int nrows = source.size()/(ncolumns);

    //sanitize header
    if((header != 0) && (header != 1))
    {
        cout << "ERROR: incorrect value for 'header' (must be 0 or 1)" << endl;
        exit(EXIT_FAILURE);
    }

    if(header == 1)
    {
        nrows = (source.size() - ncolumns)/ncolumns; //gives nrows w/o header
    }

    dest.resize(nrows);

    for(j=header;j<nrows;j++) //start at 1 because first line is assumed headers
    {
        dest[j] = source[(j*ncolumns)+column];
    }

return;
}

/*----------------------------------------------------------------------------*/

void vecmax
// finds and returns indices of local maxima 
(
vector<double> &vec,
vector<int> &maxima
)

{
    int j;
    int maxindex = 0;
    int maxcount = 0;

    // test if the starting point is a max
    if(vec[1] < vec[0])
    {
        maxima.push_back(0);
    }

    for(j=1;j<vec.size()-1;j++)
    {
        if((vec[j-1] < vec[j]) && (vec[j+1] < vec[j])) //if local maximum
        {
            maxima.push_back(j);
        }
    }

    //don't test the endpoint due to unknown end time

return;
}

/*----------------------------------------------------------------------------*/

string makeheader
// creates file column headers from list of frequencies
(
vector<double> freqvec
)

{
    int j;
    stringstream output;

    output << "time\t";

    for(j=0;j<freqvec.size();j++)
    {
        output << setw(8)<< freqvec[j] << "\t";
    }


return output.str();
}

/*----------------------------------------------------------------------------*/

void printconversions
// prints unit conversions to stderr log file
(
double ringrad,
double ztrapfreq,
double rtrapfreq,
double g11,
double g11conv,
int ouselect,
double omegaou,
double xmax,
double xmin,
double tmax
)

{
    double length, time, energy;


    // do oscillator unit conversions
    length = sqrt(hbar/(rbmass*omegaou)); // in SI
    time   = 1.0/omegaou; // in SI
    energy = hbar*omegaou; // in SI


    cerr << "Radial trap frequency: " << rtrapfreq << "Hz" << endl;
    cerr << "z trap frequency: " << ztrapfreq << "Hz" << endl;

    cerr << "Length scale oscillator unit: " << length*(1e6) << " microns" 
         << endl;
    cerr << "Time scale oscillator unit: " << time*(1e3) << " milliseconds" 
         << endl;
    cerr << "Energy scale oscillator unit: " << energy << " Joules" << endl;
    cerr << "Energy scale oscillator unit: " << (energy*(1e9))/(kboltz) 
         << " nanokelvin" << endl;


    cerr << "3D interaction parameter g3D: " << g11conv << endl;
    cerr << "Effective 1D interaction parameter g1D: " << g11 << endl;

    cerr << "Ring radius: " << ringrad << " microns = " << ringrad/(length*1e6) 
         << " osc. units" << endl;

return;
}

/*----------------------------------------------------------------------------*/

void setomegaou
(
double &omegaou,
double ztrapfreq,
double rtrapfreq,
int ouselect,
double omega1,
double husimisigma,
double husimigdepth,
int husimitype
)

{
    if(ouselect == 0)
    {
        switch(husimitype)
        {
            case 1:
                omegaou = omega1;
                cerr << "Osc. units based on harmonic trap freq." << endl;
                break;
            case 2:
                omegaou = sqrt(husimigdepth)/husimisigma; //harmonic approx.
                cerr << "Osc. units based on central harmonic approximation "
                     << "to central Gaussian" << endl;
                break;
            default:
                cerr << "Cannot set oscillator frequency!" << endl;
                break;
        }
    }
    else if(husimitype == 0)
    {
        switch(ouselect)
        {
            case 1:
                omegaou = 2.0*PI*rtrapfreq;
                cerr << "Osc. units based on radial trap frequency" << endl;
                break;
            case 2:
                omegaou = 2.0*PI*ztrapfreq;
                cerr << "Osc. units based on z-plane trap frequency" << endl;
                break;
            default:
                cerr << "Cannot set oscillator frequency!" << endl;
                break;
        }
    }

    cerr << "Oscillator frequency: " << omegaou << endl;

return;
}

/*----------------------------------------------------------------------------*/

void setoscunits
//convert quantities to osc units
(
 double omegaou,
 double &xupper,
 double &xlower,
 double &tmax,
 double &tstep
)

{
    double beta = sqrt(hbar/(rbmass*omegaou));
    cerr << "beta: " << beta << endl;

    cerr << "xupper: " << xupper << endl;
    xupper /= omegaou;
    cerr << "xupper: " << xupper << endl;
    cerr << "xlower: " << xlower << endl;
    xlower = 0.0; //temporary!
    cerr << "xlower: " << xlower << endl;
    cerr << "tmax: " << tmax << endl;
    tmax   *= omegaou;
    cerr << "tmax: " << tmax << endl;
    cerr << "tstep: " << tstep << endl;
    tstep  *= omegaou;
    cerr << "tstep: " << tstep << endl;
    

return;
}

/*----------------------------------------------------------------------------*/

vector<double> threemodeoptimise
// gradient descent function,  only optimises on a,delta,a0
(// can't pass by reference here
 double a0,     //constant background
 double a,      //average modulation amplitude
 double delta,  //mode splitting
 double l,      //angular momentum quantum number
 double omega,  //rotation rate
 double alpha,  //imprinting potential amplitude
 double res,    //search resolution
 double gii     //interaction strength
)

{
//cout << "threemodeoptimise() omega: " << omega << endl;
    vector<double> parameters;
    double a0_old, a_old, delta_old;
    double graddelta = 1.0, grada = 1.0, grada0 = 1.0;

    int loopcount = 0;

//    a0 = 0.997927;

cout << "Optimising 3MM to within " << globtol << endl;

    while((abs(grada0) > globtol) || (abs(grada) > globtol) || (abs(graddelta) > globtol)) 
    { //this condition continues the loop as long as one parameter is still improving

        loopcount++;   //keep track of iterations
        
        //set old params (if gradient descent does nothing, then exit loop
        a0_old = a0;
        a_old = a;
        delta_old = delta;
        
        // calculate derivatives using finite differences
        grada0 = -(0.5*(threemodeHexpect(a0+res,a,delta,l,omega,alpha,gii)-threemodeHexpect(a0-res,a,delta,l,omega,alpha,gii)));
        grada = -(0.5*(threemodeHexpect(a0,a+res,delta,l,omega,alpha,gii)-threemodeHexpect(a0,a-res,delta,l,omega,alpha,gii)));
        graddelta = -(0.5*(threemodeHexpect(a0,a,delta+res,l,omega,alpha,gii)-threemodeHexpect(a0,a,delta-res,l,omega,alpha,gii)));
        
//        cout << "grada0: " << grada0 << "\t" << "grada: " << grada << "\t" << "graddelta: " << graddelta << endl;

        // set new parameter values from gradient
        a0 += grada0;
        a  += grada;
        delta += graddelta;

//        cout << "a0: " << setprecision(15) << a0 << "\t" << "a: " << a << "\t" << "delta: " << delta << endl;
        // infinite loop avoidance
        if(loopcount > 1e7)
        {
            cout << "3MM search reached 10^7 steps, aborting" << endl;
            break;
        }
    }
/*
for(int j=0; j<3; j++)
{
    grada = 1.0;
    grada0 = 1.0;
    graddelta = 1.0;

    while(abs(grada) > globtol)
    {
        a_old = a;
        grada = -(0.5*(threemodeHexpect(a0,a+res,delta,l,omega,alpha)-threemodeHexpect(a0,a-res,delta,l,omega,alpha)));
        a += grada;
    }


    while(abs(grada0) > globtol)
    {
        a0_old = a0;
        grada0 = -(0.5*(threemodeHexpect(a0+res,a,delta,l,omega,alpha)-threemodeHexpect(a0-res,a,delta,l,omega,alpha)));
        a0 += grada0;
    }

    while(abs(graddelta) > globtol)
    {
        delta_old = delta;
        graddelta = -(0.5*(threemodeHexpect(a0,a,delta+res,l,omega,alpha)-threemodeHexpect(a0,a,delta-res,l,omega,alpha)));
        delta += graddelta;
    }

}
*/
cout << "grada0: " << grada0 << endl;
cout << "grada: " << grada << endl;
cout << "graddelta: " << graddelta << endl;




    //report number of iterations required
    cout << "3MM optimisation steps: " << loopcount << endl;
    cout << "Optimal a0: " << a0 << endl;
    cout << "Optimal a : " << a  << endl;
    cout << "Optimal delta: " << delta << endl;
    cout << "delta/a ratio in units of pi: " << delta/(a*PI) << endl;


    //when done, return the optimal parameters in a vector
    parameters.resize(6);
    parameters[0] = a0;
    parameters[1] = a;
    parameters[2] = delta;
    parameters[3] = l;
    parameters[4] = omega;
    parameters[5] = alpha;

return parameters;
}

/*----------------------------------------------------------------------------*/

double threemodeHexpect
// hamiltonian expectation value function for three-mode model optimisation
// ignores potential constant V0 (since it's constant)
(
 double a0,     //constant background
 double a,      //average modulation amplitude
 double delta,  //mode splitting
 double l,      //angular momentum quantum number
 double omega,  //rotation rate
 double alpha,   //imprinting potential amplitude
 double gii     // self-interaction strength
)

{
    // Interacting result
    return (1.0/(a0*a0 + 2.0*pow(a,2.0) + 2.0*pow(delta,2.0)))
            *(0.5*pow(l,2.0)*pow(2.0*PI*0.1,2.0)*(2.0*pow(a,2.0) + 2.0*pow(delta,2.0)) 
              + 4.0*l*omega*a*delta - 2.0*alpha*a0*a 
              + gii*(8.0*pow(a0,2.0)*pow(a,2.0) + pow(a,4.0) + pow(delta,4.0) 
              - 2.0*pow(a,2.0)*pow(delta,2.0)));
    
/*
 //Non-interacting variational result
    return (1.0/(a0*a0 + 2.0*a*a + 2.0*delta*delta))
            *(0.5*l*l*pow(2.0*PI*0.1,2.0)*(2.0*a*a + 2.0*delta*delta) + 4.0*l*omega*a*delta - 2.0*alpha*a0*a);
*/
    
/*
    //non-rotating custom variational
    return (1.0/(a0*a0 + 2.0*a*a + 2.0*delta*delta))
            *(0.5*l*l*pow(2.0*PI*0.1,2.0)*(2.0*a*a + 2.0*delta*delta) - 2.0*alpha*a0*a);
*/

/* 
 //Bromley 2MM result
    return (1.0/(a0*a0 + a*a))
            *(0.5*l*l*a*a - a0*a*alpha*sqrt(2.0));
*/

/* 
 //Original variational result
    return (1.0/(a0*a0 + 2.0*a*a + 2.0*delta*delta))
            *(0.5*l*l*(2.0*a*a + 2.0*delta*delta) + 4.0*l*omega*a*delta - 4.0*alpha*a0*a);
*/

}

/*----------------------------------------------------------------------------*/

double threemodeHexpect_manual
// expectation value of Hamiltonian for wavefunction constructed using
// given parameters
//
// 
// TESTING ONLY FOR NOW, SOME VALUES HARDCODED
// IF VALID, REWORK THIS FUNCTION
(
 double a0,     //constant background
 double a,      //average modulation amplitude
 double delta,  //mode splitting
 double l,      //angular momentum quantum number
 double omega,  //rotation rate
 double alpha   //imprinting potential amplitude
)

{
    vector< complex<double> > wavefunction;
    complex<double> i = {0.0,1.0};
    complex<double> energy = {0.0,0.0};
    int j;

    int Nx = 105;
    double xmin = 0.0;
    double xmax = 9.90476; 
    double xstep = 0.0952381;
    double g = 0.0;

    wavefunction.resize(Nx);
    for(j=0;j<Nx;j++)
    {
        wavefunction[j] = a0 + (a + delta)*exp( i*l*2.0*PI*((double) j/(double) Nx))
                             + (a - delta)*exp(-i*l*2.0*PI*((double) j/(double) Nx));
    }

    //ENERGY CALC - TAKEN FROM calcenergy() WITH POTENTIAL TERM REMOVED

	// energy sensitive to cyclic boundary condition
    energy += conj(wavefunction[0])*xstep*
              (
              -0.5*((wavefunction[1] - 2.0*wavefunction[0] + wavefunction[Nx-1])/(pow(xstep,2.0)))
                - (0.5*g*pow(abs(wavefunction[0]),2.0)*wavefunction[0])
                - i*omega*((xmax-xmin)/(2.0*PI))*((wavefunction[1] - wavefunction[Nx-1])/(2.0*xstep))
              );

    energy += conj(wavefunction[Nx-1])*xstep*
              (
              -0.5*((wavefunction[0]-2.0*wavefunction[Nx-1]+wavefunction[Nx-2])/(pow(xstep,2.0)))
              - (0.5*g*pow(abs(wavefunction[Nx-1]),2.0)*wavefunction[Nx-1])
              - i*omega*((xmax-xmin)/(2.0*PI))*((wavefunction[0] - wavefunction[Nx-2])/(2.0*xstep))
              );
    	// rest of the array
	for(j=1;j<(Nx-1);j++)
	{
    energy += conj(wavefunction[j])*xstep*
              (
              -0.5*((wavefunction[j+1]-2.0*wavefunction[j]+wavefunction[j-1])/(pow(xstep,2.0)))
              - (0.5*g*pow(abs(wavefunction[j]),2.0)*wavefunction[j]) // -ve as potential[] includes interactions
              - i*omega*((xmax-xmin)/(2.0*PI))*((wavefunction[j+1] - wavefunction[j-1])/(2.0*xstep))
              );
    }

return abs(energy);
}

/*----------------------------------------------------------------------------*/

vector< complex<double> > calcthreemode
// returns a three-mode model state
(
 vector<double> &x,
 double a0,
 double a,
 double delta,
 double l
)

{
    int j;
    int Nx = x.size();
    complex<double> i = {0.0,1.0};
    vector< complex<double> > state;
    double pnorm = 0.0;

    state.resize(Nx);

    for(j=0;j<Nx;j++)
    {
        state[j] = a0 + (a+delta)*exp(i*l*((double) j/(double) Nx)) 
                      + (a-delta)*exp(-i*l*((double) j/(double) Nx));
        pnorm += pow(abs(state[j]),2.0);
    }

    if(pnorm - 1.0 < globtol)
    {
        cout << "Normalising three-mode state" << endl;
        for(j=0;j<Nx;j++)
        {
            state[j] *= sqrt(1.0/pnorm);
        }
    }

return state;
}

/*----------------------------------------------------------------------------*/

vector<double> strobo_rot_measurement
// determines rotation rate at stroboscopic times based on g=0 analytics
(
 vector<double> &phasevec,
 vector<double> &timevec,
 double tstep,
 double tmax,
 double xmax,
 double xmin,
 double omega,
 int imm
)

{
    vector<double> results;
    double oneosc;
    double radius;
    double numosc;
    int j;
    int index;

    //determine stroboscopic times from analytics
    /*
     * Here we use:
     * t = (n \pi \hbar)/(E-E0)
     *   = (n \pi) / (l^2/2R^2)
     * with n=1 to give the oscilation period
     */

    //calculate the ring radius, and therefore the oscillation time
    radius = (xmax-xmin)/(2.0*PI);
    oneosc = PI/((imm*imm)/(2.0*radius*radius));

    //determine the number of full oscillations in the sim, size results vector
    numosc = floor(tmax/oneosc);
    results.resize(numosc);
    timevec.resize(numosc);


    for(j=0;j<numosc;j++)
    {
        //todo: include both odd and even measurements (should have 2n+1 total measurements,
        //      where n is the number of oscillations
        results[j] = phasevec[(int) (((j+1)*oneosc)/tstep)] - phasevec[(int) ((j*oneosc)/tstep)];
        timevec[j] = (j+1)*oneosc;
    }

return results;
}
