/*
	BECring
	Brad Mommers 2018-2019

	This program uses the Sherman-Morrison method to solve the nonlinear
	Schrodinger equation for a one- or two-component BEC confined to a 1D ring.
	Initial state is a Gaussian that can be initially relaxed with imaginary
	time evolution, in addition the damped GPE can be simulated by evolving
	forward one frame in imaginary time every m frames of real time evolution.

	Example compile:
	$ g++ BECring.cpp -o becring -lm -std=c++11 -Wall -Wextra

	Example usage:
	$ ./becring input.txt output.txt outdir

	Input: from specified input file. 
			One argument per line, all lines beginning with '#' will be ignored.
			Blank lines also ignored.
			The pre-packaged input file contains comment lines before each
			parameter giving a brief description and the variable name.

	Output:
	File "<outfilename>" containing stats, optional snapshots of wavefunction
	and probability density for animation in the frames/ (real time) and
	iframes/ (imaginary time) directories.
*/

#include<cstdlib>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<complex>
#include<vector>
#include "bec-headers.h"
#include "shermanmorrison.h"

#define PI 3.14159265358979323846264338

using namespace std;

const complex<double> i = {0.0, 1.0};

/********************************MAIN FUNCTION*********************************/
int main(int argc, char *argv[]){

// Variable declarations

	int Nx, Nt, Nit, tc, ani, ierr, ite, modulus, numco, inbetween = 0;
    int impot, imm, irotate, rrotate, nfreq, threemode;
    int husimitype;
    int imprinttype, imptimeprofile;
    int ouselect, gdim;

    double ringrad, ztrapfreq, rtrapfreq;
    double imprintamp, impoffset, impsigma;
	double a11, a22, a12, g11, g22, g12;
	double k01, k02;
	double sigma01, sigma02;
	double x01, x02;
	double x01ou, x02ou;
	double omega1, omega2, potamp;
	double xlower, xupper, xstep;
	double xlowerou, xupperou, xstepou;
	double tstep, tmax, time;
	double tstepou, tmaxou, timeou;
    double fweight;
    double dfreq;
    double omegaimag, omegareal, tgamma;
    double omegaimagou, omegarealou;
    double tmphase, modprop, modebias;
    double threshold;
    double imtcutoff, sum;
    double husimisigma, husimigdepth;
    double unitlength, unittime, unitenergy, omegaou;
    double g11conv, g22conv, g12conv;
    double fouriernorm;
    double res; // three-mode optimisation resolution
    double omegatemp;
    double bgmode;

	complex<double> zalpha, zialpha, zbeta, zibeta, zgamma, zigamma;
    complex<double> fouriertransform;
    complex<double> densityfouriertransform;

	//for IO
	string infilename;
	string outfilename;
	string framefilename;
    string outdir;
    string temp;
    string statefilename;
    string outstatefname;

    vector<string> inputlist;

	vector<double> statvec1, statvec2, aggstatvec;
	vector<double> probdensity1, probdensity2;
    vector<double> fourierlist; // list of frequencies for single-freq FT
	vector<vector<double> * > wfprintlist;
    vector<double> phasevec,gradvec;
    vector<double> threemode_params;
	
	//vectors
	vector< complex<double> > zoldpsi1, zadiag1, ziadiag1, zaupper, ziaupper;
	vector< complex<double> > zalower, zialower, zrrhs1, zxsoln1;
	vector< complex<double> > zuvec1, ziuvec1, zvvec1, zivvec1, zpsi01;
	vector< complex<double> > zoldpsi2, zadiag2, ziadiag2, zrrhs2, zxsoln2;
	vector< complex<double> > zuvec2, ziuvec2, zvvec2, zivvec2, zpsi02;
    vector< complex<double> > modeamplist;
    vector< complex<double> > modeampmatrix;
    vector<double> phaseftlist;
	vector<double> x, potential1, potential2, printpotential;
    vector<double> densitytime;
    vector<double> statsmatrix, istatsmatrix, phasextmatrix, modeampabsmatrix;


	//loop variables
	int j, k, t, rt=0;


/****************************INITIALISATION************************************/

	// check correct number of arguments supplied
	if(argc != 4){
		printf("Usage: %s <\"input filename\">"
						"<\"output filename\">\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	// read arguments
	infilename = argv[1];
	temp = argv[2];
    outdir = argv[3]; //folder restructure
    
    outfilename = outdir + "/" + temp;

    
    cout << "Loading input file" << endl;
	// retrieve input parameters from file
	ifstream inputf;
	inputf.open(infilename);
	loadinput(infilename,inputlist);


	// assign values retrieved from file
	xlower = stof(inputlist[0]);
	xupper = stof(inputlist[1]);
	Nx = (int) stof(inputlist[2]);
    tmax = stof(inputlist[3]);
	tstep = stof(inputlist[4]);
	Nt = (int) stof(inputlist[5]);
	Nit = (int) stof(inputlist[6]);
    imtcutoff = stof(inputlist[7]);
    outstatefname = inputlist[8];
	numco = stoi(inputlist[9]);
    gdim = stoi(inputlist[10]);
	a11 = stof(inputlist[11]);
	a22 = stof(inputlist[12]);
	a12 = stof(inputlist[13]);
    statefilename = inputlist[14];
	k01 = stof(inputlist[15]);
	k02 = stof(inputlist[16]);
	sigma01 = stof(inputlist[17]);
	sigma02 = stof(inputlist[18]);
	x01 = stof(inputlist[19]);
	x02 = stof(inputlist[20]);
    husimitype = stoi(inputlist[21]);
    husimisigma = stof(inputlist[22]);
    husimigdepth = stof(inputlist[23]);
	omega1 = stof(inputlist[24]);
	omega2 = stof(inputlist[25]);
	potamp = stof(inputlist[26]);
	ani = stoi(inputlist[27]);
	ite = stoi(inputlist[28]);
	modulus = stoi(inputlist[29]);
    ringrad = stof(inputlist[30]);
    ztrapfreq = stof(inputlist[31]);
    rtrapfreq = stof(inputlist[32]);
    ouselect = stoi(inputlist[33]);
    impot = stoi(inputlist[34]);
   	imm = stoi(inputlist[35]);
	fweight = stof(inputlist[36]);
    imprintamp = stof(inputlist[37]);
    imprinttype = stoi(inputlist[38]);
    imptimeprofile = stoi(inputlist[39]);
    impoffset = stof(inputlist[40]);
    impsigma = stof(inputlist[41]);
    omegaimag = stof(inputlist[42]);
    omegareal = stof(inputlist[43]);
    tgamma = stof(inputlist[44]);
    irotate = stoi(inputlist[45]);
    rrotate = stoi(inputlist[46]);
    threemode = stoi(inputlist[47]);
    tmphase = stof(inputlist[48]);
    bgmode = stof(inputlist[49]);
    modprop = stof(inputlist[50]);
    modebias = stof(inputlist[51]);
    res = stof(inputlist[52]);
    threshold = stof(inputlist[53]);
    nfreq = (int) stof(inputlist[54]);


    if(nfreq > 0)
    {
        cout << "Loading " << nfreq << " frequencies" << endl;
        fourierlist.resize(nfreq);
        for(j=0;j<nfreq;j++)
        {
            fourierlist[j] = stof(inputlist[55+j]);
        }
    } else if(nfreq == 0) //if nfreq = 0; dump the rest of inputlist as freqs
    {
        cout << "Loading all frequencies" << endl;
        fourierlist.resize( (int) inputlist.size() - 55);
        for(j=55;j< (int) inputlist.size(); j++)
        {
            fourierlist[j-55] = stof(inputlist[j]);
        }
        cout << fourierlist.size() << " frequencies loaded" << endl;
        nfreq = (int) fourierlist.size(); //needed as later code depends on
                                          // nfreq having correct value
    }


    cout << "Input file loaded" << endl;

	// declare array sizes
    phaseftlist.resize(nfreq);
    modeamplist.resize(nfreq);


//    phasevec.resize((nfreq+1)*Nt);//extra column for time values

	potential1.resize(Nx);
	potential2.resize(Nx);
    printpotential.resize(Nx);
	x.resize(Nx);

	zaupper.resize(Nx);
	ziaupper.resize(Nx);
	zalower.resize(Nx);
	zialower.resize(Nx);

	zoldpsi1.resize(Nx);
	zadiag1.resize(Nx);
	ziadiag1.resize(Nx);
	zrrhs1.resize(Nx);
	zxsoln1.resize(Nx);
	zuvec1.resize(Nx);
	ziuvec1.resize(Nx);
	zvvec1.resize(Nx);
	zivvec1.resize(Nx);
	zpsi01.resize(Nx);
	probdensity1.resize(Nx);

	zoldpsi2.resize(Nx);
	zadiag2.resize(Nx);
	ziadiag2.resize(Nx);
	zrrhs2.resize(Nx);
	zxsoln2.resize(Nx);
	zuvec2.resize(Nx);
	ziuvec2.resize(Nx);
	zvvec2.resize(Nx);
	zivvec2.resize(Nx);
	zpsi02.resize(Nx);
	probdensity2.resize(Nx);

    threemode_params.resize(6);

	statvec1.resize(10);
	if(numco == 2)
	{
		statvec2.resize(9);
		aggstatvec.resize(17);
	} else if (numco == 1)
	{
		aggstatvec.resize(10);
	}

/******************************UNITS*******************************************/


    switch(gdim)
    {
        //todo: add other components to conversions
        case 1:
            // given g_1D, want to also print g_3D
            g11conv = pow(a11,2.0)*pow(25.0/18.0,2.0)*(3.0/(4.0*pow(PI,3.0)*ztrapfreq*rtrapfreq));
            break;
        case 3:
            // given g_3D, want to print both g and use g_1D for calculations
            g11conv = a11; // this is g_3D
            a11 = (18.0/25.0)*sqrt(g11conv*((4.0*pow(PI,3.0)*ztrapfreq*rtrapfreq)/3.0));
            break;
    }

    //set omega for converting to osc. units
    setomegaou(omegaou,ztrapfreq,rtrapfreq,ouselect,omega1,husimisigma,husimigdepth,husimitype);

    //convert space and time to osc. units
    //print unit conversions to filk
    printconversions(ringrad,ztrapfreq,rtrapfreq,a11,g11conv,ouselect,omegaou,xupper,xlower,tmax);

    // omegaou does not correctly account for Husimi driving case yet!

/***************************ARRAY SETUP****************************************/

	// set up spatial grid
	setxgrid(x, xupper, xlower, Nx);
	xstep = (xupper-xlower)/((double) Nx);
    
    // convert rotation values to multiples of PI
    omegareal *= PI;
    omegaimag *= PI;
    omegatemp = omegaimag;
    
    // disable rotation based on switch values
    // this method avoids edge cases where
    // omega____ is set but the switch is off
    // such as in energy calculation (which
    // uses omega____ directly)
    if(irotate == 0)
    {
        omegaimag = 0.0;
    }

    if(rrotate == 0)
    {
        omegareal = 0.0;
    }


/*
    cout << "omegareal: " << omegareal << endl;
    cout << "omegaimag: " << omegaimag << endl;
*/    

    // if tstep=0, determine its value from Crank-Nicolson stability criteria
    // see NR eq 20.2.6
    // if tmax, tstep, Nt all defined, re-calculate tstep to match the others
    if(tstep == 0.0)
    {
        cout << "Calculating stable timestep size for xstep = " << xstep << endl;
        tstep = 0.1*pow(xstep,2.0);
        cout << "tstep set to: " << tstep << endl;
    }

    // if Nt=0, determine it from final time and timestep
    if(Nt == 0)
    {
        Nt = (int) ceil(tmax/tstep);
        cout << "Nt set to: " << Nt << endl;
    }

    //now that Nt is defined, we size modeampmatrix appropriately
    modeampmatrix.resize(nfreq*Nt);
    modeampabsmatrix.resize(nfreq*Nt);
    phasextmatrix.resize(nfreq*Nt);
    densitytime.resize(Nx*Nt);
    statsmatrix.resize((aggstatvec.size())*Nt);
    istatsmatrix.resize((aggstatvec.size())*Nit);

	// calculate alpha, beta, gamma factors for Sherman-Morrison
	//real time
	zalpha = i*(tstep/(4.0*pow(xstep,2.0)));
	zbeta = i*(tstep/2.0);
    zgamma = -(tstep/(8.0*PI*xstep))*(xupper-xlower)*omegareal; //rotation rate not in units of PI
	//imaginary time
	zialpha = (tstep/(4.0*pow(xstep,2.0)));
	zibeta = (tstep/2.0);
    zigamma = i*(tstep/(8.0*PI*xstep))*(xupper-xlower)*omegaimag; //rotation rate not in units of PI

	// set initial waveform
    if(statefilename.compare("none") != 0) //check if initial state file supplied
    {
        cout << "Loading state from file" << endl;
        loadstatefromfile(statefilename,zoldpsi1);
    }
    else if(threemode > 0) //one-component only for now
    {
        cout << "Setting initial state according to 3-mode model" << endl;
        threemode_params = initialisethreemode(zoldpsi1,x,tmphase,bgmode,modprop,modebias,imm,omegaimag,imprintamp,res,a11);
        for(j=0;j<Nx;j++)
        {
            zpsi01[j] = zoldpsi1[j]; //copy zoldpsi to avoid re-finding optimal state
        }
        //initialisethreemode(zpsi01  ,x,tmphase,modprop,modebias,imm,omegaimag,imprintamp,res);
    }
    else
    {
        cout << "Setting initial state as Gaussian" << endl;
        setinitialGaussian(zoldpsi1, x, x01, sigma01, xstep, xupper,xlower,Nx);
	    setinitialGaussian(zpsi01, x, x01, sigma01, xstep, xupper, xlower, Nx);
	    if(numco == 2)
	    {
	    	setinitialGaussian(zoldpsi2, x, x02, sigma02, xstep, xupper,xlower,Nx);
	    	setinitialGaussian(zpsi02, x, x02, sigma02, xstep, xupper, xlower, Nx);
	    } else if(numco == 1)
	    {
	    	for(j=0;j<Nx;j++)
	    	{
	    		zoldpsi2[j] = {0.0,0.0};
	    		zxsoln2[j] = {0.0,0.0};
	    		potential2[j] = 0.0;
	    	}
	    }
    }

	// initialise potentials, tridag arrays
	// calculate g values from scattering lengths
	g11 = a11;//*2.0*PI; //simplifying interactions, may restore to scattering lengths later
	if(numco == 2)
	{
		g22 = a22;//*2.0*PI
		g12 = a12;//*2.0*PI
	} else if(numco == 1)
	{
		g22 = 0.0;
		g12 = 0.0;
	}

// put a switch here to pick an external potential?
	updatepotential(x,potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,g11,g12,omega1,omega2,0.0,potamp,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,1.0);
	updatepotential(x,potential2,zoldpsi2,zoldpsi1,impot,imm,fweight,rt,g22,g12,omega1,omega2,0.0,potamp,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,1.0);
    // potential without interactions for plotting
	updatepotential(x,printpotential,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,0.0,0.0,omega1,omega2,0.0,1.0,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,0.0);



	//prepare vectors for Sherman-Morrison
	// real time
	prepSMvecs(zadiag1, zaupper, zalower, zuvec1, zvvec1, potential1, zalpha,zbeta,zgamma);
	// imaginary time
	prepSMvecs(ziadiag1,ziaupper,zialower,ziuvec1,zivvec1,potential1,zialpha,zibeta,zigamma);
	// same for second component
	if(numco == 2)
	{
		prepSMvecs(zadiag2, zaupper, zalower, zuvec2, zvvec2, potential2, zalpha, zbeta, zgamma);
		prepSMvecs(ziadiag2, ziaupper, zialower, ziuvec2, zivvec2, potential2, zialpha, zibeta, zigamma);
	}

/*****************************SET UP TIME LOOP*********************************/

	// set up files for stats output

    
    // phase extraction in separate file
    /*
    */


    // mode evolution in separate file
    /*
    ofstream modeamp;
    modeamp.open(outdir + "/" + "modes.dat");
    modeamp << "time";

    ofstream modeampabs;
    modeampabs.open(outdir + "/" + "modes-abs.dat");
    modeampabs << "time";

    modeampabs << endl;
    */

/***************************CONTROL LOOP***************************************/
//big loop - first iteration is imaginary time, second is realtime


for(rt=0;rt<2;rt++)
{

	//check time switch var has a valid value
	if((rt != 0) && (rt != 1))
	{
        // throw an exception instead of printing
		printf("Time switch (rt) error! Exiting...\n");
		exit(EXIT_FAILURE);
	}

    switch(rt)
    {
        case 0:
            cout << "Setting up imaginary time evolution" << endl;
            break;
        case 1:
            cout << "Setting up real time evolution" << endl;
            break;
    }

    //first: export result of imaginary time evo to file, if specified
    if((rt == 1) && (outstatefname.compare("none") != 0))
    {
        cout << "Saving state to file" << endl;
        savestatetofile(outstatefname,zoldpsi1);
    }

    // add kick when moving to real time
    if((rt == 1) && (k01 > 0.0))
    {
        cout << "Adding momentum kick" << endl;
        addmodekick(zoldpsi1,x,x01,k01,tmphase,Nx);
    }

    // set real time to zero (so stats etc for first frame are correct)
	t = 0;
    // update potential for real time
	updatepotential(x,potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,g11,g12,omega1,omega2,t*tstep,potamp,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,1.0);
    updatepotential(x,printpotential,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,0.0,0.0,omega1,omega2,t*tstep,1.0,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,0.0);


	// set initial waveform for real time overlap calculations
	if(rt == 1)
    {
		for(j=0;j<Nx;j++)
		{
			zpsi01[j] = zoldpsi1[j];

			if(numco == 2)
			{
				zpsi02[j] = zoldpsi2[j];
			}
		}
	}
	
    if(rt==1)
	{
        //output density to densitytime matrix
        for(j=0;j<Nx;j++)
        {
            densitytime[j] = pow(abs(zpsi01[j]),2.0);
        }
    }
    
	// calc initial stats & print to file
	//calc stats per component
	if (numco == 1)
	{
		calcstats(zoldpsi1, zpsi01, potential1, statvec1, x, xstep, t, tstep, 1, omegaou);
        if (rt == 0)
        {
		    statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,potential1,x,xstep, g11,omegaimag);
        } else {
		    statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,potential1,x,xstep, g11,omegareal);
        }
        for(j=0;j<(int) aggstatvec.size(); j++)
		{
			aggstatvec[j] = statvec1[j];
		}
	}
	else if(numco == 2)
	{
		calcstats(zoldpsi1, zpsi01, potential1,statvec1, x, xstep, t, tstep, 1, omegaou);
		statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,zoldpsi2,potential1,x,xstep,g11,g12);
		calcstats(zoldpsi2, zpsi02, potential2,statvec2, x, xstep, t, tstep, 1, omegaou);
		statvec2[statvec2.size()-1] = calcenergy(zoldpsi2,zoldpsi1,potential2,x,xstep,g22,g12);
	
	
		//aggregate into single vector -- put this into a function?
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

	}


    //phase extraction and mode evolution (real time only)
    if( (nfreq > 0) && (rt == 1))
    {
        fouriernorm = 0.0;
        for(j=0;j<nfreq;j++)
        {
            fouriertransform = fixedmodefourierint(zoldpsi1,x,fourierlist[j]);
            densityfouriertransform = fixedmodedensityfourierint(zoldpsi1,x,fourierlist[j]);
            // phase is taken from density FT, complex mode amp is from WF FT
            phaseftlist[j] = complexphase(densityfouriertransform);
            modeamplist[j] = fouriertransform;
            modeampmatrix[j*Nt] = fouriertransform;
            fouriernorm += pow(abs(fouriertransform),2.0);

            //print the t=0 mode amplitudes for the background and imprinted modes
            if((t == 0) && ((fourierlist[j] == imm) || (fourierlist[j] == -imm) || (fourierlist[j] == 0)))
            {
                cout << setw(12) << setprecision(12) << "t=0 Fourier amplitude, n=" << fourierlist[j] << ": " << abs(fouriertransform) << endl;
            }
        }
        cout << "Fourier decomp norm at t=0: " << fouriernorm << endl;



    }


	//print to file
    // REPLACED WITH POPULATING statsvecmatrix,istatsvecmatrix
	if(rt == 0)
	{
		//statfileprint(istats,aggstatvec);
        for(j=0;j<(int) aggstatvec.size();j++)
        {
            //put aggstatvec in as a row of istatsvec
            istatsmatrix[j] = aggstatvec[j];
        }

	} 
	else if(rt ==1)
	{
		//statfileprint(stats,aggstatvec);
        for(j=0;j<(int) aggstatvec.size();j++)
        {
            //put aggstatvec in as a row of istatsvec
            statsmatrix[j] = aggstatvec[j];
        }
        
	}

	// print initial wavefunction to file
	if((rt == 0) && (ani ==1))
	{
		ofstream ioutput;
		twocompfileprint(ioutput, outdir + "/" + "iframes/iframe00000.dat",x,zoldpsi1,
						 zoldpsi2, printpotential, potential2);
	}
	else if((rt == 1) && (ani == 1))
	{
		ofstream output;
		twocompfileprint(output, outdir + "/" + "frames/frame00000.dat", x, zoldpsi1, zoldpsi2,
						printpotential, printpotential);
	}

	//recalculate diagonal & u,v vectors when moving to real time
	if(rt == 1)
	{
		prepSMvecs(zadiag1, zaupper, zalower, zuvec1, zvvec1, potential1, zalpha, zbeta, zgamma);
			
		if(numco == 2)
		{
			prepSMvecs(zadiag2, zaupper, zalower, zuvec2,zvvec2, potential2, zalpha, zbeta, zgamma);
		}
	}
	
	
	// determine max number of frames for imaginary, real, real + damped time
	if((rt == 0) && (threemode == 0))
	{
		tc = Nit;//+1;
        cout << "Entering imaginary time" << endl;
	}
    else if((rt == 0) && (threemode != 0))
    {
        tc = 0;
        cout << "Using three-mode model, bypassing imag. time" << endl;
    }
	else if((rt == 1) && (ite == 1))
	{
		tc = Nt + (Nt/modulus);
        cout << "Entering damped real time" << endl;
	}
	else if((rt == 1) && (ite == 0))
	{
		tc = Nt;
        cout << "Entering real time" << endl;
	}
	else
	{
        //throw an exception here instead of printing & exiting
		printf("ite switch error. Exiting...\n");
		exit(EXIT_FAILURE);
	}

/******************************TIME LOOP***************************************/

	for(t=1;t<tc;t++)
	{
//cout << "start timestep " << t << endl;
		//if it is the modulus'th frame, do imaginary time
		if((ite == 1) && ((t%modulus) == 0) && (rt == 1))
		{
            cout << "In-between frame at t=" << t << endl;
			inbetween = 1;
		}
		// need a separate call to calcpotential() to populate printpotential vector
		// calc potential, rhs, SM vecs and do optimised split-operator type time evolution
		if((rt == 0) || (inbetween == 1))
        {
            updatepotential(x, potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt, g11, g12, omega1, omega2, 0.0, potamp,husimisigma,husimigdepth, husimitype, imprintamp, impoffset, impsigma, imprinttype, imptimeprofile, omegaou, 1.0);
            
            updatepotential(x,printpotential,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,0.0,0.0,omega1,omega2,0.0,1.0,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,0.0);
            
            ierr = osoevo(x,potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,g11,g12,omega1,omega2,0.0,potamp,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,1.0, ziadiag1, ziaupper, zialower, ziuvec1, zivvec1, zialpha, zibeta, zigamma, zrrhs1, zoldpsi1, Nx, zxsoln1,Nx,tstep);
        }
        else if (inbetween == 1)
        {
            updatepotential(x, potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt, g11, g12, omega1, omega2, (t-1)*tstep, potamp,husimisigma,husimigdepth, husimitype, imprintamp, impoffset, impsigma, imprinttype, imptimeprofile, omegaou, 1.0);
            
            updatepotential(x,printpotential,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,0.0,0.0,omega1,omega2,(t-1)*tstep,1.0,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,0.0);

            ierr = osoevo(x,potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,g11,g12,omega1,omega2,(t-1)*tstep,potamp,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,1.0, zadiag1, zaupper, zalower, zuvec1, zvvec1, zalpha, zbeta, zgamma, zrrhs1, zoldpsi1, Nx, zxsoln1,Nx,tstep);
        }
        else if (rt == 1)
        {
            updatepotential(x, potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt, g11, g12, omega1, omega2, t*tstep, potamp,husimisigma,husimigdepth, husimitype, imprintamp, impoffset, impsigma, imprinttype, imptimeprofile, omegaou, 1.0);
            
            updatepotential(x,printpotential,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,0.0,0.0,omega1,omega2,t*tstep,potamp/abs(potamp),husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,0.0);
            
            ierr = osoevo(x,potential1,zoldpsi1,zoldpsi2,impot,imm,fweight,rt,g11,g12,omega1,omega2,t*tstep,potamp,husimisigma,husimigdepth,husimitype,imprintamp,impoffset,impsigma,imprinttype,imptimeprofile,omegaou,1.0, zadiag1, zaupper, zalower, zuvec1, zvvec1, zalpha, zbeta, zgamma, zrrhs1, zoldpsi1, Nx, zxsoln1,Nx,tstep);
        }



		// calculate stats, print waveform & potential to file
		if(rt == 0) 
		{
			if(numco == 1)
			{
				calcstats(zxsoln1, zpsi01, potential1, statvec1, x, xstep, t, tstep, 1, omegaou);
				statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,potential1,x,xstep, g11,omegaimag);
			}
			else if(numco == 2)
			{
				calcstats(zxsoln1, zpsi01, potential1, statvec1, x, xstep, t, tstep, 1, omegaou);
				statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,zoldpsi2,potential1,x,xstep,g11,g12);
				calcstats(zxsoln2, zpsi02, potential2, statvec2, x, xstep, t, tstep, 1, omegaou);
				statvec2[statvec2.size()-1] = calcenergy(zoldpsi2,zoldpsi1,potential2,x,xstep,g22,g12);
			}		
		} else if((rt == 1) && (inbetween == 1))
		{
			if(numco == 1)
			{
				calcstats(zxsoln1,zpsi01,potential1,statvec1,x,xstep,t,tstep,1,modulus, omegaou);
				statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,potential1,x,xstep, g11,omegareal);
			}
			else if(numco == 2)
			{
				calcstats(zxsoln1,zpsi01,potential1,statvec1,x,xstep,t,tstep,1,modulus, omegaou);
				statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,zoldpsi2,potential1,x,xstep,g11,g12);
				calcstats(zxsoln2,zpsi02,potential2,statvec2,x,xstep,t,tstep,1,modulus, omegaou);
				statvec2[statvec2.size()-1] = calcenergy(zoldpsi2,zoldpsi1,potential2,x,xstep,g22,g12);
			}
		} else if((rt == 1) && (inbetween == 0)) //real time evolution doesn't require renormalisation
		{
			if(numco == 1)
			{
                if(ite != 0)
                {
                    calcstats(zxsoln1,zpsi01,potential1,statvec1,x,xstep,t,tstep,0,modulus, omegaou);
                    statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,potential1,x,xstep, g11,omegaimag);
                }
                else
                {
                    calcstats(zxsoln1,zpsi01,potential1,statvec1,x,xstep,t,tstep,0, omegaou);
                    statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,potential1,x,xstep, g11,omegaimag);
                }

			}
			else if(numco == 2)
			{
				calcstats(zxsoln1,zpsi01,potential1,statvec1,x,xstep,t,tstep,0, omegaou);
				statvec1[statvec1.size()-1] = calcenergy(zoldpsi1,zoldpsi2,potential1,x,xstep,g11,g12);
				calcstats(zxsoln2,zpsi02,potential2,statvec2,x,xstep,t,tstep,0, omegaou);
				statvec2[statvec2.size()-1] = calcenergy(zoldpsi2,zoldpsi1,potential2,x,xstep,g22,g12);
			}
		}


		//aggregate into single vector
		if(numco == 2)
		{
			for(j=0;j<(int) aggstatvec.size();j++)
			{
				if(j<(int) statvec1.size())
				{
					aggstatvec.at(j) = statvec1.at(j);
				}
				else if(j >= (int) statvec1.size())
				{
					aggstatvec.at(j) = statvec2.at(j-(statvec1.size()-2));
				}

			}
		} else
		{
//            cout << "aggstatvec.size() = " << aggstatvec.size() << endl;
//            cout << "statvec1.size() = " << statvec1.size() << endl;
			for(j=0;j < (int) aggstatvec.size(); j++)
			{
				aggstatvec[j] = statvec1[j];
			}
		}

		if(rt==1)
		{
            // add a row to densitytime matrix
//            cout << "writing densitytime for t=" << t << endl;
            for(j=0;j<Nx;j++)
            {
                densitytime[(Nx*t)+j] = pow(abs(zxsoln1[j]),2.0);
            }
        }
        
        //phase extraction (real time only)
        if( (nfreq > 0) && (rt == 1))
        {
            for(j=0;j<nfreq;j++)
            {
                fouriertransform = fixedmodefourierint(zxsoln1,x,fourierlist[j]);
                densityfouriertransform = fixedmodedensityfourierint(zxsoln1,x,fourierlist[j]);
                phaseftlist[j] = complexphase(densityfouriertransform);
                //modeamplist[j] = fouriertransform;
                modeampmatrix[(j*Nt)+t] = fouriertransform; 
                phasextmatrix[(t*nfreq)+j] = phaseftlist[j];
            }
        }

        // populate istatsmatrix,statsmatrix
		if(rt == 0)
		{
			//statfileprint(istats,aggstatvec);
            for(j=0;j<(int) aggstatvec.size();j++)
            {
                istatsmatrix[(t*((int) aggstatvec.size()))+j] = aggstatvec[j];
            }
		} 
		else if((rt == 1) && (ite == 1) && (inbetween == 0))
        {
            for(j=0;j<(int) aggstatvec.size();j++)
            {
                statsmatrix[((t-t/modulus)*((int) aggstatvec.size()))+j] = aggstatvec[j];
            }
        }
		else if((rt == 1) && (ite == 0))
		{
			//statfileprint(stats,aggstatvec);
            for(j=0;j<(int) aggstatvec.size();j++)
            {
                statsmatrix[(t*((int) aggstatvec.size()))+j] = aggstatvec[j];
            }
		}

		// printing waveform
		if(ani == 1) // faster to do this check first
		{
			if(rt == 0)
			{
				stringstream framenumpad;
				framenumpad << setw(5) << setfill('0') << t;
				framefilename = outdir + "/" + "iframes/iframe" + framenumpad.str() + ".dat";
				ofstream ioutput;
				twocompfileprint(ioutput, framefilename, x,zxsoln1, zxsoln2, 
								printpotential, potential2);
			}
			else if((rt == 1) && (inbetween == 0))
			{
				ostringstream framenumpad;
				framenumpad << setw(5) << setfill('0') << t;
				framefilename = outdir + "/" + "frames/frame" + framenumpad.str() + ".dat";
				ofstream output;
				twocompfileprint(output, framefilename, x,zxsoln1, zxsoln2, 
								printpotential, potential2);
			}
		}
        framefilename.clear();
//printf("completed stat output\n");


    if(rt == 0) //imaginary time convergence test
        //needs to be done while zsoln is new WF and zoldpsi is previous WF
    {
        if(abs(1.0-calcoverlap(zoldpsi1,zxsoln1,Nx,xstep)) < imtcutoff)
        {
            cout << "State converged in imaginary time within " << t << " steps, moving to real time" << endl;
            istatsmatrix.resize(aggstatvec.size()*t);
            t = Nit;
        }
    }

	//update zoldpsi, potential, zadiag (which is potential-dependent)
//cout << "About to update zoldpsi, potential, zadiag" << endl;
	for(j=0;j<Nx;j++)
	{
		zoldpsi1[j] = zxsoln1[j];
		zoldpsi2[j] = zxsoln2[j];
	}


	// reset rt if in real time and doing one frame of imaginary time evo.
	if(inbetween == 1)
	{
		inbetween = 0;
	}

//cout << "End timestep " << t << endl;
	} // time loop ends here

/***************************END TIME LOOP**************************************/

    switch(rt)
    {
        case 0:
            cout << "Max number of imaginary time frames reached" << endl;
            break;
        case 1:
            cout << "Real time evolution completed" << endl;
            break;
    }

} //rt loop ends here
/**************************END CONTROL LOOP************************************/

    //dump statsmatrix,istatsmatrix to file
    cout << "dumping stats to file" << endl;
	ofstream stats;
	stats.open(outfilename);
	stats << "#frame\t\t time \t\t norm \t\t <x1> \t\t <x^2> \t\t sqrt(var1) \t "
                "overlap1 \t <V1> \t\t <Lz> \t\t Energy1 \t\t norm2 \t\t\t <x2> \t\t "
                "sqrt(var2) \t\t overlap2 \t\t <V2> \t\t Energy2" << endl;
    for(k=0;k<Nt;k++)
    {
        for(j=0;j<(int) aggstatvec.size()-1;j++)
        {
            stats << setw(8) << statsmatrix[(k*(int) aggstatvec.size())+j] << "\t";
        }
        stats << setw(8) << statsmatrix[(k*(int) aggstatvec.size())+((int) aggstatvec.size()-1)] << endl;
    }

	ofstream istats;
	istats.open(outfilename + ".imag");
	//print headers
	istats << "#frame\t time \t norm \t\t <x1> \t\t <x^2> \t\t sqrt(var1) \t "
                "overlap1 \t\t <V1> \t\t <Lz> \t\t Energy1 \t\t norm2 \t\t\t <x2> \t\t "
                "sqrt(var2) \t\t overlap2 \t\t <V2> \t\t Energy2" << endl;
    for(k=0;k<(istatsmatrix.size()/aggstatvec.size());k++)
    {
        for(j=0;j<(int) aggstatvec.size()-1;j++)
        {
            istats << setw(8) << istatsmatrix[(k*(int) aggstatvec.size())+j] << "\t";
        }
        istats << setw(8) << istatsmatrix[(k*(int) aggstatvec.size())+((int) aggstatvec.size()-1)] << endl;
    }


    //dump densitytime matrix to file
    cout << "dumping densitytime to file" << endl;
    ofstream spacetimefile;
    spacetimefile.open(outdir + "/" + "densitytime.dat");
//    spacetimefile << setw(12);
    
    for(k=0;k<Nt;k++)
    {
        for(j=0;j<Nx-1;j++)
        {
            spacetimefile << setw(8) << densitytime[(k*Nx)+j] << "\t";
        }
        spacetimefile << setw(8) << densitytime[(k*Nx)+(Nx-1)] << endl;
    }
    
    ofstream spacetimefilelabelled;
    spacetimefilelabelled.open(outdir + "/" + "densitytime-labelled.dat");
    //first print x values as first row
    spacetimefilelabelled << "time" << "\t" ;
    for(j=0;j<Nx-1;j++)
    {
        spacetimefilelabelled << x[j] << "\t";
    }
    spacetimefilelabelled << x[Nx-1] << endl;
            
    // each row k is one timestep
    for(k=0;k<Nt-1;k++)
    {
        // print time in first column
        spacetimefilelabelled << (k)*tstep << "\t";
        
        // print density in remaining columns
        for(j=0;j<Nx-1;j++)
        {
            spacetimefilelabelled << densitytime[(k*Nx)+j] << "\t" ;
        }
    // Do the last k loop iteration without trailing newline
        spacetimefilelabelled << densitytime[(k*Nx)+(Nx-1)] << endl;
    }
    // print time in first column
    spacetimefilelabelled << (Nt-1)*tstep << "\t";
    for(j=0;j<Nx-1;j++)
    {
        spacetimefilelabelled << densitytime[((Nt-1)*Nx)+j] << "\t";
    }
        spacetimefilelabelled << densitytime[((Nt-1)*Nx)+(Nx-1)];

    //dump modeampmatrix to file in multiple formats
    //
    // todo: output headers to file! (requires reworking indices
    // current format: frequencies across a row, each row is one timestep
    
    ofstream modeampmatrixfile;

    //first: complex notation, row-wise
    modeampmatrixfile.open(outdir + "/" + "modeampmatrix-complex.dat");
    for(k=0;k<nfreq+1;k++)
    {// loop over frequencies (rows)
        for(j=0;j<Nt+1;j++)
        {// loop over timesteps (columns)
            if(k==0)
            {
                modeampmatrixfile << setw(8) << (j-1)*tstep << "\t";
            }
            else if(j==0)
            {
                modeampmatrixfile << setw(8) << fourierlist[k-1] << "\t";
            }
            else
            {
                    modeampmatrixfile << setw(8) <<  modeampmatrix[((k-1)*Nt)+(j-1)] << "\t";
            }
        }
                    modeampmatrixfile << "\n";
    }

    modeampmatrixfile.close();
    modeampmatrixfile.open(outdir + "/" + "modes.dat");

    //write headers
    modeampmatrixfile << "time";
    for(j=0;j<nfreq;j++)
    {
            modeampmatrixfile << setw(8) << "\t" << "Re(m=" << fourierlist[j]  << ")"
                    << "\t" << "Im(m=" << fourierlist[j] << ")"; 
            // adds a header for each frequency
    }
    modeampmatrixfile << endl; //newline, flush buffer

    for(k=0;k<Nt;k++)
    {// loop over timesteps (rows)
        for(j=0;j<nfreq;j++)
        {
            //write
            if(j==0)
            {
                modeampmatrixfile << setw(8) << setprecision(8) << k*tstep << "\t";
            }
            else
            {
                modeampmatrixfile << setw(8) << setprecision(8) <<
                    real(modeampmatrix[((j-1)*Nt)+k]) << "\t" <<
                    imag(modeampmatrix[((j-1)*Nt)+k]) << "\t";
                //modeampmatrixfile << setw(16) << abs(modeampmatrix[((j-1)*Nt)+k]) << "\t";
            }
        }
        modeampmatrixfile << setw(8) << setprecision(8) << 
               real(modeampmatrix[((nfreq-1)*Nt)+k]) <<
               imag(modeampmatrix[((nfreq-1)*Nt)+k]) << endl;
    }

    modeampmatrixfile.close();
    modeampmatrixfile.open(outdir + "/" + "modeampmatrix-realimag.dat");

    //second: real, imag columns
    for(k=0;k<nfreq+1;k++)
    {// loop over frequencies (rows)
        for(j=0;j<Nt+1;j++)
        {// loop over timesteps (columns)
            if(k==0)
            {
                modeampmatrixfile << (j-1)*tstep << "\t" << (j-1)*tstep << "\t";
            }
            else if(j==0)
            {
                modeampmatrixfile << fourierlist[k-1] << "\t" << fourierlist[k-1] << "\t";
            }
            else
            {
                    modeampmatrixfile << real(modeampmatrix[((k-1)*Nt)+(j-1)]) << "\t" << imag(modeampmatrix[((k-1)*Nt)+(j-1)]) << "\t";
            }
        }
                    modeampmatrixfile << "\n";
    }

    modeampmatrixfile.close();
    modeampmatrixfile.open(outdir + "/" + "modeampmatrix-real.dat");

    //third: real only
    for(k=0;k<nfreq+1;k++)
    {// loop over frequencies (rows)
        for(j=0;j<Nt+1;j++)
        {// loop over timesteps (columns)
            if(k==0)
            {
                modeampmatrixfile << (j-1)*tstep << "\t";
            }
            else if(j==0)
            {
                modeampmatrixfile << fourierlist[k-1] << "\t";
            }
            else
            {
                    modeampmatrixfile << real(modeampmatrix[((k-1)*Nt)+(j-1)]) << "\t";
            }
        }
                    modeampmatrixfile << "\n";
                    //modeampmatrixfile << real(modeampmatrix[((k-1)*Nt)+(Nt-1)]) << "\n";
    }

    modeampmatrixfile.close();
    modeampmatrixfile.open(outdir + "/" + "modeampmatrix-abs.dat");

    //fourth: magnitude
    for(k=0;k<nfreq+1;k++)
    {// loop over frequencies (rows)
        for(j=0;j<Nt+1;j++)
        {// loop over timesteps (columns)
            if(k==0)
            {
                modeampmatrixfile << (j-1)*tstep << "\t";
            }
            else if(j==0)
            {
                modeampmatrixfile << fourierlist[k-1] << "\t";
            }
            else
            {
                    modeampmatrixfile << abs(modeampmatrix[((k-1)*Nt)+(j-1)]) << "\t";
            }
        }
                    modeampmatrixfile << "\n";
                    //modeampmatrixfile << abs(modeampmatrix[((k-1)*Nt)+(Nt-1)]) << "\n";
    }

    //Absolute value
    modeampmatrixfile.close();
    modeampmatrixfile.open(outdir + "/" + "modes-abs.dat");

    modeampmatrixfile << "time";
    for(j=0;j<nfreq;j++)
    {
            modeampmatrixfile << "\t" << fourierlist[j]; 
            // adds a header for each frequency
    }
    modeampmatrixfile << endl; //newline, flush buffer

    for(k=0;k<Nt;k++)
    { // rows are timesteps
        for(j=0;j<nfreq;j++)
        { // columns are frequencies
            if(j==0)
            {   // write time to first column
                modeampmatrixfile << setw(8) << k*tstep << "\t";
            }
            else
            {
                modeampmatrixfile << setw(8) << abs(modeampmatrix[((j-1)*Nt)+k]) << "\t";
            }
        }
        modeampmatrixfile << "\n";
        //modeampmatrixfile << setw(8) << abs(modeampmatrix[((nfreq-1)*Nt)+k]) << endl;
    }

    modeampmatrixfile.close();


    ofstream phasext;
    phasext.open(outdir + "/" + "phasext.dat");
    
    phasext << setw(8) << "time"; // deliberately left hanging
    for(j=0;j<nfreq;j++)
    {
            phasext << setw(8) <<  "\t" << fourierlist[j]; // adds a header for each frequency
    }
    phasext << endl; //newline, flush buffer

    for(k=0;k<Nt;k++)
    {
        for(j=0;j<nfreq;j++)
        {
            if(j==0)
            {   // write time to first column
                phasext << setw(8) << k*tstep << "\t";
            }
            else
            {
                phasext << setw(8) << phasextmatrix[(k*nfreq)+(j-1)] << "\t";
            }
        }
        phasext << setw(8) << phasextmatrix[(k*nfreq)+(nfreq-1)] << endl;
    }

    //close phase file explicitly before opening it
    phasext.close();

    cout << "Now doing phase extraction" << endl;
    //load phase data, prepare to extract rotation using the gradient
    ifstream phasedata;
    loadmulticoldata(outdir + "/" + "phasext.dat",phasevec);
    gradvec.resize(phasevec.size());
    //extract rotation rate from phase shift
    extractrotation(phasevec,gradvec,fourierlist,tstep,threshold,Nt,nfreq,imm);
    //print rotation rates to file
    string columnhead = makeheader(fourierlist);
    multicolfp(outdir + "/" + "rotation.dat",gradvec,nfreq+1,columnhead);
    

    //load mode amplitude data to find stroboscopic time and rotation rate
    ifstream modedata;
    vector<double> modevec;
    vector<double> modesingle;
    vector<double> timevec;
    loadmulticoldata(outdir + "/" + "modes-abs.dat",modevec);
    col2vec(modevec,modesingle,nfreq+1,1,1);
    col2vec(modevec,timevec,nfreq+1,0,1);


    //determine final rotation rate result
    cout << "Taking stroboscopic rotation measurements" << endl;
    vector<double> rotationsingle;
    vector<double> phasesingle;
    vector<double> rotationmeasurement;
    col2vec(gradvec,rotationsingle,nfreq+1,1,0);
    col2vec(phasevec,phasesingle,nfreq+1,21,0);//include header in phasesingle so it can be extracted


    for(j=0;j<phasesingle.size();j++)
    {
        cout << phasesingle[j] << endl;
    }



    //rotationmeasurement = strobo_rot_measurement(phasesingle,timevec,tstep,tmax,xupper,xlower,omegareal,imm);
    rotationmeasurement = strobo_rot_measurement(phasesingle,timevec,tstep,tmax,xupper,xlower,g11,threemode_params[0],imm);


    ofstream measurement;
    measurement.open(outdir + "/rotation-measurement.dat");
    measurement << "time\t\trotation" << endl;

    for(j=0;j<rotationmeasurement.size();j++)
    {
        measurement << setw(8) 
            << timevec[j] << "\t"
            << rotationmeasurement[j] << endl;
    }

    /*
    vector<int> maxindices;
    vecmax(modesingle,maxindices);
cerr << "vecmax size: " << maxindices.size() << endl;

    ofstream measurement;
    measurement.open(outdir + "/rotation-measurement.dat");
    measurement << "time\t\trotation" << endl;

    double measurementresult;

    //correct rotation measurement calculation
    for(j=2;j<maxindices.size();j++)
    {
        if(phasesingle[maxindices[j]] < phasesingle[maxindices[j-2]])
        {
            measurementresult = ((phasesingle[maxindices[j]]+(2.0*PI))-phasesingle[maxindices[j-2]])
                                /((timevec[maxindices[j]]-timevec[maxindices[j-2]]));
        }
        else
        {
            measurementresult = (phasesingle[maxindices[j]]-phasesingle[maxindices[j-2]])
                               /((timevec[maxindices[j]]-timevec[maxindices[j-2]]));
        }

        rotationmeasurement.push_back(measurementresult);
    }



    for(j=0;j<rotationmeasurement.size();j++)
    {
        measurement << setw(8) 
            << timevec[maxindices[j]] << "\t"
            << rotationmeasurement[j]/(PI*(double) imm) << endl;
                                    // ^factor of PI is important
    }

    */


/*
 * Also need to implement: Fisher info
 * Need to investigate: Lz(t=0) vs Omgea; Omega vs quanta of ang. mom.
 */

    cout << "Rotation measurements written to file" << endl;

    cout << "Program executed succesfully." << endl;

	// done
	exit(EXIT_SUCCESS);
}

