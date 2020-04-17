/*
	LAST UPDATED: 12/13/2015

	Particle swarm optimization for measured aerosol optical depths at 6 wavelengths:
		440.699, 555.803, 674.855, 778.221, 869.738, 1019.721

	This prgram processes -aod files created by aod-sorter.  These six wavelengths will be fitted to a bimodal distribution.
*/

// c++ includes
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <new>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstdio>

using namespace std;		// required for c++ file routines

#define channels 6 		    // number of wavelength channels to use for optimization
#define particles 50		// number of particles in the swarm
#define PI 3.1415192654		// the constant PI
#define CONV_FACT 1.0E-4	// convert microns to cm
#define CONV_FACT1 1.0E-3	// convert nanometers to microns
#define CONV_FACT2 1.0E-7	// convert nanometers to cm
#define c1 1.8 				// acceleration coefficient 1
#define c2 1.8 				// acceleration coefficient 2
#define gen_max 3000 		// define the maximum number of generations to go through
#define N_tmax 100			// number of iterations to check change N0, a, b that satisfy accuracy goal
#define num_scouts 1000		// number of scouts to generate
#define pop_max 21			// maximum amount of populations (populations start at 1)

// required definition of matrix at top of source code
typedef vector<vector<double> > matrix; // define matrix object

// generates a random value between low and high numbers using log10 space
double random(double low, double high);
// calculates aod for each particle in the swarm -- call with calc_aod(swarm[x], mie[y], aodc[x][y])
void calc_aod(vector<double> &part, vector<double> &Qext, double &tauc, int distribution);
// initialize mie extinction matrix by filling it with values from files
void initQ(vector<double> &Qext, string fname);
// read in and initialize the mie files
void init_mie(matrix &mie);
// calculate the fitness value for each particle in the swarm
void calc_fitness(double aodm[], matrix &aodc, double f[]);
// initialize the swarm by passing by passing range limits for parameters
void init_swarm(matrix &swarm, double param_min[], double param_max[]);
// generates a random particle using log10 space for parameter ranges - called by init_swarm()
void init_scout(matrix &scouts, double param_min[], double param_max[]);
// initialize a particle for the swarm or scout population
void init_particle(vector<double> &part, double param_min[], double param_max[]);
// find global fitness minimum
int find_global_min(double f[], double Pg[], matrix &swarm);
// find local fitness minimum for each particle
void find_local_min(double f[], double f0[], matrix &Pl, matrix &swarm);
// adjust the swarm based off the calculated
void update_swarm(matrix &swarm, matrix &Pl, double Pg[], double param_min[], double param_max[]);
// find scout with minimum fitness function
int find_scout_min(double fj[], matrix &scouts);

int main()
{
	int i,j, N_t=0, distribution, population=1, fitness_goal; 	// generic indices for loops
	int i_g, i_j, gen=0, a, size, linesize=200; 				// indices for global minimum and generation number
	float day;
	string mystr, line; 				// holds typed in filenames
	char *fname1,*fname2,dataline[500]; // holds input and output filenames and datalines from input file
	double aodm[channels]; 				// holds measured aerosol optical depths (AODs) for a particular date & time
	double f0[particles]; 				// initial fitness vector
	double f[particles]; 				// current fitness vector
	matrix swarm(particles,vector<double>(3,0.0)); // matrix of particles each with 3 dimensions N0, a, b
	matrix aodc(particles,vector<double>(channels,0.0)); // matrix of calculated AODs for each particle at each wavelength
	matrix Pl(particles,vector<double>(3,0.0)); // matrix of local minimum adjustments
	double Pg[3]; 						// vector for global minimum adjustments for all particles
	double Pg0[3] = {1.0E9,0.2,0.5}; 	// initial global minima for checking stop conditions
	double Pj[3];						// vector for global minimum for scout particles
	matrix scouts(num_scouts,vector<double>(3,0.0)); // matrix for randomly generated particle so we don't get stuck in a local minimum - use vector object to be compatible with swarm for calls to init_particle
	matrix aodj(num_scouts, vector<double>(channels,0.0)); //vector for extra particle's aods for global minimum searching improvement
	double fj[num_scouts]; 				// fitness for extra randomly generated particle
	matrix mie(channels,vector<double>(9961,0.0)); // matrix of extinction coefficients for each wavelength at each radius
	vector<double> param(3,0.0); 		// holds initial set of N0, alpha, and beta in order to get initial fitness functions

	float dtime;						// define dtime as an array for each point

	double param_min[3], param_max[3];	//minima and maxima for N0, alpha, beta
	fstream aodfile,outfile;

	srand(time(NULL)); 					// seed the random number algorithm
	init_mie(mie);						// initialize mie files

	cout << "Enter measured AOD file to process (e.g., gYYYYMMDD-aod.csv): ";
	getline(cin,mystr);					// read in input
	fname1 = new char [mystr.size()+1];	// c_string method
	strcpy(fname1,mystr.c_str());

	cout << "Enter output file name for best number distribution parameters N0, alpha, beta (e.g., YYYYMMDD-distribution.csv): ";
	getline(cin,mystr);					// read in input
	fname2 = new char [mystr.size()+1];	// c_string method
	strcpy(fname2,mystr.c_str());

	cout << "Enter the distribution you want (1=Junge,  2=Gamma,  3=Log-Normal): ";
	cin >> distribution;				// read in input

	aodfile.open(fname1);						// open input file of aod values
	if (!aodfile.is_open()){					// check to make sure aodfile is open
		cout << "\nUnable to open file " << fname1 << "\n";
		return 0;								// if not, quit program
	}
	outfile.open(fname2,ios::app | ios::out);	// open output file for number distribution parameters N0, alpha, beta
	if (!outfile.is_open()){					// check to make sure outfile is open
		cout << "\nUnable to open file " << fname2 << "\n";
		return 0;								// if not, quit program
	}

	a=0;										// initialize read in index
	while(getline(aodfile,mystr)){				// get data from input file
		if(a>0){								// read in the rest of the file
				cout << "a is: " << a << endl;
				stringstream(mystr) >> dtime >> aodm[0] >> aodm[1] >> aodm[2] >> aodm[3] >> aodm[4] >> aodm[5];		// read in measured aod values
				cout << mystr << endl;

				cout << "dtime is: " << dtime << endl;
				outfile << "For time: " << dtime << endl;
				for(i=0;i<channels;i++){		// output measured aod values
						cout << "aodm["<<i<<"] = " << aodm[i] << endl;
				}

//		param[0]=1.0E10; 			// initial value for N0
		param_min[0]=1.0E6;
		param_max[0]=1.0E12;
//		param[1]=0.1; 				// initial value for alpha (a)
		param_min[1]=0.01;
		param_max[1]=1.0;
//		param[2]=0.5; 				// initial value for beta (b)
		param_min[2]=0.01;
		param_max[2]=1.0;

		if(distribution==1){
			fitness_goal=0.1;		// set fitness goal for Junge Distribution
		}
		else if(distribution==2){
			fitness_goal=1.0E-5;	// set fitness goal for Gamma Distribution
		}
		else if(distribution==3){
			fitness_goal=1.0E-1;	// set fitness goal for Log-Normal Distribution
		}


		do{
			init_swarm(swarm, param_min, param_max);	// initialize swarm
			for(i=0;i<particles;i++){
				for(j=0;j<3;j++){ 						// initilaize local minima
					Pl[i][j]=swarm[i][j]; 				// initialize parameters for local minimum
				}
			}

// do initial aod calculations to get initial fitness values
			cout << "initializing aods..." << endl;					// for testing

			for(i=0;i<particles;i++){ 								// cycle through particles
				for(j=0;j<channels;j++){ 							// cycle through wavelengths
					calc_aod(param,mie[j],aodc[i][j],distribution); // all initial aods calculated with first guess N0, a, b
				}
			}

	cout << "calculating initial fitness values..." << endl;	// for testing

			calc_fitness(aodm, aodc, f0); 						// find initial fitness values
			for(i=0;i<particles;i++){
				//outfile << "Initial Fitness f0["<<i<<"]: " << f0[i] << endl;
				//cout << "Initial Fitness f0["<<i<<"]: " << f0[i] << endl;
			}
			init_swarm(swarm,param_min,param_max); 				// initialize the swarm

			while(gen<gen_max){ 								// let the swarm evolve
				cout << gen << "\t" << N_t << "\t" << endl;
				for(i=0;i<particles;i++){ 						// cycle through particles
					for(j=0;j<channels;j++){ 					// cycle through wavelengths
						calc_aod(swarm[i],mie[j],aodc[i][j],distribution); // aods calculated with individual particles' N0, a, b
					}
				}
				calc_fitness(aodm, aodc, f); 					// find fitness values
				find_local_min(f, f0, Pl, swarm);				// find local minima
				i_g = find_global_min(f, Pg, swarm); 			// find global minimum
				for(i=0;i<3;i++){
					Pg[i]=swarm[i_g][i];						// redefine global minimum parameters
				}
				if(gen<10){
					init_scout(scouts, param_min, param_max); 	// generate scout particle for improving global minimum searching
					for(i=0;i<num_scouts;i++){					// cycle through particles
						for(j=0;j<channels;j++){ 				// cycle through wavelengths
							calc_aod(scouts[i],mie[j],aodj[i][j],distribution); // aods calculated with individual particles' N0, a, b
						}
					}
					calc_fitness(aodm,aodj,fj); 				// fitness for extra scout particle
					i_j = find_scout_min(fj,scouts);			// integer of best scout
					if(fj[i_j]<f[i_g]){
						fj[i_j]=f[i_g];							// if best scout's fitness is better than global best, replace with scout
						cout << "Scouts made a difference" << endl;
					}
				}
				update_swarm(swarm, Pl, Pg, param_min, param_max);	// update swarm based off equation

				if (abs(Pg[0]-Pg0[0])/Pg0[0] < 1.0E-7 && abs(Pg[1]-Pg0[1])/Pg0[1] < 1.0E-7 && abs(Pg[2]-Pg0[2])/Pg0[2] < 1.0E-7){
					N_t+=1; 			// if no significant change update number of trials for no change
				}
				else {
					N_t=0; 				// if significant change reset number of trials
				}
				if(N_t>N_tmax){ 		// stop condition reached?
					cout << "N0: " << (Pg[0]-Pg0[0])/Pg0[0] << "  a: " << (Pg[1]-Pg0[1])/Pg0[1] << "  b: " << (Pg[2]-Pg0[2])/Pg0[2] << endl;
					break; 				// if stop condition met, break out of evolution loop
				}
				else{					// if stop condition not met, update global minimum and continue
					for(i=0;i<3;i++){ 	// cycle through parameters N0, a, b
						Pg0[i]=Pg[i]; 	// to update global min
					}
				}
				gen++;
			} // end of evolution loop
			outfile << "N0: " << Pg[0] << " a: " << Pg[1] << " b: "<< Pg[2] << endl;
			outfile << "Generations: " << gen << endl;
			outfile << "Fitness: " << f[i_g] << endl;
			cout << "Generation # " << gen << ", N0 = " << Pg[0] << ", alpha = " << Pg[1] << ", beta = " << Pg[2] << endl;
			for(i=0;i<channels;i++){
					cout << "Measured: " << aodm[i] << "\t" << "  Calculated: " << aodc[i_g][i] << endl;
					outfile << "Measured: " << aodm[i] << "\t" << "  Calculated: " << aodc[i_g][i] << endl;
			}
			outfile << endl;
			cout << "Fitness = " << f[i_g] << endl;
			cout << "Information written to: " << fname2 << endl;
				N_t=0;				// if fitness goal not reached, reset N_t
				population++;		// if fitness goal not reached, increase population number
				gen=0;				// if fitness goal not reached, create new generation
			cout << "Population is: " << population << endl << endl;
			outfile << "Population is: " << population << endl;
		} while(f[i_g]>fitness_goal && population<pop_max);		// end of do while loop in populations feature
		if(population==pop_max){	// if population max is reached, quit program
			cout << "Did not converge after " << population << " populations" <<  endl;
			outfile << "Did not converge after " << population << " populations" << endl;
			for(i=0;i<channels;i++){
				outfile << "Measured: " << aodm[i] << "\t\t" << "Calculated: " << aodc[i_g][i] << endl;
			}
			outfile << endl << endl;
			population=1;			// reset population
			a++;					// go to next line
		}
		if(f[i_g]<fitness_goal){	// if fitness goal reached
			cout << "Fitness goal reached at: " << f[i_g] << endl;
			outfile << "Fitness goal reached at: " << f[i_g] << endl;
			a++;					// go to next line in input file
			outfile << "a is: " << a << endl << endl << endl;
			for(i=0;i<channels;i++){
				outfile << "Measured: " << aodm[i] << "\t\t" << "Calculated: " << aodc[i_g][i] << endl;
			}
		}
	}								// end of if loop that reads in each line of input file
a++;								// go to next line
}									// end of while(getline,mystr) loop
	cout << "closing file " << fname1 << endl;
	aodfile.close();				// close input file
	outfile << endl;
	outfile.close();				// close output file
	return 0;
}									// end of main program

double random(double low, double high)	// generate random numbers
{
	double result,range;
	range = log10(high/low); 		// log(high)-log(low)=log(high/low)
	result = log10(low) + (rand() % 1000000000)/1.0E9*range;
	return pow(10,result);
}

void calc_aod(vector<double> &part, vector<double> &Qext, double &tauc, int distribution)
{
	int i=0;
	tauc=0;
	while(i<9962){ // cycle through all radii
		if(distribution==1){
			tauc += pow((i+40)*CONV_FACT2,2.0)*Qext[i]*pow((i+40.5)*CONV_FACT2,-part[1]); // Junge distribution
		}
		else if(distribution==2){
			tauc += pow((i+40)*CONV_FACT2,2.0)*Qext[i]*pow((i+40.5)*CONV_FACT2,part[2])*exp(-part[1]*(i+40.5)*CONV_FACT2); //Gamma distribution
		}
		else if(distribution==3){
			tauc += pow((i+40)*CONV_FACT2,2.0)*Qext[i]*exp((-pow(log10(((i+40.5)*CONV_FACT1)/part[1]),2.0)/(2*pow(part[2],2.0)))); //Log-Normal Distribution
		}
		else{
			cout << "Unable to process distribution" << endl;
			return;
		}
		i++;
	}
		if(distribution==1){									// if Junge distribution
			tauc *= 2.0*PI*5.0E-8*part[0]*part[2]; 				// multiply by the constants from the sum: 2*pi, delta r, N0, beta
		}
		else if(distribution==2){								// if Gamma distribution
			tauc *= 2.0*PI*5.0E-8*part[0]; 						// multiply by the constants from the sum: 2*pi, delta r, N0
		}
		else if(distribution==3){								// use with log-normal distribution
			tauc *= sqrt(2.0*PI)*5.0E-8*part[0]*1/(part[2]); 	// multiply by the constants from the sum: 2*pi, delta r, N0, beta, 1/sqrt(2*pi)
		}
	return;
}


void init_mie(matrix &mie)	// fills the matrix mie[wavelength][radius]
{
	string fname;
	fname="mie-data/mie440699nm.csv";		// wavelength is 440.343, couldn't add period in file name
	initQ(mie[0],fname);
	fname="mie-data/mie555803nm.csv";		// wavelength is 674.855, couldn't add period in file name
	initQ(mie[1],fname);
	fname="mie-data/mie674855nm.csv";		// wavelength is 869.738, couldn't add period in file name
	initQ(mie[2],fname);
	fname="mie-data/mie778221nm.csv";		// wavelength is 1019.721, couldn't add period in file name
	initQ(mie[3],fname);
	fname="mie-data/mie869738nm.csv";		// wavelength is 869.738, couldn't add period in file name
	initQ(mie[4],fname);
	fname="mie-data/mie1019721nm.csv";		// wavelength is 1019.721, couldn't add period in file name
	initQ(mie[5],fname);
	return;
}

void initQ(vector<double> &Qext, string fname)	// fills an individual row (i.e., wavelength) of mie[wavelength][radius]
{
	int i;
	double r, Qsca, Qabs, Gsca;	// data in Mie files but not needed
	char *fname1,dataline[100];
	string mystr;
	fstream miefile;
	fname1 = new char [fname.size()+1];	/* c_string method */
	strcpy(fname1,fname.c_str());

	miefile.open(fname1);
	if (miefile.is_open())
	{
		i=0;
		getline(miefile,mystr); // throw away header line
		while(i<9962)
		{
			getline(miefile,mystr);
			strcpy(dataline,mystr.c_str());
			sscanf(dataline,"%lf,%lf,%lf,%lf,%lf",&r,&Qext[i],&Qsca,&Qabs,&Gsca); // fills the Qext[] array ignoring other values
			i++;
		}
		miefile.close();
		cout << "Processed mie file " << fname << "\n";
	}
	else cout << "\nUnable to open file " << fname << "\n";

	return;
}

void calc_fitness(double aodm[], matrix &aodc, double f[])
{
	int i=0,j=0,rows;
	double a;
	rows=aodc.size();
	for(i=0;i<rows;i++){ 						// cycle through all particles
		a=0.0;
		for(j=0;j<channels;j++){ 				// cycle through all wavelengths
			a += pow(aodm[j]-aodc[i][j],2.0); 	// sum of square deviations
		}
		f[i]=sqrt(a/6.0); 						// square root of sum over number of wavelengths
	}
	return;
}

void init_swarm(matrix &swarm, double param_min[], double param_max[])
{
	int i=0;
	int j;
	while(i<particles){ 								// cycle through all particles
		init_particle(swarm[i], param_min, param_max); 	// generate parameters
		i++;
	}
	return;
}

void init_particle(vector<double> &part, double param_min[], double param_max[])
{
	part[0]=random(param_min[0],param_max[0]); // set N0 parameter
	part[1]=random(param_min[1],param_max[1]); // set alpha (a) parameter
	part[2]=random(param_min[2],param_max[2]); // set beta (b) parameter
	return;
}

int find_global_min(double f[], double Pg[], matrix &swarm)
{
	int i;
	int i_g=0;
	for(i=0;i<particles;i++){ 			// cycle through all particles
		if (f[i]<f[i_g]) { 				// look for lowest fitness value
			i_g = i; 					// set i_g to index of lowest fitness value
		}
	}
	return i_g; 						// return index of lowest fitness value
}

void find_local_min(double f[], double f0[], matrix &Pl, matrix &swarm)
{
	int i,j;
	for(i=0;i<particles;i++){ 			// cycle through all particles
		if(f[i]<f0[i]){ 				// test if current fitness is less than previous lowest fitness
			f0[i]=f[i]; 				// if so, reset the lowest fitness
			for(j=0;j<3;j++){ 			// reset the local minimum
				Pl[i][j]=swarm[i][j];
			}
			j=0;
		}
	}
	return;
}

void update_swarm(matrix &swarm, matrix &Pl, double Pg[], double param_min[], double param_max[])
/*
Implements Eq (10) from Yuan's article:
X_i(t+1)=X_i(t)+c1*r1*[P_i(t)-X_i(t)]+c2*r2*[P_g(t)-X_i(t)]

Not sure if r1 & r2 should only change for each update or if they should change for each particle or for each particle's parameters
*/
{
	int i,j; 								// particle and parameter indices
	double r1,r2; 							// random numbers in interval [0,1]
	i=0;
	while(i<particles){
		j=0;
			r1 = (rand() % 1000000000)/1.0E9; // each particle gets a new random number
			r2 = (rand() % 1000000000)/1.0E9;	// different random numbers for r1 and r2
		while(j<3){
			swarm[i][j] = (1.0-c1*r1-c2*r2)*swarm[i][j]+c1*r1*Pl[i][j]+c2*r2*Pg[j];
			if (swarm[i][j] < param_min[j]){ // check to see if below lower bound
//				swarm[i][j] = (param_min[j]+param_max[j])/2.0; //if so, set back to average
//				swarm[i][j] = param_max[j]; //torroidal boundary conditions
//				swarm[i][j] = param[j]; 	//set it back to initial conditions
				init_particle(swarm[i], param_min, param_max);	//initializes new random particle
			}
			else if (swarm[i][j] > param_max[j]){ // check to see if above upper bound
//				swarm[i][j] = (param_max[j]+param_max[j])/2.0; //if so, set to average
//				swarm[i][j] = param_min[j];	//torroidal boundary conditions
//				swarm[i][j] = param[j]; 	//set it back to initial conditions
				init_particle(swarm[i], param_min, param_max);	//initialized new random particle
			}
			j++;
		}
		i++;
	}
	return;
}

void init_scout(matrix &scouts, double param_min[], double param_max[])
{
	int i=0;
	int j;
	while(i<num_scouts){ 								// cycle through all particles
		init_particle(scouts[i], param_min, param_max); // generate parameters
		i++;
	}
	return;
}

int find_scout_min(double fj[], matrix &scouts)
{
	int i;
	int i_j=0;
	for(i=0;i<particles;i++){ 	// cycle through all particles
		if (fj[i]<fj[i_j]) { 	// look for lowest fitness value
			i_j = i; 			// set i_g to index of lowest fitness value
		}
	}
	return i_j; 				// return index of lowest fitness value
}
