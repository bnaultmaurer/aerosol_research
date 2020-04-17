// Program to create artificial aerosol optical depths at given N0, alpha, Beta
// For use in testing particle swarm optimization program
// 13 wavelengths:412, 440, 463, 479, 500, 520, 556, 610, 675, 750, 778, 870, and 1020 nm
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

using namespace std;	// required for c++ file routines

#define wavelengths 13
#define CONV_FACT2 1.0E-7	/* convert nanometers to cm */
#define PI 3.1415192654		// the constant PI


typedef vector<vector<double> > matrix; // define matrix object

void calc_aod(vector<double> &part, vector<double> &Qext, double &tauc);
void init_mie(matrix &mie);
// fills the matrix mie[wavelength][radius]
void initQ(vector<double> &Qext, string fname);


int main()
{
	char *fname2;
	string mystr;
	fstream outfile;
	int i, time;
	double aod[wavelengths];
	matrix mie(wavelengths,vector<double>(9961,0.0)); // matrix of extinction coefficients for each wavelength at each radius
	vector<double> part(3,0.0);	//holding three elements, No, a, b
	struct tm now;
	float daynum;
	
	cout << "This program generates artifical AOD values, to be used to test particle swarm optimization" << endl << endl;

	init_mie(mie);

	cout << "Enter time (e.g. hhmm): ";
	cin >> time;

	daynum=(time/100+(time%100)/60.0)/24.0;

	cout << "Enter N0 value: ";
	cin >> part[0];

	cout << "Enter a value: ";
	cin >> part[1];

	cout << "Enter b value: ";
	cin >> part[2];

	cout << "Enter output file name for AOD (e.g., sample.txt): ";
	cin >> mystr;
	fname2 = new char [mystr.size()+1];	/* c_string method */
	strcpy(fname2,mystr.c_str());
	
	outfile.open(fname2,ios::out);	// open output file for AOd at 13 wavelengths 
	if (!outfile.is_open()){
		cout << "\nUnable to open file " << fname2 << "\n";
		return 0;
	}

//	daynum=(now.tm_hour+now.tm_min/60.0+now.tm_sec/3600.0)/24.0;
	outfile << daynum << " ";


	for( i=0; i<wavelengths; i++){
		calc_aod(part,mie[i],aod[i]);
		outfile << aod[i] << " ";
	}
	outfile << endl;

	
	return 0;
}	

void calc_aod(vector<double> &part, vector<double> &Qext, double &tauc)
{
	int i=0;
	tauc=0;
	while(i<9962){ // cycle through all radii
//		tauc += pow((i+40)*CONV_FACT2,2.0)*Qext[i]*pow((i+40.5)*CONV_FACT2,-part[1]); // Junge distribution
//		tauc += pow((i+40)*CONV_FACT2,2.0)*Qext[i]*pow((i+40.5)*CONV_FACT2,part[2])*exp(-part[1]*(i+40.5)*CONV_FACT2); //Gamma distribution
		tauc += pow((i+40)*CONV_FACT2,2.0)*Qext[i]*exp((-log10(((i+40.5)*CONV_FACT2)/part[1]),2.0)/(2*pow(part[2],2.0))); //Log-Normal Distribution
		i++;
	}
	//use with Junge distribution
//	tauc *= 2.0*PI*5.0E-8*part[0]*part[2]; // multiply by the constants from the sum: 2*pi, delta r, N0, beta
	//Use with Gamma distribution
//	tauc *= 2.0*PI*5.0E-8*part[0]; // multiply by the constants from the sum: 2*pi, delta r, N0
	//use with log-normal distribution
	tauc *= 2.0*PI*5.0E-8*part[0]*1/(part[2]*sqrt(2*PI)); // multiply by the constants from the sum: 2*pi, delta r, N0, beta, 1/sqrt(2*pi)
	return;
	return;
}

void init_mie(matrix &mie)
// fills the matrix mie[wavelength][radius]
{
	string fname;
	fname="mie-data/mie412nm.csv";
	initQ(mie[0],fname);
	fname="mie-data/mie440nm.csv";
	initQ(mie[1],fname);
	fname="mie-data/mie463nm.csv";
	initQ(mie[2],fname);
	fname="mie-data/mie479nm.csv";
	initQ(mie[3],fname);
	fname="mie-data/mie500nm.csv";
	initQ(mie[4],fname);
	fname="mie-data/mie520nm.csv";
	initQ(mie[5],fname);
	fname="mie-data/mie556nm.csv";
	initQ(mie[6],fname);
	fname="mie-data/mie610nm.csv";
	initQ(mie[7],fname);
	fname="mie-data/mie675nm.csv";
	initQ(mie[8],fname);
	fname="mie-data/mie750nm.csv";
	initQ(mie[9],fname);
	fname="mie-data/mie778nm.csv";
	initQ(mie[10],fname);
	fname="mie-data/mie870nm.csv";
	initQ(mie[11],fname);
	fname="mie-data/mie1020nm.csv";
	initQ(mie[12],fname);
	return;
}

void initQ(vector<double> &Qext, string fname)
// fills an individual row (i.e., wavelength) of mie[wavelength][radius]
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


