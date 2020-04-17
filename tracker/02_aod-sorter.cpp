/*
	LAST UPDATED 12/13/2015

	This program enters a file of sorted count values (gYYYYMMDD-sorted.csv, sorted by tracker-data-sorter.cpp)
	It opens the files, and takes the wavelength's counts at a given time, and outputs aod values for
	that time.

	This program creates a file gYYYYMMDD-aod.csv.  This file has the aod's in columns according to
	wavelength, and rows of each wavelength's aod for the given time.  This can be read into pso-tracker6.cpp
	to calculate size distributions.

	Six wavelengths: 440.699, 555.803, 674.855, 778.221, 869.738, 1019.721 (in nm)
*/

#include <fstream>		// header providing file stream classes
#include <sstream>		// header providing string stream classes
#include <iostream>		// header that defines the standard input/output stream objects
#include <string>		// introduces string types, character traits and qa set of converting functions
#include <errno.h>		// C header that defines errno (last error number)
#include <dirent.h>		// C header that defines dirent structure in main program
#include <sys/stat.h>	// ---------------------------------------------------------------------
#include <sys/types.h>	// ---------------------------------------------------------------------
#include <iomanip>		// provides parametric parameters, see http://www.cplusplus.com/reference/iomanip/
#include <cstring>		// manipulates C strings and arrays
#include <cstdlib>		// ---------------------------------------------------------------------
#include <math.h>		// declares a set of function to compute common mathematical operations and transformations
#include <new>
#include <vector>

using namespace std;
// required definition of matrix at top of source code
typedef vector<vector<float> > matrix;

void calc_aod(string input, double *wavelength, double *irr, double dtime, double *aod, int a, int linesize, float pressure, int i);	//function that does all aod calculations
void convert_time(string line, double &dtime, int linesize, int i);
void convert_counts(double *irr, string line, int linesize);
void pso(double *dtime, double *aod, double wavelength, int linesize, int b);

const int size = 6;				// defines size of array for the 6 wavelengths
int Co2 = 403;						// initial C02 density in the atm in ppm (last updated 7-31-2015)
int a, b;							// integer to be used as counter in to seperate wavelength's to their files
float PI = 3.14159265359;			// defining Pi
float latitude = 45.5807*PI/180;	// latitude of SJU (radians)
float stdpressure = 1013.25;		// standard pressure at sea level (millibars)
float height = 1200/3.2808;			// height above sea level of SJU (meters)
float A = 6.0221367E23;				// Avagadro's number
float Ns = 2.546899E19;				// molecular density of air
float pressure;						// initialize pressure as the user inputs
double wavelength[size] = {440.699, 555.803, 674.855, 778.221, 869.738, 1019.721};


int main()
{
	int i, j, linesize=6, pos;						// linesize now is the amount of wavelengths we are using
	double aod[linesize], irr[linesize];			// initalize wavelength, used in functions as definition of wavelength (368.0)
	double dtime;
	ifstream infile;								// set file to read info from
	ofstream outfile, pso;							// set name for output file
	string time, filepath, line, file, output, input, outputaod, outputpso;		// initalizes strings for directory creation
	char *fname1, *fname2;

	matrix print(7,vector<float>(linesize,0.0));	// matrix for printing to output file
	matrix counts(6,vector<float>(linesize,0.0));	// matrix for printing to output file

	cout << "Enter file to process (gYYYYMMDD-sorted.csv): ";
	getline(cin, input);							// get string
	fname1 = new char [input.size()+1];				// change string to char variable
	strcpy(fname1, input.c_str());					// change input to a c string
	pos = input.find("-");							// find end of file name
	output = input.substr(0,pos) + "-aod.csv";		// define output name as input name with -aod added

	cout << "Enter pressure for this day (mb): ";	// prompts user to input a pressure for given day
	cin >> pressure;								// defines pressure as what the user inputted

	infile.open(input.c_str());						// set input file name
	outfile.open(output.c_str()/*, ios::app*/);		// set output file name and create file // for all runs
	pso.open(outputpso.c_str());					// open output file for pso program

	if(!infile.is_open()){							// CHECK, to see if file is open
		cout << "Unable to open input file\n"; 		// CHECK, only displays if output files not opened
		return 0;
	}

	if(!outfile.is_open()){							// CHECK, to see if file is open
		cout << "Unable to open output file\n"; 	// CHECK, only displays if output files not opened
		return 0;
	}

	outfile << "dtime,440.699,555.803,674.855,778.221,869.738,1019.721" << endl;

	a=0;
	i=0;
	while(getline(infile,line)){					// while the file is open...
		if(a>0){									// skips header lines
		convert_time(line, dtime, linesize, i);		// extract time and return dtime
		outfile << dtime;							// print dtime to the file
		convert_counts(irr,line,linesize);			// extract count values and return irradience values
cout << "dtime in main is: " << dtime << endl;		// FOR TESTING ONLY
		calc_aod(input, wavelength, irr, dtime, aod, a, linesize, pressure, i);	// take irr values and calculate aod for each wavelength
		for(i=0;i<linesize;i++){					// for each aod value
			outfile << "\t" << aod[i];				// tab, then print aod value to output file
		}
		outfile << endl;							// go to next line

cout << endl << endl;								// FOR TESTING ONLY
		i++;										// go to next integer
		}				// end of if(a>1) statement
	a++;				// go to next line
	}					// end of while loop
	i=0;
	j=1;

	infile.close();		// close input file
	outfile.close();	// close output file
	pso.close();		// close pso output file
	cout << "Data saved into: " << output << endl << endl;
	cout << "Files processed successfully. Now go have fun with all those AOD values!!" << endl << endl;
	return 0;			// end of main program
}

void convert_time(string line, double &dtime, int linesize, int i){

	int pos1, pos2, pos3, hour, min, time;
	string stime, shour, smin;
	pos1 = line.find_first_of(",");							// find position at the end of the time
	stime = line.substr(0,pos1);							// write the time to appropriate time array
	pos2 = stime.find_first_of(":");						// find position of fir semi colon
	shour = stime.substr(0,pos2);							// take out shour
	smin = stime.substr(pos2+1,2);							// take out substr smin
	stime = shour.append(smin);								// write time in HHMM format
	time = atoi(stime.c_str());								// convert from string to int
	dtime = time/100 + (double)(time%100)/60;				// converts time value to decimal time
	return;
}



void convert_counts(double *irr, string line, int linesize){
	int i;
	int counts[linesize], pos1[linesize], pos2[linesize], pos3[linesize];
	string header, scounts[linesize];
	float conversion[size] = {0.0565244072, 0.0487698567, 0.063958183, 0.0320913325, 0.0292381904, 0.0482105824}; // calibration constants for spectrometer for 440.699, 674.855, 869.738, 1019.721
			line.erase(0,9);											// get rid of time in the line
	for(i=0;i<linesize;i++){											// for each wavelength's counts
			pos2[i] = line.find(",");									// find comma
			scounts[i] = line.substr(0,pos2[i]);						// extract counts values
			counts[i] = atoi(scounts[i].c_str());						// convert from string to int
			irr[i] = ((double)counts[i]*(double)conversion[i]);	// convert counts to irr using calibration constant
			line.erase(0,pos2[i]+1);									// erase counts values past comma

		}
	return;
}

void calc_aod(string input, double *wavelength, double *irr, double dtime, double *aod, int a, int linesize, float pressure, int i)		//calculates slew of data for each measurement, including year, month, day, time, checks for leap year, day number, decimal time, the m variable (optical path length), m' (secant of solar zenih angle), declination angle, RSE (Sun-Earth radius), aerosol optical depth, and the rayleigh factor.
{

// defines year, month, day, time
	int year, month, day, daynumber, time, pos1[linesize], pos2[linesize];
	string yearstr, monthstr, daystr, timestr, header;
	bool leap;

// converts timestr to integer time
	yearstr=input.substr(1,4);			// define year as string
	monthstr=input.substr(5,2);			// define month as string
	daystr=input.substr(7,2);			// define day as string
	year=atoi(yearstr.c_str());			// convert year from string to integer
	month=atoi(monthstr.c_str());		// convert month from string to integer
	day=atoi(daystr.c_str());			// convert day from string to integer

// leap year checker
	if((year%400)==0) {leap=1;}			// if year is divisible by 400 with no remainder, it is a leap year
	else if ((year%100)==0) {leap=0;}	// if year is divisible by 100 with no remainder, it is not a leap year
	else if ((year%4)==0) {leap=1;}		// if year is divisible by 4 with no remainder, it is a leap year
	else {leap=0;}						// otherwise read false, it is not a leap year

// calculates day number
	if(month== 1)						// if it is the first month of the year, it is simply the daynumber
		{daynumber= day;}				// set the daynumber as the day
	else if(month== 2)					// if it is the second month
		{daynumber= day+31;}			// add 31 days to the date
	else if(month== 3)					// if it is the third month
		{daynumber= day+59+leap;}		// add 59 days if it is a leap year
	else if(month== 4)
		{daynumber= day+90+leap;}
	else if(month== 5)
		{daynumber= day+120+leap;}
	else if(month== 6)
		{daynumber= day+151+leap;}
	else if(month== 7)
		{daynumber= day+181+leap;}
	else if(month== 8)
		{daynumber= day+212+leap;}
	else if(month== 9)
		{daynumber= day+243+leap;}
	else if(month== 10)
		{daynumber= day+273+leap;}
	else if(month== 11)
		{daynumber= day+304+leap;}
	else if(month== 12)
		{daynumber= day+334+leap;
		cout << daynumber << endl;}
	else		//CHECK
		cout << "Month value out of range" << endl;		//CHECK, only prints if error occurs with month value

for(i=0;i<linesize;i++){
//calculates hourangle(h) to use to find m'
	float hourangle;
	hourangle = (15*PI/180)*(12-dtime);					//defines hourangle from decimal time

//calculates declination
	float declination;
	declination = -23.4*PI/180 * (float)cos ( 2*PI*(daynumber+10)/365 );	//calculates declination angle
//cout << "declination is: " << declination << endl;

//calculates mprime
	float mprime;
	mprime = 1/(sin(latitude)*sin(declination) + cos(latitude)*cos(declination)*cos(hourangle));  //calculates m' for m calculation

//calculates m
	float m;
	m = mprime - 0.0018167*(mprime-1) - 0.002875*pow(mprime-1,2) - 0.0008083*pow(mprime-1,3);	//calculates optical path length

//calculates RSE (variable D, in AU)
	int delta, leap1;		//(delta = years since 1949), (leap1 = leap year indicator)
	float JD, ganom;		//(JD = Julian Day), (ganom = ellipse anomaly), (rse = sun earth radius (AU))
	double rse;
	delta = year-1949;		//years since 1949
	leap1 = delta/4;		//leap year indicator
	JD = 2432916.5+delta*365+leap1+daynumber+(dtime+6)/24.0;	//Julian day
//cout << "JD is: " << JD << endl;
	ganom=(357.528+0.9856474*(JD-2451545))*PI/180;				//ellipse anomaly
	rse=1.00014-0.01671*cos(ganom)-0.00014*cos(2*ganom);		//sun earth radius

//	irradience at the top of Earth's Atmosphere (also called air mass zero or AM0)
	double Lnot[size] = {1756.38, 1873.97, 1510.8, 1190.51, 970.088, 705.021};

//find Rayleigh factor
	long double fn2, fo2, fair, index, crosssection, multiplier, g0, g, ma, rayleigh;		//(fn2=depolarization of nitrogen) (fo2=depolarization of o2) (fair=dep. of air) (index=refractive index of air, specific to n2 and o2 levels) (crossection=scattering crossection per molecule) (multiplier=multiplier for C02) (g0=gravity at sea level) (g=gravity at height and location of measurement) (ma=mean molecular weight of dry air) (rayleigh=calculated rayleigh factor)
	wavelength[i] = wavelength[i]/1000;						//converts wavelength from nanometers to micrometers for these calculations
	multiplier = 1+(0.54*(Co2/1000000-0.0003));			// multiplier for C02
	index = (multiplier*(8060.51 +(2480990/(132.274-pow(wavelength[i],-2))) + (17455.7/(39.32957-pow(wavelength[i],-2)))))*pow(10,-8)+1;	//refractive index of air, specific to n2 and o2 levels
	fn2 = 1.034+(3.17E-4*(1/pow(wavelength[i],2)));		// depolarization of Nitrogen
	fo2 =1.096+(1.385E-3*(1/pow(wavelength[i],2)))+1.448E-4*(1/pow(wavelength[i],4));		// depolarization of Oxygen
	fair =((78.084*fn2)+(20.946*fo2)+(0.934*1.00)+(Co2*.0001*1.15))/(78.084+20.946+0.934+(Co2*.0001));		// depolarization of air
	crosssection = ((24*pow(PI,3)*pow(pow(index,2)-1,2))/(pow(wavelength[i]/10000,4)*pow(Ns,2)*pow(pow(index,2)+2,2)))*fair;	// crossection per molecule
	g0 = 980.6160*(1-0.0026373*cos(2*latitude)+0.0000059*pow(cos(2*latitude),2));		//gravity at sea level
	g = g0 - ((3.085462E-4+2.27E-7*cos(2*latitude))*height)+((7.254E-11+1.0E-13*cos(2*latitude))*pow(height,2))-(1.517E-17+6E-20*cos(2*latitude)*pow(height, 3));		//gravity at measurement's location
	ma = (15.0556*Co2/1000000)+28.9595;					// mean molecular weight of dry air
	rayleigh = crosssection * ((pressure*A)/(ma*g));	// calculated rayleigh factor
	aod[i] = ((log(Lnot[i] * pow(rse,2)/irr[i]))/(double)mprime)-(double)rayleigh;		// aod with new rayleigh factor
cout << "aod["<<i<<"] in calc_aod is: " << aod[i] << endl;
	wavelength[i] = wavelength[i]*1000;						// convert wavelength back to nanometers
	}
	return;		//end calc-data function
}
