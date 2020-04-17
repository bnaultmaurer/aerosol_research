/*
	This program enters a folder of irradience (.irr) files, opens the files, and sorts
	through the files and takes out 16 specific wavelengths (368nm, 412nm, 440nm, 463nm,
	479nm, 495nm, 500nm, 520nm, 556nm, 610nm, 673nm, 675nm, 750nm, 778nm, 862nm, 870nm).
	It then creates 3 folders (single-data, graphs, time) and creates and inserts
	files into two of them.

		In "single-data" it takes the irradiences and sorts the data points in terms of
	wavelength.  So in the 368nm.dat file we have all the 368nm measurements throughout
	the day. The information is sorted as:  wavelength,	Irradience,	Decimal time,	Aerosol Optical Depth

		The graphs folder is only created so we don't have to make it later when we
	make graphs using Mathematica

		The time folder sorts the data for each measurement.  So instead of having each
	wavelength's value of all measurements, we have each measurement's value of all wavelengths.
	The information is sorted as:	Daynumber and time		Each wavelength's aerosol optical depth.
	The files in this folder should be used with the program aerosol-num-dens2.c

	IMPORTANT.  THIS PROGRAM WILL APPEND (ADD TO) EXISTING INFORMATION.  SO RUNNING IT MULTIPLE
	TIMES WITHOUT DELETING PREVIOUS INFORMATION WILL SIMPLY MAKE FILES LARGER


*/

#include <fstream>		//header providing file stream classes
#include <sstream>		//header providing string stream classes
#include <iostream>		//header that defines the standard input/output stream objects
#include <string>		//introduces string types, character traits and qa set of converting functions
#include <errno.h>		//C header that defines errno (last error number)
#include <dirent.h>		//---------------------------------------------------------------------
#include <unistd.h>		//---------------------------------------------------------------------
#include <sys/stat.h>	//---------------------------------------------------------------------
#include <sys/types.h>	//---------------------------------------------------------------------
#include <iomanip>		//provides parametric parameters, see http://www.cplusplus.com/reference/iomanip/
#include <cstring>		//manipulates C strings and arrays
#include <cstdlib>		//---------------------------------------------------------------------
#include <math.h>		//declares a ser of function to compute common mathematical operations and transformations

using namespace std;

void calc_data(string filepath, ofstream outfile[], long double &wavelength, float &irr);		//defines calc-data funtion, written below
void write_outputfile_single(string filepath, string fnames[], ofstream outfile[]);		//defines the write-output-single function, written below
void write_outputfile_daytime(string &filepath, string &fnameout);	//defines write_outputfile_daytime funtion, written below

const int size = 16;		//defines size of array for the 16 wavelengths
int Co2 = 402;				//initial C02 density in the atm in ppm (last updated in 2014)
float PI = 3.14159265359;	//defining Pi for later use
float latitude = 45.5807*PI/180;		//latitude of SJU in radians
float stdpressure = 1013.25;		//standard pressure at sea level (millibars)
int a;						//integer to be used as counter in write_outputfile_single to seperate wavelength's to their files
float height = 1200/3.2808;		//-------------------------------------------------------------
float A = 6.0221367E23;			//Avagadro's number
float Ns = 2.546899E19;		//molecular density of air
float pressure;

int main()
{
	cout << endl <<"***************************************************************";
	cout << endl <<"***************************************************************";
	cout << endl << "WARNING.  THIS PROGRAM WILL APPEND (ADD TO) EXISTING INFORMATION.  SO RUNNING IT MULTIPLE TIMES WITHOUT DELETING PREVIOUS INFORMATION WILL SIMPLY MAKE CREATED FILES LARGER";			//warning about appending nature of program
	cout << endl <<"***************************************************************";
	cout << endl <<"***************************************************************" << endl << endl;

	int i, num, time;		//initiate integers i=___, num=____, time=____;
	long double wavelength;		//initiate wavelength, used in many sub programs as definition of wavelength (368.0)
	ifstream fin;		//----------------------------------------------------------------------------
	ofstream outfile[size];		//--------------------------------------------------------------------
	string dir, filepath, newdir, newdirtime, fnameout;		//initiates strings for directory creation
	string fnames[size] = {"368nm.dat", "412nm.dat", "440nm.dat", "463nm.dat", "479nm.dat", "495nm.dat", "500nm.dat", "520nm.dat", "556nm.dat","610nm.dat", "673nm.dat", "675nm.dat", "750nm.dat", "778nm.dat", "862nm.dat", "870nm.dat"};	//creates
	DIR *dp;	//DIR function, acesses the directory and uses pointer dp to check if it is open
	struct dirent *dirp;	//dirent function, points to dirp to access directory files
	struct stat filestat;	//stat function, stores info in filestat to transfer to new file

	cout << "dir to get files of: " << flush;	//directs user to input directory they want to sort
 	getline( cin, dir );  	//defines dir as the directory to get files out of


	dp = opendir( dir.c_str() );	//opens the folder
 	if (dp == NULL)	//if nothing in the folder, quit program
		{
		cout << "Error(" << errno << ") opening " << dir << endl;
		return errno;
		}
	newdir = dir;		// changes variable to make new directory (folder graphs)
	newdir.append("/graphs");		//define new directory's name (graphs)
	mkdir(newdir.c_str(), 0711);		//makes new folder with new name (for graphs, nothing in them now)

	newdir = dir;		//changes variable to make new directory (folder time)
	newdir.append("/time");		//defines new directory's name (time)
	mkdir(newdir.c_str(), 0711);		//makes a new directory called "time"

	newdir = dir;		//changes variable to make new directory (folder single-data)
	newdir.append("/single-data");		//define new directory's name (single-data)
	mkdir(newdir.c_str(), 0711);		//makes new directory named "single-data"


	for (i=0; i<size; i++)	//loop to create/open output files
		{
		outfile[i].open((newdir + "/" + fnames[i]).c_str(), ios::app);		//creates/opens output files
		if (!outfile[i].is_open())	//CHECK, to see if file is open
			{
			cout << "Unable to open output file\n"; //CHECK, only displays if output files not opened
			return 0;
			}
		}

	while ((dirp = readdir( dp )))	//while we are reading the directory (dirp = directory pointer) (readdir = read directory)
		{
    		filepath = dir + "/" + dirp->d_name;		//string that creates path to place the file in the correct directory
		cout << "Accessing " << filepath << endl;	//prints to screen what files are being accessed

// If the file is a directory (or is in some way invalid) we'll skip it
    		if (stat( filepath.c_str(), &filestat )) continue;		//
    		if (S_ISDIR( filestat.st_mode ))         continue;
		write_outputfile_single(filepath, fnames, outfile);		//writes information to single-data files (funtion written below)
		write_outputfile_daytime(filepath, fnameout);		//writes information to time files (funtion written below)
		}

	for (i=0; i<size; i++)	//loop to close all the output files
		{
		outfile[i].close();		//closes all output files
		}

	return 0;		//end of main program
}



void write_outputfile_single(string filepath, string fnames[], ofstream outfile[])		//opens files in user inputed directory, ignores the first two lines, then writes new .dat file with defined wavelength's info from the .irr file, sorted into wavelength specific files
{
	string line, time;		//defines strings (line = line from file that we are reading) (time = string of decimal time extracted from filename)
	ifstream infile;		//name of input file stream
	int pos3, i;			//defines integers (pos3 = used to extract string for time) (i = integer value used for loops throughout program)
	float irr;				//defines float value (irr = value of the irradience of the specified wavelength)
	long double wavelength;

	pos3 = filepath.find_last_of("-");		//position of file name to get time
	time = filepath.substr(pos3+1, 4);		//defines string time from file
	infile.open (filepath.c_str(), ios::app);		//opens the input files
	if (!infile.is_open())				//makes sure the input file is open
		{
		cout << "Unable to open input file\n";
		return;
		}
	i=0;											//starts the removal process below
	while ( getline (infile,line))		//this reads the files and removes the first two lines of the file
		{
		if (i>1)
			{
			stringstream(line) >> wavelength >> irr;	//goes through and finds all the individual variables we want
			if (wavelength==368.0)	//if this wavelength, write data out
				{
				a=0;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==412.0)
				{
				a=1;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==440.0)
				{
				a=2;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==463.0)
				{
				a=3;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==479.0)
				{
				a=4;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==495.0)
				{
				a=5;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==500.0)
				{
				a=6;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==520.0)
				{
				a=7;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==556.0)
				{
				a=8;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==610.0)
				{
				a=9;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==673.0)
				{
				a=10;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==675.0)
				{
				a=11;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==750.0)
				{
				a=12;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==778.0)
				{
				a=13;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==862.0)
				{
				a=14;
				calc_data(filepath, outfile, wavelength, irr);
				}
			if (wavelength==870.0)
				{
				a=15;
				calc_data(filepath, outfile, wavelength, irr);
				}
			}
		i++;
		}
	infile.close();	//close the input file
	return;
}

void calc_data( string filepath, ofstream outfile[], long double &wavelength, float &irr)
{

//defines year, month, day, time
	int pos1, pos2, pos3, year, month, day, daynumber, time;
	bool leap;
	string yearstr, monthstr, daystr,timestr;
	pos1 = filepath.find_first_of("-");			//defining pos1 as the date string (char)
	pos2 = filepath.find_first_of("-", pos1+1);		//defining pos2 as the month string (char)
	pos3 = filepath.find_last_of("-");		//position of file name to get time
	timestr = filepath.substr(pos3+1, 4);		//defines string time from file
	yearstr = filepath.substr(0, pos1);
	monthstr = filepath.substr(pos1+1, pos2-pos1-1);
	daystr = filepath.substr(pos2+1);
	year = atoi(yearstr.c_str());
	month = atoi(monthstr.c_str());
	day = atoi(daystr.c_str());
	time = atoi(timestr.c_str());

//leap year checker
	if((year%400)==0) {leap=1;}
	else if ((year%100)==0) {leap=0;}
	else if ((year%4)==0) {leap=1;}
	else {leap=0;}

//calculates day number
	if(month== 1)
		{daynumber= day;}
	else if(month== 2)
		{daynumber= day+31;}
	else if(month== 3)
		{daynumber= day+59+leap;}
	else if(month== 4)
		{daynumber= day+90+leap;}
	else if(month== 5)
		{daynumber= day+120+leap;;}
	else if(month== 6)
		{daynumber= day+151+leap;}
	else if(month== 7)
		{daynumber= day+181+leap;;}
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
	else
		cout << "Month value out of range" << endl;
	cout << "Wavelength is: " << wavelength << endl;
	cout << "Daynumber is: " << daynumber << endl;

//calculates decimal time
	float dtime;
	cout << "time is: " << time << endl;
	dtime = time/100 + (float)(time%100)/60;
	cout << "Decimal time is: " << dtime << endl;

//calculates hourangle(h) to use to find m'
	float hourangle;
	hourangle = (15*PI/180)*(12-dtime);
	cout << "hourangle is: " << hourangle << endl;

//calculates declination
	float declination;
	declination = -23.4*PI/180 * (float)cos ( 2*PI*(daynumber+10)/365 );
	cout << "Declination is: " << declination << endl;

//calculates mprime
	float mprime;
	mprime = 1/(sin(latitude)*sin(declination) + cos(latitude)*cos(declination)*cos(hourangle));
	cout << "latitude is: " << latitude << endl;
	cout << "mprime is: " << mprime << endl;

//calculates m
	float m;
	m = mprime - 0.0018167*(mprime-1) - 0.002875*pow(mprime-1,2) - 0.0008083*pow(mprime-1,3);
	cout << "m is: " << m << endl;

//calculates RSE (D)
	int delta, leap1;	// years since 1949, leap year indicator
	float JD, ganom, rse;	// Julian Day, ellipse anomaly, sun earth radius
	delta = year-1949;
	leap1 = delta/4;
	JD = 2432916.5+delta*365+leap1+daynumber+(dtime+6)/24.0;
	ganom=(357.528+0.9856474*(JD-2451545))*PI/180;
	rse=1.00014-0.01671*cos(ganom)-0.00014*cos(2*ganom);
	cout << "D is: " << rse << endl;

//defines pressure for each data set
	float pressure;
cout << "daynumber is: '" << daynumber << "'" << endl;
cout << "Pressure is: " << pressure << endl;
	if (daynumber == 148)
		pressure = 983;
	else if (daynumber == 154)
		pressure = 979;
	else if (daynumber == 155)
		pressure = 977;
	else if (daynumber == 161)
		pressure = 980;
	else if (daynumber == 164)
		pressure = 981;
	else if (daynumber == 171)
		pressure = 979;
	else if (daynumber == 175)
		pressure = 978.5;
	else if (daynumber == 190)
		pressure = 981;
	else if (daynumber == 191)
		pressure = 982;
	else if (daynumber == 224)
		pressure = 983;
	else if (daynumber == 225)
		pressure = 983;
	else
		{
			cout << "No pressure for this day" << endl << endl;
			return;
		}
	cout << "pressure for this day: " << pressure << endl;

//defines the arrays for beta and L0
	float beta[size] = { 0.51154, 0.31922, 0.2429, 0.19704, 0.17123, 0.1497, 0.1435, 0.1222, 0.0937, 0.06372, 0.042791, 0.042285, 0.02759, 0.023792, 0.0162614, 0.015619 };
	float Lnot[size] = { 1.1405, 1.803, 1.769, 2.008, 2.0425, 1.992, 1.9135, 1.8295, 1.859, 1.7235, 1.521, 1.51, 1.266, 1.1915, 1.9949, 0.9747 };

//calculates aod
	float aod1;
	aod1 = (log(Lnot[a] * pow(rse,2)/irr) - beta[a]*m*(pressure/stdpressure))/mprime;
	cout << "aod is: " << aod1 << endl;


//find the new Rayleigh factor
	long double fn2, fo2, fair, index, crosssection, multiplier, g0, g, ma, rayleigh;
	wavelength = wavelength/1000;
	cout << "wavelength is(micrometers): " << wavelength << endl;
	multiplier = 1+(0.54*(Co2/1000000-0.0003));
	index = (multiplier*(8060.51 +(2480990/(132.274-pow(wavelength,-2))) + (17455.7/(39.32957-pow(wavelength,-2)))))*pow(10,-8)+1;
	fn2 = 1.034+(3.17E-4*(1/pow(wavelength,2)));
	fo2 =1.096+(1.385E-3*(1/pow(wavelength,2)))+1.448E-4*(1/pow(wavelength,4));
	fair =((78.084*fn2)+(20.946*fo2)+(0.934*1.00)+(Co2*.0001*1.15))/(78.084+20.946+0.934+(Co2*.0001));
	crosssection = ((24*pow(PI,3)*pow(pow(index,2)-1,2))/(pow(wavelength/10000,4)*pow(Ns,2)*pow(pow(index,2)+2,2)))*fair;
		cout << "Cross section is: " << crosssection<< endl;
	g0 = 980.6160*(1-0.0026373*cos(2*latitude)+0.0000059*pow(cos(2*latitude),2));
	g = g0 - ((3.085462E-4+2.27E-7*cos(2*latitude))*height)+((7.254E-11+1.0E-13*cos(2*latitude))*pow(height,2))-(1.517E-17+6E-20*cos(2*latitude)*pow(height, 3));
		cout << "g is: " << g << endl;
	ma = (15.0556*Co2/1000000)+28.9595;
		cout << "ma is: " << ma << endl;
	rayleigh = crosssection * ((pressure*A)/(ma*g));
	cout << "Rayleigh factor is: " << rayleigh << endl;

//with new Rayleigh factor
	float aod;
	aod = ((log(Lnot[a] * pow(rse,2)/irr))/mprime)-rayleigh;
	cout << "AOD is: " << aod << endl;

//writes everything to the output file
	wavelength = wavelength*1000;
	outfile[a] <<  (int)wavelength << "\t" << scientific << irr << "\t" <<  fixed << dtime << /*"\t"<< aod1  <<*/ "\t" << aod <<endl;
	cout << "Processed file " << filepath << " for wavelength " << (int)wavelength << "nm" << endl << endl;	//let the user know it processed the file


	return;
}

void write_outputfile_daytime(string &filepath, string &fnameout)		//takes the files, opens them, deletes the first two lines, then writes new .dat file with info in the .irr file
{
	string line, line1, want, keep, timestr;
	ifstream infile;
	ofstream outfile;
	int pos1, pos2, i, daynumber, time;
	float f1, f2;
	double dtime, daytime;

	cout << "filepath is: " << filepath << endl;
	pos1 = filepath.find_last_of("/");	//this chunk allows me to change the name of output file
	keep = filepath.substr(0,pos1+1);
	timestr = filepath.substr(22,4);	//find time
	time = atoi(timestr.c_str());		//time is int
	cout << "timestr is: " << timestr << endl;
	pos2 = keep.find_first_of("/");		//finding the position to insert "dat/"
	keep.insert(pos2, "/time");
	fnameout = keep + timestr + ".csv";

//find the values from filepath for daynumber calc
	int year, month, day;
	string yearstr, monthstr, daystr;
	yearstr = filepath.substr(0,4);		//find year
	monthstr = filepath.substr(5,2);	//find month
	daystr = filepath.substr(8,9);		//find day
	year = atoi(yearstr.c_str());		//year is int
	month = atoi(monthstr.c_str());		//month is int
	day = atoi(daystr.c_str());			//day is int
	cout << "date is: " << year << " " << month << " " << day << " " << time << endl;

//leap year checker
	bool leap;
	if((year%400)==0) {leap=1;}
	else if ((year%100)==0) {leap=0;}
	else if ((year%4)==0) {leap=1;}
	else {leap=0;}

//calculates day number
	if(month== 1)
		{daynumber= day;}
	else if(month== 2)
		{daynumber= day+31;}
	else if(month== 3)
		{daynumber= day+59+leap;}
	else if(month== 4)
		{daynumber= day+90+leap;}
	else if(month== 5)
		{daynumber= day+120+leap;;}
	else if(month== 6)
		{daynumber= day+151+leap;}
	else if(month== 7)
		{daynumber= day+181+leap;;}
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
	else
		cout << "Month value out of range" << endl;
	cout << "Daynumber is: " << daynumber << endl;


//finds the daytime value to insert at the front of file
	cout << "time is: " << time << endl;
	dtime = time / 100 + (float)(time%100)/60;
	dtime = dtime/24;
	daytime = daynumber + dtime;
	cout << "daytime is: " << daytime << endl;

	infile.open (filepath.c_str(), ios::app);		//opens the input files
	outfile.open (fnameout.c_str(), ios::app);		//opens the output files

	if (!infile.is_open())				//makes sure the input file is open
		{
		cout << "Unable to open input file\n";
		return;
		}

	if (!outfile.is_open())				//makes sure the output file is open
		{
		cout << "Unable to open output file\n";
		return;
		}

	outfile << daytime << ", ";
	i=0;								//starts the removal process below
	while ( getline (infile,line))		//this reads the files and removes the first two lines of the file
		{
	if (i>1)
		{
		stringstream(line) >> f1 >> f2;
		if (f1==412.0 || f1==440.0 || f1==463.0 || f1==479.0 || f1==500.0 || f1==520.0 || f1==556.0 || f1==610.0 || f1==675.0 || f1==750.0 || f1==778.0 || f1==870.0 )
			{
		outfile << scientific << f2 << ", ";
			}
		}
		i++;
		}

	cout << endl;
	outfile.close();
	infile.close();
	return;
}
