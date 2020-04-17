/*
	LAST UPDATED: 12/13/2015

		This program sorts through a data file from the Kipp & Zonen PGS-100 spectrometer, in a .csv or .txt
	format and it only outputs data that meets a certain criteria.  Right now this criteria is that the wavelength
	closest to the oxygen absorption line, 781.855nm, meets a minimum count of 20,000.  If the data at that time
	meets the requirement (is above 20,000 counts) then the count values for the wavelengths that we desire are
	written to a file with "-sorted.csv" appended to the end of the original file name.

	Six wavelengths are: 440.699, 555.803, 674.855, 778.221, 869.738, 1019.721 (all in nm)

	Output files are in the format of :
		HH:MM:SS, wavelength 1, wavelength 2, ....
		HH:MM:SS, wavelength 1, wavelentgh 2, ....

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

typedef vector<vector<float> > matrix;		// required definition of matrix at top of source code

void extract_counts(string line, int linesize, int *counts);	// takes out count values for measurements
void extract_time(string line, string *time, int linesize, string *stime);	// takes out time in format HH:MM:SS and puts them in array stime

int a;														// initialize integer to keep track of lines in/from file

int main()
{
	int i, j, linesize=300, pos, check;						// set linesize too big so that we can redefine it later
	int counts[linesize];
	string ptime[linesize];
	double wavelength, *aod, *irr;
	string *time, stime[linesize];
	ifstream infile;										// set name to read info from
	ofstream outfile;										// set name to output info to
	string filepath, line, file, output, input;				// initalizes strings for directory creation
	char *fname1, *fname2;									// variables to be used to define the input/output files
	matrix print(7,vector<float>(linesize,0.0));			// matrix for printing to output file

	cout << "Enter file to process (gYYYYMMDD.txt): ";
	getline(cin, input);									// get string
	fname1 = new char [input.size()+1];						// change string to char variable
	strcpy(fname1, input.c_str());							// change input to a c string
	pos = input.find(".");									// find end of file name
	output = input.substr(0,pos) + "-sorted.csv";			// define output name as input name with -sorted added

	infile.open(input.c_str());								// set input file name
	outfile.open(output.c_str()/*, ios::app*/);				// set output file name and create file // for all runs
	outfile << "time, 440.699nm, 555.803nm, 674.855nm, 778.221nm, 869.738nm, 1019.721nm" << endl;

	if(!infile.is_open()){									// CHECK, to see if file is open
		cout << "Unable to open output file\n"; 			// CHECK, only displays if output files not opened
		return 0;
	}

	if(!outfile.is_open()){									// CHECK, to see if file is open
		cout << "Unable to open output file\n"; 			// CHECK, only displays if output files not opened
		return 0;
	}

	a=0;													// initialize integer to keep track of lines
	while(getline(infile,line)){
		if(a==0){
			pos=0;											// initialize pos
			i=0;											// initialize integer
			do{
				pos = line.find(",",pos+1);					// find position of comma
				i++;										// increase integer for next comma
			}while(pos>0);									// do while there are still commas in front
			linesize=i;										// redefine linesize when run out of commas
		}
		if(a==2){											// if time line, read in times (for testing)
			extract_time(line, time, linesize, stime);
		}
		if(a>7){											// start after header lines
		stringstream(line) >> wavelength >> line;			// read in wavelength string and counts string
			for(i=1;i<linesize;i++){
				if(wavelength==781.855){					// print reference variable for time check
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[0][i] = counts[i];			// fill in print matrix with counts values
					}
				}
				if(wavelength==440.699){
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[1][i] = counts[i];			// fill in print matrix with counts values
					}
				}
				if(wavelength==555.803){
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[2][i] = counts[i];			// fill in print matrix with counts values
					}
				}
				if(wavelength==674.855){
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[3][i] = counts[i];			// fill in print matrix with counts values
					}
				}
				if(wavelength==440.699){
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[4][i] = counts[i];			// fill in print matrix with counts values
					}
				}
				if(wavelength==778.221){
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[5][i] = counts[i];			// fill in print matrix with counts values
					}
				}
				if(wavelength==1019.721){
					extract_counts(line, linesize, counts);	// read in counts values and seperate them into vector
					for(i=1;i<linesize;i++){
						print[6][i] = counts[i];			// fill in print matrix with counts values
					}
				}
			}		// end of for(i++...) loop
		}			// end of if(a>7) loop
		a++;		// increase a for next iteration
		}			// end of while loop

// running the check and printing out good measurements
	for(i=0;i<linesize;i++){								// for each time
//		cout << "This is the count value that is checked: " << print[0][i] << endl;	 // for testing
		if(print[0][i]>20000){								// if the minimum count check is reached
			outfile << stime[i];							// write time
			for(j=1;j<7;j++){								// for each wavelength
				outfile << "," << print[j][i];				// write count value
			}
			outfile << endl;								// go to next line
		}
	}														// end of the for loop for check
	cout << endl << "File successfully processed.  Information written to " << output << endl << endl;
	return 0;												// return for end of main program
}															// end of main program

void extract_counts(string line, int linesize, int *counts){

	int i;													// set integer
	int pos1[linesize], pos2[linesize], pos3[linesize];		// define variables
	string header, scounts[linesize];

	for(i=0;i<linesize;i++){								// find comma
		pos2[i] = line.find(",");							// find position of comma
		scounts[i] = line.substr(0,pos2[i]);				// extract counts values
		counts[i] = atoi(scounts[i].c_str());				// convert from string to int
		line.erase(0,pos2[i]+1);							// erase counts values past comma
		}
	return;
}

void extract_time(string line, string *time, int linesize, string *stime){
	int i, f;
	int pos1[linesize];
	string header, shour[linesize], smin[linesize];

	for(i=0;i<linesize;i++){
		if(i==0){
			header = line.substr(0,5);						// take out the "time,"
		}
		else if(i>0){
			pos1[i] = line.find(",")+1;						// find the next comma
			stime[i] = line.substr(pos1[i],7);				// write the to appropriate time array
			line.erase(0,pos1[i]);							// erase that part of the string
		}
	}
	return;
}
