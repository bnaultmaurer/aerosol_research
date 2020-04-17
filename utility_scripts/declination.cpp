/* this file gives us the declination variable which we will use to find the m' variable 
First it takes input from the user in the form of YYYY-MMM-DD and then it takes that string, seperates it into strings yearstr, monthstr, and daystr, then uses those values to decide if it is a leap year, then gives the value of the daynumber for any given date.  Then it takes that value and plugs the daynumber value into the formula given to find the solar declination*/

#include <cstdlib>
#include <math.h>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <iomanip>

using namespace std;

bool leap_year(int year);	//will check for leap year
void extract (string date, int &year, int &month, int &day);	//the function that will extract the information from the given YYYY-MM-DD
int convert_date(bool leap, int year, int month, int day);	
void find_declination(int daynumber,  float &declination);
float PI = 3.1415;

int main()
{
	bool leap;
	int year, month, day, daynumber;
	string date;
	float declination;
	bool leap_year(int year);

	cout << endl << "Enter the date (YYYY-MM-DD): ";
	cin >> date;

	extract(date, year, month, day);

	leap = leap_year(year);

	cout << "\nLeap year? " << boolalpha << leap << noboolalpha << endl;
	daynumber = convert_date(leap, year, month, day);

	cout << date << " is the day number " << daynumber << endl;
	find_declination(daynumber, declination);

	cout << "The declination is: " << declination << endl << endl;

	return 0;	
}	

bool leap_year(int year)	//will check for leap year
{
	if((year%400)==0) {return true;}
	else if ((year%100)==0) {return false;}
	else if ((year%4)==0) {return true;}
	else {return false;}
}

void extract(string date, int &year, int &month, int &day)	//will extract the year, month, and day from a string of form YYYY-MM-DD. Will then alter the year, month, and day integer variables using the address of (&). 
{
	int pos1, pos2;
	string yearstr, monthstr, daystr;
	pos1 = date.find_first_of("-");			//defining pos1 as the date string (char)
	pos2 = date.find_first_of("-", pos1+1);		//defining pos2 as the month string (char)
	yearstr = date.substr(0, pos1);
	monthstr = date.substr(pos1+1, pos2-pos1-1);
	daystr = date.substr(pos2+1);
	year = atoi(yearstr.c_str());
	month = atoi(monthstr.c_str());
	day = atoi(daystr.c_str());
	return;
}

int convert_date(bool leap, int year, int month, int day)
/* Calculate day number */
{
	switch (month) {
	case 1:
		return day;
	case 2:
		return day+31;
	case 3:
		return day+59+leap;
	case 4:
		return day+90+leap;
	case 5:
		return day+120+leap;
	case 6:
		return day+151+leap;
	case 7:
		return day+181+leap;
	case 8:
		return day+212+leap;
	case 9:
		return day+243+leap;
	case 10:
		return day+273+leap;
	case 11:
		return day+304+leap;
	case 12:
		return day+334+leap;
	default:
		cout << "Month value out of range" << endl;
		return 0;
	}
}

void find_declination(int daynumber, float &declination)
{
	declination = -23.4*PI/180 * (float)cos ( 2*PI*(daynumber+10)/365 );

	return;
}
