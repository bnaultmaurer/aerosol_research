/*

This script will take an input and calculate the day number and the decimal time for that input.  Then it will put the two together to find the time at that daynumber.

ex. for 2014-06-20 (day 171) and time 0930 (9.75), the "daytime" is 171.406

*/

#include <iostream>
#include <cstring>
#include <cstdlib>
using namespace std;

bool leap_year(int year);	//will check for leap year
void extract_assign(string date, int &year, int &month, int &day);	//the function that will extract the information from the given YYYY-MM-DD
int convert_date(bool leap, int year, int month, int day);
double calc_daytime(double &dtime, int time, int daynumber, double &daytime);


int main()
{
	bool leap;
	int year, month, day, daynumber, time;
	double dtime, daytime;
	string date;
	bool leap_year(int year);
	cout << "Enter the date (YYYY-MM-DD): ";
	cin >> date;
	cout << "Enter the time (ex. 0945): ";
	cin >> time;
	extract_assign(date, year, month, day);
	leap = leap_year(year);
	cout << "\nLeap year? " << boolalpha << leap << noboolalpha << endl;
	daynumber = convert_date(leap, year, month, day);
	cout << date << " is the day number " << daynumber << endl;
	calc_daytime(dtime, time, daynumber, daytime);
	cout << "daytime is: " << daytime << endl << endl;

	return 0;
}

bool leap_year(int year)	//will check for leap year
{
	if((year%400)==0) {return true;}
	else if ((year%100)==0) {return false;}
	else if ((year%4)==0) {return true;}
	else {return false;}
}

void extract_assign(string date, int &year, int &month, int &day)	//will extract the year, month, and day from a string of form YYYY-MM-DD. Will then alter the year, month, and day integer variables using the address of (&).
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

double calc_daytime(double &dtime, int time, int daynumber, double &daytime)
{
	dtime = time / 100 + (float)(time%100)/60;
	dtime = dtime/24;
	daytime = daynumber + dtime;

	return dtime;
}
