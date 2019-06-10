#include<iostream>
#include<ctime>
#include <string>
using namespace std;

int main()
{

	time_t now = time(0);
	tm *ltm = localtime(&now);
	string Mydata;
	Mydata = to_string(ltm->tm_mday);
	Mydata += "-";
	Mydata += to_string(1+ltm->tm_mon);
	Mydata += "-";
	Mydata += to_string(1900+ltm->tm_year);
	Mydata += ".txt";
	cout<<Mydata<<endl;

	return 0;
}

