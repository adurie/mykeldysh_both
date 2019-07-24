#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

int main(){
	ifstream infile("22-7-2019-Fe_DOS.txt");
	string line;
	double DOS = 0;
	double E = 0;
	/* double step = 0.001; */
	double step = 0.0026;
	/* double step = 0.0001; */
	double a, b;
	/* while (E < 22) */ 
	while (E < 16) 
	/* while (E < 11) */ 
	/* while (E < 11.308) */ 
	/* while (E < 4.8762) */ 
	/* while (E < 5.80353) */ 
	/* while (E < 4.8745) */ 
	/* while (E < 2) */ 
	{
		getline(infile, line);
		istringstream iss(line);
		if (!(iss >> a >> b)) {break;}
		b*=8*M_PI*M_PI;
		DOS+=b;
		E = DOS*step;
	}
	cout<<a<<endl;
	cout<<E<<endl;
	return 0;
}
