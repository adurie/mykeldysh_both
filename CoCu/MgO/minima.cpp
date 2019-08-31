#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
using namespace std;

int main(){
	ifstream infile("Cu_fermi.txt");
	string line;
	double a, b, c;
	while (getline(infile, line)){
		istringstream iss(line);
		if (!(iss >> a >> b >> c)) {break;}
		if (b == c){
			if (a > M_PI)
				a -= 2*M_PI;
			if (b > M_PI)
				b -= 2*M_PI;
			if (b < -M_PI)
				b += 2*M_PI;
			if (a > 0)
				b = -b;
			cout<<b<<" "<<a<<endl;
		}
	}
	return 0;
}
