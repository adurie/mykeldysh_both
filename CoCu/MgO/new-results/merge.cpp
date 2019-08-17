#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char **argv){
	ofstream Myfile;
	string Mydata;
	Mydata = "grand_out.txt";
	Myfile.open(Mydata.c_str(), ios::trunc);
	for (int i = 1; i < argc; i++){
		ifstream infile(argv[i]);
		string line;
		Myfile<<0.06 - (i-1)*0.01<<" ";
		while (getline(infile, line)){
			istringstream iss(line);
			double a, b;
			if (!(iss >> a >> b))
				break;
			Myfile<<b<<" ";
		}
		Myfile<<endl;
	}
	Myfile.close();
	return 0;
}
