#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
/* #include "TB.h" */
#include "TB_altered.h"
/* #include "TBdynamic.h" */

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
//calculates the bandstructure of fcc Cu

int main(){

	/* cout<<"Name the data file\n"; */
	string Mydata;
	/* getline(cin, Mydata); */
	ofstream Myfile;	
	/* Mydata += ".txt"; */
	Mydata += "bandido.txt";
	Myfile.open( Mydata.c_str(),ios::trunc );

	Vector3d d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, d_9;
	Vector3d d_10, d_11, d_12, d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	/* double a = 6.692; */
	double a = 0.5;
	double ryd = 13.6058;
	
	//position vectors of nearest neighbours in fcc
	d_1 << a, a, 0;
	d_2 << -a, -a, 0;
	d_3 << a, 0, a;
	d_4 << -a, 0, -a;
	d_5 << 0, a, a;
	d_6 << 0, -a, -a;
	d_7 << -a, a, 0;
	d_8 << a, -a, 0;
	d_9 << -a, 0, a;
	d_10 << a, 0, -a;
	d_11 << 0, -a, a;
	d_12 << 0, a, -a;

	//position vectors of next nearest neighbours
	d_13 << 2*a, 0, 0;
	d_14 << -2*a, 0, 0;
	d_15 << 0, 2*a, 0;
	d_16 << 0, -2*a, 0;
	d_17 << 0, 0, 2*a;
	d_18 << 0, 0, -2*a;

	//initialise onsite and hopping matrices for each nn
	Matrix<dcomp, 9, 9> u;
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	Matrix<dcomp, 9, 9> t_10, t_11, t_12, t_13, t_14, t_15, t_16;
	Matrix<dcomp, 9, 9> t_17, t_18;
	u = TB(2, 0, 0, 8, d_1);
	t_1 = TB(2, 1, 0, 8, d_1);
	t_2 = TB(2, 1, 0, 8, d_2);
	t_3 = TB(2, 1, 0, 8, d_3);
	t_4 = TB(2, 1, 0, 8, d_4);
	t_5 = TB(2, 1, 0, 8, d_5);
	t_6 = TB(2, 1, 0, 8, d_6);
	t_7 = TB(2, 1, 0, 8, d_7);
	t_8 = TB(2, 1, 0, 8, d_8);
	t_9 = TB(2, 1, 0, 8, d_9);
	t_10 = TB(2, 1, 0, 8, d_10);
	t_11 = TB(2, 1, 0, 8, d_11);
	t_12 = TB(2, 1, 0, 8, d_12);

	t_13 = TB(2, 1, 1, 8, d_13);
	t_14 = TB(2, 1, 1, 8, d_14);
	t_15 = TB(2, 1, 1, 8, d_15);
	t_16 = TB(2, 1, 1, 8, d_16);
	t_17 = TB(2, 1, 1, 8, d_17);
	t_18 = TB(2, 1, 1, 8, d_18);
	/* Myfile<<"P X Y"<<endl; */
	Myfile<<"X Y"<<endl;

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> E;
	Matrix<dcomp, 9, 9> I = Matrix<dcomp, 9, 9>::Identity();

	/* Matrix<dcomp, 18, 18> E; */
	/* Matrix<dcomp, 9, 9> u_11, u_12; */
	/* Matrix<dcomp, 18, 18> U, T; */
	/* Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero(); */

	int nt = 251;
	VectorXd vec1(nt), vec2(nt), vec3(nt), vec4(nt), vec5(nt), vec6(nt), vec7(nt), vec8(nt), vec9(nt);
	double k_x, k_y, k_z, pi;
	Vector3d K;
	for (int k = 0; k < 251; k++)
	{

		//This to create E(k_y) presently at neck
		k_x = 2.532374;
		k_z = 2.532374;
		k_y = 2*M_PI*k/250.;

		/* //This for conventional bandstructure */
		/* if (k < 101){ */
		/* 	pi = 2.*M_PI*k/100.; */
		/* 	k_x = pi; */
		/* 	k_y = 0; */
		/* 	k_z = 0; */
		/* } */
		/* if ((k > 100) && (k < 151)){ */
		/* 	pi = M_PI*(k-100)/50.; */
		/* 	k_x = 2.*M_PI; */
		/* 	k_y = pi; */
		/* 	k_z = 0; */
		/* } */	
		/* if ((k > 150) && (k < 201)){ */
		/* 	pi = M_PI*(k-150)/50.; */
		/* 	k_x = 2.*M_PI - pi; */
		/* 	k_y = M_PI; */
		/* 	k_z = pi; */
		/* } */
		/* if ((k > 200) && (k < 251)){ */
		/* 	pi = M_PI*(k-200)/50.; */
		/* 	k_x = M_PI-pi; */
		/* 	k_y = M_PI-pi; */
		/* 	k_z = M_PI-pi; */
		/* } */

		K(0) = k_x;
		K(1) = k_y;
		K(2) = k_z;

		//fully diagonalised Hamiltonian
		/* E = u -I*(0.2 + 0.57553) + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K)) */
		/* E = u -I*(0.11 + 0.57553) + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K)) */
		E = u -I*(0.29 + 0.57553) + t_1*exp(i*d_1.dot(K))+ t_2*exp(i*d_2.dot(K))+ t_3*exp(i*d_3.dot(K))
			+ t_4*exp(i*d_4.dot(K)) + t_5*exp(i*d_5.dot(K)) + t_6*exp(i*d_6.dot(K))
				+ t_7*exp(i*d_7.dot(K)) + t_8*exp(i*d_8.dot(K))
				+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K))
				+ t_11*exp(i*d_11.dot(K)) + t_12*exp(i*d_12.dot(K))
				+ t_13*exp(i*d_13.dot(K)) + t_14*exp(i*d_14.dot(K))
				+ t_15*exp(i*d_15.dot(K)) + t_16*exp(i*d_16.dot(K))
				+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
		SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> es;
		es.compute(E);
		Matrix<double, 9, 1> O;
		O = es.eigenvalues();

		/* //Hamiltonian firstly diagonalised in-plane */
		/* u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + */ 
		/* 	t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K)); */
		/* u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K)); */
		/* U << u_11, u_12, u_12.adjoint(), u_11; */
		/* T << t_15, zero, u_12.adjoint(), t_16; */
		/* E = U + T*exp(i*d_15.dot(K)) + T.adjoint()*exp(i*d_16.dot(K)); */
		/* SelfAdjointEigenSolver<Matrix<dcomp, 18, 18>> es; */
		/* es.compute(E); */
		/* Matrix<double, 18, 1> O; */
		/* O = es.eigenvalues(); */

		/* Myfile<<"A"<<" "<<k<<" "<<O(0)<<endl; */
		/* Myfile<<"B"<<" "<<k<<" "<<O(1)<<endl; */
		/* Myfile<<"C"<<" "<<k<<" "<<O(2)<<endl; */
		/* Myfile<<"D"<<" "<<k<<" "<<O(3)<<endl; */
		/* Myfile<<"E"<<" "<<k<<" "<<O(4)<<endl; */
		/* Myfile<<"F"<<" "<<k<<" "<<O(5)<<endl; */
		/* Myfile<<"G"<<" "<<k<<" "<<O(6)<<endl; */
		/* Myfile<<"H"<<" "<<k<<" "<<O(7)<<endl; */
		/* Myfile<<"I"<<" "<<k<<" "<<O(8)<<endl; */

		vec1(k) = O(0);
		vec2(k) = O(1);
		vec3(k) = O(2);
		vec4(k) = O(3);
		vec5(k) = O(4);
		vec6(k) = O(5);
		vec7(k) = O(6);
		vec8(k) = O(7);
		vec9(k) = O(8);

		/* Myfile<<"J"<<" "<<k<<" "<<O(9)<<endl; */
		/* Myfile<<"K"<<" "<<k<<" "<<O(10)<<endl; */
		/* Myfile<<"L"<<" "<<k<<" "<<O(11)<<endl; */
		/* Myfile<<"M"<<" "<<k<<" "<<O(12)<<endl; */
		/* Myfile<<"N"<<" "<<k<<" "<<O(13)<<endl; */
		/* Myfile<<"O"<<" "<<k<<" "<<O(14)<<endl; */
		/* Myfile<<"P"<<" "<<k<<" "<<O(15)<<endl; */
		/* Myfile<<"Q"<<" "<<k<<" "<<O(16)<<endl; */
		/* Myfile<<"R"<<" "<<k<<" "<<O(17)<<endl; */

	}

	for (int k = 0; k < vec1.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec1(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec2.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec2(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec3.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec3(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec4.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec4(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec5.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec5(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec6.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec6(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec7.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec7(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec8.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec8(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec9.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec9(k)<<endl;

	Myfile.close();
	return 0;
}
