#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <ctime>
#include "CoCuCo.h"
#include <nag.h>
#include <nagc05.h>
#include <nag_stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif
  static double NAG_CALL f(double x, Nag_Comm *comm);
#ifdef __cplusplus
}
#endif 

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<complex<double>, 9, 9> M9;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vec3;
typedef vector<Matrix<dcomp, 9, 9>, aligned_allocator<Matrix<dcomp, 9, 9>>> vM;
//calculates the bandstructure of fcc Cu

typedef struct
	{
		double *k_y;
		double *k_z;
		vec3 *Cu_pos;
	        vM *copper;
	}
variables;

M9 H(const vec3 &pos, const vM &U, const double x, const double y, const double z){
	Vector3d K;
	dcomp i = -1;
	i = sqrt(i);
	K << x, y, z;
	M9 result;
	result.fill(0.);
	Vector3d tmp_vec;
	for (int k = 0; k < pos.size(); k++){
		tmp_vec = pos[k];
		result = result + U[k]*exp(i*tmp_vec.dot(K));
	}
	return result;
}

M9 InPlaneH(const vec3 &pos, const Vector3d &basis, const vM &U, const double x, const double y, const double z){
	Vector3d K;
	dcomp i = -1;
	i = sqrt(i);
	K << x, y, z;
	M9 result;
	result.fill(0.);
	Vector3d tmp_vec;
	for (int k = 0; k < pos.size(); k++){
		tmp_vec = pos[k] - basis;
		result = result + U[k]*exp(i*tmp_vec.dot(K));
	}
	return result;
}

static double NAG_CALL f(double k_x, Nag_Comm *comm)
{
	variables * send = (variables *) comm->p;
	Matrix<dcomp, 9, 9> E;
	E = H(*send->Cu_pos, *send->copper, k_x, *send->k_y, *send->k_z);// - fermi*I;
	SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> es;
	es.compute(E, EigenvaluesOnly);
	Matrix<double, 9, 1> O;
	O = es.eigenvalues();
	return O(5);
}

int pass(variables * send, double &k_x){
	double a, b;
	double eta, eps;
	Integer exit_status = 0;
	NagError fail;
	Nag_Comm comm;
	comm.p = send;
	k_x = 0;

	INIT_FAIL(fail);

		/* For communication with user-supplied functions: */

	/* a = 2.5; */
	/* b = 3.; */
	double h = 0.03;
	eps = 1e-07;
	eta = 0.0;
	c05auc(&k_x, h, eps, eta, f, &a, &b, &comm, &fail);
	/* c05ayc(a, b, eps, eta, f, &k_x, &comm, &fail); */
	int test;
	if (fail.code != NE_NOERROR) {
		printf("%s\n", fail.message);
		test = 1;
	}
	else 
		test = 0;
	/* 	if (fail.code == NE_TOO_SMALL || fail.code == NE_PROBABLE_POLE) */
	/* 		printf("Final point = %12.5f\n", k_x); */
	/* 	exit_status = 1; */
	/* } */

	return test;
}

int main(){
	vector<double> Cu1, Cu2;
	Cu1.reserve(10); Cu2.reserve(10);
	Cu1 = param(2,1); Cu2 = param(2,2);
	M9 Cu;
	Cu = U(2,0);

	//lattice vectors for the whole system;
	Vector3d lat_vec1, lat_vec2, Cu_lat_oop;
	lat_vec1 << 0.5, 0, 0.5;
	lat_vec2 << 0.5, 0, -0.5;
	Cu_lat_oop << 0.5, 0.5, 0;

	//Distance information for n.n and n.n.n
	double Cu_nn_dist, Cu_nnn_dist, Cu_nnnn_dist;
	Cu_nn_dist = M_SQRT2/2.;
	Cu_nnn_dist = 1.;
	Cu_nnnn_dist = 0;//this tells the code to ignore

	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	Vector3d tmp_vec;
	vM copper;
	vec3 Cu_pos;
	Matrix<dcomp, 9, 9> tmp_mat;
	double distance;
	for (int i1 = -2; i1 < 3; i1++){
		for (int i2 = -2; i2 < 3; i2++){
			for (int i3 = -2; i3 < 3; i3++){
				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*Cu_lat_oop;
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Cu_pos.emplace_back(tmp_vec);
					copper.emplace_back(Cu);
				}
				else if (distance < Cu_nn_dist + 1e-3){
					Cu_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Cu1, x, y, z);
					copper.emplace_back(tmp_mat);
				}
				else if (distance < Cu_nnn_dist + 1e-3){
					Cu_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Cu2, x, y, z);
					copper.emplace_back(tmp_mat);
				}
			}
		}
	}

	time_t now = time(0);
	tm *ltm = localtime(&now);
	string Mydata;
	Mydata = to_string(ltm->tm_mday);
	Mydata += "-";
	Mydata += to_string(1+ltm->tm_mon);
	Mydata += "-";
	Mydata += to_string(1900+ltm->tm_year);
	Mydata += "-Cu_fermi.txt";

	ofstream Myfile;	

	Myfile.open( Mydata.c_str(),ios::trunc );

	/* double ryd = 13.6058; */
	double ryd = 1.;
	/* double fermi = 0.5074;//Fermi level for Au */
	/* double fermi = 0.57553; */
	/* double fermi = 0.7466;//Fermi level for Fe */
	/* double fermi = 0.3191;//Fermi level for MgO */
	double fermi = 0;//note revisions to AuMgOFe.h

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> E;
	Matrix<dcomp, 9, 9> I = Matrix<dcomp, 9, 9>::Identity();
	Matrix<dcomp, 18, 18> Ibig = Matrix<dcomp, 18, 18>::Identity();
	M9 ins_11, ins_12, ins_21, ins_22;
	Matrix<dcomp, 18, 18> ins, E2;

	double k_x, k_y, k_z;
	variables send;
	send.copper = &copper;
	send.Cu_pos = &Cu_pos;

	int size = 500;
	int test;
	for (int k = 0; k < size+1; k++){
		for (int l = 0; l < size+1; l++){
			k_z = 4.*M_PI*k/(1.*size) - 2.*M_PI;
			k_y = 4.*M_PI*l/(1.*size) - 2.*M_PI;
			send.k_y = &k_y;
			send.k_z = &k_z;
			test = pass(&send, k_x);
			if (test == 0)
				Myfile<<k_x<<" "<<k_y<<" "<<k_z<<endl;
		}
	}

	Myfile.close();
	return 0;
}
