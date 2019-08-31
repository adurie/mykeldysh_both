#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <ctime>
#include "AuMgOFe.h"
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
		vec3 *Au_pos;
	        vM *gold;
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
	E = H(*send->Au_pos, *send->gold, k_x, *send->k_y, *send->k_z);// - fermi*I;
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
	//This block creates the SK tight binding Hamiltonians for onsite, first 
	//second and third neighbours for Au in fcc, MgO in halite and Fe in bcc
	vector<double> Au1, Au2, Fe_u1, Fe_u2, Fe_u3, Fe_d1, Fe_d2, Fe_d3, Mg1, Mg2, Mg3, O1, O2, O3;
	Au1.reserve(10); Au2.reserve(10); Fe_u1.reserve(10); Fe_u2.reserve(10); Fe_u3.reserve(10);
	Fe_d1.reserve(10); Fe_d2.reserve(10); Fe_d3.reserve(10);
	Mg1.reserve(10); Mg2.reserve(10); Mg3.reserve(10); O1.reserve(10); O2.reserve(10); O3.reserve(10);
	Au1 = param(1,1); Au2 = param(1,2); Fe_u1 = param(2,1); Fe_u2 = param(2,2); Fe_u3 = param(2,3);
	Fe_d1 = param(3,1); Fe_d2 = param(3,2); Fe_d3 = param(3,3); Mg1 = param(4,1); Mg2 = param(4,2);
	Mg3 = param(4,3); O1 = param(5,1); O2 = param(5,2); O3 = param(5,3);
	M9 Au, Fe_u, Fe_d, Mg, O;
	Fe_d = U(1,0);
	Fe_u = U(1,1);
	Au = U(2,0);
	Mg = U(3,0);
	O = U(4,0);

	//lattice vectors for the whole system;
	Vector3d lat_vec1, lat_vec2, Au_lat_oop;
	lat_vec1 << 0.5, 0, 0.5;
	lat_vec2 << 0.5, 0, -0.5;
	Au_lat_oop << 0.5, 0.5, 0;
	Vector3d Fe_lat_oop;
	Fe_lat_oop << 0.5, 0.25*M_SQRT2, 0.;
	//Distance information for n.n and n.n.n
	double Fe_nn_dist, Fe_nnn_dist, Fe_nnnn_dist;
	Fe_nn_dist = 0.5*sqrt(3.)/M_SQRT2;
	Fe_nnn_dist = 0.5*M_SQRT2;
	Fe_nnnn_dist = 1.;

	//Distance information for n.n and n.n.n
	double Au_nn_dist, Au_nnn_dist, Au_nnnn_dist;
	Au_nn_dist = M_SQRT2/2.;
	Au_nnn_dist = 1.;
	Au_nnnn_dist = 0;//this tells the code to ignore

	//This section defines the basis atoms and lattice vectors for MgO 
	Vector3d MgO_bas1, MgO_bas2;
	vec3 MgO_basis;
	MgO_basis.reserve(2);//magic 2 is number of subatoms
	MgO_bas1<< 0., 0., 0.;//Mg
	MgO_bas2<< 0.5, 0., 0.;//O
	MgO_basis.emplace_back(MgO_bas1);
	MgO_basis.emplace_back(MgO_bas2);
	//This section defines the out of plane lattice vector
	Vector3d MgO_lat_oop;
	MgO_lat_oop << 0.5, 0.5, 0;//TODO question mark as to whether the cosine vectors for T are correct e.g. -1.5, 0, 0.5 
	//Distance information for n.n and n.n.n
	double MgO_nn_dist, MgO_nnn_dist, MgO_nnnn_dist;
	MgO_nn_dist = 0.5;
	MgO_nnn_dist = M_SQRT2/2.;
	MgO_nnnn_dist = 0.5*sqrt(3.);

	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	Vector3d tmp_vec;
	vM iron_up, iron_dn, gold, magnesium_11, magnesium_12, magnesium_21, oxide_11, oxide_12, oxide_21;
	vec3 Au_pos, MgO_pos_11, MgO_pos_12, MgO_pos_21, Fe_pos;
	Matrix<dcomp, 9, 9> tmp_mat;
	double distance;
	for (int i1 = -2; i1 < 3; i1++){
		for (int i2 = -2; i2 < 3; i2++){
			for (int i3 = -2; i3 < 3; i3++){
				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*MgO_lat_oop + MgO_basis[0];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_pos_11.emplace_back(tmp_vec);
					oxide_11.emplace_back(O);
					magnesium_11.emplace_back(Mg);
				}
				else if (distance < MgO_nn_dist + 1e-3){
					MgO_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(O1, x, y, z);
					oxide_11.emplace_back(tmp_mat);
					tmp_mat = eint1(Mg1, x, y, z);
					magnesium_11.emplace_back(tmp_mat);
				}
				else if (distance < MgO_nnn_dist + 1e-3){
					MgO_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(O2, x, y, z);
					oxide_11.emplace_back(tmp_mat);
					tmp_mat = eint1(Mg2, x, y, z);
					magnesium_11.emplace_back(tmp_mat);
				}
				else if (distance < MgO_nnnn_dist + 1e-3){
					MgO_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(O3, x, y, z);
					oxide_11.emplace_back(tmp_mat);
					tmp_mat = eint1(Mg3, x, y, z);
					magnesium_11.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*MgO_lat_oop + MgO_basis[1];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_pos_12.emplace_back(tmp_vec);
					oxide_12.emplace_back(O);
					magnesium_12.emplace_back(Mg);
				}
				else if (distance < MgO_nn_dist + 1e-3){
					MgO_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Mg1, x, y, z);
					oxide_12.emplace_back(tmp_mat);
					magnesium_12.emplace_back(tmp_mat);
				}
				else if (distance < MgO_nnn_dist + 1e-3){
					MgO_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Mg2, x, y, z);
					oxide_12.emplace_back(tmp_mat);
					magnesium_12.emplace_back(tmp_mat);
				}
				else if (distance < MgO_nnnn_dist + 1e-3){
					MgO_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Mg3, x, y, z);
					oxide_12.emplace_back(tmp_mat);
					magnesium_12.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*MgO_lat_oop - MgO_basis[1];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_pos_21.emplace_back(tmp_vec);
					oxide_21.emplace_back(O);
					magnesium_21.emplace_back(Mg);
				}
				else if (distance < MgO_nn_dist + 1e-3){
					MgO_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Mg1, x, y, z);
					oxide_21.emplace_back(tmp_mat);
					magnesium_21.emplace_back(tmp_mat);
				}
				else if (distance < MgO_nnn_dist + 1e-3){
					MgO_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Mg2, x, y, z);
					oxide_21.emplace_back(tmp_mat);
					magnesium_21.emplace_back(tmp_mat);
				}
				else if (distance < MgO_nnnn_dist + 1e-3){
					MgO_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Mg3, x, y, z);
					oxide_21.emplace_back(tmp_mat);
					magnesium_21.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*Au_lat_oop;
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Au_pos.emplace_back(tmp_vec);
					gold.emplace_back(Au);
				}
				else if (distance < Au_nn_dist + 1e-3){
					Au_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Au1, x, y, z);
					gold.emplace_back(tmp_mat);
				}
				else if (distance < Au_nnn_dist + 1e-3){
					Au_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Au2, x, y, z);
					gold.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*Fe_lat_oop;
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Fe_pos.emplace_back(tmp_vec);
					iron_up.emplace_back(Fe_u);
					iron_dn.emplace_back(Fe_d);
				}
				else if (distance < Fe_nn_dist + 1e-3){
					Fe_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Fe_u1, x, y, z);
					iron_up.emplace_back(tmp_mat);
					tmp_mat = eint1(Fe_d1, x, y, z);
					iron_dn.emplace_back(tmp_mat);
				}
				else if (distance < Fe_nnn_dist + 1e-3){
					Fe_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Fe_u2, x, y, z);
					iron_up.emplace_back(tmp_mat);
					tmp_mat = eint1(Fe_d2, x, y, z);
					iron_dn.emplace_back(tmp_mat);
				}
				else if (distance < Fe_nnnn_dist + 1e-3){
					Fe_pos.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(Fe_u3, x, y, z);
					iron_up.emplace_back(tmp_mat);
					tmp_mat = eint1(Fe_d3, x, y, z);
					iron_dn.emplace_back(tmp_mat);
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
	Mydata += "-Au_fermi.txt";

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
	send.gold = &gold;
	send.Au_pos = &Au_pos;

	int size = 1000;
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

	/* int size = 4000; */
	/* int j = 0; */
	/* for (int k = 0; k < size; k++){ */
	/* 	/1* for (int j = 0; j < size; j++){ *1/ */
	/* 		for (int m = 0; m < size; m++){ */

	/* 			//This to create E(k_y) presently at neck */
	/* 			k_x = M_PI*k/(1.*size); */
	/* 			k_y = M_PI*m/(1.*size); */
	/* 			k_z = -M_PI;//*j/(1.*size); */

	/* 			//fully diagonalised Hamiltonian */
	/* 			E = H(Au_pos, gold, k_x, k_y, k_z);// - fermi*I; */
	/* 			/1* E = H(Fe_pos, iron_dn, k_x, k_y, k_z) - fermi*I; *1/ */
	/* 			/1* E2 = ins - fermi*Ibig; *1/ */
	/* 			/1* E = H(Fe_pos, iron_up, k_x, k_y, k_z) - fermi*I; *1/ */
	/* 			SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> es; */
	/* 			/1* SelfAdjointEigenSolver<Matrix<dcomp, 18, 18>> es; *1/ */
	/* 			es.compute(E, EigenvaluesOnly); */
	/* 			/1* es.compute(E2); *1/ */
	/* 			Matrix<double, 9, 1> O; */
	/* 			/1* Matrix<double, 18, 1> O; *1/ */
	/* 			O = es.eigenvalues(); */
	/* 			for (int l = 0; l < O.size(); l++){ */
	/* 				if (abs(O(l)) < 1e-2){ */

	/* 					Myfile<<k_x<<" "<<k_y<<endl; */
	/* 					Myfile<<-k_x<<" "<<k_y<<endl; */
	/* 					Myfile<<k_x<<" "<<-k_y<<endl; */
	/* 					Myfile<<-k_x<<" "<<-k_y<<endl; */

	/* 					/1* Myfile<<k_x<<" "<<k_y<<" "<<k_z<<endl; *1/ */
	/* 					/1* Myfile<<-k_x<<" "<<k_y<<" "<<k_z<<endl; *1/ */
	/* 					/1* Myfile<<k_x<<" "<<-k_y<<" "<<k_z<<endl; *1/ */
	/* 					/1* Myfile<<k_x<<" "<<k_y<<" "<<-k_z<<endl; *1/ */
	/* 					/1* Myfile<<-k_x<<" "<<-k_y<<" "<<k_z<<endl; *1/ */
	/* 					/1* Myfile<<k_x<<" "<<-k_y<<" "<<-k_z<<endl; *1/ */
	/* 					/1* Myfile<<-k_x<<" "<<k_y<<" "<<-k_z<<endl; *1/ */
	/* 					/1* Myfile<<-k_x<<" "<<-k_y<<" "<<-k_z<<endl; *1/ */

	/* 				} */
	/* 			} */
	/* 		/1* } *1/ */
	/* 	} */
	/* } */

	/* Myfile.close(); */
	/* return 0; */
/* } */
