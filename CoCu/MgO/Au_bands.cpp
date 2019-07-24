#include <iostream>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <ctime>
#include "AuMgOFe.h"

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<complex<double>, 9, 9> M9;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vec3;
typedef vector<Matrix<dcomp, 9, 9>, aligned_allocator<Matrix<dcomp, 9, 9>>> vM;
//calculates the bandstructure of fcc Cu

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
	Mydata += "-MgO_DOS.txt";

	ofstream Myfile;	

	Myfile.open( Mydata.c_str(),ios::trunc );

	/* double ryd = 13.6058; */
	double ryd = 1.;
	/* double fermi = 0.5074;//Fermi level for Au */
	/* double fermi = 0.57553; */
	/* double fermi = 0.7466;//Fermi level for Fe */
	double fermi = 0.;

	dcomp i;
	i = -1.;
	i = sqrt(i);

	Matrix<dcomp, 9, 9> E;
	Matrix<dcomp, 9, 9> I = Matrix<dcomp, 9, 9>::Identity();
	M9 ins_11, ins_12, ins_21, ins_22;
	Matrix<dcomp, 18, 18> ins, E2;

	int nt = 251;
	VectorXd vec1(nt), vec2(nt), vec3(nt), vec4(nt), vec5(nt), vec6(nt), vec7(nt), vec8(nt), vec9(nt),
		 vec10(nt), vec11(nt), vec12(nt), vec13(nt), vec14(nt), vec15(nt), vec16(nt), vec17(nt), vec18(nt);
	double k_x, k_y, k_z, pi;
	Vector3d K;
	for (int k = 0; k < 251; k++)
	{

		/* //This to create E(k_y) presently at neck */
		/* k_x = 2.532374; */
		/* k_z = 2.532374; */
		/* k_y = 2*M_PI*k/250.; */

		//This for conventional bandstructure
		if (k < 101){
			pi = 2.*M_PI*k/100.;
			k_x = pi;
			k_y = 0;
			k_z = 0;
		}
		if ((k > 100) && (k < 151)){
			pi = M_PI*(k-100)/50.;
			k_x = 2.*M_PI;
			k_y = 0;
			k_z = pi;
		}	
		if ((k > 150) && (k < 201)){
			pi = M_PI*(k-150)/50.;
			k_x = 2.*M_PI - pi;
			k_y = pi;
			k_z = M_PI;
		}
		if ((k > 200) && (k < 251)){
			pi = M_PI*(k-200)/50.;
			k_x = M_PI-pi;
			k_y = M_PI-pi;
			k_z = M_PI-pi;
		}

		/* K(0) = k_x; */
		/* K(1) = k_y; */
		/* K(2) = k_z; */

		//generate in plane Hamiltonians for simple bilayers
		ins_11 = InPlaneH(MgO_pos_11,  MgO_basis[0], magnesium_11, k_x, k_y, k_z);//TODO only 2nd NN and onsite depend on atom type, catered for in 11 and 22 only
		ins_12 = InPlaneH(MgO_pos_12,  MgO_basis[1], oxide_12, k_x, k_y, k_z);//in theory this should contain
		ins_21 = InPlaneH(MgO_pos_21, -MgO_basis[1], oxide_21, k_x, k_y, k_z);//the same hoppings as magnesium
		ins_22 = InPlaneH(MgO_pos_11,  MgO_basis[0], oxide_11, k_x, k_y, k_z);

		ins.topLeftCorner(9,9) = ins_11;
		ins.topRightCorner(9,9) = ins_12;
		ins.bottomLeftCorner(9,9) = ins_21;
		ins.bottomRightCorner(9,9) = ins_22;

		//fully diagonalised Hamiltonian
		/* E = H(Au_pos, gold, k_x, k_y, k_z) - fermi*I; */
		/* E = H(Fe_pos, iron_dn, k_x, k_y, k_z) - fermi*I; */
		E2 = ins; //- fermi*I;
		/* E = H(Fe_pos, iron_up, k_x, k_y, k_z) - fermi*I; */
		/* SelfAdjointEigenSolver<Matrix<dcomp, 9, 9>> es; */
		SelfAdjointEigenSolver<Matrix<dcomp, 18, 18>> es;
		/* es.compute(E); */
		es.compute(E2);
		/* Matrix<double, 9, 1> O; */
		Matrix<double, 18, 1> O;
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

		vec10(k) = O(9);
		vec11(k) = O(10);
		vec12(k) = O(11);
		vec13(k) = O(12);
		vec14(k) = O(13);
		vec15(k) = O(14);
		vec16(k) = O(15);
		vec17(k) = O(16);
		vec18(k) = O(17);

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

	Myfile<<endl;
	for (int k = 0; k < vec10.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec10(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec11.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec11(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec12.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec12(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec13.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec13(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec14.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec14(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec15.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec15(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec16.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec16(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec17.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec17(k)<<endl;
	Myfile<<endl;
	for (int k = 0; k < vec18.size(); k++)
		Myfile<<k/250.<<" "<<ryd*vec18(k)<<endl;

	Myfile.close();
	return 0;
}
