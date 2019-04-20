#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/src/Core/util/MKL_support.h>
#include <vector>
#include "CoCuCo.h"
#define EIGEN_USE_MKL_ALL
// important note - presently the code assumes fcc only - integration
// in addition it is assumed SK for neighbour terms are the same for each spin

using namespace Eigen;
using namespace std;
typedef complex<double> dcomp;
typedef Matrix<complex<double>, 9, 9> M9;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vec3;
typedef vector<Matrix<dcomp, 9, 9>, aligned_allocator<Matrix<dcomp, 9, 9>>> vM;
typedef Matrix<dcomp, 18, 18> ddmat;
typedef Matrix<dcomp, 36, 36> dddmat;

//calculates g at the interfaces
ddmat gs(const ddmat &OM, const ddmat &T)
{
	ddmat zero = ddmat::Zero();
	ddmat Tinv;
	Tinv = T.inverse();
	Matrix<dcomp, 36, 36> X,O;
	X.topLeftCorner(18, 18) = zero;
	X.topRightCorner(18, 18) = Tinv;
	X.bottomLeftCorner(18, 18) = -T.adjoint();
	X.bottomRightCorner(18, 18) = OM*Tinv;
	ComplexEigenSolver<Matrix<dcomp, 36, 36>> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	ddmat b = O.topRightCorner(18, 18);
	ddmat d = O.bottomRightCorner(18, 18);
	ddmat GR;
	GR = b*d.inverse();
	return GR;
}

dcomp fermi(const dcomp arg, const double Ef){
	const double k = 8.617e-5/13.6058;
	const double T = 300;
	double kT = k*T;
	return 1./(1.+exp((arg-Ef)/kT));
}

M9 InPlaneH(const vec3 &pos, const Vector3d &basis, const vM &U, const double x, const double z){
	double max_dist = 1. + 1e-4; //TODO make general - this caters for 2nd NN only
	double distance;
	Vector3d K;
	dcomp i = -1;
	i = sqrt(i);
	K << x, 0., z;
	M9 result;
	result.fill(0.);
	Vector3d tmp_vec;
	for (int k = 0; k < pos.size(); k++){
		tmp_vec = pos[k] - basis;
		if (abs(tmp_vec(1)) < 1e-5){
			distance = 0;
			for (int l = 0; l < 3; l++)
				distance += tmp_vec(l)*tmp_vec(l);
			distance = sqrt(distance);
			if (distance <= max_dist)
				result = result + U[k]*exp(i*tmp_vec.dot(K));
		}
	}
	return result;
}

vector<double> f(const double x, const double z, const double a, const dcomp E, const double Ef, const int N,
	       	const double theta, const int myswitch, const double V,	const Vector3d &t, const vec3 &pos, const vec3 &basis, 
		const vM &copper, const vM &cobalt_up, const vM &cobalt_dn, const vM &cob_cop_up, const vM &cob_cop_dn) {
// ...NM|ins|FM(0)|NM(n)|FM(theta)...
	dcomp i = -1;
	i = sqrt(i);

	M9 U, U12, U21;
	U = InPlaneH(pos, basis[0], copper, x, z);
	U12 = InPlaneH(pos, basis[1], copper, x, z);
	U21 = InPlaneH(pos, -basis[1], copper, x, z);
	ddmat NM;
	NM.topLeftCorner(9,9) = U;
	NM.topRightCorner(9,9) = U12;
	NM.bottomLeftCorner(9,9) = U21;
	NM.bottomRightCorner(9,9) = U;

	U = InPlaneH(pos, t + basis[0], copper, x, z);
	U12 = InPlaneH(pos, t + basis[1], copper, x, z);
	U21 = InPlaneH(pos, t - basis[1], copper, x, z);
	ddmat NM_T;
	NM_T.topLeftCorner(9,9) = U;
	NM_T.topRightCorner(9,9) = U12;
	NM_T.bottomLeftCorner(9,9) = U21;
	NM_T.bottomRightCorner(9,9) = U;

	U = InPlaneH(pos, basis[0], cobalt_up, x, z);
	U12 = InPlaneH(pos, basis[1], cobalt_up, x, z);
	U21 = InPlaneH(pos, -basis[1], cobalt_up, x, z);
	ddmat FM_up;
	FM_up.topLeftCorner(9,9) = U;
	FM_up.topRightCorner(9,9) = U12;
	FM_up.bottomLeftCorner(9,9) = U21;
	FM_up.bottomRightCorner(9,9) = U;

	U = InPlaneH(pos, t + basis[0], cobalt_up, x, z);
	U12 = InPlaneH(pos, t + basis[1], cobalt_up, x, z);
	U21 = InPlaneH(pos, t - basis[1], cobalt_up, x, z);
	ddmat FM_T_up;
	FM_T_up.topLeftCorner(9,9) = U;
	FM_T_up.topRightCorner(9,9) = U12;
	FM_T_up.bottomLeftCorner(9,9) = U21;
	FM_T_up.bottomRightCorner(9,9) = U;

	U = InPlaneH(pos, basis[0], cobalt_dn, x, z);
	U12 = InPlaneH(pos, basis[1], cobalt_dn, x, z);
	U21 = InPlaneH(pos, -basis[1], cobalt_dn, x, z);
	ddmat FM_dn;
	FM_dn.topLeftCorner(9,9) = U;
	FM_dn.topRightCorner(9,9) = U12;
	FM_dn.bottomLeftCorner(9,9) = U21;
	FM_dn.bottomRightCorner(9,9) = U;

	//TODO in this scheme so far hopping for spin up is the same as spin down
	/* U = InPlaneH(pos, t + basis[0], cobalt_dn, x, z); */
	/* U12 = InPlaneH(pos, t + basis[1], cobalt_dn, x, z); */
	/* U21 = InPlaneH(pos, t - basis[1], cobalt_dn, x, z); */
	/* ddmat FM_T_dn; */
	/* FM_T_dn.topLeftCorner(9,9) = U; */
	/* FM_T_dn.topRightCorner(9,9) = U12; */
	/* FM_T_dn.bottomLeftCorner(9,9) = U21; */
	/* FM_T_dn.bottomRightCorner(9,9) = U; */

	//TODO in this scheme so far hopping for spin up is the same as spin down
	//TODO there is a very high chance this will need to be edited so only
	//differing atoms have gmean version! WE MAY NEED DIFFERING LAYERS FOR EACH INTERFACE!!
	U = InPlaneH(pos, basis[0], cob_cop_dn, x, z);
	U12 = InPlaneH(pos, basis[1], cob_cop_dn, x, z);
	U21 = InPlaneH(pos, -basis[1], cobalt_dn, x, z);
	ddmat FM_NM_T;
	FM_NM_T.topLeftCorner(9,9) = U;
	FM_NM_T.topRightCorner(9,9) = U12;
	FM_NM_T.bottomLeftCorner(9,9) = U21;
	FM_NM_T.bottomRightCorner(9,9) = U;
	U21 = InPlaneH(pos, -basis[1], copper, x, z);
	ddmat NM_FM_T;
	NM_FM_T.topLeftCorner(9,9) = U;
	NM_FM_T.topRightCorner(9,9) = U12;
	NM_FM_T.bottomLeftCorner(9,9) = U21;
	NM_FM_T.bottomRightCorner(9,9) = U;

	const double t = 0.5;
	const double nab = -0.175;
	Matrix2cd T, NM1, NM2, FM1, FM2, ins, I, S;
	//same hopping throughout for simplicity
	T << t,0.,0.,t;

	NM1 << -2.5, 0., 0., -2.5;
	NM2 = NM1;
	ins << 5.0, 0., 0., 5.0;
	FM1 << -2.425 + nab, 0., 0., -2.425 - nab;
	FM2 = FM1; 
	I << 1.,0.,0.,1.;
	FM1 = FM1 - I*V;
	S << cos(theta/2.),sin(theta/2.),-sin(theta/2.),cos(theta/2.);
	FM2 = FM2 - I*V;
	FM2 = S.inverse()*FM2*S;
//apply the bias to the RHS
	NM2 = NM2 - I*V;
	/* ins = ins - I*V; */
	/* ins = ins + I*V; */

	Matrix2cd OMV2=E*I-FM2;

	Matrix4cd X2,O2;
	X2 << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV2(0,0)/t,OMV2(0,1)/t,
		0,	-t,	OMV2(1,0)/t,OMV2(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces2;
	ces2.compute(X2);
	O2 = ces2.eigenvectors();
	Matrix2cd b2 = O2.topRightCorner(2,2);
	Matrix2cd d2 = O2.bottomRightCorner(2,2);
	Matrix2cd GR = b2*d2.inverse();
	Matrix2cd GR_dagg = GR.adjoint();
	Matrix2cd T_dagg = T.adjoint();

	Matrix2cd OM = E*I;

	/* Matrix2cd OMV1=E*I-NM1; */
	Matrix2cd OMV1=E*I-FM1;

	Matrix4cd X,O;
	X << 	0,	0,	1/t,	0,
		0,	0,	0,	1/t,
		-t,	0,	OMV1(0,0)/t,OMV1(0,1)/t,
		0,	-t,	OMV1(1,0)/t,OMV1(1,1)/t;
	ComplexEigenSolver<Matrix4cd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	Matrix2cd b = O.topRightCorner(2,2);
	Matrix2cd d = O.bottomRightCorner(2,2);
	Matrix2cd GL = b*d.inverse();
	/* cout<<"energy = "<<E<<endl<<endl; */
	/* cout<<"LHS SGF"<<endl; */
	/* cout<<GL<<endl<<endl; */
	/* cout<<"RHS SGF"<<endl; */
	/* cout<<GR<<endl<<endl; */

	Matrix2cd Pauli, xPau, yPau, zPau;
	xPau << 0.,1.,1.,0.;
	yPau << 0.,-i,i,0.;
	zPau << 1.,0.,0.,-1.;

	Pauli = yPau;
	double spincurrent;
	Matrix2cd A,B,TOT,GM;
//lim is thickness of layer 2
	const int lim = 1;
//build thickness of layer 2 to lim layers
	for (int it=0; it < lim; ++it){
		/* if (lim > 1) */
		/* 	ins = ins - I*(V*it/(lim*1.-1)); */
		GL = (OM - ins -T*GL*T).inverse();
	}
//lim2 is thickness of layer 3
	const int lim2 = 10;
//build thickness of layer 3 to lim2 layers
	for (int it=0; it < lim2; ++it){

		GL = (OM - FM1 -T*GL*T).inverse();
	}
//adlayer layer 2 from layer 1 to spacer thickness, N
	vector<double> result;
	result.reserve(N);
	/* if (myswitch == 1) */
	/* 	cout<<GR<<endl<<endl<<GL<<endl<<endl; */
	Matrix2cd tmp1, tmp2;
	for (int it=0; it < N; ++it){
		/* tmp1 = GL*T; */
		/* tmp2 = T_dagg*tmp1; */
		/* tmp1 = GR*tmp2; */
		/* tmp2 = I - tmp1; */
		/* A = tmp2.inverse(); */
		A = (I-GR*T_dagg*GL*T).inverse();
		/* tmp2 = GL.adjoint(); */
		/* tmp1 = tmp2*T; */
		/* tmp2 = T_dagg*tmp1; */
		/* tmp1 = GR_dagg*tmp2; */
		/* tmp2 = I - tmp1; */
		/* B = tmp2.inverse(); */
		B = (I-GR_dagg*T_dagg*GL.adjoint()*T).inverse();
		/* B = (I-T_dagg*GL*T*GR).inverse(); */
		if (myswitch == 0){
			TOT = (GL*T*A*B*GR_dagg*T_dagg-A*B+0.5*(A+B))*Pauli;
			spincurrent = (1./(2.*M_PI))*real(TOT.trace()*(fermi(E,Ef)-fermi(E,Ef-V)));
		}
		if (myswitch == 1){
			TOT = (B.adjoint()-A)*Pauli;
			spincurrent = .5*imag(TOT.trace());
		}
		result.emplace_back(spincurrent);
		GL = (OM - NM2 -T*GL*T).inverse();
		/* cout<<"N = "<<it<<endl; */
		/* cout<<"LHS SGF"<<endl; */
		/* cout<<GL<<endl; */
	}
	/* for (int ii = 0; ii < N; ii++) */
	/* 	result_tot[ii] = result1[ii] + result2[ii]; */
	return result;
	/* return result_tot; */
}

vector<double> int_theta(const double x, const double z, const double a, const dcomp E,
	       	const double Ef, const int N, const int myswitch, const double V, const vec3 &pos, const vec3 &basis, 
		const vM &copper, const vM &cobalt_up, const vM &cobalt_dn, const vM &cob_cop_up, const vM &cob_cop_dn) {
	vector<double> result;
	vector<double> integrate;
	result.reserve(N);
	integrate.reserve(N);
	for (int i = 0; i < N; i++)
		result[i] = 0.;
	double theta;

	const int n = 10;
	/* const int n = 1; */
	for (int k=0; k<n+1; k++) {
		theta = k*M_PI/n;
		integrate = f(x, z, a, E, Ef, N, theta, myswitch, V, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += M_PI*(0.5/n)*integrate[i];
			else 
				result[i] += (M_PI/n)*integrate[i];
		}
	}	
	return result;
}

vector<double> int_energy(const double x, const double z, const double a, const double Ef, const int N, const double V, const vec3 &pos, const vec3 &basis, 
		const vM &copper, const vM &cobalt_up, const vM &cobalt_dn, const vM &cob_cop_up, const vM &cob_cop_dn) {
	vector<double> result;
	vector<double> integrate;
	result.reserve(N);
	integrate.reserve(N);
	for (int i = 0; i < N; i++)
		result[i] = 0.;

	double E;
	dcomp E_send;
	dcomp im = -1;
	im = sqrt(im);
	double end = 0.1;
	double start = -0.4;

	/* gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000); */
	/* double error, ans; */
	/* gsl_function F; */
	/* F.function = &pass; */
	/* vector<double> args; */
	/* args.emplace_back(x); */
	/* args.emplace_back(z); */
	/* args.emplace_back(a); */
	/* args.emplace_back(Ef); */
	/* args.emplace_back(N); */
	/* vector<vector<double>> hold_all; */
	/* hold_all.emplace_back(result); */
	/* hold_all.emplace_back(args); */
	/* F.params = &hold_all; */
	/* gsl_integration_qags (&F, start, end, 0, 1e-2, 1000, w, &ans, &error); */
	/* gsl_integration_workspace_free (w); */

	const int n = 1000;
	double factor = (end - start)/(n*1.);
	for (int k=0; k<n+1; k++) {
		E = start + k*(end-start)/(n*1.);
		E_send = E + 1e-6*im;
		integrate = int_theta(x, z, a, E_send, Ef, N, 0, V, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += 0.5*factor*integrate[i];
			else 
				result[i] += factor*integrate[i];
		}
	}	

	return result;
}

vector<double> switching(const double x, const double z, const double a, const double Ef, const int N, const vec3 &pos, const vec3 &basis, 
		const vM &copper, const vM &cobalt_up, const vM &cobalt_dn, const vM &cob_cop_up, const vM &cob_cop_dn) {
	vector<double> result1, result2, integrate;
	result1.reserve(N);
	/* double V = 0.0; */
	double V = 0.3;
	result1 = int_energy(x, z, a, Ef, N, V, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
	integrate.reserve(N);
	result2.reserve(N);
	for (int l = 0; l < N; l++)
		result2[l] = 0.;
	dcomp i;
	i = -1.;
	i = sqrt(i);
	dcomp E = 0.;
	const double k = 8.617e-5/13.6058;//TODO this exist above as well..
	const double T = 300;// need to reconcile...
	double kT = k*T;
	for (int j=0; j!=15; j++){
		E = Ef + (2.*j + 1.)*kT*M_PI*i;
		integrate = int_theta(x, z, a, E, Ef, N, 1, V, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		for (int l = 0; l < N; l++)
			result2[l] += kT*integrate[l]; 
		E = Ef - V + (2.*j + 1.)*kT*M_PI*i;
		integrate = int_theta(x, z, a, E, Ef, N, 1, V, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		for (int l = 0; l < N; l++)
			result2[l] += kT*integrate[l]; 
	}
	vector<double> total;
	total.reserve(N);
	for (int l = 0; l < N; l++)
		total[l] = result1[l] + result2[l];
	/* return result2; */
	return total;
}

vector<double> int_kpoints(const double a, const double Ef, const int N){
	vector<double> result;
	vector<double> integrate;
	result.reserve(N);
	integrate.reserve(N);
	double x, z;
	for (int i = 0; i < N; i++)
		result[i] = 0.;

	int n = 20;
	int counter = 0;
	double factor = 8./(n*n);
	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = M_PI*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					counter++;
					z = M_PI*l/n;
					/* integrate = int_theta(x, z, 1, Ef, Ef, N); */
					/* integrate = int_energy(x, z, 1, Ef, N); */
					/* integrate = switching(x, z, 1, Ef, N); */
					for (int i = 0; i < N; i++){
						if ((k==1) && (l==1))
							result[i] += factor*0.5*integrate[i];
						else if (k==l)
							result[i] += factor*0.5*integrate[i];
						else
							result[i] += factor*integrate[i];
					}
					if (counter%4000 == 0)
						cout<<"     "<<counter*800./(n*n*1. + n)<<"% completed"<<endl;
				}
			}
		}
	}
	cout<<"     100% completed"<<endl;
	return result;
}

int main() 
{
	//number of atomic planes
	// plot output of spincurrent against energy
	string Mydata;
	ofstream Myfile;	
	/* Mydata = "AB.txt"; */
	Mydata = "sc_fixed_k.txt";
	/* Mydata = "sc.txt"; */
	Myfile.open( Mydata.c_str(),ios::trunc );
	const double Ef = 0.0;
	/* const double Ef = -0.3; */

	//This block creates the SK tight binding Hamiltonians for onsite, first 
	//and second neighbours for Co and Cu in fcc
	vector<double> Co1, Co2, Cu1, Cu2;
	Co1.reserve(10); Co2.reserve(10); Cu1.reserve(10); Cu2.reserve(10);
	Co1 = param(1,1); Co2 = param(1,2); Cu1 = param(2,1); Cu2 = param(2,2);
	M9 Co_u, Co_d, Cu;
	Co_d = U(1,0);
	Co_u = U(1,1);
	Cu = U(2,0);
	vector<double> CoCu1, CoCu2;
	CoCu1.reserve(10); CoCu2.reserve(10);
	double tmp;
	//This loop creates the geometric means used at the interfaces between elements
	for (int k = 0; k < 10; k++){
		tmp = gmean(Co1[k], Cu1[k]);
		CoCu1.emplace_back(tmp);
		tmp = gmean(Co2[k], Cu2[k]);
		CoCu2.emplace_back(tmp);
	}
	//This section defines the basis atoms 
	Vector3d bas1, bas2, tmp_vec;
	vec3 basis;
	basis.reserve(2);//magic 2 is number of subatoms
	bas1<< 0., 0., 0.;
	bas2<< 0.5, 0.5, 0.;
	Vector3d t; //TODO distance between principle layers
	t << 0., 1., 0.;
	basis.emplace_back(bas1);
	basis.emplace_back(bas2);
	//This section generates the Hamiltonians from SK parameters and NN positions
	double xx, yy;
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	vM cobalt_up, cobalt_dn, copper, cob_cop_up, cob_cop_dn;
	vec3 pos;
	pos.reserve(19);
	cobalt_up.reserve(19); cobalt_dn.reserve(19); copper.reserve(19); cob_cop_up.reserve(19);
	cob_cop_dn.reserve(19);
	tmp_vec << 0., 0., 0.;
	pos.emplace_back(tmp_vec);
	cobalt_up.emplace_back(Co_u);
	cobalt_dn.emplace_back(Co_d);
	copper.emplace_back(Cu);
	cob_cop_up.emplace_back(Cu); // In theory this will not be used
	cob_cop_dn.emplace_back(Cu); // In theory this will not be used
	//magic 19 above is num onsite + num nn + num nnn = 1 + 12 + 6
	Matrix<dcomp, 9, 9> tmp_mat;
	//This for 1st neighbours
	for (int k = -1; k < 2; k += 2){//TODO this needs to be automated
		for (int l = -1; l < 2; l += 2){
			xx = 0.5*k;
			yy = 0.5*l;
			tmp_vec << xx, yy, 0;
			pos.emplace_back(tmp_vec);
			x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			tmp_mat = eint1(Co1, x, y, z);
			cobalt_up.emplace_back(tmp_mat);
			cobalt_dn.emplace_back(tmp_mat);
			tmp_mat = eint1(Cu1, x, y, z);
			copper.emplace_back(tmp_mat);
			tmp_mat = eint1(CoCu1, x, y, z);
			cob_cop_up.emplace_back(tmp_mat);
			cob_cop_dn.emplace_back(tmp_mat);

			tmp_vec << 0, xx, yy;
			pos.emplace_back(tmp_vec);
			x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			tmp_mat = eint1(Co1, x, y, z);
			cobalt_up.emplace_back(tmp_mat);
			cobalt_dn.emplace_back(tmp_mat);
			tmp_mat = eint1(Cu1, x, y, z);
			copper.emplace_back(tmp_mat);
			tmp_mat = eint1(CoCu1, x, y, z);
			cob_cop_up.emplace_back(tmp_mat);
			cob_cop_dn.emplace_back(tmp_mat);

			tmp_vec << xx, 0, yy;
			pos.emplace_back(tmp_vec);
			x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
			tmp_mat = eint1(Co1, x, y, z);
			cobalt_up.emplace_back(tmp_mat);
			cobalt_dn.emplace_back(tmp_mat);
			tmp_mat = eint1(Cu1, x, y, z);
			copper.emplace_back(tmp_mat);
			tmp_mat = eint1(CoCu1, x, y, z);
			cob_cop_up.emplace_back(tmp_mat);
			cob_cop_dn.emplace_back(tmp_mat);
		}
	}

	//This for 2nd neighbours
	for (int k = -1; k < 2; k += 2){
		xx = 1.*k;
		tmp_vec << xx, 0, 0;
		pos.emplace_back(tmp_vec);
		x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		tmp_mat = eint1(Co2, x, y, z);
		cobalt_up.emplace_back(tmp_mat);
		cobalt_dn.emplace_back(tmp_mat);
		tmp_mat = eint1(Cu2, x, y, z);
		copper.emplace_back(tmp_mat);
		tmp_mat = eint1(CoCu2, x, y, z);
		cob_cop_up.emplace_back(tmp_mat);
		cob_cop_dn.emplace_back(tmp_mat);

		tmp_vec << 0, xx, 0;
		pos.emplace_back(tmp_vec);
		x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		tmp_mat = eint1(Co2, x, y, z);
		cobalt_up.emplace_back(tmp_mat);
		cobalt_dn.emplace_back(tmp_mat);
		tmp_mat = eint1(Cu2, x, y, z);
		copper.emplace_back(tmp_mat);
		tmp_mat = eint1(CoCu2, x, y, z);
		cob_cop_up.emplace_back(tmp_mat);
		cob_cop_dn.emplace_back(tmp_mat);

		tmp_vec << 0, 0, xx;
		pos.emplace_back(tmp_vec);
		x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
		tmp_mat = eint1(Co2, x, y, z);
		cobalt_up.emplace_back(tmp_mat);
		cobalt_dn.emplace_back(tmp_mat);
		tmp_mat = eint1(Cu2, x, y, z);
		copper.emplace_back(tmp_mat);
		tmp_mat = eint1(CoCu2, x, y, z);
		cob_cop_up.emplace_back(tmp_mat);
		cob_cop_dn.emplace_back(tmp_mat);
	}

	/* M9 U, U12, U21; */
	/* U = InPlaneH(pos, basis[0], copper, 0.8, 2.2); */
	/* U12 = InPlaneH(pos, basis[1], copper, 0.8, 2.2); */
	/* U21 = InPlaneH(pos, -basis[1], copper, 0.8, 2.2); */
	/* ddmat UU; */
	/* UU.topLeftCorner(9,9) = U; */
	/* UU.topRightCorner(9,9) = U12; */
	/* UU.bottomLeftCorner(9,9) = U21; */
	/* UU.bottomRightCorner(9,9) = U; */
	/* /1* cout<<UU.real()<<endl<<endl; *1/ */
	/* cout<<UU.real()<<endl<<endl; */

	// number of spacer layers
	int N = 11;
	vector<double> answer;
	answer.reserve(N);
	/* answer = int_theta(0, 0, 1,  0.1, Ef, N); */
	/* answer = switching(0, 0, 1, Ef, N, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn); */
	/* answer = int_energy(0, 0, 1, Ef, N); */
	/* answer = int_kpoints(1, Ef, N); */
	/* answer = f(0, 0, 1, Ef, Ef, i, 0); */
	/* for (int i = 1; i < N; i++){ */
	/* 	Myfile<<scientific<<i<<" "<<answer[i]<<endl; */
	/* } */
	return 0;
}
