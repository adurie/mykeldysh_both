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

	M9 NM_ii, NM_12, NM_21, FM_up_ii, FM_up_12, FM_up_21, 
	   FM_dn_ii, FM_dn_12, FM_dn_21, NM_T_ii, NM_T_12, NM_T_21,
	   FM_T_ii, FM_T_12, FM_T_21, FM_NM_12, FM_NM_21, 
	   FM_NM_T_ii, FM_NM_T_12, FM_NM_T_21;
	//generate in plane Hamiltonians for simple bilayers
	NM_ii = InPlaneH(pos,  basis[0], copper, x, z);
	NM_12 = InPlaneH(pos,  basis[1], copper, x, z);
	NM_21 = InPlaneH(pos, -basis[1], copper, x, z);
	FM_up_ii = InPlaneH(pos,  basis[0], cobalt_up, x, z);
	FM_up_12 = InPlaneH(pos,  basis[1], cobalt_up, x, z);
	FM_up_21 = InPlaneH(pos, -basis[1], cobalt_up, x, z);
	FM_dn_ii = InPlaneH(pos,  basis[0], cobalt_dn, x, z);
	FM_dn_12 = InPlaneH(pos,  basis[1], cobalt_dn, x, z);
	FM_dn_21 = InPlaneH(pos, -basis[1], cobalt_dn, x, z);
	//generate hopping between simple bilayers of the same type
	NM_T_ii = InPlaneH(pos, t + basis[0], copper, x, z);
	NM_T_12 = InPlaneH(pos, t + basis[1], copper, x, z);
	NM_T_21 = InPlaneH(pos, t - basis[1], copper, x, z);
	//TODO in this scheme so far hopping for spin up is the same as spin down
	FM_T_ii = InPlaneH(pos, t + basis[0], cobalt_up, x, z);
	FM_T_12 = InPlaneH(pos, t + basis[1], cobalt_up, x, z);
	FM_T_21 = InPlaneH(pos, t - basis[1], cobalt_up, x, z);
	//additional off diagonal Hamiltonians needed for bilayers
	//made of different atom types
	//TODO in this scheme so far hopping for spin up is the same as spin down
	FM_NM_12 = InPlaneH(pos,  basis[1], cob_cop_dn, x, z);
	FM_NM_21 = InPlaneH(pos, -basis[1], cob_cop_dn, x, z);
	FM_NM_T_ii = InPlaneH(pos, t + basis[0], cob_cop_dn, x, z);
	FM_NM_T_12 = InPlaneH(pos, t + basis[1], cob_cop_dn, x, z);
	FM_NM_T_21 = InPlaneH(pos, t - basis[1], cob_cop_dn, x, z);

	ddmat FM_up, FM_dn, NM, NM_T, FM_T, FM_NM_T, odd_l1_up, odd_l1_dn, odd_l1_T1, odd_l1_T2;

	NM.topLeftCorner(9,9) = NM_ii;
	NM.topRightCorner(9,9) = NM_12;
	NM.bottomLeftCorner(9,9) = NM_21;
	NM.bottomRightCorner(9,9) = NM_ii;

	NM_T.topLeftCorner(9,9) = NM_T_ii;
	NM_T.topRightCorner(9,9) = NM_T_12;
	NM_T.bottomLeftCorner(9,9) = NM_T_21;
	NM_T.bottomRightCorner(9,9) = NM_T_ii;

	FM_up.topLeftCorner(9,9) = FM_up_ii;
	FM_up.topRightCorner(9,9) = FM_up_12;
	FM_up.bottomLeftCorner(9,9) = FM_up_21;
	FM_up.bottomRightCorner(9,9) = FM_up_ii;

	FM_dn.topLeftCorner(9,9) = FM_dn_ii;
	FM_dn.topRightCorner(9,9) = FM_dn_12;
	FM_dn.bottomLeftCorner(9,9) = FM_dn_21;
	FM_dn.bottomRightCorner(9,9) = FM_dn_ii;

	FM_T.topLeftCorner(9,9) = FM_T_ii;
	FM_T.topRightCorner(9,9) = FM_T_12;
	FM_T.bottomLeftCorner(9,9) = FM_T_21;
	FM_T.bottomRightCorner(9,9) = FM_T_ii;

	FM_NM_T.topLeftCorner(9,9) = FM_NM_T_ii;
	FM_NM_T.topRightCorner(9,9) = FM_NM_T_12; 
	FM_NM_T.bottomLeftCorner(9,9) = FM_NM_T_21;
	FM_NM_T.bottomRightCorner(9,9) = FM_NM_T_ii;

	//TODO in this scheme so far hopping for spin up is the same as spin down
	odd_l1_up.topLeftCorner(9,9) = FM_up_ii;
	odd_l1_up.topRightCorner(9,9) = FM_NM_12;
	odd_l1_up.bottomLeftCorner(9,9) = FM_NM_21;
	odd_l1_up.bottomRightCorner(9,9) = NM_ii;

	//TODO in this scheme so far hopping for spin up is the same as spin down
	odd_l1_dn.topLeftCorner(9,9) = FM_dn_ii;
	odd_l1_dn.topRightCorner(9,9) = FM_NM_12;
	odd_l1_dn.bottomLeftCorner(9,9) = FM_NM_21;
	odd_l1_dn.bottomRightCorner(9,9) = NM_ii;

	odd_l1_T1.topLeftCorner(9,9) = FM_T_ii;
	odd_l1_T1.topRightCorner(9,9) = FM_NM_T_12;
	odd_l1_T1.bottomLeftCorner(9,9) = FM_T_21;
	odd_l1_T1.bottomRightCorner(9,9) = FM_NM_T_ii;

	odd_l1_T2.topLeftCorner(9,9) = FM_NM_T_ii;
	odd_l1_T2.topRightCorner(9,9) = FM_NM_T_12;
	odd_l1_T2.bottomLeftCorner(9,9) = NM_T_21;
	odd_l1_T2.bottomRightCorner(9,9) = NM_T_ii;

	ddmat NM_T_dagg;
	NM_T_dagg = NM_T.adjoint();
	ddmat I = ddmat::Identity();
	dddmat S;
	ddmat S11, S12;
	S11 = cos(theta/2.)*I;
	S12 = sin(theta/2.)*I;
	S.topLeftCorner(18,18) = S11;
	S.topRightCorner(18,18) = S12;
	S.bottomLeftCorner(18,18) = -S12;
	S.bottomRightCorner(18,18) = S11;

	ddmat OMup=E*I-(FM_up - V*I);
	ddmat OMdn=E*I-(FM_dn - V*I);
	ddmat OM = E*I;

	ddmat FM_T_dagg = FM_T.adjoint();
	ddmat GR_up = gs(OMup, FM_T_dagg);
	ddmat GR_dn = gs(OMdn, FM_T_dagg);
	ddmat FM_NM_T_dagg = FM_NM_T.adjoint();
	
	/* //this for trilayer */
	/* ddmat GL_up_even = gs(OMup, FM_T); */
	/* ddmat GL_dn_even = gs(OMdn, FM_T); */

	//this below block for 5 layer
	ddmat GL_up_even = gs(OM - NM, NM_T);
	ddmat GL_dn_even = GL_up_even;
//lim is thickness of layer 2
	const int lim = 1;
	ddmat ins;
	M9 Ismall = M9::Identity();
	ins.fill(0.);
	ins.bottomRightCorner(9,9) = 5.*Ismall;
//build thickness of layer 2 to lim layers
//add one layer of artificial insulater (first layer is NM) TODO
	for (int it=0; it < lim; ++it){
		/* if (lim > 1) */
		/* 	ins = ins - I*(V*it/(lim*1.-1)); */
		//TODO in one band model this was one layer.. figure out how to do drop properly in bilayer
		//TODO perhaps V isn't being treated correctly throughout this section
		GL_up_even = (OM - (NM + ins) -NM_T_dagg*GL_up_even*NM_T).inverse();
		GL_dn_even = (OM - (NM + ins) -NM_T_dagg*GL_dn_even*NM_T).inverse();
	}
//lim2 is thickness of layer 3
	const int lim2 = 5;
//build thickness of layer 3 to lim2 layers
//add 5 bilayers i.e. 10 layers of FM
	GL_up_even = (OM - (FM_up - V*I) -FM_NM_T_dagg*GL_up_even*FM_NM_T).inverse();
	GL_dn_even = (OM - (FM_dn - V*I) -FM_NM_T_dagg*GL_dn_even*FM_NM_T).inverse();
	for (int it=0; it < lim2 - 1; ++it){
		GL_up_even = (OM - (FM_up - V*I) -FM_T_dagg*GL_up_even*FM_T).inverse();
		GL_dn_even = (OM - (FM_dn - V*I) -FM_T_dagg*GL_dn_even*FM_T).inverse();
	}

	ddmat GL_up_odd = GL_up_even;
	ddmat GL_dn_odd = GL_dn_even;
	ddmat odd_l1_T1_dagg = odd_l1_T1.adjoint();
	ddmat odd_l1_T2_dagg = odd_l1_T2.adjoint();
	//adlayer one bilayer onto LHS G_even to ensure gmean is correct
	//this means 2 layers are on before we begin!
	GL_up_even = (OM - (NM - V*I) -FM_NM_T_dagg*GL_up_even*FM_NM_T).inverse();
	GL_dn_even = (OM - (NM - V*I) -FM_NM_T_dagg*GL_dn_even*FM_NM_T).inverse();
	//adlayer one bilayer of CoCu onto LHS G for odd layers, then adlayer a 
	//further bilayer of Cu to ensure gmean is correct. This means 3 layers are on before we begin!
	GL_up_odd = (OM - (odd_l1_up - V*I) -odd_l1_T1_dagg*GL_up_odd*odd_l1_T1).inverse();
	GL_dn_odd = (OM - (odd_l1_dn - V*I) -odd_l1_T1_dagg*GL_dn_odd*odd_l1_T1).inverse();
	GL_up_odd = (OM - (NM - V*I) -odd_l1_T2_dagg*GL_up_odd*odd_l1_T2).inverse();
	GL_dn_odd = (OM - (NM - V*I) -odd_l1_T2_dagg*GL_dn_odd*odd_l1_T2).inverse();

	dddmat GR, GL_even, GL_odd, GR_dagg;
	GR.fill(0.);
	GL_even.fill(0.);
	GL_odd.fill(0.);
	GR.topLeftCorner(18,18) = GR_up;
	GR.bottomRightCorner(18,18) = GR_dn;
	GR = S.inverse()*GR*S;

	dddmat Ibig = dddmat::Identity();
	dddmat Tmean, Tmeandagg;
	Tmean.fill(0.);
	Tmean.topLeftCorner(18,18) = FM_NM_T;
	Tmean.bottomRightCorner(18,18) = FM_NM_T;
	Tmeandagg.topLeftCorner(18,18) = FM_NM_T_dagg;
	Tmeandagg.bottomRightCorner(18,18) = FM_NM_T_dagg;
	dddmat OMbig = E*Ibig;
	dddmat NMbig;
	NMbig.fill(0.);
	NMbig.topLeftCorner(18,18) = NM;
	NMbig.bottomRightCorner(18,18) = NM;
	//adlayer one bilayer onto RHS G to ensure gmean is correct
	//this means 2 layers are on before we begin!
	GR = (OMbig - (NMbig - V*Ibig)-Tmean*GR*Tmeandagg).inverse();
	GR_dagg = GR.adjoint();

	dddmat Pauli;//This is the y Pauli sigma Matrix
	Pauli.fill(0.);
	Pauli.topRightCorner(18,18) = -i*I;
	Pauli.bottomLeftCorner(18,18) = i*I;

	double spincurrent_even, spincurrent_odd;
	dddmat A_even, A_odd, B_even, B_odd, TOT_even, TOT_odd;
	dddmat T, Tdagg;
	T.fill(0.);
	T.topLeftCorner(18,18) = NM_T;
	T.bottomRightCorner(18,18) = NM_T;
	Tdagg = T.adjoint();
	dddmat GR_T_dagg, GR_dagg_T_dagg;
	GR_T_dagg = GR*Tdagg;
	GR_dagg_T_dagg = GR_dagg*Tdagg;
	dddmat tmp1, tmp2;
	//TODO at the moment, this is only accurate from N = 4...
	//because of gmean behaviour. See questions.txt
//adlayer layer 2 from layer 1 to spacer thickness, N
	vector<double> result;
	result.reserve(N);
	for (int it=0; it < N/2; ++it){
		GL_even.topLeftCorner(18,18) = GL_up_even;
		GL_even.bottomRightCorner(18,18) = GL_dn_even;
		GL_odd.topLeftCorner(18,18) = GL_up_odd;
		GL_odd.bottomRightCorner(18,18) = GL_dn_odd;
		A_even = (Ibig-GR_T_dagg*GL_even*T).inverse();
		B_even = (Ibig-GR_dagg_T_dagg*GL_even.adjoint()*T).inverse();
		A_odd = (Ibig-GR_T_dagg*GL_odd*T).inverse();
		B_odd = (Ibig-GR_dagg_T_dagg*GL_odd.adjoint()*T).inverse();
		if (myswitch == 0){
			tmp1 = B_even*GR_dagg_T_dagg;
			tmp2 = A_even*tmp1;
			tmp1 = T*tmp2;
			tmp2 = GL_even*tmp1;
			TOT_even = (tmp2-A_even*B_even+0.5*(A_even+B_even))*Pauli;
			tmp1 = B_odd*GR_dagg_T_dagg;
			tmp2 = A_odd*tmp1;
			tmp1 = T*tmp2;
			tmp2 = GL_odd*tmp1;
			TOT_odd = (tmp2-A_odd*B_odd+0.5*(A_odd+B_odd))*Pauli;
			spincurrent_even = (1./(2.*M_PI))*real(TOT_even.trace()*(fermi(E,Ef)-fermi(E,Ef-V)));
			spincurrent_odd = (1./(2.*M_PI))*real(TOT_odd.trace()*(fermi(E,Ef)-fermi(E,Ef-V)));
		}
		if (myswitch == 1){
			TOT_even = (B_even.adjoint()-A_even)*Pauli;
			TOT_odd = (B_odd.adjoint()-A_odd)*Pauli;
			spincurrent_even = .5*imag(TOT_even.trace());
			spincurrent_odd = .5*imag(TOT_odd.trace());
		}
		result.emplace_back(spincurrent_even);
		result.emplace_back(spincurrent_odd);
		GL_up_even = (OM - (NM - V*I) -NM_T_dagg*GL_up_even*NM_T).inverse();
		GL_dn_even = (OM - (NM - V*I) -NM_T_dagg*GL_dn_even*NM_T).inverse();
		GL_up_odd = (OM - (NM - V*I) -NM_T_dagg*GL_up_odd*NM_T).inverse();
		GL_dn_odd = (OM - (NM - V*I) -NM_T_dagg*GL_dn_odd*NM_T).inverse();
	}
	return result;
}

vector<double> int_theta(const double x, const double z, const double a, const dcomp E,
	       	const double Ef, const int N, const int myswitch, const double V, const Vector3d &t, const vec3 &pos, const vec3 &basis, 
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
		integrate = f(x, z, a, E, Ef, N, theta, myswitch, V, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += M_PI*(0.5/n)*integrate[i];
			else 
				result[i] += (M_PI/n)*integrate[i];
		}
	}	
	return result;
}

vector<double> int_energy(const double x, const double z, const double a, const double Ef, const int N, const double V, const Vector3d &t, 
		const vec3 &pos, const vec3 &basis, const vM &copper, const vM &cobalt_up, const vM &cobalt_dn, const vM &cob_cop_up, const vM &cob_cop_dn) {
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
	//TODO integration doesn't seem as expected... seem to require end > than present
	double end = Ef + 0.1;
	double start = Ef - V - 0.1;

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
		integrate = int_theta(x, z, a, E_send, Ef, N, 0, V, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += 0.5*factor*integrate[i];
			else 
				result[i] += factor*integrate[i];
		}
	}	

	return result;
}

vector<double> switching(const double x, const double z, const double a, const double Ef, const int N, const Vector3d &t, const vec3 &pos, const vec3 &basis, 
		const vM &copper, const vM &cobalt_up, const vM &cobalt_dn, const vM &cob_cop_up, const vM &cob_cop_dn) {
	vector<double> result1, result2, integrate;
	result1.reserve(N);
	double V = 0.0;
	/* double V = 0.3; */
	if (abs(V) > 1e-4)
		result1 = int_energy(x, z, a, Ef, N, V, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
	else {
		for (int l = 0; l < N; l++)
			result1[l] = 0.;
	}
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
		integrate = int_theta(x, z, a, E, Ef, N, 1, V, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
		if (abs(V) < 1e-5){
			for (int l = 0; l < N; l++)
				result2[l] += 2.*kT*integrate[l]; 
		}
		else {
			for (int l = 0; l < N; l++)
				result2[l] += kT*integrate[l]; 
			E = Ef - V + (2.*j + 1.)*kT*M_PI*i;
			integrate = int_theta(x, z, a, E, Ef, N, 1, V, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
			for (int l = 0; l < N; l++)
				result2[l] += kT*integrate[l]; 
		}
	}
	vector<double> total;
	total.reserve(N);
	for (int l = 0; l < N; l++)
		total[l] = result1[l] + result2[l];
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
	/* Mydata = "SC_fixed_k_no_V.txt"; */
	Mydata = "SC_fixed_k.txt";
	/* Mydata = "sc.txt"; */
	Myfile.open( Mydata.c_str(),ios::trunc );
	const double Ef = 0.57553;
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
	for (int k = -1; k < 2; k += 2){//TODO this needs to be automated - use lattice vectors
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
	int N = 30;
	vector<double> answer;
	answer.reserve(N);
	/* answer = int_theta(0, 0, 1,  0.1, Ef, N); */
	answer = switching(0, 0, 1, Ef, N, t, pos, basis, copper, cobalt_up, cobalt_dn, cob_cop_up, cob_cop_dn);
	/* answer = int_energy(0, 0, 1, Ef, N); */
	/* answer = int_kpoints(1, Ef, N); */
	/* answer = f(0, 0, 1, Ef, Ef, i, 0); */
	//magic 4 below due to this being the number of Cu planes before spincurrent is calculated
	for (int i = 0; i < N; i++){
		Myfile<<scientific<<i+4<<" "<<answer[i]<<endl;
	}
	return 0;
}
