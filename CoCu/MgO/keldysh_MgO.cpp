#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/src/Core/util/MKL_support.h>
/* #include "vector_integration.h" */
/* #include <gsl/gsl_integration.h> */
#include <nag.h>
#include <nagd01.h>
#include <nag_stdlib.h>
#include <vector>
#include "AuMgOFe.h"
#include <ctime>
#include "/home/alex/INTEL/impi/2019.1.144/intel64/include/mpi.h"
#include <iomanip>
#define EIGEN_DONT_PARALLELIZE
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

typedef struct
	{
		int lim;
		int lim2;
		int N;
		double Ef;
		double x;
		double z;
		double V;
		double kT;
		Vector3d *t;
		vec3 *basis;
		vec3 *pos;
		vM *cobalt_up;
		vM *cobalt_dn;
	        vM *copper;
	       	vM *cob_cop_up; 
		vM *cob_cop_dn;
		vM *ins_copper;
		vM *INS;
		vM *cob_ins_up;
		vM *cob_ins_dn;
		ddmat *NM;
		ddmat *NM_T;
		ddmat *FM_up;
		ddmat *FM_dn;
		ddmat *FM_T;
		ddmat *FM_NM_T;
		ddmat *odd_l1_up;
		ddmat *odd_l1_dn;
		ddmat *odd_l1_T1;
		ddmat *odd_l1_T2;
		ddmat *ins;
		ddmat *ins_T;
		ddmat *ins_FM_T;
		ddmat *ins_NM_T;

	}
variables;

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

dcomp fermi(const dcomp arg, const double Ef, const double kT){
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

vector<double> f_vec(const double theta, const dcomp E, variables * send, const ddmat &GL_UP_e, const ddmat &GL_DN_e, 
		const ddmat &GL_UP_o, const ddmat &GL_DN_o, const dddmat &Gr) {
	dcomp i = -1;
	i = sqrt(i);
	double V = send->V;
	double kT = send->kT;
	double Ef = send->Ef;
	int N = send->N;
	ddmat NM = *send->NM;
	ddmat NM_T = *send->NM_T;
	ddmat FM_T = *send->FM_T;
	ddmat FM_NM_T = *send->FM_NM_T;

	ddmat GL_up_even, GL_dn_even, GL_up_odd, GL_dn_odd;
	GL_up_even = GL_UP_e;
	GL_dn_even = GL_DN_e;
	GL_up_odd = GL_UP_o;
	GL_dn_odd = GL_DN_o;
	dddmat GR;
	ddmat NM_T_dagg;
	NM_T_dagg = NM_T.adjoint();
	ddmat FM_NM_T_dagg;
	FM_NM_T_dagg = FM_NM_T.adjoint();
	ddmat I = ddmat::Identity();
	dddmat S;
	ddmat S11, S12;
	S11 = cos(theta/2.)*I;
	S12 = sin(theta/2.)*I;
	S.topLeftCorner(18,18) = S11;
	S.topRightCorner(18,18) = S12;
	S.bottomLeftCorner(18,18) = -S12;
	S.bottomRightCorner(18,18) = S11;

	dddmat GL_even, GL_odd, GR_dagg;
	GL_even.fill(0.);
	GL_odd.fill(0.);
	GR = S.inverse()*Gr*S;

	ddmat OM = E*I;
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
	/* Pauli.topRightCorner(18,18) = I; */
	/* Pauli.bottomLeftCorner(18,18) = I; */

	/* dddmat Pauli = Ibig; */

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
	vector<double> result;
	result.reserve(N);
	//TODO at the moment, this is only accurate from N = 4...
	//because of gmean behaviour. See questions.txt
//adlayer layer 2 from layer 1 to spacer thickness, N
	for (int it=0; it < N/2; ++it){
		GL_even.topLeftCorner(18,18) = GL_up_even;
		GL_even.bottomRightCorner(18,18) = GL_dn_even;
		GL_odd.topLeftCorner(18,18) = GL_up_odd;
		GL_odd.bottomRightCorner(18,18) = GL_dn_odd;
		A_even = (Ibig-GR_T_dagg*GL_even*T).inverse();
		B_even = (Ibig-GR_dagg_T_dagg*GL_even.adjoint()*T).inverse();
		A_odd = (Ibig-GR_T_dagg*GL_odd*T).inverse();
		B_odd = (Ibig-GR_dagg_T_dagg*GL_odd.adjoint()*T).inverse();
		TOT_even = (B_even.adjoint()-A_even)*Pauli;
		TOT_odd = (B_odd.adjoint()-A_odd)*Pauli;
		spincurrent_even = .25*imag(TOT_even.trace());
		spincurrent_odd = .25*imag(TOT_odd.trace());
		result.emplace_back(spincurrent_even);
		result.emplace_back(spincurrent_odd);
		GL_up_even = (OM - (NM - V*I) -NM_T_dagg*GL_up_even*NM_T).inverse();
		GL_dn_even = (OM - (NM - V*I) -NM_T_dagg*GL_dn_even*NM_T).inverse();
		GL_up_odd = (OM - (NM - V*I) -NM_T_dagg*GL_up_odd*NM_T).inverse();
		GL_dn_odd = (OM - (NM - V*I) -NM_T_dagg*GL_dn_odd*NM_T).inverse();
	}
	return result;
}

void f(double * spincurrent, const double theta, const dcomp E, int k, Integer * needi, variables * send, const ddmat &GL_UP_even, 
		const ddmat &GL_DN_even, const ddmat &GL_UP_odd, const ddmat &GL_DN_odd, const dddmat &Gr) {
// ...NM|ins|FM(0)|NM(n)|FM(theta)...
	dcomp i = -1;
	i = sqrt(i);
	double V = send->V;
	double kT = send->kT;
	double Ef = send->Ef;
	int N = send->N;
	ddmat NM = *send->NM;
	ddmat NM_T = *send->NM_T;
	ddmat FM_T = *send->FM_T;
	ddmat FM_NM_T = *send->FM_NM_T;

	ddmat GL_up_even, GL_dn_even, GL_up_odd, GL_dn_odd;
	GL_up_even = GL_UP_even;
	GL_dn_even = GL_DN_even;
	GL_up_odd = GL_UP_odd;
	GL_dn_odd = GL_DN_odd;
	dddmat GR;
	ddmat NM_T_dagg;
	NM_T_dagg = NM_T.adjoint();
	ddmat FM_NM_T_dagg;
	FM_NM_T_dagg = FM_NM_T.adjoint();
	ddmat I = ddmat::Identity();
	dddmat S;
	ddmat S11, S12;
	S11 = cos(theta/2.)*I;
	S12 = sin(theta/2.)*I;
	S.topLeftCorner(18,18) = S11;
	S.topRightCorner(18,18) = S12;
	S.bottomLeftCorner(18,18) = -S12;
	S.bottomRightCorner(18,18) = S11;

	dddmat GL_even, GL_odd, GR_dagg;
	GL_even.fill(0.);
	GL_odd.fill(0.);
	GR = S.inverse()*Gr*S;

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

	ddmat OM = E*I;
	dddmat A, B, TOT;
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
	int kl;
	for (kl = N - 1; kl >= 0; kl--){
		if (needi[kl] == 1)
			break;
	}
	kl++;
	if (kl % 2 == 1)
		kl++;

	for (int it=0; it < kl/2; ++it){
		if (needi[2*it] == 1){
			GL_even.topLeftCorner(18,18) = GL_up_even;
			GL_even.bottomRightCorner(18,18) = GL_dn_even;
			A = (Ibig-GR_T_dagg*GL_even*T).inverse();
			B = (Ibig-GR_dagg_T_dagg*GL_even.adjoint()*T).inverse();
			tmp1 = B*GR_dagg_T_dagg;
			tmp2 = A*tmp1;
			tmp1 = T*tmp2;
			tmp2 = GL_even*tmp1;
			TOT = (tmp2-A*B+0.5*(A+B))*Pauli;
			spincurrent[2*it] = (1./(4.*M_PI))*real(TOT.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT)));
		}

		if (needi[2*it + 1] == 1){
			GL_odd.topLeftCorner(18,18) = GL_up_odd;
			GL_odd.bottomRightCorner(18,18) = GL_dn_odd;
			A= (Ibig-GR_T_dagg*GL_odd*T).inverse();
			B= (Ibig-GR_dagg_T_dagg*GL_odd.adjoint()*T).inverse();
			tmp1 = B*GR_dagg_T_dagg;
			tmp2 = A*tmp1;
			tmp1 = T*tmp2;
			tmp2 = GL_odd*tmp1;
			TOT = (tmp2-A*B+0.5*(A+B))*Pauli;
			spincurrent[2*it + 1] = (1./(4.*M_PI))*real(TOT.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT)));
		}

		if (it != kl/2 - 1){//saves wasting unused results
			GL_up_even = (OM - (NM - V*I) -NM_T_dagg*GL_up_even*NM_T).inverse();
			GL_dn_even = (OM - (NM - V*I) -NM_T_dagg*GL_dn_even*NM_T).inverse();
			GL_up_odd = (OM - (NM - V*I) -NM_T_dagg*GL_up_odd*NM_T).inverse();
			GL_dn_odd = (OM - (NM - V*I) -NM_T_dagg*GL_dn_odd*NM_T).inverse();
		}
	}
}

vector<double> int_theta(const dcomp E, variables * send) {
	vector<double> result;
	vector<double> integrate;
	ddmat FM_up = *send->FM_up;
	ddmat FM_dn = *send->FM_dn;
	ddmat NM = *send->NM;
	ddmat NM_T = *send->NM_T;
	ddmat FM_T = *send->FM_T;
	ddmat FM_NM_T = *send->FM_NM_T;
	ddmat odd_l1_up = *send->odd_l1_up;
	ddmat odd_l1_dn = *send->odd_l1_dn;
	ddmat odd_l1_T1 = *send->odd_l1_T1;
	ddmat odd_l1_T2 = *send->odd_l1_T2;
	ddmat ins = *send->ins;
	ddmat ins_T = *send->ins_T;
	ddmat ins_FM_T = *send->ins_FM_T;
	ddmat ins_NM_T = *send->ins_NM_T;
	/* cout<<NM<<endl<<endl; */
	double V = send->V;

	ddmat I = ddmat::Identity();
	ddmat OMup=E*I-(FM_up - V*I);
	ddmat OMdn=E*I-(FM_dn - V*I);
	ddmat OM = E*I;

	ddmat FM_T_dagg = FM_T.adjoint();
	ddmat GR_up = gs(OMup, FM_T_dagg);
	ddmat GR_dn = gs(OMdn, FM_T_dagg);
	ddmat FM_NM_T_dagg = FM_NM_T.adjoint();
	ddmat NM_T_dagg = NM_T.adjoint();
	
	/* //this for trilayer */
	/* ddmat GL_up = gs(OMup, FM_T); */
	/* ddmat GL_dn = gs(OMdn, FM_T); */

	//this below block for 5 layer
	ddmat GL_up_even = gs(OM - NM, NM_T);
	ddmat GL_dn_even = GL_up_even;
	//this below block for 5 layer
//lim is thickness of layer 2
	const int lim = send->lim;
	const int lim2 = send->lim2;
	ddmat ins_T_dagg = ins_T.adjoint();
	M9 Ismall = M9::Identity();
	ddmat ins_NM_T_dagg = ins_NM_T.adjoint();
//build thickness of layer 2 to lim layers
	for (int it=0; it < lim; ++it){//TODO the diagonal elements need to be shifted by the same amount in halites
		ins.topLeftCorner(9,9) = ins.topLeftCorner(9,9) - Ismall*(V*2.*it/(lim*2.));//TODO changed this so that full bias isn't on last layer so if lim = 1
		ins.bottomRightCorner(9,9) = ins.bottomRightCorner(9,9) - Ismall*(V*(2.*it + 1)/(lim*2.));// then shift = 0, 1/2 (so 1 falls on next layer)
		if (it == 0){
			GL_up_even = (OM - ins -ins_T_dagg*GL_up_even*ins_T).inverse();
			GL_dn_even = (OM - ins -ins_T_dagg*GL_dn_even*ins_T).inverse();
		}
		else {
			GL_up_even = (OM - ins -ins_NM_T_dagg*GL_up_even*ins_NM_T).inverse();
			GL_dn_even = (OM - ins -ins_NM_T_dagg*GL_dn_even*ins_NM_T).inverse();
		}
	}
//lim2 is thickness of layer 3
//build thickness of layer 3 to lim2 layers
//add 10 bilayers i.e. 20 layers of FM
	GL_up_even = (OM - (FM_up - V*I) -ins_FM_T.adjoint()*GL_up_even*ins_FM_T).inverse();
	GL_dn_even = (OM - (FM_dn - V*I) -ins_FM_T.adjoint()*GL_dn_even*ins_FM_T).inverse();
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
	dddmat GR;
	GR.fill(0.);
	GR.topLeftCorner(18,18) = GR_up;
	GR.bottomRightCorner(18,18) = GR_dn;
	int N = send->N;
	result.reserve(N);
	integrate.reserve(N);
	for (int i = 0; i < N; i++){
		result[i] = 0.;
		integrate[i] = 0.;
	}
	double theta;

	const int n = 10;
	/* const int n = 1; */
	for (int k=0; k<n+1; k++) {
		theta = k*M_PI/n;
		integrate = f_vec(theta, E, send, GL_up_even, GL_dn_even, GL_up_odd, GL_dn_odd, GR);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += M_PI*(0.5/n)*integrate[i];
			else 
				result[i] += (M_PI/n)*integrate[i];
		}
	}	
	return result;
}

void int_theta_E(const dcomp E, int k, double * fm, Integer * needi, variables * send) {
	double theta;
	int N = send->N;
	double result[N];
	ddmat FM_up = *send->FM_up;
	ddmat FM_dn = *send->FM_dn;
	ddmat NM = *send->NM;
	ddmat NM_T = *send->NM_T;
	ddmat FM_T = *send->FM_T;
	ddmat FM_NM_T = *send->FM_NM_T;
	ddmat odd_l1_up = *send->odd_l1_up;
	ddmat odd_l1_dn = *send->odd_l1_dn;
	ddmat odd_l1_T1 = *send->odd_l1_T1;
	ddmat odd_l1_T2 = *send->odd_l1_T2;
	ddmat ins = *send->ins;
	ddmat ins_T = *send->ins_T;
	ddmat ins_FM_T = *send->ins_FM_T;
	ddmat ins_NM_T = *send->ins_NM_T;
	/* cout<<NM<<endl<<endl; */
	double V = send->V;

	ddmat I = ddmat::Identity();
	ddmat OMup=E*I-(FM_up - V*I);
	ddmat OMdn=E*I-(FM_dn - V*I);
	ddmat OM = E*I;

	ddmat FM_T_dagg = FM_T.adjoint();
	ddmat GR_up = gs(OMup, FM_T_dagg);
	ddmat GR_dn = gs(OMdn, FM_T_dagg);
	ddmat FM_NM_T_dagg = FM_NM_T.adjoint();
	ddmat NM_T_dagg = NM_T.adjoint();
	
	/* //this for trilayer */
	/* ddmat GL_up = gs(OMup, FM_T); */
	/* ddmat GL_dn = gs(OMdn, FM_T); */

	//this below block for 5 layer
	ddmat GL_up_even = gs(OM - NM, NM_T);
	ddmat GL_dn_even = GL_up_even;
	//this below block for 5 layer
//lim is thickness of layer 2
	const int lim = send->lim;
	double lim2 = send->lim2;
	ddmat ins_T_dagg = ins_T.adjoint();
	M9 Ismall = M9::Identity();
	ddmat ins_NM_T_dagg = ins_NM_T.adjoint();
//build thickness of layer 2 to lim layers
	for (int it=0; it < lim; ++it){//TODO the diagonal elements need to be shifted by the same amount in halites
		ins.topLeftCorner(9,9) = ins.topLeftCorner(9,9) - Ismall*(V*2.*it/(lim*2.));//TODO changed this so that full bias isn't on last layer so if lim = 1
		ins.bottomRightCorner(9,9) = ins.bottomRightCorner(9,9) - Ismall*(V*(2.*it + 1)/(lim*2.));// then shift = 0, 1/2 (so 1 falls on next layer)
		if (it == 0){
			GL_up_even = (OM - ins -ins_T_dagg*GL_up_even*ins_T).inverse();
			GL_dn_even = (OM - ins -ins_T_dagg*GL_dn_even*ins_T).inverse();
		}
		else {
			GL_up_even = (OM - ins -ins_NM_T_dagg*GL_up_even*ins_NM_T).inverse();
			GL_dn_even = (OM - ins -ins_NM_T_dagg*GL_dn_even*ins_NM_T).inverse();
		}
	}
//lim2 is thickness of layer 3
//build thickness of layer 3 to lim2 layers
	GL_up_even = (OM - (FM_up - V*I) -ins_FM_T.adjoint()*GL_up_even*ins_FM_T).inverse();
	GL_dn_even = (OM - (FM_dn - V*I) -ins_FM_T.adjoint()*GL_dn_even*ins_FM_T).inverse();
	for (int it=0; it < lim2 - 1; ++it){
		GL_up_even = (OM - (FM_up - V*I) -FM_T_dagg*GL_up_even*FM_T).inverse();
		GL_dn_even = (OM - (FM_dn - V*I) -FM_T_dagg*GL_dn_even*FM_T).inverse();
	}

	ddmat odd_l1_T1_dagg = odd_l1_T1.adjoint();
	ddmat odd_l1_T2_dagg = odd_l1_T2.adjoint();
	ddmat GL_up_odd = GL_up_even;
	ddmat GL_dn_odd = GL_dn_even;
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
	dddmat GR;
	GR.fill(0.);
	GR.topLeftCorner(18,18) = GR_up;
	GR.bottomRightCorner(18,18) = GR_dn;

	//initialise integration so that they don't all get summed over every E!
	for (int ll = 0; ll < N; ll++){
		if (needi[ll] == 1)
			fm[ll + k] = 0.;
	}
	const int n = 10;
	/* const int n = 1; */
	for (int kk=0; kk<n+1; kk++) {
		theta = kk*M_PI/n;
		f(result, theta, E, k, needi, send, GL_up_even, GL_dn_even, GL_up_odd, GL_dn_odd, GR);
		if ((kk==0)||(kk==n)){
			for (int ll = 0; ll < N; ll++){
				if (needi[ll] == 1)
					fm[ll + k] += M_PI*(0.5/n)*result[ll];
			}
		}
		else {
			for (int ll = 0; ll < N; ll++){
				if (needi[ll] == 1)
					fm[ll + k] += (M_PI/n)*result[ll];
			}
		}
	}	
	/* NAG_FREE(result); */
}

void pass(const double E[], Integer nx, Integer ldfm, double * fm, Integer * needi, variables * send) {
	dcomp E_send;
	dcomp im = -1;
	im = sqrt(im);
	int ksd;
	for (int k = 0; k < nx; k++){
		/* cout<<setprecision(8)<<E[k]<<endl; */
		E_send = E[k] + 1e-6*im;//TODO Andrey has 1e-8 here
		ksd = k*ldfm;
		int_theta_E(E_send, ksd, fm, needi, send);
	}
}

vector<double> int_energy(variables * send) {
	Integer irevcm, lcmax, lcmin, lcom, ldfm, ldfmrq,
       		lenx, lenxrq, licmax, licmin, licom, liopts, lopts, ni, nx,
       		sdfm, sdfmrq, sid;
	  /* Arrays */
	char cvalue[17];
	double *com = 0, *dinest = 0, *errest = 0, *fm = 0, *opts = 0, *x = 0;
	Integer *icom = 0, *iopts = 0, *needi = 0;

	  /* NAG types */
	Nag_VariableType optype;
	NagError fail;

	  /* Setup phase. */
	  /* Set problem parameters. */
	ni = send->N;
	double Ef = send->Ef;
	double left = Ef;
	double right = Ef - send->V;

	double b = max(left,right) + 0.04375;
	double a = min(left,right) - 0.04375;
	
	liopts = 100;
	lopts = 100;
	if (!(opts = NAG_ALLOC((lopts), double)) || !(iopts = NAG_ALLOC((liopts), Integer))){
		cout<<"Allocation failure"<<endl;
		exit(EXIT_FAILURE);
	}

	INIT_FAIL(fail);
	/* Initialize option arrays using nag_quad_opt_set (d01zkc). */
	nag_quad_opt_set("Initialize = nag_quad_1d_gen_vec_multi_rcomm", iopts, liopts, opts, lopts, &fail);
	if (fail.code != NE_NOERROR) {
		cout<<"Error from nag_quad_opt_set (d01zkc)."<<endl<<fail.message<<endl;
		exit(EXIT_FAILURE);
	}
	nag_quad_opt_set("Quadrature Rule = gk15", iopts, liopts, opts, lopts, &fail);
	/* nag_quad_opt_set("Quadrature Rule = gk21", iopts, liopts, opts, lopts, &fail); */
	/* nag_quad_opt_set("Quadrature Rule = gk31", iopts, liopts, opts, lopts, &fail); */
	/* nag_quad_opt_set("Quadrature Rule = gk41", iopts, liopts, opts, lopts, &fail); */
	/* nag_quad_opt_set("Quadrature Rule = gk51", iopts, liopts, opts, lopts, &fail); */
	/* nag_quad_opt_set("Quadrature Rule = gk61", iopts, liopts, opts, lopts, &fail); */
	nag_quad_opt_set("Absolute Tolerance = 1.0e-6", iopts, liopts, opts, lopts, &fail);
	nag_quad_opt_set("Relative Tolerance = 1.0e-6", iopts, liopts, opts, lopts, &fail);

	/* Determine required array dimensions for
	 * nag_quad_1d_gen_vec_multi_rcomm (d01rac) using
	 * nag_quad_1d_gen_vec_multi_dimreq (d01rcc).
	 */
	nag_quad_1d_gen_vec_multi_dimreq(ni, &lenxrq, &ldfmrq, &sdfmrq, &licmin, &licmax, &lcmin, &lcmax,
                                   iopts, opts, &fail);
	if (fail.code != NE_NOERROR) {
		cout<<"Error from nag_quad_1d_gen_vec_multi_dimreq (d01rcc)."<<endl<<fail.message<<endl;
		exit(EXIT_FAILURE);
	}
	ldfm = ldfmrq;
	sdfm = sdfmrq;
	lenx = lenxrq;
	licom = licmax;
	lcom = lcmax;

	/* Allocate remaining arrays. */
	if (!(x = NAG_ALLOC((lenx), double)) ||	!(needi = NAG_ALLOC((ni), Integer)) || !(fm = NAG_ALLOC((ldfm) * (sdfm), double)) ||
		!(dinest = NAG_ALLOC((ni), double)) || !(errest = NAG_ALLOC((ni), double)) ||
	       	!(com = NAG_ALLOC((lcom), double)) || !(icom = NAG_ALLOC((licom), Integer))){
		cout<<"Allocation failure"<<endl;
		exit(EXIT_FAILURE);
	}

	/* Solve phase. */
	INIT_FAIL(fail);
	/* Set initial irevcm. */
	irevcm = 1;
	while (irevcm) {
		/* nag_quad_1d_gen_vec_multi_rcomm (d01rac).
		 * One-dimensional quadrature, adaptive, vectorized, multi-integral,
		 * reverse communication.
		 */
		nag_quad_1d_gen_vec_multi_rcomm(&irevcm, ni, a, b, &sid, needi, x, lenx, &nx, fm, ldfm,
                                    dinest, errest, iopts, opts, icom, licom, com, lcom, &fail);
		switch (irevcm) {
			case 11:
				/* Initial returns.
				 * These will occur during the non-adaptive phase.
				 * All values must be supplied.
				 * dinest and errest do not contain approximations over the complete
			 	 * interval at this stage.
				 */
				pass(x, nx, ldfm, fm, needi, send);
				break;
			case 12:
				/* Intermediate returns.
				 * These will occur during the adaptive phase.
				 * All requested values must be supplied.
				 * dinest and errest contain approximations over the complete
				 * interval at this stage.
				 */
				pass(x, nx, ldfm, fm, needi, send);
				break;
		}
	}
	if (fail.code != NE_NOERROR)
		cout<<"For x = "<<send->x<<", z = "<<send->z<<" "<<fail.message<<endl<<endl;

	vector<double> result;
	result.reserve(ni);
	for (int kk = 0; kk < ni; kk++)
		result.emplace_back(dinest[kk]);
	NAG_FREE(com);
	NAG_FREE(dinest);
	NAG_FREE(errest);
	NAG_FREE(fm);
	NAG_FREE(opts);
	NAG_FREE(x);
	NAG_FREE(icom);
	NAG_FREE(iopts);
	NAG_FREE(needi);
	return result;
}

vector<double> switching(variables * send) {
	double x = send->x;
	double z = send->z;
	vec3 basis = *send->basis;
	vec3 pos = *send->pos;
	Vector3d t = *send->t;
	vM cobalt_up = *send->cobalt_up;
	vM cobalt_dn = *send->cobalt_dn;
        vM copper = *send->copper;
       	vM cob_cop_up = *send->cob_cop_up; 
	vM cob_cop_dn = *send->cob_cop_dn;
	vM ins_copper = *send->ins_copper;
	vM INS = *send->INS;
	vM cob_ins_up = *send->cob_ins_up;
	vM cob_ins_dn = *send->cob_ins_dn;
	M9 ins_ii, ins_12, ins_21, ins_T_ii, ins_T_12, ins_T_21, ins_NM_12,
	   ins_NM_21, ins_NM_T_ii, ins_NM_T_12, ins_NM_T_21, ins_FM_12,
	   ins_FM_21, ins_FM_T_ii, ins_FM_T_12, ins_FM_T_21;

	//generate in plane Hamiltonians for simple bilayers
	ins_ii = InPlaneH(pos,  basis[0], INS, x, z);
	ins_12 = InPlaneH(pos,  basis[1], INS, x, z);
	ins_21 = InPlaneH(pos, -basis[1], INS, x, z);
	//generate hopping between simple bilayers of the same type
	ins_T_ii = InPlaneH(pos, t + basis[0], INS, x, z);
	ins_T_12 = InPlaneH(pos, t + basis[1], INS, x, z);
	ins_T_21 = InPlaneH(pos, t - basis[1], INS, x, z);
	//additional off diagonal Hamiltonians needed for bilayers
	//made of different atom types
	ins_NM_12 = InPlaneH(pos,  basis[1], ins_copper, x, z);
	ins_NM_21 = InPlaneH(pos, -basis[1], ins_copper, x, z);
	ins_NM_T_ii = InPlaneH(pos, t + basis[0], ins_copper, x, z);
	ins_NM_T_12 = InPlaneH(pos, t + basis[1], ins_copper, x, z);
	ins_NM_T_21 = InPlaneH(pos, t - basis[1], ins_copper, x, z);
	//TODO in this scheme so far hopping for spin up is the same as spin down
	ins_FM_12 = InPlaneH(pos,  basis[1], cob_ins_up, x, z);
	ins_FM_21 = InPlaneH(pos, -basis[1], cob_ins_up, x, z);
	ins_FM_T_ii = InPlaneH(pos, t + basis[0], cob_ins_up, x, z);
	ins_FM_T_12 = InPlaneH(pos, t + basis[1], cob_ins_up, x, z);
	ins_FM_T_21 = InPlaneH(pos, t - basis[1], cob_ins_up, x, z);
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
	ddmat ins, ins_T, ins_NM, ins_NM_T, ins_FM, ins_FM_T;

	ins.topLeftCorner(9,9) = ins_ii;
	ins.topRightCorner(9,9) = ins_12;
	ins.bottomLeftCorner(9,9) = ins_21;
	ins.bottomRightCorner(9,9) = ins_ii;
	ins_T.topLeftCorner(9,9) = ins_T_ii;
	ins_T.topRightCorner(9,9) = ins_T_12;
	ins_T.bottomLeftCorner(9,9) = ins_T_21;
	ins_T.bottomRightCorner(9,9) = ins_T_ii;
	ins_NM_T.topLeftCorner(9,9) = ins_NM_T_ii;
	ins_NM_T.topRightCorner(9,9) = ins_NM_T_12; 
	ins_NM_T.bottomLeftCorner(9,9) = ins_NM_T_21;
	ins_NM_T.bottomRightCorner(9,9) = ins_NM_T_ii;
	ins_FM_T.topLeftCorner(9,9) = ins_FM_T_ii;
	ins_FM_T.topRightCorner(9,9) = ins_FM_T_12; 
	ins_FM_T.bottomLeftCorner(9,9) = ins_FM_T_21;
	ins_FM_T.bottomRightCorner(9,9) = ins_FM_T_ii;

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

	send->NM = &NM;
	send->NM_T = &NM_T;
	send->FM_up = &FM_up;
	send->FM_dn = &FM_dn;
	send->FM_T = &FM_T;
	send->FM_NM_T = &FM_NM_T;
	send->odd_l1_up = &odd_l1_up;
	send->odd_l1_dn = &odd_l1_dn;
	send->odd_l1_T1 = &odd_l1_T1;
	send->odd_l1_T2 = &odd_l1_T2;
	send->ins = &ins;
	send->ins_T = &ins_T;
	send->ins_NM_T = &ins_NM_T;
	send->ins_FM_T = &ins_FM_T;

	vector<double> result1, result2, integrate;
	int N = send->N;
	double V = send->V;
	result1.reserve(N);
	if (abs(V) > 1e-9)
		result1 = int_energy(send);
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
	double kT = send->kT;
	for (int j=0; j!=15; j++){
		E = send->Ef + (2.*j + 1.)*kT*M_PI*i;
		integrate = int_theta(E, send);
		if (abs(V) < 1e-9){
			for (int l = 0; l < N; l++)
				result2[l] += 2.*kT*integrate[l]; 
		}
		else {
			for (int l = 0; l < N; l++)
				result2[l] += kT*integrate[l]; 
			E = send->Ef - V + (2.*j + 1.)*kT*M_PI*i;
			integrate = int_theta(E, send);
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

int main() 
{
	// plot output of spincurrent against energy
	const double Ef = 0.57553;
	// number of spacer layers
	int N = 10;
	// set bias
	double V = 0.05;
	//set number of insulator principle layers
	int lim = 2;//remember in halite structures the two basis atoms are on the same plane
	//set number of RH FM bilayers
	int lim2 = 10;
	
	const double k = 8.617e-5/13.6058;//boltzmann constant (in Ryds)
	const double T = 300;//set the temperature

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
	vector<double> AuMg1, AuMg2, AuMg3, AuO1, AuO2, AuO3, MgFe_u1, MgFe_u2, MgFe_u3, MgFe_d1, MgFe_d2,
		MgFe_d3, FeAu_u1, FeAu_u2, FeAu_u3, FeAu_d1, FeAu_d2, FeAu_d3, Mg_O2;
	AuMg1.reserve(10); AuMg2.reserve(10); AuMg3.reserve(10); AuO1.reserve(10); AuO2.reserve(10);
	AuO3.reserve(10); MgFe_u1.reserve(10); MgFe_u2.reserve(10); MgFe_u3.reserve(10); MgFe_d1.reserve(10); 
	MgFe_d2.reserve(10); MgFe_d3.reserve(10); FeAu_u1.reserve(10); FeAu_u2.reserve(10); FeAu_u3.reserve(10);
	FeAu_d1.reserve(10); FeAu_d2.reserve(10); FeAu_d3.reserve(10); Mg_O2.reserve(10);
	vector<double> zilch;
	zilch.reserve(10);
	for (int zz = 0; zz < 10; zz++)
		zilch.emplace_back(0.);
	double tmp;
	//This loop creates the geometric means used at the interfaces between elements
	for (int k = 0; k < 10; k++){
		tmp = gmean(Au1[k], Mg1[k]);
		AuMg1.emplace_back(tmp);
		tmp = gmean(Au2[k], Mg2[k]);
		AuMg2.emplace_back(tmp);
		tmp = gmean(zilch[k], Mg3[k]);
		AuMg3.emplace_back(tmp);
		tmp = gmean(Au1[k], O1[k]);
		AuO1.emplace_back(tmp);
		tmp = gmean(Au2[k], O2[k]);
		AuO2.emplace_back(tmp);
		tmp = gmean(zilch[k], O3[k]);
		AuO3.emplace_back(tmp);
		tmp = gmean(Mg1[k], Fe_u1[k]);
		MgFe_u1.emplace_back(tmp);
		tmp = gmean(Mg2[k], Fe_u2[k]);
		MgFe_u2.emplace_back(tmp);
		tmp = gmean(Mg3[k], Fe_u3[k]);
		MgFe_u3.emplace_back(tmp);
		tmp = gmean(Mg1[k], Fe_d1[k]);
		MgFe_d1.emplace_back(tmp);
		tmp = gmean(Mg2[k], Fe_d2[k]);
		MgFe_d2.emplace_back(tmp);
		tmp = gmean(Mg3[k], Fe_d3[k]);
		MgFe_d3.emplace_back(tmp);
		tmp = gmean(Au1[k], Fe_u1[k]);
		FeAu_u1.emplace_back(tmp);
		tmp = gmean(Au2[k], Fe_u2[k]);
		FeAu_u2.emplace_back(tmp);
		tmp = gmean(zilch[k], Fe_u3[k]);
		FeAu_u3.emplace_back(tmp);
		tmp = gmean(Au1[k], Fe_d1[k]);
		FeAu_d1.emplace_back(tmp);
		tmp = gmean(Au2[k], Fe_d2[k]);
		FeAu_d2.emplace_back(tmp);
		tmp = gmean(zilch[k], Fe_d3[k]);
		FeAu_d3.emplace_back(tmp);
		tmp = gmean(Mg2[k], O2[k]);
		Mg_O2.emplace_back(tmp);
	}

	//in-plane lattice vectors for the whole system;
	Vector3d lat_vec1, lat_vec2;
	lat_vec1 << 0.5, 0, 0.5;
	lat_vec2 << 0.5, 0, -0.5;

	//This section defines the basis atoms and lattice vectors for Au 
	Vector3d Au_bas1, Au_bas2;
	vec3 Au_basis;
	Au_basis.reserve(2);//magic 2 is number of subatoms
	Au_bas1<< 0., 0., 0.;
	Au_bas2<< 0.5, 0.5, 0.;
	Au_basis.emplace_back(Au_bas1);
	Au_basis.emplace_back(Au_bas2);
	Vector3d Au_lat_oop;
	//This section defines the out of plane lattice vector
	Au_lat_oop << 0, 1, 0;
	//Distance information for n.n and n.n.n
	double Au_nn_dist, Au_nnn_dist, Au_nnnn_dist;
	Au_nn_dist = M_SQRT2/2.;
	Au_nnn_dist = 1.;
	Au_nnnn_dist = 0;//this tells the code to ignore

	//This section defines the basis atoms and lattice vectors for MgO 
	Vector3d MgO_bas1, MgO_bas2, tmp_vec;
	vec3 MgO_basis;
	MgO_basis.reserve(2);//magic 2 is number of subatoms
	MgO_bas1<< 0., 0., 0.;//Mg
	MgO_bas2<< 0.5, 0., 0.;//O
	MgO_basis.emplace_back(MgO_bas1);
	MgO_basis.emplace_back(MgO_bas2);
	//This section defines the out of plane lattice vector
	Vector3d MgO_lat_oop;
	MgO_lat_oop << 0.5, 0.5, 0;
	//Distance information for n.n and n.n.n
	double MgO_nn_dist, MgO_nnn_dist, MgO_nnnn_dist;
	MgO_nn_dist = 0.5;
	MgO_nnn_dist = M_SQRT2/2.;
	MgO_nnnn_dist = 0.5*sqrt(3.);

	//This section defines the basis atoms and lattice vectors for Fe 
	Vector3d Fe_bas1, Fe_bas2;
	vec3 Fe_basis;
	Fe_basis.reserve(2);//magic 2 is number of subatoms
	Fe_bas1<< 0., 0., 0.;
	Fe_bas2<< 0.5, 0.25*M_SQRT2, 0.;
	Fe_basis.emplace_back(Fe_bas1);
	Fe_basis.emplace_back(Fe_bas2);
	//This section defines the out of plane lattice vector
	Vector3d Fe_lat_oop;
	Fe_lat_oop << 0., 0.5*M_SQRT2, 0.;
	//Distance information for n.n and n.n.n
	double Fe_nn_dist, Fe_nnn_dist, Fe_nnnn_dist;
	Fe_nn_dist = 0.5*sqrt(3.)/M_SQRT2;
	Fe_nnn_dist = 0.5*M_SQRT2;
	Fe_nnnn_dist = 1.;

	Vector3d lat_Au_MgO, lat_MgO_Fe, lat_Au_Fe, lat_Fe_Au;
	lat_Au_MgO << 0, 0.605 + 0.5, 0;// distance taken as 2.47 A from dx.doi.org/10.1021/jz4022975 | J. Phys. Chem. Lett. 2014, 5, 131−137
	lat_MgO_Fe << 0.5, 0.505, 0;// distance taken as 2.06 A from J. Phys. D: Appl. Phys. 42 (2009) 015003 (5pp) doi:10.1088/0022-3727/42/1/015003
	lat_Au_Fe << 0, 0.446 + 0.5, 0;// distance taken as 1.82 A from Surface Science 370 (1997)293-310
	lat_Fe_Au << 0, 0.446 + 0.25*M_SQRT2, 0;
	/* ins_met_lat_oop1 << 0., 1., 0.; */
	/* ins_met_lat_oop2 << 0., .5, 0.; */
	/* ins_met_lat_oop << 0., 0.6393, 0.;// this from LiCl paper detailing distance between LiCl and Co/Cu. */

	double Au_Fe_nn, Au_Fe_nnn, Au_Fe_nnnn, Au_MgO_nn, Au_MgO_nnn, Au_MgO_nnnn, MgO_Fe_nn, MgO_Fe_nnn, MgO_Fe_nnnn;
	Au_Fe_nn = 0.670012;
	Au_Fe_nnn = 0.946;//this is the larger of the two, as Fe to Au basis 1 is 0.7996 TODO is this right..?
	Au_Fe_nnnn = 1.181066;//this is the larger of the two, as the smaller is 1.06737; TODO is this right..?
	Au_MgO_nn = 0.605; // Au to O;
	Au_MgO_nnn = 0.784873;//Au to Mg
	Au_MgO_nnnn = 0.930605;//Au to O;
	MgO_Fe_nn = 0.505;//O to Fe
	MgO_Fe_nnn = 0.710652;//Mg to Fe
	MgO_Fe_nnnn = 0.868922;//O to Fe, but need to include the fact that basis 2 Fe is 0.8585534 from Mg!

	//This section generates the Hamiltonians from SK parameters and NN positions
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	vM iron_up, iron_dn, gold, iron_gold_up, iron_gold_dn, gold_MgO, magnesium, oxide, iron_up_MgO, iron_dn_MgO,
	   gold_iron_up, gold_iron_dn;
	vec3 Au_pos, MgO_pos, Fe_pos, Au_MgO_pos, MgO_Fe_pos, Au_Fe_pos, Fe_Au_pos;
	Au_pos.reserve(19); Fe_pos.reserve(19); MgO_pos.reserve(19); Au_MgO_pos.reserve(19); MgO_Fe_pos.reserve(19);
	Au_Fe_pos.reserve(19); Fe_Au_pos.reserve(19);
	iron_up.reserve(19); iron_dn.reserve(19); gold.reserve(19); iron_gold_up.reserve(19); iron_gold_dn.reserve(19);
	gold_MgO.reserve(19); magnesium.reserve(19); oxide.reserve(19); iron_up_MgO.reserve(19); iron_dn_MgO.reserve(19);
	gold_iron_up.reserve(19); gold_iron_dn.reserve(19);
	//magic 19 above is num onsite + num nn + num nnn = 1 + 12 + 6
	Matrix<dcomp, 9, 9> tmp_mat;
	double distance;
	for (int i1 = -1; i1 < 2; i1++){
		for (int i2 = -1; i2 < 2; i2++){
			for (int i3 = -1; i3 < 2; i3++){
				for (int i4 = 0; i4 < Au_basis.size(); i4++){
					tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*Au_lat_oop + Au_basis[i4];
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

					tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*MgO_lat_oop + MgO_basis[i4];
					distance = 0;
					for (int l = 0; l < 3; l++)
						distance += tmp_vec(l)*tmp_vec(l);
					distance = sqrt(distance);
					if (distance < 1e-5){
						MgO_pos.emplace_back(tmp_vec);
						oxide.emplace_back(O);
						magnesium.emplace_back(Mg);
					}
					else if (distance < MgO_nn_dist + 1e-3){
						MgO_pos.emplace_back(tmp_vec);
						x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						tmp_mat = eint1(Mg1, x, y, z);
						oxide.emplace_back(tmp_mat);
						magnesium.emplace_back(tmp_mat);
					}
					else if (distance < MgO_nnn_dist + 1e-3){
						MgO_pos.emplace_back(tmp_vec);
						x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						tmp_mat = eint1(Mg_O2, x, y, z);
						oxide.emplace_back(tmp_mat);
						magnesium.emplace_back(tmp_mat);
					}
					else if (distance < MgO_nnnn_dist + 1e-3){
						MgO_pos.emplace_back(tmp_vec);
						x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
						tmp_mat = eint1(Mg3, x, y, z);
						oxide.emplace_back(tmp_mat);
						magnesium.emplace_back(tmp_mat);
					}

					tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*Fe_lat_oop + Fe_basis[i4];
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

					for (int i5 = 0; i5 < Au_basis.size(); i5++){
						tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Au_MgO + MgO_basis[i4] - Au_basis[i5];
						distance = 0;
						for (int l = 0; l < 3; l++)
							distance += tmp_vec(l)*tmp_vec(l);
						distance = sqrt(distance);
						if (distance < 1e-5){
							cout<<"There probably shouldn't be anything here: Au-MgO"<<endl;
						}
						else if (distance < Au_MgO_nn + 1e-3){
							Au_MgO_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(AuO1, x, y, z);
							gold_MgO.emplace_back(tmp_mat);
						}
						else if (distance < Au_MgO_nnn + 1e-3){
							Au_MgO_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(AuMg2, x, y, z);
							gold_MgO.emplace_back(tmp_mat);
						}
						else if (distance < Au_MgO_nnnn + 1e-3){
							Au_MgO_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(AuO3, x, y, z);
							gold_MgO.emplace_back(tmp_mat);
						}

						tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_MgO_Fe + Fe_basis[i4] - MgO_basis[i5];
						distance = 0;
						for (int l = 0; l < 3; l++)
							distance += tmp_vec(l)*tmp_vec(l);
						distance = sqrt(distance);
						if (distance < 1e-5){
							cout<<"There probably shouldn't be anything here: MgO-Fe"<<endl;
						}
						else if (distance < MgO_Fe_nn + 1e-3){
							MgO_Fe_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(MgFe_u1, x, y, z);
							iron_up_MgO.emplace_back(tmp_mat);
							tmp_mat = eint1(MgFe_d1, x, y, z);
							iron_dn_MgO.emplace_back(tmp_mat);
						}
						else if (distance < MgO_Fe_nnn + 1e-3){
							MgO_Fe_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(MgFe_u1, x, y, z);
							iron_up_MgO.emplace_back(tmp_mat);
							tmp_mat = eint1(MgFe_d1, x, y, z);
							iron_dn_MgO.emplace_back(tmp_mat);
						}
						else if (distance < MgO_Fe_nnnn + 1e-3){
							MgO_Fe_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(MgFe_u1, x, y, z);
							iron_up_MgO.emplace_back(tmp_mat);
							tmp_mat = eint1(MgFe_d1, x, y, z);
							iron_dn_MgO.emplace_back(tmp_mat);
						}

						tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Fe_Au + Au_basis[i4] - Fe_basis[i5];
						distance = 0;
						for (int l = 0; l < 3; l++)
							distance += tmp_vec(l)*tmp_vec(l);
						distance = sqrt(distance);
						if (distance < 1e-5){
							cout<<"There probably shouldn't be anything here: Fe-Au"<<endl;
						}
						else if (distance < Au_Fe_nn + 1e-3){
							Fe_Au_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(FeAu_u1, x, y, z);
							iron_gold_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d1, x, y, z);
							iron_gold_dn.emplace_back(tmp_mat);
						}
						else if (distance < Au_Fe_nnn + 1e-3){
							Fe_Au_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(FeAu_u2, x, y, z);
							iron_gold_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d2, x, y, z);
							iron_gold_dn.emplace_back(tmp_mat);
						}
						else if (distance < Au_Fe_nnnn + 1e-3){
							Fe_Au_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(FeAu_u3, x, y, z);
							iron_gold_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d3, x, y, z);
							iron_gold_dn.emplace_back(tmp_mat);
						}

						tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Au_Fe + Fe_basis[i4] - Au_basis[i5];
						distance = 0;
						for (int l = 0; l < 3; l++)
							distance += tmp_vec(l)*tmp_vec(l);
						distance = sqrt(distance);
						if (distance < 1e-5){
							cout<<"There probably shouldn't be anything here: Au-Fe"<<endl;
						}
						else if (distance < Au_Fe_nn + 1e-3){
							Au_Fe_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(FeAu_u1, x, y, z);
							gold_iron_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d1, x, y, z);
							gold_iron_dn.emplace_back(tmp_mat);
						}
						else if (distance < Au_Fe_nnn + 1e-3){
							Au_Fe_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(FeAu_u2, x, y, z);
							gold_iron_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d2, x, y, z);
							gold_iron_dn.emplace_back(tmp_mat);
						}
						else if (distance < Au_Fe_nnnn + 1e-3){
							Au_Fe_pos.emplace_back(tmp_vec);
							x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
							tmp_mat = eint1(FeAu_u3, x, y, z);
							gold_iron_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d3, x, y, z);
							gold_iron_dn.emplace_back(tmp_mat);
						}
					}
				}
			}
		}
	}

      vector<string> atname;
      atname.reserve(4);
      atname.emplace_back("Fe");
      atname.emplace_back("Au");
      atname.emplace_back("Mg");
      atname.emplace_back("O");

//     whole cluster
      ofstream Myfile2;	
      string Mydata2 = "pos0.dat";
      Myfile2.open( Mydata2.c_str(),ios::trunc );

      int idum0=0;
      int cluster = 2 + lim + lim2 + N/2 + 2;
      vector<int>itype;
      for (int iii = 0; iii < 2; iii++)
	      itype.emplace_back(1);
      for (int iii = 0; iii < lim; iii++)
	      itype.emplace_back(2);
      for (int iii = 0; iii < lim2; iii++)
	      itype.emplace_back(0);
      for (int iii = 0; iii < N/2; iii++)
	      itype.emplace_back(1);
      for (int iii = 0; iii < 2; iii++)
	      itype.emplace_back(0);

      Vector3d rr;
      Vector3d counter;
      counter << 0, 0, 0;
      for (int ilay = 0; ilay < cluster; ilay++){
	for (int i3 = -5; i3 <= 3; i3++){
          for (int i1=-5; i1 <= 3; i1++){
            for (int i2= 0; i2 < Fe_basis.size(); i2++){
		    if ((ilay > 1) && (ilay < 2 + lim))
              		rr= i1*lat_vec1+i3*lat_vec2 + MgO_basis[i2] + counter;
		    else if (ilay < 2)
              		rr= i1*lat_vec1+i3*lat_vec2 + Au_basis[i2] + counter;
		    else if ((ilay > 1 + lim) && (ilay < 2 + lim + lim2))
              		rr= i1*lat_vec1+i3*lat_vec2 + Fe_basis[i2] + counter;
		    else if ((ilay > 1 + lim + lim2) && (ilay < 2 + lim + lim2 + N/2))
              		rr= i1*lat_vec1+i3*lat_vec2 + Au_basis[i2] + counter;
		    else
              		rr= i1*lat_vec1+i3*lat_vec2 + Fe_basis[i2] + counter;
	      if ((abs(rr(0)) < 1.3001) && (abs(rr(2)) < 1.3001))
       	         idum0++;
	    }
      	  }
	}
	if (ilay < 1)
		counter = counter + Au_lat_oop;
	else if (ilay == 1)
		counter = counter + lat_Au_MgO; 
	else if ((ilay > 1) && (ilay < 1 + lim))
		counter = counter + MgO_lat_oop;
	else if (ilay == 1 + lim)
		counter = counter + lat_MgO_Fe; 
	else if ((ilay > 1 + lim) && (ilay < 1 + lim + lim2))
		counter = counter + Fe_lat_oop; 
	else if (ilay == 1 + lim + lim2)
		counter = counter + lat_Fe_Au; 
	else if ((ilay > 1 + lim + lim2) && (ilay < 1 + lim + lim2 + N/2))
		counter = counter + Au_lat_oop; 
	else if (ilay == 1 + lim + lim2 + N/2)
		counter = counter + lat_Au_Fe; 
	else
		counter = counter + Fe_lat_oop; 
      }

      Myfile2<<idum0<<endl<<"foo"<<endl;

      counter << 0, 0, 0;
      for (int ilay = 0; ilay < cluster; ilay++){
	for (int i3 = -5; i3 <= 3; i3++){
          for (int i1=-5; i1 <= 3; i1++){
            for (int i2= 0; i2 < Fe_basis.size(); i2++){
		    if ((ilay > 1) && (ilay < 2 + lim))
              		rr= i1*lat_vec1+i3*lat_vec2 + MgO_basis[i2] + counter;
		    else if (ilay < 2)
              		rr= i1*lat_vec1+i3*lat_vec2 + Au_basis[i2] + counter;
		    else if ((ilay > 1 + lim) && (ilay < 2 + lim + lim2))
              		rr= i1*lat_vec1+i3*lat_vec2 + Fe_basis[i2] + counter;
		    else if ((ilay > 1 + lim + lim2) && (ilay < 2 + lim + lim2 + N/2))
              		rr= i1*lat_vec1+i3*lat_vec2 + Au_basis[i2] + counter;
		    else
              		rr= i1*lat_vec1+i3*lat_vec2 + Fe_basis[i2] + counter;
	      if ((abs(rr(0)) < 1.3001) && (abs(rr(2)) < 1.3001)){
		      if ((itype[ilay] == 2) && (i2 == 1))
               		Myfile2<<atname[itype[ilay]+1]<<" "<<4*rr.transpose()<<endl;
		      else
               		Myfile2<<atname[itype[ilay]]<<" "<<4*rr.transpose()<<endl;
	      }
	    }
      	  }
	}
	if (ilay < 1)
		counter = counter + Au_lat_oop;
	else if (ilay == 1)
		counter = counter + lat_Au_MgO; 
	else if ((ilay > 1) && (ilay < 1 + lim))
		counter = counter + MgO_lat_oop;
	else if (ilay == 1 + lim)
		counter = counter + lat_MgO_Fe; 
	else if ((ilay > 1 + lim) && (ilay < 1 + lim + lim2))
		counter = counter + Fe_lat_oop; 
	else if (ilay == 1 + lim + lim2)
		counter = counter + lat_Fe_Au; 
	else if ((ilay > 1 + lim + lim2) && (ilay < 1 + lim + lim2 + N/2))
		counter = counter + Au_lat_oop; 
	else if (ilay == 1 + lim + lim2 + N/2)
		counter = counter + lat_Au_Fe; 
	else
		counter = counter + Fe_lat_oop; 
      }
      Myfile2.close();

	/* double kT = k*T; */
	/* //set up the variables to send */
	variables send;
      send.t = &Au_lat_oop;
	/* send.kT = kT; */
	/* send.Ef = Ef; */
	/* send.x = 2.532374; //for now! */
	/* send.z = 2.532374; //for now! */
	/* send.pos = &pos; */
	/* send.basis = &basis; */
	/* send.copper = &copper; */
	/* send.cobalt_up = &cobalt_up; */
	/* send.cobalt_dn = &cobalt_dn; */
	/* send.cob_cop_up = &cob_cop_up; */
	/* send.cob_cop_dn = &cob_cop_dn; */
	/* send.ins_copper = &ins_copper; */
	/* send.INS = &INS; */
	/* send.cob_ins_up = &cob_ins_up; */
	/* send.cob_ins_dn = &cob_ins_dn; */
	/* send.V = V; */
	/* send.lim = lim; */
	/* send.lim2 = lim2; */
	/* send.N = N; */

	vector<double> answer;
	answer.reserve(N);

	time_t now = time(0);
	tm *ltm = localtime(&now);
	string Mydata;
	Mydata = to_string(ltm->tm_mday);
	Mydata += "-";
	Mydata += to_string(1+ltm->tm_mon);
	Mydata += "-";
	Mydata += to_string(1900+ltm->tm_year);

	ofstream Myfile;	

	if (abs(V) < 1e-4)
		Mydata += "-Keldysh_V0.txt";
	else
		Mydata += "-Keldysh_V.txt";
	/* answer = switching(&send); */

	/* vector<double> result; */
	/* vector<double> integrate; */
	/* result.reserve(N); */
	/* integrate.reserve(N); */
	/* for (int i = 0; i < N; i++) */
	/* 	result[i] = 0.; */
	/* /1* int n = 350;//set the number of k-points along x axis *1/ */
	/* int n = 2;//set the number of k-points along x axis */
	/* int counter = 0; */
	/* double factor = 2./(n*n); */
	/* int sumk = n*(n + 1)/2; */
	/* int p, q; */
	/* int product1; */
	/* int myid, numprocs; */
	/* time_t timer; */
	/* int kk, l, i; */
	/* int start_time = time(&timer); */
	/* MPI_Init(NULL,NULL); */
	/* MPI_Comm_size(MPI_COMM_WORLD, &numprocs); */
	/* MPI_Comm_rank(MPI_COMM_WORLD, &myid); */
	/* for (kk = 2*myid + 1; kk<2*sumk; kk+=2*numprocs){ */
	/* 	for (l = 0; l < n; l++){ */
	/* 		product1 = (2*n - l)*(l + 1); */
	/* 		if ( kk < product1){ */
	/* 			p = 2*(kk/(2*n - l)) + 1; */
	/* 			q = (l*(l + 1) + kk)%(2*n); */
	/* 			break; */
	/* 		} */
	/* 	} */
	/* 	send.x = M_PI*(p + q)/(2*n);//this translates the grid for fcc 1st BZ */
	/* 	send.z = M_PI*(p - q)/(2*n);//this translates the grid for fcc 1st BZ */
	/* 	integrate = switching(&send); */
	/* 	for (i = 0; i < N; i++){ */
	/* 		if ((p==1) && (q==1)) */
	/* 			result[i] += factor*0.5*integrate[i]; */
	/* 		else if (p==q) */
	/* 			result[i] += factor*0.5*integrate[i]; */
	/* 		else */
	/* 			result[i] += factor*integrate[i]; */
	/* 	} */
	/* 	counter++; */
	/* } */
	/* cout<<"process "<<myid<<" took "<<time(&timer)-start_time<<"s"<<endl; */
	/* cout<<"process "<<myid<<" performed "<<counter<<" computations"<<endl; */
	/* for (i = 0; i < N; i++) */
	/* 	MPI_Reduce(&result[i], &answer[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); */
	/* if (myid == 0){ */

		//magic 4 below due to this being the number of Cu planes before spincurrent is calculated
		Myfile.open( Mydata.c_str(),ios::trunc );
		for (int i = 0; i < N; i++)
			Myfile<<scientific<<i+4<<" "<<answer[i]<<endl;

	/* } */
	/* MPI_Finalize(); */

	return 0;
}
