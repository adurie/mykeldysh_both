#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/src/Core/util/MKL_support.h>
/* #include "vector_integration.h" */
#include <gsl/gsl_integration.h>
#include <vector>
#include "CoCuCo.h"
#include <ctime>
#include "/home/alex/INTEL/impi/2019.1.144/intel64/include/mpi.h"
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
typedef vector<Matrix<dcomp, 18, 18>, aligned_allocator<Matrix<dcomp, 18, 18>>> vMM;
typedef vector<Matrix<dcomp, 36, 36>, aligned_allocator<Matrix<dcomp, 36, 36>>> vMMM;
typedef Matrix<dcomp, 18, 18> ddmat;
typedef Matrix<dcomp, 36, 36> dddmat;

typedef struct
	{
		int N;
		int E_N;
		double Ef;
		double x;
		double z;
		double V;
		double kT;
		vec3 *basis;
		vec3 *pos;
		Vector3d *t;
		vM *cobalt_up;
		vM *cobalt_dn;
	        vM *copper;
	       	vM *cob_cop_up; 
		vM *cob_cop_dn;
		vector<int> *N_old;
		vMM *GL_U;
		vMM *GL_D;
		vMMM *GR;
		ddmat *NM_T;
		ddmat *NM;
		ddmat *FM_NM_T;
		vector<double> *E_list;
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

double f(const double theta, const dcomp E, variables * send, const int myswitch, vector<double> &result) {
// ...NM|ins|FM(0)|NM(n)|FM(theta)...
	dcomp i = -1;
	i = sqrt(i);
	double x = send->x;
	double z = send->z;
	double V = send->V;
	int N = send->N;
	double kT = send->kT;
	vec3 basis = *send->basis;
	double Ef = send->Ef;
	vec3 pos = *send->pos;
	Vector3d t = *send->t;
	vM cobalt_up = *send->cobalt_up;
	vM cobalt_dn = *send->cobalt_dn;
        vM copper = *send->copper;
       	vM cob_cop_up = *send->cob_cop_up; 
	vM cob_cop_dn = *send->cob_cop_dn;
	int E_N = send->E_N;

	/* cout<<"x = "<<x<<" z = "<<z<<" V = "<<V<<" N = "<<N<<" Ef = "<<Ef<<endl; */

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
	const int lim = 10;
	ddmat ins;
	M9 Ismall = M9::Identity();
	ins.fill(0.);
	ins.topLeftCorner(9,9) = 500.*Ismall;
	ins.bottomRightCorner(9,9) = 500.*Ismall;
//build thickness of layer 2 to lim layers
//add ten bilayers of artificial insulater 
	for (int it=0; it < lim; ++it){
		if (lim > 1){
			ins.topLeftCorner(9,9) = ins.topLeftCorner(9,9) - Ismall*(V*2.*it/(lim*2.-1));
			ins.bottomRightCorner(9,9) = ins.bottomRightCorner(9,9) - Ismall*(V*(2.*it + 1)/(lim*2.-1));
		}
		GL_up_even = (OM - (NM + ins) -NM_T_dagg*GL_up_even*NM_T).inverse();
		GL_dn_even = (OM - (NM + ins) -NM_T_dagg*GL_dn_even*NM_T).inverse();
	}
//lim2 is thickness of layer 3
	const int lim2 = 10;
//build thickness of layer 3 to lim2 layers
//add 10 bilayers i.e. 20 layers of FM
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
	if ((myswitch == 1) || (E_N % 2 == 1)){
		GL_up_odd = (OM - (odd_l1_up - V*I) -odd_l1_T1_dagg*GL_up_odd*odd_l1_T1).inverse();
		GL_dn_odd = (OM - (odd_l1_dn - V*I) -odd_l1_T1_dagg*GL_dn_odd*odd_l1_T1).inverse();
		GL_up_odd = (OM - (NM - V*I) -odd_l1_T2_dagg*GL_up_odd*odd_l1_T2).inverse();
		GL_dn_odd = (OM - (NM - V*I) -odd_l1_T2_dagg*GL_dn_odd*odd_l1_T2).inverse();
	}

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
	Tmeandagg.fill(0.);
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
	double spincurrent = 0;
	if (myswitch == 0){
		if (E_N % 2 == 0){
			for (int kk = 0; kk < E_N/2; kk++){
				GL_up_even = (OM - (NM - V*I) -NM_T_dagg*GL_up_even*NM_T).inverse();
				GL_dn_even = (OM - (NM - V*I) -NM_T_dagg*GL_dn_even*NM_T).inverse();
			}
			GL_even.topLeftCorner(18,18) = GL_up_even;
			GL_even.bottomRightCorner(18,18) = GL_dn_even;
			A_even = (Ibig-GR_T_dagg*GL_even*T).inverse();
			B_even = (Ibig-GR_dagg_T_dagg*GL_even.adjoint()*T).inverse();
			tmp1 = B_even*GR_dagg_T_dagg;
			tmp2 = A_even*tmp1;
			tmp1 = T*tmp2;
			tmp2 = GL_even*tmp1;
			TOT_even = (tmp2-A_even*B_even+0.5*(A_even+B_even))*Pauli;
			spincurrent = (-1./(4.*M_PI))*real(TOT_even.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT)));
		}
		else if (E_N % 2 == 1){
			for (int kk = 0; kk < E_N/2; kk++){
				GL_up_odd = (OM - (NM - V*I) -NM_T_dagg*GL_up_odd*NM_T).inverse();
				GL_dn_odd = (OM - (NM - V*I) -NM_T_dagg*GL_dn_odd*NM_T).inverse();
			}
			GL_odd.topLeftCorner(18,18) = GL_up_odd;
			GL_odd.bottomRightCorner(18,18) = GL_dn_odd;
			A_odd = (Ibig-GR_T_dagg*GL_odd*T).inverse();
			B_odd = (Ibig-GR_dagg_T_dagg*GL_odd.adjoint()*T).inverse();
			tmp1 = B_odd*GR_dagg_T_dagg;
			tmp2 = A_odd*tmp1;
			tmp1 = T*tmp2;
			tmp2 = GL_odd*tmp1;
			TOT_odd = (tmp2-A_odd*B_odd+0.5*(A_odd+B_odd))*Pauli;
			spincurrent = (-1./(4.*M_PI))*real(TOT_odd.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT)));
		}
		else
			cout<<"error here"<<endl;
	}
	if (myswitch == 1){
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
			spincurrent_even = -.25*imag(TOT_even.trace());
			spincurrent_odd = -.25*imag(TOT_odd.trace());
			result.emplace_back(spincurrent_even);
			result.emplace_back(spincurrent_odd);
			GL_up_even = (OM - (NM - V*I) -NM_T_dagg*GL_up_even*NM_T).inverse();
			GL_dn_even = (OM - (NM - V*I) -NM_T_dagg*GL_dn_even*NM_T).inverse();
			GL_up_odd = (OM - (NM - V*I) -NM_T_dagg*GL_up_odd*NM_T).inverse();
			GL_dn_odd = (OM - (NM - V*I) -NM_T_dagg*GL_dn_odd*NM_T).inverse();
		}
	}
	return spincurrent;
}

/* double mem(const double theta, int N_old, variables * send, ddmat &GL_up, ddmat &GL_dn, dddmat &GR, dcomp E, ddmat &NM_T, ddmat &NM, ddmat &FM_NM_T){ */
/* 	dcomp i = -1; */
/* 	i = sqrt(i); */
/* 	double x = send->x; */
/* 	double z = send->z; */
/* 	double V = send->V; */
/* 	int N = send->N; */
/* 	double kT = send->kT; */
/* 	double Ef = send->Ef; */
/* 	int E_N = send->E_N; */
/* 	ddmat NM_T_dagg;//TODO */
/* 	NM_T_dagg = NM_T.adjoint();//TODO */
/* 	ddmat I = ddmat::Identity();//TODO */
/* 	dddmat S; */
/* 	ddmat S11, S12; */
/* 	S11 = cos(theta/2.)*I; */
/* 	S12 = sin(theta/2.)*I; */
/* 	S.topLeftCorner(18,18) = S11; */
/* 	S.topRightCorner(18,18) = S12; */
/* 	S.bottomLeftCorner(18,18) = -S12; */
/* 	S.bottomRightCorner(18,18) = S11; */
/* 	ddmat OM = E*I; */
/* 	GR = S.inverse()*GR*S; */
/* 	dddmat Ibig = dddmat::Identity();//TODO */
/* 	dddmat OMbig = E*Ibig; */
/* 	dddmat NMbig;//TODO */
/* 	NMbig.fill(0.);//TODO */
/* 	NMbig.topLeftCorner(18,18) = NM;//TODO */
/* 	NMbig.bottomRightCorner(18,18) = NM;//TODO */
/* 	dddmat Tmean, Tmeandagg;//TODO */
/* 	Tmean.fill(0.);//TODO */
/* 	Tmeandagg.fill(0.);//TODO */
/* 	ddmat FM_T_dagg = FM_T.adjoint();//TODO */
/* 	Tmean.topLeftCorner(18,18) = FM_NM_T;//TODO */
/* 	Tmean.bottomRightCorner(18,18) = FM_NM_T;//TODO */
/* 	Tmeandagg.topLeftCorner(18,18) = FM_NM_T_dagg;//TODO */
/* 	Tmeandagg.bottomRightCorner(18,18) = FM_NM_T_dagg;//TODO */
/* 	//adlayer one bilayer onto RHS G to ensure gmean is correct */
/* 	//this means 2 layers are on before we begin! */
/* 	GR = (OMbig - (NMbig - V*Ibig)-Tmean*GR*Tmeandagg).inverse(); */
/* 	GR_dagg = GR.adjoint(); */

	/* dddmat Pauli;//This is the y Pauli sigma Matrix TODO */
	/* Pauli.fill(0.);//TODO */
	/* Pauli.topRightCorner(18,18) = -i*I;//TODO */
	/* Pauli.bottomLeftCorner(18,18) = i*I;//TODO */

	/* double spincurrent; */
	/* dddmat A, B, TOT; */
	/* dddmat T, Tdagg; */
	/* T.fill(0.);//TODO */
	/* T.topLeftCorner(18,18) = NM_T;//TODO */
	/* T.bottomRightCorner(18,18) = NM_T;//TODO */
	/* Tdagg = T.adjoint();//TODO */
	/* dddmat GR_T_dagg, GR_dagg_T_dagg; */
	/* GR_T_dagg = GR*Tdagg; */
	/* GR_dagg_T_dagg = GR_dagg*Tdagg; */
	/* dddmat tmp1, tmp2; */
	/* int N_diff = E_N - N_old; */
	/* for (int kk = 0; kk < N_old/2; kk++){ */
	/* 	GL_up = (OM - (NM - V*I) -NM_T_dagg*GL*NM_T).inverse(); */
	/* 	GL_dn = (OM - (NM - V*I) -NM_T_dagg*GL*NM_T).inverse(); */
	/* } */
	/* GL.topLeftCorner(18,18) = GL_up; */
	/* GL.bottomRightCorner(18,18) = GL_dn; */
	/* A = (Ibig-GR_T_dagg*GL*T).inverse(); */
	/* B = (Ibig-GR_dagg_T_dagg*GL.adjoint()*T).inverse(); */
	/* tmp1 = B*GR_dagg_T_dagg; */
	/* tmp2 = A*tmp1; */
	/* tmp1 = T*tmp2; */
	/* tmp2 = GL*tmp1; */
	/* TOT = (tmp2-A*B+0.5*(A+B))*Pauli; */
	/* spincurrent = (-1./(4.*M_PI))*real(TOT.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT))); */
	/* return spincurrent; */
/* } */

double int_theta_E(double E, void * params) {
	variables * send = (variables *) params;
	dcomp E_send;
	dcomp im = -1;
	im = sqrt(im);
	E_send = E + 1e-6*im;
	double result = 0.;
	double integrate;
	double theta;
	vector<double> dummy;
	vector<double> E_list;
	E_list = *send->E_list;
	int cond = 0;
	for (int kk = 0; kk < E_list.size(); kk++){
		if (abs(E-E_list[kk]) < 1e-10){
			cond = 1;
			break;
		}
	}

	const int n = 10;
	for (int k=0; k<n+1; k++) {
		theta = k*M_PI/n;
		if ((k==0) && (cond == 0)){
			integrate = f(theta, E_send, send, 0, dummy);
			E_list.emplace_back(E);
		}
		else {
			integrate = mem();
			/* integrate = f(theta, E_send, send, 0, dummy); */
			/* if (k==0) */
			/* 	cout<<"test"<<endl; */
		}
		if ((k==0)||(k==n))
			result += M_PI*(0.5/n)*integrate;
		else 
			result += (M_PI/n)*integrate;
	}	
	*send->E_list = E_list;
	return result;
}

vector<double> int_theta(const dcomp E, variables * send, const int myswitch) {
	vector<double> result;
	vector<double> integrate;
	int N = send->N;
	result.reserve(N);
	integrate.reserve(N);
	for (int i = 0; i < N; i++){
		result[i] = 0.;
		integrate[i] = 0.;
	}
	double theta;
	double dummy;

	const int n = 10;
	/* const int n = 1; */
	for (int k=0; k<n+1; k++) {
		theta = k*M_PI/n;
		integrate.clear();
		dummy = f(theta, E, send, myswitch, integrate);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += M_PI*(0.5/n)*integrate[i];
			else 
				result[i] += (M_PI/n)*integrate[i];
		}
	}	
	return result;
}

vector<double> int_energy(variables * send) {
	vector<double> result;
	int N = send->N;
	double Ef = send->Ef;

	result.reserve(N);

	double end = Ef + 0.1;
	double start = Ef - send->V - 0.1;
	double dresult;

	vector<int> N_old;
	vMM GL_U;
	vMM GL_D;
	vMMM GR;
	vector<double> E_list;
	send->GL_U = &GL_U;
	send->GL_D = &GL_D;
	send->GR = &GR;
	send->N_old = &N_old;
	send->E_list = &E_list;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &int_theta_E;
	double tol = 1e2;
	int max_it = 10;
	double error;
	/* int key = 1; */
	int status;
	for (int i = 0; i < N; i++){
		send->E_N = i;
		F.params = send;
		status = gsl_integration_qags(&F, start, end, 0, tol, max_it, w, &dresult, &error);
		/* vector<double> test = *send->E_list; */
		/* for (int kk = 0; kk < test.size(); kk++) */
		/* 	cout<<test[kk]<<endl; */
		cout<<status<<endl;
		result.emplace_back(dresult);
	}
	gsl_integration_workspace_free (w);

	/* double E; */
	/* const int n = 2000; */
	/* double factor = (end - start)/(n*1.); */
	/* for (int i = 0; i < N; i++){ */
	/* 	double tmp_result = 0; */
	/* 	for (int k=0; k<n+1; k++) { */
	/* 		E = start + k*(end-start)/(n*1.); */
	/* 		send->E_N = i; */
	/* 		dresult = pass(E, send); */
	/* 		if ((k==0)||(k==n)) */
	/* 			tmp_result += 0.5*factor*dresult; */
	/* 		else */ 
	/* 			tmp_result += factor*dresult; */
	/* 	} */	
	/* 	result.emplace_back(tmp_result); */
	/* } */

	return result;
}

vector<double> switching(variables * send) {
	vector<double> result1, result2, integrate;
	int N = send->N;
	double V = send->V;
	result1.reserve(N);
	if (abs(V) > 1e-4)
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
		integrate = int_theta(E, send, 1);
		if (abs(V) < 1e-5){
			for (int l = 0; l < N; l++)
				result2[l] += 2.*kT*integrate[l]; 
		}
		else {
			for (int l = 0; l < N; l++)
				result2[l] += kT*integrate[l]; 
			E = send->Ef - V + (2.*j + 1.)*kT*M_PI*i;
			integrate = int_theta(E, send, 1);
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
	//number of atomic planes
	// plot output of spincurrent against energy
	const double Ef = 0.57553;

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
	//This section defines the primitive lattice vectors
	Vector3d lat1, lat2, lat3;
	lat1 = bas2;
	lat2 << 0, 0.5, 0.5;
	lat3 << 0.5, 0, 0.5;
	//Distance information for n.n and n.n.n
	double nn_dist, nnn_dist;
	nn_dist = M_SQRT2/2.;
	nnn_dist = 1.;
	//This section generates the Hamiltonians from SK parameters and NN positions
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
	//magic 19 above is num onsite + num nn + num nnn = 1 + 12 + 6
	Matrix<dcomp, 9, 9> tmp_mat;
	double distance;
	for (int i1 = -1; i1 < 2; i1++){
		for (int i2 = -1; i2 < 2; i2++){
			for (int i3 = -1; i3 < 2; i3++){
				tmp_vec = i1*lat1 + i2*lat2 + i3*lat3;
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					pos.emplace_back(tmp_vec);
					cobalt_up.emplace_back(Co_u);
					cobalt_dn.emplace_back(Co_d);
					copper.emplace_back(Cu);
					cob_cop_up.emplace_back(Cu); // In theory this will not be used
					cob_cop_dn.emplace_back(Cu); // In theory this will not be used
				}
				else if (distance < nn_dist + 1e-3){
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
				else if (distance < nnn_dist + 1e-3){
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
			}
		}
	}

	// number of spacer layers
	int N = 10;
	// set bias
	/* double V = 0.0; */
	double V = 0.3;

	const double k = 8.617e-5/13.6058;//boltzmann constant (in Ryds)
	const double T = 300;//set the temperature
	double kT = k*T;
	//set up the variables to send
	variables send;
	send.kT = kT;
	send.N = N;
	send.Ef = Ef;
	send.x = 2.532374; //for now!
	send.z = 2.532374; //for now!
	send.t = &t;
	send.pos = &pos;
	send.basis = &basis;
	send.copper = &copper;
	send.cobalt_up = &cobalt_up;
	send.cobalt_dn = &cobalt_dn;
	send.cob_cop_up = &cob_cop_up;
	send.cob_cop_dn = &cob_cop_dn;
	send.V = V;

	vector<double> answer;
	answer.reserve(N);
	answer = switching(&send);

	string Mydata;
	ofstream Myfile;	
	if (abs(V) < 1e-4)
		Mydata = "Keldysh_V0.txt";
	else
		Mydata = "Keldysh_V.txt";

	/* vector<double> result; */
	/* vector<double> integrate; */
	/* result.reserve(N); */
	/* integrate.reserve(N); */
	/* for (int i = 0; i < N; i++) */
	/* 	result[i] = 0.; */
	/* int n = 350;//set the number of k-points along x axis */
	/* /1* int n = 2;//set the number of k-points along x axis *1/ */
	/* int counter = 0; */
	/* double factor = 2./(n*n); */
	/* int sumk = n*(n + 1)/2; */
	/* int p, q; */
	/* int product1; */
	/* int myid, numprocs; */
	/* time_t timer; */
	/* int k, l, i; */
	/* int start_time = time(&timer); */
	/* MPI_Init(NULL,NULL); */
	/* MPI_Comm_size(MPI_COMM_WORLD, &numprocs); */
	/* MPI_Comm_rank(MPI_COMM_WORLD, &myid); */
	/* for (k = 2*myid + 1; k<2*sumk; k+=2*numprocs){ */
	/* 	for (l = 0; l < n; l++){ */
	/* 		product1 = (2*n - l)*(l + 1); */
	/* 		if ( k < product1){ */
	/* 			p = 2*(k/(2*n - l)) + 1; */
	/* 			q = (l*(l + 1) + k)%(2*n); */
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
