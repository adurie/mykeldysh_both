#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/src/Core/util/MKL_support.h>
#include <nag.h>
#include <nagd01.h>
#include <nag_stdlib.h>
#include <vector>
#include <unordered_map>
#include "AuMgOFe.h"
#include <ctime>
/* #include "/home/alex/INTEL/impi/2019.1.144/intel64/include/mpi.h" */
#include <iomanip>
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_USE_MKL_ALL
//this now has modified onsite Fe pots when next to MgO, and MgO - atom gmeans
//TODO very important, the above note means that a different system will require editing

using namespace Eigen;
using namespace std;
typedef complex<double> dcomp;
typedef Matrix<complex<double>, 9, 9> M9;
typedef vector<Vector3d, aligned_allocator<Vector3d>> vec3;
typedef vector<Matrix<dcomp, 9, 9>, aligned_allocator<Matrix<dcomp, 9, 9>>> vM;
typedef vector<vec3> vvec3;
typedef vector<MatrixXcd, aligned_allocator<MatrixXcd>> vMXd;

typedef struct
	{
		int numbas;
		int numlay;
		int N;
		int it;
		double Ef;
		double x;
		double z;
		double V;
		double kT;
		Vector3d *b1;
		Vector3d *b2;
		unordered_map<double, unordered_map<double, vector<double>>> *all_the_data;
		vector<bool> *isMag;
		vector<vector<string>> *attype;
		vector<vector<string>> *oddattype;
		vector<int> *thick;
		unordered_map<string, Vector3d> *lat_oop;
		vec3 *odd_oop;
		unordered_map<string, vector<vvec3>> *pos;
		unordered_map<string, vector<vvec3>> *odd_pos;
		vvec3 *basis;
		vvec3 *odd_basis;
		unordered_map<string, vector<vector<vM>>> *atom_up;
		unordered_map<string, vector<vector<vM>>> *atom_dn;
		unordered_map<string, vector<vector<vM>>> *odd_atom_up;
		unordered_map<string, vector<vector<vM>>> *odd_atom_dn;
		vMXd *HUu;
		vMXd *HUd;
		vMXd *HTu;
		vMXd *HTd;
		vMXd *odd_Hu;
		vMXd *odd_Hd;
	}
variables;

static double NAG_CALL phi1(double y, Nag_Comm *comm)
{
	return 0.0;
}

static double NAG_CALL phi2(double y, Nag_Comm *comm)
{
	return y;
}

void recip(const Vector3d &a1, const Vector3d &a2, Vector3d &b1, Vector3d &b2){
	//construct the third (normalised) lattice vector
      Vector3d d;
      d = a1.cross(a2);
      d.normalize();
//     -----------------------------------------------------------------
//     now determine the reciprocal lattice vectors
      Vector3d an1, an2;
      an1=a1;
      an2=a2;
      b1 = an2.cross(d);
      b2 = an1.cross(d);
      double a1b1, a2b2;
      a1b1 = a1.dot(b1);
      a2b2 = a2.dot(b2);
      b1 = 2*M_PI*b1/a1b1;
      b2 = 2*M_PI*b2/a2b2;
//     now determine whether b2 or -b2 has the smallest angle with b1
      double b1b1, b2b2, b1b2;
      b1b1 = b1.dot(b1);
      b2b2 = b2.dot(b2);
      b1b2 = b1.dot(b2);
      if (b1b2 < 0)
        b2=-b2;
      return;
}

//calculates g at the interfaces
MatrixXcd gs(const MatrixXcd &OM, const MatrixXcd &T, const int numbas)
{
	int n = numbas*9;
	int n2 = n*2;
	MatrixXcd zero = MatrixXcd::Zero(n,n);
	MatrixXcd Tinv;
	Tinv = T.inverse();
	MatrixXcd X(n2,n2), O(n2,n2);
	X.topLeftCorner(n, n) = zero;
	X.topRightCorner(n, n) = Tinv;
	X.bottomLeftCorner(n, n) = -T.adjoint();
	X.bottomRightCorner(n, n) = OM*Tinv;
	ComplexEigenSolver<MatrixXcd> ces;
	ces.compute(X);
	O = ces.eigenvectors();
	MatrixXcd b = O.topRightCorner(n, n);
	MatrixXcd d = O.bottomRightCorner(n, n);
	MatrixXcd GR;
	GR = b*d.inverse();
	return GR;
}

dcomp fermi(const dcomp arg, const double Ef, const double kT){
	return 1./(1.+exp((arg-Ef)/kT));
}

M9 InPlaneH(const vec3 &pos, const Vector3d &basis, const vM &U, const double x, const double z){
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
			/* cout<<tmp_vec.transpose()<<endl; */
			result = result + U[k]*exp(i*tmp_vec.dot(K));
		}
	}
	return result;
}

 vector<double> f_vec(const double theta, const dcomp E, variables * send, const MatrixXcd &GR_UP_e, const MatrixXcd &GR_DN_e,  
 		const MatrixXcd &GR_UP_o, const MatrixXcd &GR_DN_o, const MatrixXcd &Gl_up, const MatrixXcd &Gl_dn) { 
 	dcomp i = -1; 
 	i = sqrt(i); 
 	double V = send->V; 
 	double kT = send->kT; 
 	double Ef = send->Ef; 
 	int N = send->N; 
	int numbas = send->numbas;
	int n = numbas*9;
	int n2 = n*2;
	vMXd HUu = *send->HUu;
	vMXd HTu = *send->HTu;
	vMXd HTd = *send->HTd;
 	MatrixXcd NM = HUu[3]; 
 	MatrixXcd NM_T = HTu[6]; 

 	MatrixXcd GR_up_even, GR_dn_even, GR_up_odd, GR_dn_odd; 
 	GR_up_even = GR_UP_e; 
 	GR_dn_even = GR_DN_e; 
 	GR_up_odd = GR_UP_o; 
 	GR_dn_odd = GR_DN_o; 
 	MatrixXcd GL_up, GL_dn; 
	GL_up = Gl_up;
	GL_dn = Gl_dn;
 	MatrixXcd GL(n2,n2); 
 	MatrixXcd NM_T_dagg; 
 	NM_T_dagg = NM_T.adjoint(); 
 	MatrixXcd I = MatrixXcd::Identity(n,n); 
 	MatrixXcd S(n2,n2); 
 	MatrixXcd S11, S12; 
 	S11 = cos(theta/2.)*I; 
 	S12 = sin(theta/2.)*I; 
 	S.topLeftCorner(n,n) = S11; 
 	S.topRightCorner(n,n) = S12; 
 	S.bottomLeftCorner(n,n) = -S12; 
 	S.bottomRightCorner(n,n) = S11; 

 	MatrixXcd OM = E*I; 
 	MatrixXcd GR_even(n2,n2), GR_odd(n2,n2), GR_dagg_even(n2,n2), GR_dagg_odd(n2,n2); 
 	GR_even.fill(0.); 
 	GR_odd.fill(0.); 
 	GR_even.topLeftCorner(n,n) = GR_up_even; 
 	GR_even.bottomRightCorner(n,n) = GR_dn_even; 
 	GR_odd.topLeftCorner(n,n) = GR_up_odd; 
 	GR_odd.bottomRightCorner(n,n) = GR_dn_odd; 

 	GL.fill(0.); 
 	GR_even = S.inverse()*GR_even*S; 
 	GR_odd = S.inverse()*GR_odd*S; 

 	MatrixXcd Ibig = MatrixXcd::Identity(n2,n2); 
 	MatrixXcd OMbig = E*Ibig; 
 	MatrixXcd NMbig(n2,n2); 
 	NMbig.fill(0.); 
 	NMbig.topLeftCorner(n,n) = NM; 
 	NMbig.bottomRightCorner(n,n) = NM; 
 	GR_dagg_even = GR_even.adjoint(); 
 	GR_dagg_odd = GR_odd.adjoint(); 

	 MatrixXcd Pauli(n2,n2);//This is the y Pauli sigma Matrix 
	 Pauli.fill(0.); 
	 Pauli.topRightCorner(n,n) = -i*I; 
	 Pauli.bottomLeftCorner(n,n) = i*I; 
	 /* Pauli.topRightCorner(n,n) = I; */ 
	 /* Pauli.bottomLeftCorner(n,n) = I; */ 

 	/* MatrixXcd Pauli = Ibig; */ 

 	double spincurrent_even, spincurrent_odd; 
 	MatrixXcd A_even, A_odd, B_even, B_odd, TOT_even, TOT_odd; 
 	MatrixXcd T(n2,n2), Tdagg(n2,n2); 
 	T.fill(0.); 
 	T.topLeftCorner(n,n) = NM_T; 
 	T.bottomRightCorner(n,n) = NM_T; 
 	Tdagg = T.adjoint(); 
 	MatrixXcd GR_T_dagg_even, GR_dagg_T_dagg_even; 
 	MatrixXcd GR_T_dagg_odd, GR_dagg_T_dagg_odd; 
 	GR_T_dagg_even = GR_even*Tdagg; 
 	GR_T_dagg_odd = GR_odd*Tdagg; 
 	GR_dagg_T_dagg_even = GR_dagg_even*Tdagg; 
 	GR_dagg_T_dagg_odd = GR_dagg_odd*Tdagg; 
 	MatrixXcd tmp1, tmp2; 
 	vector<double> result; 
 	result.reserve(N); 
 //adlayer layer 2 from layer 1 to spacer thickness, N 
 	for (int it=0; it < N/2; ++it){ 
 		GL.topLeftCorner(n,n) = GL_up; 
 		GL.bottomRightCorner(n,n) = GL_dn; 
 		A_even = (Ibig-GR_T_dagg_even*GL*T).inverse(); 
 		B_even = (Ibig-GR_dagg_T_dagg_even*GL.adjoint()*T).inverse(); 
 		A_odd = (Ibig-GR_T_dagg_odd*GL*T).inverse(); 
 		B_odd = (Ibig-GR_dagg_T_dagg_odd*GL.adjoint()*T).inverse(); 
 		TOT_even = (B_even.adjoint()-A_even)*Pauli; 
 		TOT_odd = (B_odd.adjoint()-A_odd)*Pauli; 
 		spincurrent_even = .25*imag(TOT_even.trace()); 
 		spincurrent_odd = .25*imag(TOT_odd.trace()); 
 		result.emplace_back(spincurrent_even); 
 		result.emplace_back(spincurrent_odd); 
 		GL_up = (OM - (NM - V*I) -NM_T_dagg*GL_up*NM_T).inverse(); 
 		GL_dn = (OM - (NM - V*I) -NM_T_dagg*GL_dn*NM_T).inverse(); 
 	} 
 	return result; 
 } 

 void f(double * spincurrent, const double theta, const dcomp E, int k, Integer * needi, variables * send, const MatrixXcd &GR_UP_even,  
 		const MatrixXcd &GR_DN_even, const MatrixXcd &GR_UP_odd, const MatrixXcd &GR_DN_odd, const MatrixXcd &Gl_up, const MatrixXcd &Gl_dn) { 
 // ...NM|ins|FM(0)|NM(n)|FM(theta)... 
 	dcomp i = -1; 
 	i = sqrt(i); 
	int numbas = send->numbas;
	int n = numbas*9;
	int n2 = 2*n;
 	double V = send->V; 
 	double kT = send->kT; 
 	double Ef = send->Ef; 
 	int N = send->N; 
	vMXd HUu = *send->HUu;
	vMXd HTu = *send->HTu;
	vMXd HTd = *send->HTd;
 	MatrixXcd NM = HUu[3]; 
 	MatrixXcd NM_T = HTu[6]; 

 	MatrixXcd GR_up_even, GR_dn_even, GR_up_odd, GR_dn_odd; 
 	GR_up_even = GR_UP_even; 
 	GR_dn_even = GR_DN_even; 
 	GR_up_odd = GR_UP_odd; 
 	GR_dn_odd = GR_DN_odd; 
	MatrixXcd GR_odd(n2,n2), GR_even(n2,n2);
 	MatrixXcd GL_up, GL_dn; 
	GL_up = Gl_up;
	GL_dn = Gl_dn;
 	MatrixXcd NM_T_dagg; 
 	NM_T_dagg = NM_T.adjoint(); 
 	MatrixXcd I = MatrixXcd::Identity(n,n); 
 	MatrixXcd S(n2,n2); 
 	MatrixXcd S11, S12; 
 	S11 = cos(theta/2.)*I; 
 	S12 = sin(theta/2.)*I; 
 	S.topLeftCorner(n,n) = S11; 
 	S.topRightCorner(n,n) = S12; 
 	S.bottomLeftCorner(n,n) = -S12; 
 	S.bottomRightCorner(n,n) = S11; 

 	MatrixXcd OM = E*I; 
 	GR_even.fill(0.); 
 	GR_even.topLeftCorner(n,n) = GR_up_even; 
 	GR_even.bottomRightCorner(n,n) = GR_dn_even; 
 	GR_odd.fill(0.); 
 	GR_odd.topLeftCorner(n,n) = GR_up_odd; 
 	GR_odd.bottomRightCorner(n,n) = GR_dn_odd; 

 	MatrixXcd GL(n2,n2), GR_dagg_even(n2,n2), GR_dagg_odd(n2,n2); 
 	GL.fill(0.); 
 	GR_even = S.inverse()*GR_even*S; 
 	GR_odd = S.inverse()*GR_odd*S; 

 	MatrixXcd Ibig = MatrixXcd::Identity(n2,n2); 
 	MatrixXcd OMbig = E*Ibig; 
 	MatrixXcd NMbig(n2,n2); 
 	NMbig.fill(0.); 
 	NMbig.topLeftCorner(n,n) = NM; 
 	NMbig.bottomRightCorner(n,n) = NM; 
 	GR_dagg_even = GR_even.adjoint(); 
 	GR_dagg_odd = GR_odd.adjoint(); 

 	MatrixXcd Pauli(n2,n2);//This is the y Pauli sigma Matrix 
 	Pauli.fill(0.); 
 	Pauli.topRightCorner(n,n) = -i*I; 
 	Pauli.bottomLeftCorner(n,n) = i*I; 

 	MatrixXcd A, B, TOT; 
 	MatrixXcd T(n2,n2), Tdagg(n2,n2); 
 	T.fill(0.); 
 	T.topLeftCorner(n,n) = NM_T; 
 	T.bottomRightCorner(n,n) = NM_T; 
 	Tdagg = T.adjoint(); 
 	MatrixXcd GR_T_dagg_even, GR_dagg_T_dagg_even; 
 	MatrixXcd GR_T_dagg_odd, GR_dagg_T_dagg_odd; 
 	GR_T_dagg_even = GR_even*Tdagg; 
 	GR_T_dagg_odd = GR_odd*Tdagg; 
 	GR_dagg_T_dagg_even = GR_dagg_even*Tdagg; 
 	GR_dagg_T_dagg_odd = GR_dagg_odd*Tdagg; 
 	MatrixXcd tmp1, tmp2; 
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
 			GL.topLeftCorner(n,n) = GL_up; 
 			GL.bottomRightCorner(n,n) = GL_dn; 
 			A = (Ibig-GR_T_dagg_even*GL*T).inverse(); 
 			B = (Ibig-GR_dagg_T_dagg_even*GL.adjoint()*T).inverse(); 
 			tmp1 = B*GR_dagg_T_dagg_even; 
 			tmp2 = A*tmp1; 
 			tmp1 = T*tmp2; 
 			tmp2 = GL*tmp1; 
 			TOT = (tmp2-A*B+0.5*(A+B))*Pauli; 
 			spincurrent[2*it] = (1./(4.*M_PI))*real(TOT.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT))); 
 		} 

 		if (needi[2*it + 1] == 1){ 
 			GL.topLeftCorner(n,n) = GL_up; 
 			GL.bottomRightCorner(n,n) = GL_dn; 
 			A= (Ibig-GR_T_dagg_odd*GL*T).inverse(); 
 			B= (Ibig-GR_dagg_T_dagg_odd*GL.adjoint()*T).inverse(); 
 			tmp1 = B*GR_dagg_T_dagg_odd; 
 			tmp2 = A*tmp1; 
 			tmp1 = T*tmp2; 
 			tmp2 = GL*tmp1; 
 			TOT = (tmp2-A*B+0.5*(A+B))*Pauli; 
 			spincurrent[2*it + 1] = (1./(4.*M_PI))*real(TOT.trace()*(fermi(E,Ef,kT)-fermi(E,Ef-V,kT))); 
 		} 

 		if (it != kl/2 - 1){//saves wasting unused results 
 			GL_up = (OM - (NM - V*I) -NM_T_dagg*GL_up*NM_T).inverse(); 
 			GL_dn = (OM - (NM - V*I) -NM_T_dagg*GL_dn*NM_T).inverse(); 
 		} 
 	} 
 } 

 vector<double> int_theta(const dcomp E, variables * send) { 
 	vector<double> result; 
 	vector<double> integrate; 
	int numbas = send->numbas;
	int n = numbas*9;
	vMXd HUu = *send->HUu;
	vMXd HUd = *send->HUd;
	vMXd HTu = *send->HTu;
	vMXd HTd = *send->HTd;
	vMXd odd_Hu = *send->odd_Hu;
	vMXd odd_Hd = *send->odd_Hd;
 	MatrixXcd NM = HUu[0]; 
 	MatrixXcd NM_T = HTu[0]; 
 	MatrixXcd FM_up = HUu[2]; 
 	MatrixXcd FM_dn = HUd[2]; 
 	MatrixXcd FM_up_T = HTu[4]; 
 	MatrixXcd FM_dn_T = HTd[4]; 
 	MatrixXcd FM_NM_up_T = HTu[5]; 
 	MatrixXcd FM_NM_dn_T = HTd[5]; 
 	MatrixXcd odd_l1_up = odd_Hu[0]; 
 	MatrixXcd odd_l1_dn = odd_Hd[0]; 
 	MatrixXcd odd_l1_up_T1 = odd_Hu[1]; 
 	MatrixXcd odd_l1_dn_T1 = odd_Hd[1]; 
 	MatrixXcd odd_l1_up_T2 = odd_Hu[2]; 
 	MatrixXcd odd_l1_dn_T2 = odd_Hd[2]; 
 	MatrixXcd ins = HUu[1]; 
 	MatrixXcd ins_T = HTu[2]; 
 	MatrixXcd ins_FM_up_T = HTu[3]; 
 	MatrixXcd ins_FM_dn_T = HTd[3]; 
 	MatrixXcd ins_NM_T = HTu[1]; 
 	MatrixXcd NM_FM_up_T = HTu[7]; 
 	MatrixXcd NM_FM_dn_T = HTd[7]; 
 	MatrixXcd NM_FM_up_T_dagg; 
 	MatrixXcd NM_FM_dn_T_dagg; 
 	NM_FM_up_T_dagg = NM_FM_up_T.adjoint(); 
 	NM_FM_dn_T_dagg = NM_FM_dn_T.adjoint(); 
 	double V = send->V; 

 	MatrixXcd I = MatrixXcd::Identity(n,n); 
 	MatrixXcd OMup=E*I-(FM_up - V*I); 
 	MatrixXcd OMdn=E*I-(FM_dn - V*I); 
 	MatrixXcd OM = E*I; 

 	MatrixXcd FM_up_T_dagg = FM_up_T.adjoint(); 
 	MatrixXcd FM_dn_T_dagg = FM_dn_T.adjoint(); 
 	MatrixXcd GR_up_even = gs(OMup, FM_up_T_dagg, numbas); 
 	MatrixXcd GR_dn_even = gs(OMdn, FM_dn_T_dagg, numbas); 
 	MatrixXcd FM_NM_up_T_dagg = FM_NM_up_T.adjoint(); 
 	MatrixXcd FM_NM_dn_T_dagg = FM_NM_dn_T.adjoint(); 
 	MatrixXcd NM_T_dagg = NM_T.adjoint(); 
	
 	/* //this for trilayer */ 
 	/* MatrixXcd GL_up = gs(OMup, FM_T, numbas); */ 
 	/* MatrixXcd GL_dn = gs(OMdn, FM_T, numbas); */ 

 	//this below block for 5 layer 
 	MatrixXcd GL_up = gs(OM - NM, NM_T, numbas); 
 	MatrixXcd GL_dn = GL_up; 
 	//this below block for 5 layer 
 //lim is thickness of layer 2 
 	MatrixXcd ins_T_dagg = ins_T.adjoint(); 
 	MatrixXcd ins_NM_T_dagg = ins_NM_T.adjoint(); 
	vector<int> thick = *send->thick;
	int lim = thick[0];
	int lim2 = thick[1];
 //build thickness of layer 2 to lim layers 
 	for (int it=0; it < lim; ++it){//TODO the diagonal elements need to be shifted by the same amount in halites 
 		ins = ins - I*(V*it/(lim*1.));//TODO changed this so that full bias isn't on last layer so if lim = 1 then shift = 0, 1/2 (so 1 falls on next layer) 
 		if (it == 0){ 
 			GL_up = (OM - ins -ins_T_dagg*GL_up*ins_T).inverse(); 
 			GL_dn = (OM - ins -ins_T_dagg*GL_dn*ins_T).inverse(); 
 		} 
 		else { 
 			GL_up = (OM - ins -ins_NM_T_dagg*GL_up*ins_NM_T).inverse(); 
 			GL_dn = (OM - ins -ins_NM_T_dagg*GL_dn*ins_NM_T).inverse(); 
 		} 
 	} 
 //lim2 is thickness of layer 3 
 //build thickness of layer 3 to lim2 layers 
 //add 10 bilayers i.e. 20 layers of FM 
	//TODO very important -- this will need to be altered or even removed for a different system
	int csp5 = -0.310269;//Andrey's technique, modify interface pots of Fe against MgO to get correct magnetic moment
	int cd5 = 4.9467782e-3;
	MatrixXcd tmp_up(n,n), tmp_dn(n,n);
	tmp_up.fill(0.); tmp_dn.fill(0.);
	int mat_count;
	for (int kl = 0; kl < numbas; kl++){
		mat_count = 9*kl;
		tmp_up(mat_count, mat_count) = csp5;//s
		tmp_up(1+mat_count, 1+mat_count) = csp5;//p
		tmp_up(2+mat_count, 2+mat_count) = csp5;//p
		tmp_up(3+mat_count, 3+mat_count) = csp5;//p
		tmp_up(4+mat_count, 4+mat_count) = -cd5;//dt
		tmp_up(5+mat_count, 5+mat_count) = -cd5;//dt
		tmp_up(6+mat_count, 6+mat_count) = -cd5;//dt
		tmp_up(7+mat_count, 7+mat_count) = -cd5;//de
		tmp_up(8+mat_count, 8+mat_count) = -cd5;//de

		tmp_dn(mat_count, mat_count) = csp5;//s
		tmp_dn(1+mat_count, 1+mat_count) = csp5;//p
		tmp_dn(2+mat_count, 2+mat_count) = csp5;//p
		tmp_dn(3+mat_count, 3+mat_count) = csp5;//p
		tmp_dn(4+mat_count, 4+mat_count) = cd5;//dt
		tmp_dn(5+mat_count, 5+mat_count) = cd5;//dt
		tmp_dn(6+mat_count, 6+mat_count) = cd5;//dt
		tmp_dn(7+mat_count, 7+mat_count) = cd5;//de
		tmp_dn(8+mat_count, 8+mat_count) = cd5;//de
	}
 	GL_up = (OM - (FM_up + tmp_up - V*I) -ins_FM_up_T.adjoint()*GL_up*ins_FM_up_T).inverse(); 
 	GL_dn = (OM - (FM_dn + tmp_dn - V*I) -ins_FM_dn_T.adjoint()*GL_dn*ins_FM_dn_T).inverse(); 
 	for (int it=0; it < lim2 - 1; ++it){ 
 		GL_up = (OM - (FM_up - V*I) -FM_up_T_dagg*GL_up*FM_up_T).inverse(); 
 		GL_dn = (OM - (FM_dn - V*I) -FM_dn_T_dagg*GL_dn*FM_dn_T).inverse(); 
 	} 

 	MatrixXcd odd_l1_up_T1_dagg = odd_l1_up_T1.adjoint(); 
 	MatrixXcd odd_l1_dn_T1_dagg = odd_l1_dn_T1.adjoint(); 
 	MatrixXcd odd_l1_up_T2_dagg = odd_l1_up_T2.adjoint(); 
 	MatrixXcd odd_l1_dn_T2_dagg = odd_l1_dn_T2.adjoint(); 
 	MatrixXcd GR_up_odd = GR_up_even; 
 	MatrixXcd GR_dn_odd = GR_dn_even; 
 	//adlayer one bilayer onto RHS & LHS G to ensure gmean is correct 
 	//this means 2 layers are on before we begin! 
 	GL_up = (OM - (NM - V*I) -FM_NM_up_T_dagg*GL_up*FM_NM_up_T).inverse(); 
 	GL_dn = (OM - (NM - V*I) -FM_NM_dn_T_dagg*GL_dn*FM_NM_dn_T).inverse(); 
 	GR_up_even = (OM - (NM - V*I)-NM_FM_up_T*GR_up_even*NM_FM_up_T_dagg).inverse();
 	GR_dn_even = (OM - (NM - V*I)-NM_FM_dn_T*GR_dn_even*NM_FM_dn_T_dagg).inverse(); 
 	//adlayer one bilayer of CoCu onto LHS G for odd layers, then adlayer a  
 	//further bilayer of Cu to ensure gmean is correct. This means 3 layers are on before we begin! 
 	GR_up_odd = (OM - (odd_l1_up - V*I) -odd_l1_up_T2*GR_up_odd*odd_l1_up_T2_dagg).inverse(); 
 	GR_dn_odd = (OM - (odd_l1_dn - V*I) -odd_l1_dn_T2*GR_dn_odd*odd_l1_dn_T2_dagg).inverse(); 
 	GR_up_odd = (OM - (NM - V*I) -odd_l1_up_T1*GR_up_odd*odd_l1_up_T1_dagg).inverse(); 
 	GR_dn_odd = (OM - (NM - V*I) -odd_l1_dn_T1*GR_dn_odd*odd_l1_dn_T1_dagg).inverse(); 
 	//adlayer one bilayer onto LHS G_even to ensure gmean is correct 
 	//this means 2 layers are on before we begin! 
 	int N = send->N; 
 	result.reserve(N); 
 	integrate.reserve(N); 
 	for (int i = 0; i < N; i++){ 
 		result[i] = 0.; 
 		integrate[i] = 0.; 
 	} 
 	double theta; 

 	const int nn = 10; 
 	/* const int nn = 1; */ 
 	for (int k=0; k<nn+1; k++) { 
 		theta = k*M_PI/nn; 
 		integrate = f_vec(theta, E, send, GR_up_even, GR_dn_even, GR_up_odd, GR_dn_odd, GL_up, GL_dn); 
 		for (int i = 0; i < N; i++){ 
 			if ((k==0)||(k==nn)) 
 				result[i] += M_PI*(0.5/nn)*integrate[i]; 
 			else  
 				result[i] += (M_PI/nn)*integrate[i]; 
 		} 
 	} 	
 	return result; 
 } 

 void int_theta_E(const dcomp E, int k, double * fm, Integer * needi, variables * send) { 
 	double theta; 
 	int N = send->N; 
 	double result[N]; 
	int numbas = send->numbas;
	int n = numbas*9;
	vMXd HUu = *send->HUu;
	vMXd HUd = *send->HUd;
	vMXd HTu = *send->HTu;
	vMXd HTd = *send->HTd;
	vMXd odd_Hu = *send->odd_Hu;
	vMXd odd_Hd = *send->odd_Hd;
 	MatrixXcd NM = HUu[0]; 
 	MatrixXcd NM_T = HTu[0]; 
 	MatrixXcd FM_up = HUu[2]; 
 	MatrixXcd FM_dn = HUd[2]; 
 	MatrixXcd FM_up_T = HTu[4]; 
 	MatrixXcd FM_dn_T = HTd[4]; 
 	MatrixXcd FM_NM_up_T = HTu[5]; 
 	MatrixXcd FM_NM_dn_T = HTd[5]; 
 	MatrixXcd odd_l1_up = odd_Hu[0]; 
 	MatrixXcd odd_l1_dn = odd_Hd[0]; 
 	MatrixXcd odd_l1_up_T1 = odd_Hu[1]; 
 	MatrixXcd odd_l1_dn_T1 = odd_Hd[1]; 
 	MatrixXcd odd_l1_up_T2 = odd_Hu[2]; 
 	MatrixXcd odd_l1_dn_T2 = odd_Hd[2]; 
 	MatrixXcd ins = HUu[1]; 
 	MatrixXcd ins_T = HTu[2]; 
 	MatrixXcd ins_FM_up_T = HTu[3]; 
 	MatrixXcd ins_FM_dn_T = HTd[3]; 
 	MatrixXcd ins_NM_T = HTu[1]; 
 	MatrixXcd NM_FM_up_T = HTu[7]; 
 	MatrixXcd NM_FM_dn_T = HTd[7]; 
 	MatrixXcd NM_FM_up_T_dagg; 
 	MatrixXcd NM_FM_dn_T_dagg; 
 	NM_FM_up_T_dagg = NM_FM_up_T.adjoint(); 
 	NM_FM_dn_T_dagg = NM_FM_dn_T.adjoint(); 
 	double V = send->V; 

 	MatrixXcd I = MatrixXcd::Identity(n,n); 
 	MatrixXcd OMup=E*I-(FM_up - V*I); 
 	MatrixXcd OMdn=E*I-(FM_dn - V*I); 
 	MatrixXcd OM = E*I; 

 	MatrixXcd FM_up_T_dagg = FM_up_T.adjoint(); 
 	MatrixXcd FM_dn_T_dagg = FM_dn_T.adjoint(); 
 	MatrixXcd GR_up_even = gs(OMup, FM_up_T_dagg, numbas); 
 	MatrixXcd GR_dn_even = gs(OMdn, FM_dn_T_dagg, numbas); 
 	MatrixXcd FM_NM_up_T_dagg = FM_NM_up_T.adjoint(); 
 	MatrixXcd FM_NM_dn_T_dagg = FM_NM_dn_T.adjoint(); 
 	MatrixXcd NM_T_dagg = NM_T.adjoint(); 
	
 	/* //this for trilayer */ 
 	/* MatrixXcd GL_up = gs(OMup, FM_T, numbas); */ 
 	/* MatrixXcd GL_dn = gs(OMdn, FM_T, numbas); */ 

 	//this below block for 5 layer 
 	MatrixXcd GL_up = gs(OM - NM, NM_T, numbas); 
 	MatrixXcd GL_dn = GL_up; 
 	//this below block for 5 layer 
 //lim is thickness of layer 2 
	vector<int> thick = *send->thick;
	int lim = thick[0];
	int lim2 = thick[1];
 	MatrixXcd ins_T_dagg = ins_T.adjoint(); 
 	MatrixXcd ins_NM_T_dagg = ins_NM_T.adjoint(); 
 //build thickness of layer 2 to lim layers 
 	for (int it=0; it < lim; ++it){//TODO the diagonal elements need to be shifted by the same amount in halites 
 		ins = ins - I*(V*it/(lim*1.));//TODO changed this so that full bias isn't on last layer so if lim = 1 then shift = 0, 1/2 (so 1 falls on next layer) 
 		if (it == 0){ 
 			GL_up = (OM - ins -ins_T_dagg*GL_up*ins_T).inverse(); 
 			GL_dn = (OM - ins -ins_T_dagg*GL_dn*ins_T).inverse(); 
 		} 
 		else { 
 			GL_up = (OM - ins -ins_NM_T_dagg*GL_up*ins_NM_T).inverse(); 
 			GL_dn = (OM - ins -ins_NM_T_dagg*GL_dn*ins_NM_T).inverse(); 
 		} 
 	} 
 //lim2 is thickness of layer 3 
 //build thickness of layer 3 to lim2 layers 
	//TODO very important -- this will need to be altered or even removed for a different system
	int csp5 = -0.310269;//Andrey's technique, modify interface pots of Fe against MgO to get correct magnetic moment
	int cd5 = 4.9467782e-3;
	MatrixXcd tmp_up(n,n), tmp_dn(n,n);
	tmp_up.fill(0.); tmp_dn.fill(0.);
	int mat_count;
	for (int kl = 0; kl < numbas; kl++){
		mat_count = 9*kl;
		tmp_up(mat_count, mat_count) = csp5;//s
		tmp_up(1+mat_count, 1+mat_count) = csp5;//p
		tmp_up(2+mat_count, 2+mat_count) = csp5;//p
		tmp_up(3+mat_count, 3+mat_count) = csp5;//p
		tmp_up(4+mat_count, 4+mat_count) = -cd5;//dt
		tmp_up(5+mat_count, 5+mat_count) = -cd5;//dt
		tmp_up(6+mat_count, 6+mat_count) = -cd5;//dt
		tmp_up(7+mat_count, 7+mat_count) = -cd5;//de
		tmp_up(8+mat_count, 8+mat_count) = -cd5;//de

		tmp_dn(mat_count, mat_count) = csp5;//s
		tmp_dn(1+mat_count, 1+mat_count) = csp5;//p
		tmp_dn(2+mat_count, 2+mat_count) = csp5;//p
		tmp_dn(3+mat_count, 3+mat_count) = csp5;//p
		tmp_dn(4+mat_count, 4+mat_count) = cd5;//dt
		tmp_dn(5+mat_count, 5+mat_count) = cd5;//dt
		tmp_dn(6+mat_count, 6+mat_count) = cd5;//dt
		tmp_dn(7+mat_count, 7+mat_count) = cd5;//de
		tmp_dn(8+mat_count, 8+mat_count) = cd5;//de
	}
 	GL_up = (OM - (FM_up + tmp_up - V*I) -ins_FM_up_T.adjoint()*GL_up*ins_FM_up_T).inverse(); 
 	GL_dn = (OM - (FM_dn + tmp_dn - V*I) -ins_FM_dn_T.adjoint()*GL_dn*ins_FM_dn_T).inverse(); 
 	for (int it=0; it < lim2 - 1; ++it){ 
 		GL_up = (OM - (FM_up - V*I) -FM_up_T_dagg*GL_up*FM_up_T).inverse(); 
 		GL_dn = (OM - (FM_dn - V*I) -FM_dn_T_dagg*GL_dn*FM_dn_T).inverse(); 
 	} 

 	MatrixXcd odd_l1_up_T1_dagg = odd_l1_up_T1.adjoint(); 
 	MatrixXcd odd_l1_up_T2_dagg = odd_l1_up_T2.adjoint(); 
 	MatrixXcd odd_l1_dn_T1_dagg = odd_l1_dn_T1.adjoint(); 
 	MatrixXcd odd_l1_dn_T2_dagg = odd_l1_dn_T2.adjoint(); 
 	MatrixXcd GR_up_odd = GR_up_even; 
 	MatrixXcd GR_dn_odd = GR_dn_even; 
 	//adlayer one bilayer onto RHS & LHS G to ensure gmean is correct 
 	//this means 2 layers are on before we begin! 
 	GL_up = (OM - (NM - V*I) -FM_NM_up_T_dagg*GL_up*FM_NM_up_T).inverse(); 
 	GL_dn = (OM - (NM - V*I) -FM_NM_dn_T_dagg*GL_dn*FM_NM_dn_T).inverse(); 
 	GR_up_even = (OM - (NM - V*I)-NM_FM_up_T*GR_up_even*NM_FM_up_T_dagg).inverse();
 	GR_dn_even = (OM - (NM - V*I)-NM_FM_dn_T*GR_dn_even*NM_FM_dn_T_dagg).inverse(); 
 	//adlayer one bilayer of CoCu onto LHS G for odd layers, then adlayer a  
 	//further bilayer of Cu to ensure gmean is correct. This means 3 layers are on before we begin! 
 	GR_up_odd = (OM - (odd_l1_up - V*I) -odd_l1_up_T2*GR_up_odd*odd_l1_up_T2_dagg).inverse(); 
 	GR_dn_odd = (OM - (odd_l1_dn - V*I) -odd_l1_dn_T2*GR_dn_odd*odd_l1_dn_T2_dagg).inverse(); 
 	GR_up_odd = (OM - (NM - V*I) -odd_l1_up_T1*GR_up_odd*odd_l1_up_T1_dagg).inverse(); 
 	GR_dn_odd = (OM - (NM - V*I) -odd_l1_dn_T1*GR_dn_odd*odd_l1_dn_T1_dagg).inverse(); 

 	//initialise integration so that they don't all get summed over every E! 
 	for (int ll = 0; ll < N; ll++){ 
 		if (needi[ll] == 1) 
 			fm[ll + k] = 0.; 
 	} 
 	const int nn = 10; 
 	/* const int nn = 1; */ 
 	for (int kk=0; kk<nn+1; kk++) { 
 		theta = kk*M_PI/nn; 
 		f(result, theta, E, k, needi, send, GR_up_even, GR_dn_even, GR_up_odd, GR_dn_odd, GL_up, GL_dn); 
 		if ((kk==0)||(kk==nn)){ 
 			for (int ll = 0; ll < N; ll++){ 
 				if (needi[ll] == 1) 
 					fm[ll + k] += M_PI*(0.5/nn)*result[ll]; 
 			} 
 		} 
 		else { 
 			for (int ll = 0; ll < N; ll++){ 
 				if (needi[ll] == 1) 
 					fm[ll + k] += (M_PI/nn)*result[ll]; 
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
	   //Arrays 
 	char cvalue[17];
 	double *com = 0, *dinest = 0, *errest = 0, *fm = 0, *opts = 0, *x = 0;
 	Integer *icom = 0, *iopts = 0, *needi = 0;

	  // NAG types 
 	Nag_VariableType optype;
 	NagError fail;

	   //Setup phase. 
	   //Set problem parameters. 
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
	 //Initialize option arrays using nag_quad_opt_set (d01zkc). 
 	nag_quad_opt_set("Initialize = nag_quad_1d_gen_vec_multi_rcomm", iopts, liopts, opts, lopts, &fail);
 	if (fail.code != NE_NOERROR) {
 		cout<<"Error from nag_quad_opt_set (d01zkc)."<<endl<<fail.message<<endl;
 		exit(EXIT_FAILURE);
 	}
 	nag_quad_opt_set("Quadrature Rule = gk15", iopts, liopts, opts, lopts, &fail);
	 //nag_quad_opt_set("Quadrature Rule = gk21", iopts, liopts, opts, lopts, &fail); 
	 //nag_quad_opt_set("Quadrature Rule = gk31", iopts, liopts, opts, lopts, &fail); 
	 //nag_quad_opt_set("Quadrature Rule = gk41", iopts, liopts, opts, lopts, &fail); 
	 //nag_quad_opt_set("Quadrature Rule = gk51", iopts, liopts, opts, lopts, &fail); 
	 //nag_quad_opt_set("Quadrature Rule = gk61", iopts, liopts, opts, lopts, &fail); 
 	nag_quad_opt_set("Absolute Tolerance = 1.0e-6", iopts, liopts, opts, lopts, &fail);
 	nag_quad_opt_set("Relative Tolerance = 1.0e-6", iopts, liopts, opts, lopts, &fail);

	 // Determine required array dimensions for 
	 // nag_quad_1d_gen_vec_multi_rcomm (d01rac) using 
	 // nag_quad_1d_gen_vec_multi_dimreq (d01rcc). 
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

	 //Allocate remaining arrays. 
 	if (!(x = NAG_ALLOC((lenx), double)) ||	!(needi = NAG_ALLOC((ni), Integer)) || !(fm = NAG_ALLOC((ldfm) * (sdfm), double)) ||
 		!(dinest = NAG_ALLOC((ni), double)) || !(errest = NAG_ALLOC((ni), double)) ||
 	       	!(com = NAG_ALLOC((lcom), double)) || !(icom = NAG_ALLOC((licom), Integer))){
 		cout<<"Allocation failure"<<endl;
 		exit(EXIT_FAILURE);
 	}

	 //Solve phase. 
 	INIT_FAIL(fail);
	 //Set initial irevcm. 
 	irevcm = 1;
 	while (irevcm) {
		 // nag_quad_1d_gen_vec_multi_rcomm (d01rac). 
		 // One-dimensional quadrature, adaptive, vectorized, multi-integral, 
		 // reverse communication. 
		 
 		nag_quad_1d_gen_vec_multi_rcomm(&irevcm, ni, a, b, &sid, needi, x, lenx, &nx, fm, ldfm,
                                     dinest, errest, iopts, opts, icom, licom, com, lcom, &fail);
 		switch (irevcm) {
 			case 11:
				 // Initial returns. 
				 // These will occur during the non-adaptive phase. 
				 // All values must be supplied. 
				 // dinest and errest do not contain approximations over the complete 
			 	 // interval at this stage. 
				 // 
 				pass(x, nx, ldfm, fm, needi, send);
 				break;
 			case 12:
				 // Intermediate returns. 
				 // These will occur during the adaptive phase. 
				 // All requested values must be supplied. 
				 // dinest and errest contain approximations over the complete 
				 // interval at this stage. 
				  // 
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

vector<double> switching(variables * send) {//TODO we need to check that spin up/down is catered for in the Hams below
	double x = send->x;
	double z = send->z;
	int numlay = send->numlay;
	int numbas = send->numbas;
	vector<bool> isMag = *send->isMag;
	unordered_map<string, Vector3d> lat_oop = *send->lat_oop;
	vec3 odd_oop = *send->odd_oop;
	unordered_map<string, vector<vvec3>> pos = *send->pos;
	unordered_map<string, vector<vvec3>> odd_pos = *send->odd_pos;
	unordered_map<string, vector<vector<vM>>> atom_up = *send->atom_up;
	unordered_map<string, vector<vector<vM>>> atom_dn = *send->atom_dn;
	unordered_map<string, vector<vector<vM>>> odd_atom_up = *send->odd_atom_up;
	unordered_map<string, vector<vector<vM>>> odd_atom_dn = *send->odd_atom_dn;
	vector<vector<string>> attype = *send->attype;
	vector<vector<string>> oddattype = *send->oddattype;
	vvec3 basis = *send->basis;
	vvec3 odd_basis = *send->odd_basis;
	vMXd HUu, HUd, HTu, HTd, odd_Hu, odd_Hd;
	HUu.reserve(numlay);
	HUd.reserve(numlay);
	HTu.reserve(2*numlay - 1);
	HTd.reserve(2*numlay - 1);
	odd_Hu.reserve(3);
	odd_Hd.reserve(3);
	MatrixXcd Htmp_up(numbas*9, numbas*9);
	MatrixXcd Htmp_dn(numbas*9, numbas*9);

	M9 tmp_mat_up, tmp_mat_dn;
	int calc, calc2;
	string atomtype;

	//diagonalised in-plane
	//TODO maybe recheck some of the tricks below when using 4-atom basis
	for (int at = 0; at < numlay; at++){
	  for (int ii = 0; ii < numbas; ii++){
	    for (int jj = 0; jj < numbas; jj++){
	      if (attype[at][ii] == attype[at][jj])
		atomtype = attype[at][ii];
	      else
		atomtype = attype[at][jj] + attype[at][ii];
	      calc = 0;
	      calc2 = 0;
	      for (auto const& k : pos){
		if (k.first == atomtype){
		  if (ii == jj){//this section to avoid duplication of work. Checks to see if calculation has already been made
		    for (int kk = 0; kk < ii; kk++){
		      if(attype[at][kk] == attype[at][ii]){
		        tmp_mat_up = Htmp_up.block<9,9>(9*kk,9*kk);
		        tmp_mat_dn = Htmp_dn.block<9,9>(9*kk,9*kk);
			calc2 = 1;
			break;
		      }
		    }
		  }
		  if (calc2 == 0){
		    for (int kk = 0; kk < at; kk++){//this section to avoid duplication of work. Checks to see if calculation has already been made
	              if ((attype[at][ii] == attype[kk][ii]) && (attype[at][jj] == attype[kk][jj])){
		        tmp_mat_up = HUu[kk].block<9,9>(9*jj,9*ii);
		        tmp_mat_dn = HUd[kk].block<9,9>(9*jj,9*ii);
		        calc = 1;
		        break;
		      }
		    }
		    if (calc == 0){
		      tmp_mat_up = InPlaneH( pos[atomtype][ii][jj], basis[at][ii] - basis[at][jj], atom_up[atomtype][ii][jj], x, z);
		      if (isMag[at] == true)
		        tmp_mat_dn = InPlaneH( pos[atomtype][ii][jj], basis[at][ii] - basis[at][jj], atom_dn[atomtype][ii][jj], x, z);
		      else 
		        tmp_mat_dn = tmp_mat_up;
		    }
		  }
		  Htmp_up.block<9,9>(9*jj,9*ii) = tmp_mat_up;
		  Htmp_dn.block<9,9>(9*jj,9*ii) = tmp_mat_dn;
		  break;
		}
              }
	    }
	  }
	  HUu.emplace_back(Htmp_up);
	  HUd.emplace_back(Htmp_dn);
	}

	//hopping out of plane
	for (int at = 0; at < numlay; at++){
	  for (int inter = 1; inter > -1; inter--){
            if ((at == 0) && (inter != 0))//don't want to reference index < 0 or duplicate data
              continue;
	    for (int ii = 0; ii < numbas; ii++){
	      for (int jj = 0; jj < numbas; jj++){
	        if (attype[at][ii] == attype[at - inter][jj])
		  atomtype = attype[at][ii];
		else
		  atomtype = attype[at - inter][jj] + attype[at][ii];
	        calc = 0;
	        for (auto const& k : pos){
		  if (k.first == atomtype){
		    if ((ii == jj) && (inter == 0)){//this section to avoid duplication of work. Checks to see if calculation has already been made
		      for (int kk = 0; kk < ii; kk++){
		        if(attype[at][kk] == attype[at][ii]){
		          tmp_mat_up = Htmp_up.block<9,9>(9*kk,9*kk);
		          tmp_mat_dn = Htmp_dn.block<9,9>(9*kk,9*kk);
		  	  calc = 1;
			  break;
		        }
		      }
		    }
		    if (calc == 0){
		      tmp_mat_up = InPlaneH( pos[atomtype][ii][jj], lat_oop[atomtype] + basis[at][ii] - basis[at - inter][jj], atom_up[atomtype][ii][jj], x, z);
		      if ((isMag[at] == true) || (isMag[at - inter] == true))
		        tmp_mat_dn = InPlaneH( pos[atomtype][ii][jj], lat_oop[atomtype] + basis[at][ii] - basis[at - inter][jj], atom_dn[atomtype][ii][jj], x, z);
		      else 
		        tmp_mat_dn = tmp_mat_up;
		    }
		    Htmp_up.block<9,9>(9*jj,9*ii) = tmp_mat_up;
		    Htmp_dn.block<9,9>(9*jj,9*ii) = tmp_mat_dn;
		    break;
		  }
		}
              }
	    }
	    HTu.emplace_back(Htmp_up);
	    HTd.emplace_back(Htmp_dn);
	  }
	}

	/* for (int kk = 0; kk < HUu.size(); kk++) */
	/* 	cout<<HUu[kk]<<endl<<endl; */
	/* for (int kk = 0; kk < HUd.size(); kk++) */
	/* 	cout<<HUd[kk]<<endl<<endl; */

	int at;
	for (int c = 0; c < 3; c++){
	  for (int ii = 0; ii < numbas; ii++){
	    for (int jj = 0; jj < numbas; jj++){
	      if (c == 0){
		at = 1;
	        if (oddattype[at][ii] == oddattype[at][jj])
	          atomtype = oddattype[at][ii];
	        else
	          atomtype = oddattype[at][jj] + oddattype[at][ii];
	      }
	      else {
		at = c;
	        if (oddattype[at][ii] == oddattype[at - 1][jj])
	          atomtype = oddattype[at][ii];
	        else
	          atomtype = oddattype[at - 1][jj] + oddattype[at][ii];
	      }
	      for (auto const& k : odd_pos){
	        if (k.first == atomtype){
		  if (c == 0){
	            tmp_mat_up = InPlaneH( odd_pos[atomtype][ii][jj], odd_basis[at][ii] - odd_basis[at][jj], odd_atom_up[atomtype][ii][jj], x, z);
	            tmp_mat_dn = InPlaneH( odd_pos[atomtype][ii][jj], odd_basis[at][ii] - odd_basis[at][jj], odd_atom_dn[atomtype][ii][jj], x, z);
		  }
		  else {
		    tmp_mat_up = InPlaneH( odd_pos[atomtype][ii][jj], odd_oop[at-1] + odd_basis[at][ii] - odd_basis[at - 1][jj], odd_atom_up[atomtype][ii][jj], x, z);
		    tmp_mat_dn = InPlaneH( odd_pos[atomtype][ii][jj], odd_oop[at-1] + odd_basis[at][ii] - odd_basis[at - 1][jj], odd_atom_dn[atomtype][ii][jj], x, z);
		  }
		  Htmp_up.block<9,9>(9*jj,9*ii) = tmp_mat_up;
		  Htmp_dn.block<9,9>(9*jj,9*ii) = tmp_mat_dn;
		  break;
		}
              }
	    }
	  }
	  odd_Hu.emplace_back(Htmp_up);
	  odd_Hd.emplace_back(Htmp_dn);
	}

	send->HUu = &HUu;
	send->HUd = &HUd;
	send->HTu = &HTu;
	send->HTd = &HTd;
	send->odd_Hu = &odd_Hu;
	send->odd_Hd = &odd_Hd;

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

static double NAG_CALL fa(double x, double y, Nag_Comm *comm)
{
	cout<<"Loop beginning"<<endl;
	variables * send = (variables *) comm->p;
	int i = send->it;
	cout<<i<<endl;
	Vector3d b1 = *send->b1;
	Vector3d b2 = *send->b2;
	double result;
	Vector3d xk;
	cout<<i<<endl;
	unordered_map<double, unordered_map<double, vector<double>>> all_the_data = *send->all_the_data;
	cout<<i<<endl;
	if (!all_the_data.count(x)){
		vector<double> integrate;
		xk = 0.5*x*b1 + 0.5*y*b2;
		send->x = xk(0);
		send->z = xk(2);
		cout<<"in 1st if"<<endl;
		cout<<i<<endl;
		integrate = switching(send);
		cout<<i<<endl;
		all_the_data[x][y] = integrate;
		cout<<i<<endl;
		result = integrate[i];
		cout<<i<<endl;
		send->all_the_data = &all_the_data;
		cout<<i<<endl;
	}
	else if (!all_the_data[x].count(y)){
		vector<double> integrate;
		xk = 0.5*x*b1 + 0.5*y*b2;
		send->x = xk(0);
		send->z = xk(2);
		cout<<"in 2nd if"<<endl;
		cout<<i<<endl;
		integrate = switching(send);
		cout<<i<<endl;
		all_the_data[x][y] = integrate;
		cout<<i<<endl;
		result = integrate[i];
		cout<<i<<endl;
		send->all_the_data = &all_the_data;
		cout<<i<<endl;
	}
	else {
		cout<<"in 3rd if"<<endl;
		cout<<i<<endl;
		result = all_the_data[x][y][i];
		cout<<i<<endl;
	}
	cout<<"out"<<endl;
	cout<<i<<endl;
	comm->p = send;
	cout<<i<<endl;
	return result;
}

int main() 
{
	//EDITABLE INPUT DATA HERE TODO COULD LOOK AT QE/WANNIER LIKE INPUT FILE
	// plot output of spincurrent against energy
	const double Ef = 0.;//set Fermi level to zero throughout in header file
	// number of spacer layers
	const int N = 10;
	// set bias
	const double V = 0.05;
	const int numlay = 5;//number of slabs i.e. Au|MgO|Fe|Au|Fe
	//figure out a way to determine this from header file
	vector<bool> isMag;//is each layer magnetic?
	isMag.emplace_back(false);
	isMag.emplace_back(false);
	isMag.emplace_back(true);
	isMag.emplace_back(false);
	isMag.emplace_back(true);
	//thickness of intermediate layers (except spacer)
	vector<int> thick;
	//set number of insulator principle layers
	thick.emplace_back(2);//remember in halite structures the two basis atoms are on the same plane
	//set number of LH FM bilayers
	thick.emplace_back(10);
	if (thick.size() + 3 - numlay != 0){
		cerr<<"Error! Intermediate layer thickness have not all been set"<<endl;
		exit(EXIT_FAILURE);

	}
	
	const double k = 8.617e-5/13.6058;//boltzmann constant (in Ryds)
	const double T = 300;//set the temperature
	const int numbas = 2;//number of basis atoms
	//TODO in this version of the code, I still need to manually define each layer

	//in-plane lattice vectors for the whole system;
	Vector3d lat_vec1, lat_vec2;
	lat_vec1 << 0.5, 0, 0.5;
	lat_vec2 << 0.5, 0, -0.5;

	unordered_map<string, Vector3d> lat_oop;
	unordered_map<string, vector<double>> dist;
	vector<double> dist_tmp;

	vvec3 basis;
	vec3 test_vec;
	test_vec.reserve(numbas);
	basis.reserve(numlay);
	Vector3d tmp_vec;
	tmp_vec << 0, 0, 0;
	vector<string> tmp_at;
	vector<vector<string>> attype;
	vector<vector<string>> oddattype;
	string tmp_atom = "";
	for (int k = 0; k < numbas; k++){
		test_vec.emplace_back(tmp_vec);
		tmp_at.emplace_back(tmp_atom);
	}
	for (int k = 0; k < numlay; k++){
		basis.emplace_back(test_vec);
		attype.emplace_back(tmp_at);
	}
	for (int k = 0; k < 3; k++)
		oddattype.emplace_back(tmp_at);
	
	//This section defines the basis atoms and lattice vectors for Au 
	basis[0][0] << 0., 0., 0.;
	basis[0][1] << 0.5, 0.5, 0.;
	tmp_atom = species(1);//1 is Au
	attype[0][0] = tmp_atom;
	attype[0][1] = tmp_atom;
	lat_oop[tmp_atom] << 0, 1, 0;
	//Distance information for n.n and n.n.n
	dist_tmp.emplace_back(M_SQRT2/2.);
	dist_tmp.emplace_back(1.);
	dist_tmp.emplace_back(1.2247);//the hoppings are zero here
	dist[tmp_atom] = dist_tmp;
	dist_tmp.clear();

	//This section defines the basis atoms and lattice vectors for MgO 
	basis[1][0]<< 0., 0., 0.;//Mg
	basis[1][1]<< 0.5, 0., 0.;//O
	attype[1][0] = species(3);//3 is Mg
	attype[1][1] = species(4);//4 is O 
	dist_tmp.emplace_back(0.5);
	dist_tmp.emplace_back(M_SQRT2/2.);
	dist_tmp.emplace_back(0.5*sqrt(3.));
	lat_oop[species(3)] << 0.5, 0.5, 0;
	lat_oop[species(4)] << 0.5, 0.5, 0;
	tmp_atom = species(3) + species(4);
	dist[tmp_atom] = dist_tmp;
	lat_oop[tmp_atom] << 0.5, 0.5, 0;
	tmp_atom = species(4) + species(3);
	dist[tmp_atom] = dist_tmp;
	lat_oop[tmp_atom] << 0.5, 0.5, 0;
	dist[species(3)] = dist_tmp;
	dist[species(4)] = dist_tmp;
	//clear for re-use
	dist_tmp.clear();

	//This section defines the basis atoms and lattice vectors for Fe 
	basis[2][0] << 0., 0., 0.;
	basis[2][1] << 0.5, 0.25*M_SQRT2, 0.;
	tmp_atom = species(2);//2 is Fe
	attype[2][0] = tmp_atom;
	attype[2][1] = tmp_atom;
	lat_oop[tmp_atom] << 0, 0.5*M_SQRT2, 0;
	//Distance information for n.n and n.n.n
	dist_tmp.emplace_back(0.5*sqrt(3.)/M_SQRT2);
	dist_tmp.emplace_back(0.5*M_SQRT2);
	dist_tmp.emplace_back(1.);
	dist[tmp_atom] = dist_tmp;
	dist_tmp.clear();

	//then the duplicate layers
	basis[3][0] = basis[0][0];
	basis[3][1] = basis[0][1];
	tmp_atom = species(1);//1 is Au
	attype[3][0] = tmp_atom;
	attype[3][1] = tmp_atom;
	basis[4][0] = basis[2][0];
	basis[4][1] = basis[2][1];
	tmp_atom = species(2);//2 is Fe
	attype[4][0] = tmp_atom;
	attype[4][1] = tmp_atom;

	lat_oop[species(1) + species(3)] << 0, 0.605 + 0.5, 0;// Au-MgO distance taken as 2.47 A from dx.doi.org/10.1021/jz4022975 | J. Phys. Chem. Lett. 2014, 5, 131−137
	lat_oop[species(1) + species(4)] << 0, 0.605 + 0.5, 0;// distance taken as 2.47 A from dx.doi.org/10.1021/jz4022975 | J. Phys. Chem. Lett. 2014, 5, 131−137
	lat_oop[species(3) + species(2)] << 0.5, 0.505, 0;// MgO-Fe distance taken as 2.06 A from J. Phys. D: Appl. Phys. 42 (2009) 015003 (5pp) doi:10.1088/0022-3727/42/1/015003
	lat_oop[species(4) + species(2)] << 0.5, 0.505, 0;// MgO-Fe distance taken as 2.06 A from J. Phys. D: Appl. Phys. 42 (2009) 015003 (5pp) doi:10.1088/0022-3727/42/1/015003
	lat_oop[species(1) + species(2)] << 0, 0.446 + 0.5, 0;// Au-Fe distance taken as 1.82 A from Surface Science 370 (1997)293-310
	lat_oop[species(2) + species(1)] << 0, 0.446 + 0.25*M_SQRT2, 0;
	/* ins_met_lat_oop1 << 0., 1., 0.; */
	/* ins_met_lat_oop2 << 0., .5, 0.; */
	/* ins_met_lat_oop << 0., 0.6393, 0.;// this from LiCl paper detailing distance between LiCl and Co/Cu. */

	unordered_map<string, vector<int>> ntype;//keeps track of which distance is for which neighbour
	vector<int> ntype_tmp;
	int numnns = numnn();
	//O-Mg (and Mg-O)
	dist_tmp.emplace_back(0.5);
	dist_tmp.emplace_back(M_SQRT2/2.);
	dist_tmp.emplace_back(0.5*sqrt(3.));
	ntype_tmp.emplace_back(1);
	ntype_tmp.emplace_back(2);
	ntype_tmp.emplace_back(3);
	dist[species(3) + species(4)] = dist_tmp;
	dist[species(4) + species(3)] = dist_tmp;
	ntype[species(3) + species(4)] = ntype_tmp;
	ntype[species(4) + species(3)] = ntype_tmp;
	dist_tmp.clear();
	ntype_tmp.clear();
	
	//Au-Fe (and Fe-Au)
	dist_tmp.emplace_back(0.670012);//the code takes into account differing NN distances due to lattice compression
	dist_tmp.emplace_back(0.7996);//or expansion at interfaces - 2nd nn
	dist_tmp.emplace_back(0.946);// - 2nd nn
	dist_tmp.emplace_back(1.06737);//3rd nn
	dist_tmp.emplace_back(1.181066);//3rd nn
	ntype_tmp.emplace_back(1);
	ntype_tmp.emplace_back(2);
	ntype_tmp.emplace_back(2);
	ntype_tmp.emplace_back(3);
	ntype_tmp.emplace_back(3);
	dist[species(1) + species(2)] = dist_tmp;
	dist[species(2) + species(1)] = dist_tmp;
	ntype[species(1) + species(2)] = ntype_tmp;
	ntype[species(2) + species(1)] = ntype_tmp;
	dist_tmp.clear();
	ntype_tmp.clear();

	//Au-Mg, Au-O
	dist_tmp.emplace_back(0.605);//Au to O 
	dist_tmp.emplace_back(0.784873);//Au to Mg
	dist_tmp.emplace_back(0.930605);//Au to O
	ntype_tmp.emplace_back(1);
	ntype_tmp.emplace_back(2);
	ntype_tmp.emplace_back(3);
	dist[species(1) + species(3)] = dist_tmp;
	dist[species(1) + species(4)] = dist_tmp;
	ntype[species(1) + species(3)] = ntype_tmp;
	ntype[species(1) + species(4)] = ntype_tmp;
	dist_tmp.clear();
	ntype_tmp.clear();

	//Mg-Fe, O-Fe
	dist_tmp.emplace_back(0.505);//O to Fe 
	dist_tmp.emplace_back(0.710652);//Mg to Fe
	dist_tmp.emplace_back(0.8585534);//Mg to Fe 3rd nn
	dist_tmp.emplace_back(0.868922);//O to Fe 3rd nn
	ntype_tmp.emplace_back(1);
	ntype_tmp.emplace_back(2);
	ntype_tmp.emplace_back(3);
	ntype_tmp.emplace_back(3);
	dist[species(3) + species(2)] = dist_tmp;
	dist[species(4) + species(2)] = dist_tmp;
	ntype[species(3) + species(2)] = ntype_tmp;
	ntype[species(4) + species(2)] = ntype_tmp;
	dist_tmp.clear();
	ntype_tmp.clear();
	
	//odd layer stuff TODO this will need editing if placement of odd layer changes
	vvec3 odd_basis;
	vec3 odd_oop;
	//old position
	/* odd_basis.emplace_back(basis[2]);//Fe edit */
	/* odd_basis.emplace_back(basis[3]);//Au edit */
	/* odd_basis.emplace_back(basis[3]);//Au edit */
	/* odd_basis[1][1] << 0.5, 0.446, 0; */
	/* odd_oop.emplace_back(lat_oop[species(2)]); */
	/* odd_oop.emplace_back(lat_oop[species(1) + species(2)]); */
	/* odd_oop.emplace_back(lat_oop[species(1)]); */
	/* oddattype[0][0] = species(2);//edit */
	/* oddattype[0][1] = species(2);//edit */
	/* oddattype[1][0] = species(2);//edit */
	/* oddattype[1][1] = species(1);//edit */
	/* oddattype[2][0] = species(1);//edit */
	/* oddattype[2][1] = species(1);//edit */
	//new position - next to LH lead
	odd_basis.emplace_back(basis[3]);//Fe edit
	odd_basis.emplace_back(basis[4]);//Au edit
	odd_basis.emplace_back(basis[4]);//Au edit
	odd_basis[1][1] << 0.5, 0.446, 0;
	odd_oop.emplace_back(lat_oop[species(1)]);
	odd_oop.emplace_back(lat_oop[species(2) + species(1)]);
	odd_oop.emplace_back(lat_oop[species(2)]);
	oddattype[0][0] = species(1);//edit
	oddattype[0][1] = species(1);//edit
	oddattype[1][0] = species(1);//edit
	oddattype[1][1] = species(2);//edit
	oddattype[2][0] = species(2);//edit
	oddattype[2][1] = species(2);//edit
	//END OF EDITABLE DATA, UNLESS DEBUGGING OR IMPROVING, OR REIMPLEMENTING, DON'T CHANGE THE CODE BELOW

	//This block creates the SK tight binding Hamiltonians
	unordered_map<string, M9> onsite_up, onsite_dn;
	unordered_map<string, vector<vector<double>>> hop_up, hop_dn;
	int numats = numatoms();
	vector<double> temporary_vector;
	for (int jj = 1; jj < numats + 1; jj++){
		onsite_up[species(jj)] = U(jj, 0);
		onsite_dn[species(jj)] = U(jj, 1);
		for (int kk = 1; kk < numnns + 1; kk++){
			temporary_vector = param(jj, kk, 0);
			hop_up[species(jj)].emplace_back(temporary_vector);
			temporary_vector = param(jj, kk, 1);
			hop_dn[species(jj)].emplace_back(temporary_vector);
		}
	}
	temporary_vector.clear();
	vector<double> AuMg1;
	double tmp;
	double order;
	string tmpat1, tmpat2;
	int l;
	//This loop creates the geometric means used at the interfaces between elements
	for (int ii = 1; ii < numats + 1; ii++){
		for (int jj = 1; jj < numats + 1; jj++){
			tmpat1 = species(ii);
			tmpat2 = species(jj);
			tmp_atom = tmpat1 + tmpat2;
			if (dist.count(tmp_atom)){
				for (int m = 0; m < dist[tmp_atom].size(); m++){
					l = ntype[tmp_atom][m] - 1;
					for (int kk = 0; kk < 10; kk++){
						if (kk < 4)//this block gives the Harrison formula for the distance dependence of the SK potentials
							order = 2.;
						else if (kk < 7)
							order = 3.5;
						else 
							order = 5.;
						//if Au/Fe and MgO then just scale Fe/Au to the differing distances e.g. no gmean to MgO
						if ((tmpat1 == "Au") && ((tmpat2 == "Mg") || (tmpat2 == "O")))//atom names are explicit here as we only need to do this for this system
							tmp = gmean(hop_up[tmpat1][l][kk], hop_up[tmpat1][l][kk], dist[tmpat1][l], dist[tmpat1][l], dist[tmp_atom][m], order);
						else if ((tmpat2 == "Au") && ((tmpat1 == "Mg") || (tmpat1 == "O")))//atom names are explicit here as we only need to do this for this system
							tmp = gmean(hop_up[tmpat2][l][kk], hop_up[tmpat2][l][kk], dist[tmpat2][l], dist[tmpat2][l], dist[tmp_atom][m], order);
						else if ((tmpat1 == "Fe") && ((tmpat2 == "Mg") || (tmpat2 == "O")))
							tmp = gmean(hop_up[tmpat1][l][kk], hop_up[tmpat1][l][kk], dist[tmpat1][l], dist[tmpat1][l], dist[tmp_atom][m], order);
						else if ((tmpat2 == "Fe") && ((tmpat1 == "Mg") || (tmpat1 == "O")))
							tmp = gmean(hop_up[tmpat2][l][kk], hop_up[tmpat2][l][kk], dist[tmpat2][l], dist[tmpat2][l], dist[tmp_atom][m], order);
						else
							tmp = gmean(hop_up[tmpat1][l][kk], hop_up[tmpat2][l][kk], dist[tmpat1][l], dist[tmpat2][l], dist[tmp_atom][m], order);
						temporary_vector.emplace_back(tmp);
					}
					hop_up[tmp_atom].emplace_back(temporary_vector);
					temporary_vector.clear();
					for (int kk = 0; kk < 10; kk++){
						if (kk < 4)//this block gives the Harrison formula for the distance dependence of the SK potentials
							order = 2.;
						else if (kk < 7)
							order = 3.5;
						else 
							order = 5.;
						if ((tmpat1 == "Au") && ((tmpat2 == "Mg") || (tmpat2 == "O")))//atom names are explicit here as we only need to do this for this system
							tmp = gmean(hop_dn[tmpat1][l][kk], hop_dn[tmpat1][l][kk], dist[tmpat1][l], dist[tmpat1][l], dist[tmp_atom][m], order);
						else if ((tmpat2 == "Au") && ((tmpat1 == "Mg") || (tmpat1 == "O")))//atom names are explicit here as we only need to do this for this system
							tmp = gmean(hop_dn[tmpat2][l][kk], hop_dn[tmpat2][l][kk], dist[tmpat2][l], dist[tmpat2][l], dist[tmp_atom][m], order);
						else if ((tmpat1 == "Fe") && ((tmpat2 == "Mg") || (tmpat2 == "O")))
							tmp = gmean(hop_dn[tmpat1][l][kk], hop_dn[tmpat1][l][kk], dist[tmpat1][l], dist[tmpat1][l], dist[tmp_atom][m], order);
						else if ((tmpat2 == "Fe") && ((tmpat1 == "Mg") || (tmpat1 == "O")))
							tmp = gmean(hop_dn[tmpat2][l][kk], hop_dn[tmpat2][l][kk], dist[tmpat2][l], dist[tmpat2][l], dist[tmp_atom][m], order);
						else
							tmp = gmean(hop_dn[tmpat1][l][kk], hop_dn[tmpat2][l][kk], dist[tmpat1][l], dist[tmpat2][l], dist[tmp_atom][m], order);
						temporary_vector.emplace_back(tmp);
					}
					hop_dn[tmp_atom].emplace_back(temporary_vector);
					temporary_vector.clear();
				}
			}
		}
	}

	//This section generates the Hamiltonians from SK parameters and NN positions
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	vM iron_up;

	Matrix<dcomp, 9, 9> tmp_mat_up, tmp_mat_dn;
	double distance, odd_distance;
	string atomtype;
	int i1start = -2;
	int i1end = 3;
	int i2start = -2;
	int i2end = 3;
	int i3start = 0;
	int i3end = 2;
	//the below block is designed to eliminate duplicate entries per atom type 
	unordered_map<string, vector<vector<int>>> isdone;
	vector<int> isdone_tmp1;
	vector<vector<int>> isdone_tmp2;
	for (int a1 = 0; a1 < numbas; a1++){
	  for (int a2 = 0; a2 < numbas; a2++)
 	    isdone_tmp1.emplace_back(0);
	  isdone_tmp2.emplace_back(isdone_tmp1);
	}

	string oddat;
	Vector3d odd_tmp_vec;
	unordered_map<string, vector<vector<vec3>>> pos, odd_pos;
	unordered_map<string, vector<vector<vM>>> atom_up, atom_dn, odd_atom_up, odd_atom_dn;
	M9 zero;
	Vector3d vec_zero;
	vec_zero << 0, 0, 0;
	zero = M9::Zero();
	vM tmp1;
	vector<vM> tmp2;
	vec3 tmp3;
	vector<vec3> tmp4;
	tmp1.emplace_back(zero);
	tmp3.emplace_back(vec_zero);
	for (int ii = 0; ii < numbas; ii++){
		tmp2.emplace_back(tmp1);
		tmp4.emplace_back(tmp3);
	}

			    //odd block
	for (int at = 0; at < 3; at++){
	  for (int inter = 1; inter > -1; inter--){
            if ((at == 0) && (inter != 0))//don't want to reference index < 0 or duplicate data
              continue;
	    for (int ii = 0; ii < numbas; ii++){
	      for (int jj = 0; jj < numbas; jj++){
		if (oddattype[at][ii] == oddattype[at - inter][jj])
		  oddat = oddattype[at][ii];
		else
		  oddat = oddattype[at - inter][jj] + oddattype[at][ii];
	        for (int i1 = i1start; i1 < i1end; i1++){
	          for (int i2 = i2start; i2 < i2end; i2++){
			if (((at == 1) && (inter == 0)) || ((at == 1) && (inter == 1)) || ((at == 2) && (inter == 1))){
			  int i3 = 1;
			  if ((at == 1) && (inter == 0))
				  i3 = 0;
	              	  odd_tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*odd_oop[at-inter] + odd_basis[at][ii] - odd_basis[at - inter][jj];
	                  odd_distance = 0;
	                  for (int l = 0; l < 3; l++)
	                    odd_distance += odd_tmp_vec(l)*odd_tmp_vec(l);
	                  odd_distance = sqrt(odd_distance);
	                if ((inter == 0) && (odd_distance < 1e-5)){
		            if (!(odd_pos[oddat].size() > 0)){
			      for (int ki = 0; ki < numbas; ki++){//this to avoid segfault
			        odd_pos[oddat].emplace_back(tmp4);
			        odd_atom_up[oddat].emplace_back(tmp2);
			        odd_atom_dn[oddat].emplace_back(tmp2);
	                      }
	                      odd_pos[oddat][ii][jj][0] = odd_tmp_vec;
	                      odd_atom_up[oddat][ii][jj][0] = onsite_up[oddat];
	                      odd_atom_dn[oddat][ii][jj][0] = onsite_dn[oddat];
			    }
			    else {
	                      odd_pos[oddat][ii][jj].emplace_back(odd_tmp_vec);
	                      odd_atom_up[oddat][ii][jj].emplace_back(onsite_up[oddat]);
	                      odd_atom_dn[oddat][ii][jj].emplace_back(onsite_dn[oddat]);
			    }
		        }
		        else if (odd_distance > 1e-5){
	                    for (int kk = 0; kk < dist[oddat].size(); kk++){
	                      if (odd_distance < dist[oddat][kk] + 1e-3){
	                        x = odd_tmp_vec.dot(X)/sqrt(odd_tmp_vec(0)*odd_tmp_vec(0) + odd_tmp_vec(1)*odd_tmp_vec(1) + odd_tmp_vec(2)*odd_tmp_vec(2)); 
	                        y = odd_tmp_vec.dot(Y)/sqrt(odd_tmp_vec(0)*odd_tmp_vec(0) + odd_tmp_vec(1)*odd_tmp_vec(1) + odd_tmp_vec(2)*odd_tmp_vec(2)); 
	                        z = odd_tmp_vec.dot(Z)/sqrt(odd_tmp_vec(0)*odd_tmp_vec(0) + odd_tmp_vec(1)*odd_tmp_vec(1) + odd_tmp_vec(2)*odd_tmp_vec(2)); 
	                        tmp_mat_up = eint1(hop_up[oddat][kk], x, y, z);
	                        tmp_mat_dn = eint1(hop_dn[oddat][kk], x, y, z);
			        if (!(odd_pos[oddat].size() > 0)){
    			          for (int ki = 0; ki < numbas; ki++){//this to avoid segfault
    			            odd_pos[oddat].emplace_back(tmp4);
    			            odd_atom_up[oddat].emplace_back(tmp2);
    			            odd_atom_dn[oddat].emplace_back(tmp2);
    	                          }
	                          odd_pos[oddat][ii][jj][0] = odd_tmp_vec;
	                          odd_atom_up[oddat][ii][jj][0] = tmp_mat_up;
	                          odd_atom_dn[oddat][ii][jj][0] = tmp_mat_dn;
			        }
			        else {
	                          odd_pos[oddat][ii][jj].emplace_back(odd_tmp_vec);
	                          odd_atom_up[oddat][ii][jj].emplace_back(tmp_mat_up);
	                          odd_atom_dn[oddat][ii][jj].emplace_back(tmp_mat_dn);
			        }
	                        break;
			      }
	                    }
		        }
		    }
	          }
	        }
	      }
	    }
	  }
	}

	for (int at = 0; at < numlay; at++){
	  for (int inter = 1; inter > -1; inter--){
            if ((at == 0) && (inter != 0))//don't want to reference index < 0 or duplicate data
              continue;
	    for (int ii = 0; ii < numbas; ii++){
	      for (int jj = 0; jj < numbas; jj++){
		if (attype[at][ii] == attype[at - inter][jj])
		  atomtype = attype[at][ii];
		else
		  atomtype = attype[at - inter][jj] + attype[at][ii];
		//the below block is a continuation designed to eliminate duplicate entries per atom type 
		if (!atom_up.count(atomtype))//this will only happen at ii == jj == 0, if first time this atom is found
		  isdone[atomtype] = isdone_tmp2;
		if (isdone[atomtype][ii][jj] == 1)//this criteria will appear below
		  continue;//no need to break the whole loop if we plan on making more general
	        for (int i1 = i1start; i1 < i1end; i1++){
	          for (int i2 = i2start; i2 < i2end; i2++){
	            for (int i3 = i3start; i3 < i3end; i3++){

	              tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_oop[atomtype] + basis[at][ii] - basis[at - inter][jj];
	              distance = 0;
	              for (int l = 0; l < 3; l++)
	                distance += tmp_vec(l)*tmp_vec(l);
	              distance = sqrt(distance);
	              if ((inter == 0) && (distance < 1e-5)){
			if (!(pos[atomtype].size() > 0)){
			  for (int ki = 0; ki < numbas; ki++){//this to avoid segfault
			    pos[atomtype].emplace_back(tmp4);
			    atom_up[atomtype].emplace_back(tmp2);
			    atom_dn[atomtype].emplace_back(tmp2);
	                  }
	                  pos[atomtype][ii][jj][0] = tmp_vec;
	                  atom_up[atomtype][ii][jj][0] = onsite_up[atomtype];
	                  atom_dn[atomtype][ii][jj][0] = onsite_dn[atomtype];
			}
			else {
	                  pos[atomtype][ii][jj].emplace_back(tmp_vec);
	                  atom_up[atomtype][ii][jj].emplace_back(onsite_up[atomtype]);
	                  atom_dn[atomtype][ii][jj].emplace_back(onsite_dn[atomtype]);
			}
	              }
		      else if (distance > 1e-5){
	                for (int kk = 0; kk < dist[atomtype].size(); kk++){
	                  if (distance < dist[atomtype][kk] + 1e-3){
	                    x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
	                    y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
	                    z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
	                    tmp_mat_up = eint1(hop_up[atomtype][kk], x, y, z);
	                    tmp_mat_dn = eint1(hop_dn[atomtype][kk], x, y, z);
			    if (!(pos[atomtype].size() > 0)){
    			      for (int ki = 0; ki < numbas; ki++){//this to avoid segfault
    			        pos[atomtype].emplace_back(tmp4);
    			        atom_up[atomtype].emplace_back(tmp2);
    			        atom_dn[atomtype].emplace_back(tmp2);
    	                      }
	                      pos[atomtype][ii][jj][0] = tmp_vec;
	                      atom_up[atomtype][ii][jj][0] = tmp_mat_up;
	                      atom_dn[atomtype][ii][jj][0] = tmp_mat_dn;
			    }
			    else {
	                      pos[atomtype][ii][jj].emplace_back(tmp_vec);
	                      atom_up[atomtype][ii][jj].emplace_back(tmp_mat_up);
	                      atom_dn[atomtype][ii][jj].emplace_back(tmp_mat_dn);
			    }
	                    break;
			  }
	                }
	              }
	            }
	          }
	        }
		isdone[atomtype][ii][jj] = 1;//this will mean that we won't add spurious duplicate elements to pos, atom_up or atom_dn
	      }
	    }
	  }
	}

	/* M9 tmp_mat; */

	/* for (int at = 0; at < numlay; at++){ */
	/*   for (int ii = 0; ii < numbas; ii++){ */
	/*     for (int jj = 0; jj < numbas; jj++){ */
	/*       if (attype[at][ii] == attype[at][jj]) */
	/* 	atomtype = attype[at][ii]; */
	/*       else */
	/* 	atomtype = attype[at][jj] + attype[at][ii]; */
	/*       for (auto const& k : pos){ */
	/* 	if (k.first == atomtype){ */
	/*  	  cout<<k.first<<endl; */
	/* 	  if ((ii == 0) && (jj == 0)) */
	/* 	    cout<<"topLeft"<<endl; */
	/* 	  if ((ii == 0) && (jj == 1)) */
	/* 	    cout<<"bottomLeft"<<endl; */
	/* 	  if ((ii == 1) && (jj == 0)) */
	/*  	    cout<<"topRight"<<endl; */
	/* 	  if ((ii == 1) && (jj == 1)) */
	/* 	    cout<<"bottomRight"<<endl; */
	/* 	      /1* if ((ii == 1) && (jj == 1)){ *1/ */
	/* 		      /1* for (int kk = 0; kk < pos[atomtype][ii][jj].size(); kk++){ *1/ */
	/* 		      /1* if ((pos[atomtype][ii][jj][kk](0) == 0) && (pos[atomtype][ii][jj][kk](1) == 0) && (pos[atomtype][ii][jj][kk](2) == 0)) *1/ */
	/* 			      /1* cout<<atom_up[atomtype][ii][jj][kk]<<endl<<endl; *1/ */
	/* 		      /1* } *1/ */
	/* 	      /1* } *1/ */
	/* 	  tmp_mat = InPlaneH( pos[atomtype][ii][jj], basis[at][ii] - basis[at][jj], atom_up[atomtype][ii][jj], 2.75, 2.75); */
	/* 	  cout<<endl; */
	/* 	  break; */
	/* 	} */
              /* } */
	/*     } */
	/*   } */
	/* } */

	/* //hopping out of plane */
	/* for (int at = 0; at < numlay; at++){ */
	/*   for (int inter = 1; inter > -1; inter--){ */
            /* if ((at == 0) && (inter != 0))//don't want to reference index < 0 or duplicate data */
              /* continue; */
	/*     for (int ii = 0; ii < numbas; ii++){ */
	/*       for (int jj = 0; jj < numbas; jj++){ */
	/*         if (attype[at][ii] == attype[at - inter][jj]) */
	/* 	  atomtype = attype[at][ii]; */
	/* 	else */
	/* 	  atomtype = attype[at - inter][jj] + attype[at][ii]; */
	/*         for (auto const& k : pos){ */
	/* 	  if (k.first == atomtype){ */
	/*  	    cout<<k.first<<endl; */
	/* 	    if ((ii == 0) && (jj == 0)) */
	/* 	      cout<<"topLeft"<<endl; */
	/* 	    if ((ii == 0) && (jj == 1)) */
	/* 	      cout<<"bottomLeft"<<endl; */
	/* 	    if ((ii == 1) && (jj == 0)) */
	/*  	      cout<<"topRight"<<endl; */
	/* 	    if ((ii == 1) && (jj == 1)) */
	/* 	      cout<<"bottomRight"<<endl; */
	/* 	    tmp_mat = InPlaneH( pos[atomtype][ii][jj], lat_oop[atomtype] + basis[at][ii] - basis[at - inter][jj], atom_up[atomtype][ii][jj], 2.75, 2.75); */
	/* 	    cout<<endl; */
	/* 	    break; */
	/* 	  } */
	/* 	} */
	/*       } */
	/*     } */
	/*   } */
	/* } */

	//odd file
	ofstream Myfile3;	
	string Mydata3 = "odd_interface.dat";
	Myfile3.open( Mydata3.c_str(),ios::trunc );

	int idum0=0;
	vector<vector<string>> itype2;
	vvec3 new_bas2;

	for (int iii = 0; iii < 2; iii++){
		itype2.emplace_back(oddattype[0]);
		new_bas2.emplace_back(odd_basis[0]);
	}
	itype2.emplace_back(oddattype[1]);
	new_bas2.emplace_back(odd_basis[1]);
	for (int iii = 0; iii < 2; iii++){
		itype2.emplace_back(oddattype[2]);
		new_bas2.emplace_back(odd_basis[2]);
	}

	Vector3d rr2;
	Vector3d counter2;
	counter2 << 0, 0, 0;
	vector<string> store_string2;
	vec3 store_vec2;

	for (int ilay = 0; ilay < 5; ilay++){
	  for (int i3 = -5; i3 <= 3; i3++){
	    for (int i1=-5; i1 <= 3; i1++){
	      for (int i2= 0; i2 < new_bas2[ilay].size(); i2++){
	       	rr2= i1*lat_vec1+i3*lat_vec2 + new_bas2[ilay][i2] + counter2;
		if ((abs(rr2(0)) < 1.3001) && (abs(rr2(2)) < 1.3001)){
	          idum0++;
		  store_string2.emplace_back(itype2[ilay][i2]);
		  store_vec2.emplace_back(4*rr2);
		}
	      }
	    }
	  }
	  if (ilay == 4)
	    break;
	  if (itype2[ilay][0] == itype2[ilay + 1][0])
	    atomtype = itype2[ilay][0];
	  else
	    atomtype = itype2[ilay][0] + itype2[ilay + 1][0];
	  if (ilay < 2)
	  	counter2 = counter2 + odd_oop[0];
	  else if (ilay < 3)
	  	counter2 = counter2 + odd_oop[1];
	  else
	  	counter2 = counter2 + odd_oop[2];
	}

	Myfile3<<idum0<<endl<<"foo"<<endl;
	for (int kk = 0; kk < store_vec2.size(); kk++)
 	     Myfile3<<store_string2[kk]<<" "<<store_vec2[kk].transpose()<<endl;

	Myfile3.close();

//      whole cluster
	ofstream Myfile2;	
	string Mydata2 = "atoms.dat";
	Myfile2.open( Mydata2.c_str(),ios::trunc );

	idum0 = 0;
	int cluster = 2 + N/2 + 2;
	for (int ii = 0; ii < thick.size(); ii++)
		cluster += thick[ii];
	vector<vector<string>> itype;
	vvec3 new_bas;
	new_bas.reserve(cluster);

	for (int iii = 0; iii < 2; iii++){
		itype.emplace_back(attype[0]);
		new_bas.emplace_back(basis[0]);
	}
	for (int kk = 0; kk < thick.size(); kk++){
		for (int iii = 0; iii < thick[kk]; iii++){
			itype.emplace_back(attype[kk + 1]);
			new_bas.emplace_back(basis[kk + 1]);
		}
	}
	for (int iii = 0; iii < N/2; iii++){
		itype.emplace_back(attype[3]);
		new_bas.emplace_back(basis[3]);
	}
	for (int iii = 0; iii < 2; iii++){
		itype.emplace_back(attype[4]);
		new_bas.emplace_back(basis[4]);
	}

	Vector3d rr;
	Vector3d counter;
	counter << 0, 0, 0;
	vector<string> store_string;
	vec3 store_vec;

	for (int ilay = 0; ilay < cluster; ilay++){
	  for (int i3 = -5; i3 <= 3; i3++){
	    for (int i1=-5; i1 <= 3; i1++){
	      for (int i2= 0; i2 < new_bas[ilay].size(); i2++){
	       	rr= i1*lat_vec1+i3*lat_vec2 + new_bas[ilay][i2] + counter;
		if ((abs(rr(0)) < 1.3001) && (abs(rr(2)) < 1.3001)){
	          idum0++;
		  store_string.emplace_back(itype[ilay][i2]);
		  store_vec.emplace_back(4*rr);
		}
	      }
	    }
	  }
	  if (ilay == cluster - 1)
	    break;
	  if (itype[ilay][0] == itype[ilay + 1][0])
	    atomtype = itype[ilay][0];
	  else
	    atomtype = itype[ilay][0] + itype[ilay + 1][0];
	  counter = counter + lat_oop[atomtype];
	}

	Myfile2<<idum0<<endl<<"foo"<<endl;
	for (int kk = 0; kk < store_vec.size(); kk ++)
 	     Myfile2<<store_string[kk]<<" "<<store_vec[kk].transpose()<<endl;

	Myfile2.close();

	double kT = k*T;
	//set up the variables to send
	variables send;
	send.isMag = &isMag;
	send.attype = &attype;
	send.oddattype = &oddattype;
	send.numbas = numbas;
	send.numlay = numlay;
	send.kT = kT;
	send.Ef = Ef;
	/* send.x = 2.532374; //for now! Neck position of Cu */
	/* send.z = 2.532374; //for now! Neck position of Cu */
	send.x = 2.745644; //for now! Neck position of Au
	send.z = 2.745644; //for now! Neck position of Au
	send.lat_oop = &lat_oop;
	send.odd_oop = &odd_oop;
	send.pos = &pos;
	send.odd_pos = &odd_pos;
	send.basis = &basis;
	send.odd_basis = &odd_basis;
	send.atom_up = &atom_up;
	send.atom_dn = &atom_dn;
	send.odd_atom_up = &odd_atom_up;
	send.odd_atom_dn = &odd_atom_dn;
	send.V = V;
	send.thick = &thick;
	send.N = N;

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
	Myfile.open( Mydata.c_str(),ios::trunc );

	Vector3d b1, b2;
	recip(lat_vec1, lat_vec2, b1, b2);
	send.b1 = &b1;
	send.b2 = &b2;

	unordered_map<double, unordered_map<double, vector<double>>> all_the_data;
	send.all_the_data = &all_the_data;

	Integer exit_status = 0;
	Integer npts;
	double absacc, ans, ya, yb;
	ya = 0;
	yb = 1;
	absacc = 1e-1;
	Nag_Comm comm;
	NagError fail;

	INIT_FAIL(fail);

	/* For communication with user-supplied functions: */

	for (int it = 0; it < N; it++){
		send.it = it;
		comm.p = &send;
		nag_quad_2d_fin(ya, yb, phi1, phi2, fa, absacc, &ans, &npts, &comm, &fail);
		if (fail.code != NE_NOERROR) {
			cout<<"Error from nag_quad_2d_fin (d01dac). "<<fail.message<<endl;
				/* exit(EXIT_FAILURE); */
		}
		cout<<"Number of function evaluations = "<<npts<<endl;
		Myfile<<scientific<<it+4<<" "<<ans<<endl;
	}

	return 0;
}
