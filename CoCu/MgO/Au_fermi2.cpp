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
/* #include <nag.h> */
/* #include <nagd01.h> */
/* #include <nag_stdlib.h> */
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
		Vector3d *Au_lat_oop;
		Vector3d *MgO_lat_oop;
		Vector3d *Fe_lat_oop;
		Vector3d *lat_Au_MgO;
		Vector3d *lat_MgO_Fe;
		Vector3d *lat_Fe_Au;
		Vector3d *lat_Au_Fe;
		vec3 *Au_pos;
		vec3 *MgO_pos_11;
		vec3 *MgO_pos_12;
		vec3 *MgO_pos_21;
		vec3 *Fe_pos;
		vec3 *Au_MgO_pos_11;
		vec3 *Au_MgO_pos_12;
		vec3 *Au_MgO_pos_21;
		vec3 *Au_MgO_pos_22;
		vec3 *MgO_Fe_pos_11;
		vec3 *MgO_Fe_pos_12;
		vec3 *MgO_Fe_pos_21;
		vec3 *MgO_Fe_pos_22;
		vec3 *Au_Fe_pos;
		vec3 *Fe_Au_pos;
		vec3 *Fe_basis;
		vec3 *Au_basis;
		vec3 *MgO_basis;
		vM *iron_up;
		vM *iron_dn;
	        vM *gold;
		vM *magnesium_11;
		vM *magnesium_12;
		vM *magnesium_21;
		vM *oxide_11;
		vM *oxide_12;
		vM *oxide_21;
		vM *gold_iron_up;
		vM *gold_iron_dn;
		vM *iron_gold_up;
		vM *iron_gold_dn;
		vM *iron_dn_MgO_11;
		vM *iron_dn_MgO_12;
		vM *iron_dn_MgO_21;
		vM *iron_dn_MgO_22;
		vM *iron_up_MgO_11;
		vM *iron_up_MgO_12;
		vM *iron_up_MgO_21;
		vM *iron_up_MgO_22;
		vM *gold_MgO_11;
		vM *gold_MgO_12;
		vM *gold_MgO_21;
		vM *gold_MgO_22;
		ddmat *NM;
		ddmat *NM_T;
		ddmat *FM_up;
		ddmat *FM_dn;
		ddmat *FM_up_T;
		ddmat *FM_dn_T;
		ddmat *FM_NM_up_T;
		ddmat *FM_NM_dn_T;
		ddmat *NM_FM_up_T;
		ddmat *NM_FM_dn_T;
		ddmat *odd_l1_up;
		ddmat *odd_l1_dn;
		ddmat *odd_l1_up_T1;
		ddmat *odd_l1_up_T2;
		ddmat *odd_l1_dn_T1;
		ddmat *odd_l1_dn_T2;
		ddmat *ins;
		ddmat *ins_T;
		ddmat *ins_FM_up_T;
		ddmat *ins_FM_dn_T;
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
	Vector3d K;
	dcomp i = -1;
	i = sqrt(i);
	K << x, 0., z;
	M9 result;
	result.fill(0.);
	Vector3d tmp_vec;
	for (int k = 0; k < pos.size(); k++){
		tmp_vec = pos[k] - basis;
		if (abs(tmp_vec(1)) < 1e-5)
			result = result + U[k]*exp(i*tmp_vec.dot(K));
	}
	return result;
}

double int_theta(const dcomp E, variables * send) {
	double result;
	double integrate;
	ddmat FM_up = *send->FM_up;
	ddmat FM_dn = *send->FM_dn;
	ddmat NM = *send->NM;
	ddmat NM_T = *send->NM_T;
	ddmat FM_up_T = *send->FM_up_T;
	ddmat FM_dn_T = *send->FM_dn_T;
	ddmat FM_NM_up_T = *send->FM_NM_up_T;
	ddmat FM_NM_dn_T = *send->FM_NM_dn_T;
	ddmat odd_l1_up = *send->odd_l1_up;
	ddmat odd_l1_dn = *send->odd_l1_dn;
	ddmat odd_l1_up_T1 = *send->odd_l1_up_T1;
	ddmat odd_l1_dn_T1 = *send->odd_l1_dn_T1;
	ddmat odd_l1_up_T2 = *send->odd_l1_up_T2;
	ddmat odd_l1_dn_T2 = *send->odd_l1_dn_T2;
	ddmat ins = *send->ins;
	ddmat ins_T = *send->ins_T;
	ddmat ins_FM_up_T = *send->ins_FM_up_T;
	ddmat ins_FM_dn_T = *send->ins_FM_dn_T;
	ddmat ins_NM_T = *send->ins_NM_T;
	/* cout<<NM<<endl<<endl; */
	double V = send->V;

	ddmat I = ddmat::Identity();
	ddmat OMup=E*I-(FM_up - V*I);
	ddmat OMdn=E*I-(FM_dn - V*I);
	ddmat OM = E*I;

	/* cout<<ins_T.real()<<endl<<endl; */
	ddmat NM_T_dagg = NM_T.adjoint();
	/* ddmat ins_T_dagg = ins_T.adjoint(); */
	/* ddmat FM_up_T_dagg = FM_up_T.adjoint(); */
	/* ddmat FM_dn_T_dagg = FM_dn_T.adjoint(); */
	/* ddmat GR_up = gs(OMup, FM_up_T_dagg); */
	/* ddmat GR_dn = gs(OMdn, FM_dn_T_dagg); */
	ddmat GR_up = gs(OM-NM, NM_T_dagg);
	/* ddmat GR_up = gs(OM-ins, ins_T_dagg); */
	ddmat GR_dn = GR_up;
	/* ddmat FM_NM_up_T_dagg = FM_NM_up_T.adjoint(); */
	/* ddmat FM_NM_dn_T_dagg = FM_NM_dn_T.adjoint(); */
	
	ddmat GL_up = gs(OM - NM, NM_T);
	/* ddmat GL_up = gs(OM - ins, ins_T); */
	ddmat GL_dn = GL_up;
	/* ddmat GL_up = gs(OMup, FM_up_T); */
	/* ddmat GL_dn = gs(OMdn, FM_dn_T); */

	dddmat GR, GL, T;
	GR.fill(0.); GL.fill(0.); T.fill(0.);
	T.topLeftCorner(18,18) = NM_T;
	T.bottomRightCorner(18,18) = NM_T;
	/* T.topLeftCorner(18,18) = ins_T; */
	/* T.bottomRightCorner(18,18) = ins_T; */
	/* T.topLeftCorner(18,18) = FM_up_T; */
	/* T.bottomRightCorner(18,18) = FM_dn_T; */
	dddmat Tdagg = T.adjoint();
	GR.topLeftCorner(18,18) = GR_up;
	GR.bottomRightCorner(18,18) = GR_dn;
	GL.topLeftCorner(18,18) = GL_up;
	GL.bottomRightCorner(18,18) = GL_dn;
	/* cout<<GL.trace()<<" "<<GR.trace()<<endl; */
	dddmat GRinv = GR.inverse();
	dddmat GNinv = GRinv - Tdagg*GL*T;
	dddmat GN = GNinv.inverse();

	return -imag(GN.trace());
}

double switching(double e, variables * send) {//TODO we need to check that spin up/down is catered for in the Hams below
	double x = send->x;
	double z = send->z;
	vec3 Fe_basis = *send->Fe_basis;
	vec3 Au_basis = *send->Au_basis;
	vec3 MgO_basis = *send->MgO_basis;
	vec3 Au_pos = *send->Au_pos;
	vec3 MgO_pos_11 = *send->MgO_pos_11;
	vec3 MgO_pos_12 = *send->MgO_pos_12;
	vec3 MgO_pos_21 = *send->MgO_pos_21;
	vec3 Fe_pos = *send->Fe_pos;
	vec3 Au_MgO_pos_11 = *send->Au_MgO_pos_11;
	vec3 Au_MgO_pos_12 = *send->Au_MgO_pos_12;
	vec3 Au_MgO_pos_21 = *send->Au_MgO_pos_21;
	vec3 Au_MgO_pos_22 = *send->Au_MgO_pos_22;
	vec3 MgO_Fe_pos_11 = *send->MgO_Fe_pos_11;
	vec3 MgO_Fe_pos_12 = *send->MgO_Fe_pos_12;
	vec3 MgO_Fe_pos_21 = *send->MgO_Fe_pos_21;
	vec3 MgO_Fe_pos_22 = *send->MgO_Fe_pos_22;
	vec3 Au_Fe_pos = *send->Au_Fe_pos;
	vec3 Fe_Au_pos = *send->Fe_Au_pos;
	Vector3d Au_lat_oop = *send->Au_lat_oop;
	Vector3d Fe_lat_oop = *send->Fe_lat_oop;
	Vector3d MgO_lat_oop = *send->MgO_lat_oop;
	Vector3d lat_MgO_Fe = *send->lat_MgO_Fe;
	Vector3d lat_Au_MgO = *send->lat_Au_MgO;
	Vector3d lat_Fe_Au = *send->lat_Fe_Au;
	Vector3d lat_Au_Fe = *send->lat_Au_Fe;
	vM iron_up = *send->iron_up;
	vM iron_dn = *send->iron_dn;
	vM gold_iron_dn = *send->gold_iron_dn;
	vM gold_iron_up = *send->gold_iron_up;
        vM gold = *send->gold;
	vM iron_gold_dn = *send->iron_gold_dn;
	vM iron_gold_up = *send->iron_gold_up;
	vM magnesium_11 = *send->magnesium_11;
	vM magnesium_12 = *send->magnesium_12;
	vM magnesium_21 = *send->magnesium_21;
	vM oxide_11 = *send->oxide_11;
	vM oxide_12 = *send->oxide_12;
	vM oxide_21 = *send->oxide_21;
	vM iron_dn_MgO_11 = *send->iron_dn_MgO_11;
	vM iron_dn_MgO_12 = *send->iron_dn_MgO_12;
	vM iron_dn_MgO_21 = *send->iron_dn_MgO_21;
	vM iron_dn_MgO_22 = *send->iron_dn_MgO_22;
	vM iron_up_MgO_11 = *send->iron_up_MgO_11;
	vM iron_up_MgO_12 = *send->iron_up_MgO_12;
	vM iron_up_MgO_21 = *send->iron_up_MgO_21;
	vM iron_up_MgO_22 = *send->iron_up_MgO_22;
	vM gold_MgO_11 = *send->gold_MgO_11;
	vM gold_MgO_12 = *send->gold_MgO_12;
	vM gold_MgO_21 = *send->gold_MgO_21;
	vM gold_MgO_22 = *send->gold_MgO_22;
	M9 ins_11, ins_12, ins_21, ins_22, ins_T_11, ins_T_12, ins_T_21, ins_T_22, 
	   ins_NM_T_11, ins_NM_T_12, ins_NM_T_21, ins_NM_T_22, ins_FM_up_T_11, 
	   ins_FM_up_T_12, ins_FM_up_T_21, ins_FM_up_T_22, ins_FM_dn_T_11, 
	   ins_FM_dn_T_12, ins_FM_dn_T_21, ins_FM_dn_T_22;

	//generate in plane Hamiltonians for simple bilayers
	ins_11 = InPlaneH(MgO_pos_11,  MgO_basis[0], magnesium_11, x, z);//TODO only 2nd NN and onsite depend on atom type, catered for in 11 and 22 only
	ins_12 = InPlaneH(MgO_pos_12,  MgO_basis[1], oxide_12, x, z);//in theory this should contain
	ins_21 = InPlaneH(MgO_pos_21, -MgO_basis[1], oxide_21, x, z);//the same hoppings as magnesium
	ins_22 = InPlaneH(MgO_pos_11,  MgO_basis[0], oxide_11, x, z);
	//generate hopping between simple bilayers of the same type
	ins_T_11 = InPlaneH(MgO_pos_11, MgO_lat_oop + MgO_basis[0], magnesium_11, x, z);//TODO as above except no onsite terms
	ins_T_12 = InPlaneH(MgO_pos_12, MgO_lat_oop + MgO_basis[1], oxide_12, x, z);
	ins_T_21 = InPlaneH(MgO_pos_21, MgO_lat_oop - MgO_basis[1], oxide_21, x, z);
	ins_T_22 = InPlaneH(MgO_pos_11, MgO_lat_oop + MgO_basis[0], oxide_11, x, z);
	//additional off diagonal Hamiltonians needed for bilayers
	//made of different atom types
	ins_NM_T_11 = InPlaneH(Au_MgO_pos_11, lat_Au_MgO + MgO_basis[0] - Au_basis[0], gold_MgO_11, x, z);
	ins_NM_T_12 = InPlaneH(Au_MgO_pos_12, lat_Au_MgO + MgO_basis[1] - Au_basis[0], gold_MgO_12, x, z);
	ins_NM_T_21 = InPlaneH(Au_MgO_pos_21, lat_Au_MgO + MgO_basis[0] - Au_basis[1], gold_MgO_21, x, z);
	ins_NM_T_22 = InPlaneH(Au_MgO_pos_22, lat_Au_MgO + MgO_basis[1] - Au_basis[1], gold_MgO_22, x, z);

	ins_FM_up_T_11 = InPlaneH(MgO_Fe_pos_11, lat_MgO_Fe + Fe_basis[0] - MgO_basis[0], iron_up_MgO_11, x, z);
	ins_FM_up_T_12 = InPlaneH(MgO_Fe_pos_12, lat_MgO_Fe + Fe_basis[1] - MgO_basis[0], iron_up_MgO_12, x, z);
	ins_FM_up_T_21 = InPlaneH(MgO_Fe_pos_21, lat_MgO_Fe + Fe_basis[0] - MgO_basis[1], iron_up_MgO_21, x, z);
	ins_FM_up_T_22 = InPlaneH(MgO_Fe_pos_22, lat_MgO_Fe + Fe_basis[1] - MgO_basis[1], iron_up_MgO_22, x, z);

	ins_FM_dn_T_11 = InPlaneH(MgO_Fe_pos_11, lat_MgO_Fe + Fe_basis[0] - MgO_basis[0], iron_dn_MgO_11, x, z);
	ins_FM_dn_T_12 = InPlaneH(MgO_Fe_pos_12, lat_MgO_Fe + Fe_basis[1] - MgO_basis[0], iron_dn_MgO_12, x, z);
	ins_FM_dn_T_21 = InPlaneH(MgO_Fe_pos_21, lat_MgO_Fe + Fe_basis[0] - MgO_basis[1], iron_dn_MgO_21, x, z);
	ins_FM_dn_T_22 = InPlaneH(MgO_Fe_pos_22, lat_MgO_Fe + Fe_basis[1] - MgO_basis[1], iron_dn_MgO_22, x, z);

	M9 NM_ii, NM_12, NM_21, FM_up_ii, FM_up_12, FM_up_21, FM_dn_ii, FM_dn_12, FM_dn_21, 
	   NM_T_ii, NM_T_12, NM_T_21, FM_up_T_ii, FM_up_T_12, FM_up_T_21, FM_dn_T_ii, 
	   FM_dn_T_12, FM_dn_T_21, FM_NM_up_12, FM_NM_up_21, FM_NM_dn_12, FM_NM_dn_21, 
	   FM_NM_up_T_11, FM_NM_up_T_12, FM_NM_up_T_21, FM_NM_up_T_22,
	   FM_NM_dn_T_11, FM_NM_dn_T_12, FM_NM_dn_T_21, FM_NM_dn_T_22,
	   NM_FM_up_T_11, NM_FM_up_T_12, NM_FM_up_T_21, NM_FM_up_T_22,
	   NM_FM_dn_T_11, NM_FM_dn_T_12, NM_FM_dn_T_21, NM_FM_dn_T_22;
	//generate in plane Hamiltonians for simple bilayers
	NM_ii = InPlaneH( Au_pos,  Au_basis[0], gold, x, z);
	NM_12 = InPlaneH( Au_pos,  Au_basis[1], gold, x, z);
	NM_21 = InPlaneH( Au_pos, -Au_basis[1], gold, x, z);

	FM_up_ii = InPlaneH( Fe_pos,  Fe_basis[0], iron_up, x, z);
	FM_up_12 = InPlaneH( Fe_pos,  Fe_basis[1], iron_up, x, z);
	FM_up_21 = InPlaneH( Fe_pos, -Fe_basis[1], iron_up, x, z);

	FM_dn_ii = InPlaneH( Fe_pos,  Fe_basis[0], iron_dn, x, z);
	FM_dn_12 = InPlaneH( Fe_pos,  Fe_basis[1], iron_dn, x, z);
	FM_dn_21 = InPlaneH( Fe_pos, -Fe_basis[1], iron_dn, x, z);
	//generate hopping between simple bilayers of the same type
	NM_T_ii = InPlaneH( Au_pos, Au_lat_oop + Au_basis[0], gold, x, z);
	NM_T_12 = InPlaneH( Au_pos, Au_lat_oop + Au_basis[1], gold, x, z);
	NM_T_21 = InPlaneH( Au_pos, Au_lat_oop - Au_basis[1], gold, x, z);

	FM_up_T_ii = InPlaneH( Fe_pos, Fe_lat_oop + Fe_basis[0], iron_up, x, z);
	FM_up_T_12 = InPlaneH( Fe_pos, Fe_lat_oop + Fe_basis[1], iron_up, x, z);
	FM_up_T_21 = InPlaneH( Fe_pos, Fe_lat_oop - Fe_basis[1], iron_up, x, z);

	FM_dn_T_ii = InPlaneH( Fe_pos, Fe_lat_oop + Fe_basis[0], iron_dn, x, z);
	FM_dn_T_12 = InPlaneH( Fe_pos, Fe_lat_oop + Fe_basis[1], iron_dn, x, z);
	FM_dn_T_21 = InPlaneH( Fe_pos, Fe_lat_oop - Fe_basis[1], iron_dn, x, z);
	//additional off diagonal Hamiltonians needed for bilayers
	//made of different atom types
	Vector3d X;
	X << 1, 0, 0;//TODO do we need to do this above, where MgO puts NN with exp beyond 1..?
	
	FM_NM_dn_12 = InPlaneH(Fe_Au_pos, X + lat_Fe_Au - Fe_basis[1], iron_gold_dn, x, z);
	FM_NM_dn_21 = InPlaneH(Fe_Au_pos, -(X + lat_Fe_Au - Fe_basis[1]), iron_gold_dn, x, z);
	FM_NM_up_12 = InPlaneH(Fe_Au_pos, X + lat_Fe_Au - Fe_basis[1], iron_gold_up, x, z);
	FM_NM_up_21 = InPlaneH(Fe_Au_pos, -(X + lat_Fe_Au - Fe_basis[1]), iron_gold_up, x, z);

	NM_FM_up_T_11 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[0] - Au_basis[0], gold_iron_up, x, z);
	NM_FM_up_T_12 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[1] - Au_basis[0], gold_iron_up, x, z);
	NM_FM_up_T_21 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[0] - Au_basis[1], gold_iron_up, x, z);
	NM_FM_up_T_22 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[1] - Au_basis[1], gold_iron_up, x, z);

	NM_FM_dn_T_11 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[0] - Au_basis[0], gold_iron_dn, x, z);
	NM_FM_dn_T_12 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[1] - Au_basis[0], gold_iron_dn, x, z);
	NM_FM_dn_T_21 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[0] - Au_basis[1], gold_iron_dn, x, z);
	NM_FM_dn_T_22 = InPlaneH(Au_Fe_pos, lat_Au_Fe + Fe_basis[1] - Au_basis[1], gold_iron_dn, x, z);

	FM_NM_dn_T_11 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[0] - Fe_basis[0], iron_gold_dn, x, z);
	FM_NM_dn_T_12 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[1] - Fe_basis[0], iron_gold_dn, x, z);
	FM_NM_dn_T_21 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[0] - Fe_basis[1], iron_gold_dn, x, z);
	FM_NM_dn_T_22 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[1] - Fe_basis[1], iron_gold_dn, x, z);

	FM_NM_up_T_11 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[0] - Fe_basis[0], iron_gold_up, x, z);
	FM_NM_up_T_12 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[1] - Fe_basis[0], iron_gold_up, x, z);
	FM_NM_up_T_21 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[0] - Fe_basis[1], iron_gold_up, x, z);
	FM_NM_up_T_22 = InPlaneH(Fe_Au_pos, lat_Fe_Au + Au_basis[1] - Fe_basis[1], iron_gold_up, x, z);

	ddmat FM_up, FM_dn, NM, NM_T, FM_up_T, FM_dn_T, FM_NM_up_T, FM_NM_dn_T, odd_l1_up, odd_l1_dn, odd_l1_up_T1, 
	      odd_l1_dn_T1, odd_l1_up_T2, odd_l1_dn_T2;
	ddmat ins, ins_T, ins_NM_T, ins_FM_up_T, ins_FM_dn_T, NM_FM_up_T, NM_FM_dn_T;

	ins.topLeftCorner(9,9) = ins_11;
	ins.topRightCorner(9,9) = ins_12;
	ins.bottomLeftCorner(9,9) = ins_21;
	ins.bottomRightCorner(9,9) = ins_22;
	ins_T.topLeftCorner(9,9) = ins_T_11;
	ins_T.topRightCorner(9,9) = ins_T_12;
	ins_T.bottomLeftCorner(9,9) = ins_T_21;
	ins_T.bottomRightCorner(9,9) = ins_T_22;
	ins_NM_T.topLeftCorner(9,9) = ins_NM_T_11;
	ins_NM_T.topRightCorner(9,9) = ins_NM_T_12; 
	ins_NM_T.bottomLeftCorner(9,9) = ins_NM_T_21;
	ins_NM_T.bottomRightCorner(9,9) = ins_NM_T_22;
	ins_FM_up_T.topLeftCorner(9,9) = ins_FM_up_T_11;
	ins_FM_up_T.topRightCorner(9,9) = ins_FM_up_T_12; 
	ins_FM_up_T.bottomLeftCorner(9,9) = ins_FM_up_T_21;
	ins_FM_up_T.bottomRightCorner(9,9) = ins_FM_up_T_22;
	ins_FM_dn_T.topLeftCorner(9,9) = ins_FM_dn_T_11;
	ins_FM_dn_T.topRightCorner(9,9) = ins_FM_dn_T_12; 
	ins_FM_dn_T.bottomLeftCorner(9,9) = ins_FM_dn_T_21;
	ins_FM_dn_T.bottomRightCorner(9,9) = ins_FM_dn_T_22;

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

	FM_up_T.topLeftCorner(9,9) = FM_up_T_ii;
	FM_up_T.topRightCorner(9,9) = FM_up_T_12;
	FM_up_T.bottomLeftCorner(9,9) = FM_up_T_21;
	FM_up_T.bottomRightCorner(9,9) = FM_up_T_ii;

	FM_dn_T.topLeftCorner(9,9) = FM_dn_T_ii;
	FM_dn_T.topRightCorner(9,9) = FM_dn_T_12;
	FM_dn_T.bottomLeftCorner(9,9) = FM_dn_T_21;
	FM_dn_T.bottomRightCorner(9,9) = FM_dn_T_ii;

	FM_NM_up_T.topLeftCorner(9,9) = FM_NM_up_T_11;
	FM_NM_up_T.topRightCorner(9,9) = FM_NM_up_T_12; 
	FM_NM_up_T.bottomLeftCorner(9,9) = FM_NM_up_T_21;
	FM_NM_up_T.bottomRightCorner(9,9) = FM_NM_up_T_22;

	FM_NM_dn_T.topLeftCorner(9,9) = FM_NM_dn_T_11;
	FM_NM_dn_T.topRightCorner(9,9) = FM_NM_dn_T_12; 
	FM_NM_dn_T.bottomLeftCorner(9,9) = FM_NM_dn_T_21;
	FM_NM_dn_T.bottomRightCorner(9,9) = FM_NM_dn_T_22;

	NM_FM_up_T.topLeftCorner(9,9) = NM_FM_up_T_11;
	NM_FM_up_T.topRightCorner(9,9) = NM_FM_up_T_12; 
	NM_FM_up_T.bottomLeftCorner(9,9) = NM_FM_up_T_21;
	NM_FM_up_T.bottomRightCorner(9,9) = NM_FM_up_T_22;

	NM_FM_dn_T.topLeftCorner(9,9) = NM_FM_dn_T_11;
	NM_FM_dn_T.topRightCorner(9,9) = NM_FM_dn_T_12; 
	NM_FM_dn_T.bottomLeftCorner(9,9) = NM_FM_dn_T_21;
	NM_FM_dn_T.bottomRightCorner(9,9) = NM_FM_dn_T_22;

	odd_l1_up.topLeftCorner(9,9) = FM_up_ii;
	odd_l1_up.topRightCorner(9,9) = FM_NM_up_12;
	odd_l1_up.bottomLeftCorner(9,9) = FM_NM_up_21;
	odd_l1_up.bottomRightCorner(9,9) = NM_ii;

	odd_l1_dn.topLeftCorner(9,9) = FM_dn_ii;
	odd_l1_dn.topRightCorner(9,9) = FM_NM_dn_12;
	odd_l1_dn.bottomLeftCorner(9,9) = FM_NM_dn_21;
	odd_l1_dn.bottomRightCorner(9,9) = NM_ii;

	odd_l1_up_T1.topLeftCorner(9,9) = FM_up_T_ii;
	odd_l1_up_T1.topRightCorner(9,9) = FM_NM_up_T_12;
	odd_l1_up_T1.bottomLeftCorner(9,9) = FM_up_T_21;
	odd_l1_up_T1.bottomRightCorner(9,9) = FM_NM_up_T_11;//TODO I think this is right.. thinking about distances...

	odd_l1_dn_T1.topLeftCorner(9,9) = FM_dn_T_ii;
	odd_l1_dn_T1.topRightCorner(9,9) = FM_NM_dn_T_12;
	odd_l1_dn_T1.bottomLeftCorner(9,9) = FM_dn_T_21;
	odd_l1_dn_T1.bottomRightCorner(9,9) = FM_NM_dn_T_11;//TODO I think this is right.. thinking about distances...

	odd_l1_up_T2.topLeftCorner(9,9) = FM_NM_up_T_22;//TODO as above, but opposite..
	odd_l1_up_T2.topRightCorner(9,9) = FM_NM_up_T_12;
	odd_l1_up_T2.bottomLeftCorner(9,9) = NM_T_21;
	odd_l1_up_T2.bottomRightCorner(9,9) = NM_T_ii;

	odd_l1_dn_T2.topLeftCorner(9,9) = FM_NM_dn_T_22;
	odd_l1_dn_T2.topRightCorner(9,9) = FM_NM_dn_T_12;
	odd_l1_dn_T2.bottomLeftCorner(9,9) = NM_T_21;
	odd_l1_dn_T2.bottomRightCorner(9,9) = NM_T_ii;

	send->NM = &NM;
	send->NM_T = &NM_T;
	send->FM_up = &FM_up;
	send->FM_dn = &FM_dn;
	send->FM_up_T = &FM_up_T;
	send->FM_dn_T = &FM_dn_T;
	send->FM_NM_up_T = &FM_NM_up_T;
	send->FM_NM_dn_T = &FM_NM_dn_T;
	send->NM_FM_up_T = &NM_FM_up_T;
	send->NM_FM_dn_T = &NM_FM_dn_T;
	send->odd_l1_up = &odd_l1_up;
	send->odd_l1_dn = &odd_l1_dn;
	send->odd_l1_up_T1 = &odd_l1_up_T1;
	send->odd_l1_up_T2 = &odd_l1_up_T2;
	send->odd_l1_dn_T1 = &odd_l1_dn_T1;
	send->odd_l1_dn_T2 = &odd_l1_dn_T2;
	send->ins = &ins;
	send->ins_T = &ins_T;
	send->ins_NM_T = &ins_NM_T;
	send->ins_FM_up_T = &ins_FM_up_T;
	send->ins_FM_dn_T = &ins_FM_dn_T;

	/* vector<double> result1, result2, integrate; */
	double integrate;
	int N = send->N;
	double V = send->V;
	/* result1.reserve(N); */
	/* if (abs(V) > 1e-9) */
	/* 	result1 = int_energy(send); */
	/* else { */
	/* 	for (int l = 0; l < N; l++) */
	/* 		result1[l] = 0.; */
	/* } */
	/* result2.reserve(N); */
	/* for (int l = 0; l < N; l++) */
	/* 	result2[l] = 0.; */
	dcomp i;
	i = -1.;
	i = sqrt(i);
	dcomp E = e + i*1e-5;
	/* double kT = send->kT; */
	/* for (int j=0; j!=15; j++){ */
	/* 	E = send->Ef + (2.*j + 1.)*kT*M_PI*i; */
		integrate = int_theta(E, send);
		/* if (abs(V) < 1e-9){ */
		/* 	for (int l = 0; l < N; l++) */
		/* 		result2[l] += 2.*kT*integrate[l]; */ 
		/* } */
		/* else { */
		/* 	for (int l = 0; l < N; l++) */
		/* 		result2[l] += kT*integrate[l]; */ 
		/* 	E = send->Ef - V + (2.*j + 1.)*kT*M_PI*i; */
		/* 	integrate = int_theta(E, send); */
		/* 	for (int l = 0; l < N; l++) */
		/* 		result2[l] += kT*integrate[l]; */ 
		/* } */
	/* } */
	/* vector<double> total; */
	/* total.reserve(N); */
	/* for (int l = 0; l < N; l++) */
	/* 	total[l] = result1[l] + result2[l]; */
	return integrate;
	/* return total; */
}

int main() 
{
	// plot output of spincurrent against energy
	const double Ef = 0.57553;//TODO
	// number of spacer layers
	int N = 10;
	// set bias
	double V = 0.00;
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
	MgO_lat_oop << 0.5, 0.5, 0;//TODO question mark as to whether the cosine vectors for T are correct e.g. -1.5, 0, 0.5 
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
	//TODO be aware that the code takes advantage of the fact that only second NN hoppings are different between
	//Mg and O... the code will break if this changes

	//This section generates the Hamiltonians from SK parameters and NN positions
	double x, y, z;
	Vector3d X, Y, Z;
	X << 1, 0, 0;
	Y << 0, 1, 0;
	Z << 0, 0, 1;
	vM iron_up, iron_dn, gold, iron_gold_up, iron_gold_dn, gold_MgO_11, gold_MgO_12, gold_MgO_21, gold_MgO_22, 
	   magnesium_11, magnesium_12, magnesium_21, oxide_11, oxide_12, oxide_21, iron_up_MgO_11, iron_up_MgO_12, 
	   iron_up_MgO_21, iron_up_MgO_22, iron_dn_MgO_11, iron_dn_MgO_12, iron_dn_MgO_21, iron_dn_MgO_22, gold_iron_up, gold_iron_dn;
	vec3 Au_pos, MgO_pos_11, MgO_pos_12, MgO_pos_21, Fe_pos, Au_MgO_pos_11, MgO_Fe_pos_11, Au_MgO_pos_12, MgO_Fe_pos_12, Au_MgO_pos_21, 
	     MgO_Fe_pos_21, Au_MgO_pos_22, MgO_Fe_pos_22, Au_Fe_pos, Fe_Au_pos;
	Au_pos.reserve(19); Fe_pos.reserve(19);
	Au_Fe_pos.reserve(19); Fe_Au_pos.reserve(19);
	iron_up.reserve(19); iron_dn.reserve(19); gold.reserve(19); iron_gold_up.reserve(19); iron_gold_dn.reserve(19);
	gold_iron_up.reserve(19); gold_iron_dn.reserve(19);
	//magic 19 above is num onsite + num nn + num nnn = 1 + 12 + 6
	Matrix<dcomp, 9, 9> tmp_mat;
	double distance;
	for (int i1 = -2; i1 < 3; i1++){
		for (int i2 = -2; i2 < 3; i2++){
			for (int i3 = -1; i3 < 2; i3++){
				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*MgO_lat_oop + MgO_basis[0];
				//TODO this needs sorting out.. per quadrant, like the hybrids below. Complicated due to U&T...
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
				//TODO this needs sorting out.. per quadrant, like the hybrids below. Complicated due to U&T...
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
					tmp_mat = eint1(Mg_O2, x, y, z);
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
				//TODO this needs sorting out.. per quadrant, like the hybrids below. Complicated due to U&T...
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
					tmp_mat = eint1(Mg_O2, x, y, z);
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

					tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Fe_Au - Au_basis[i4] + Fe_basis[1];
					/* if ((tmp_vec(1) < -0.44) && (tmp_vec(1) > -0.45)) */
					/* 	cout<<tmp_vec.transpose()<<endl; */
					distance = 0;
					for (int l = 0; l < 3; l++)
						distance += tmp_vec(l)*tmp_vec(l);
					distance = sqrt(distance);
					if (distance < 1e-5){
						Fe_Au_pos.emplace_back(tmp_vec);
						/* cout<<"There probably shouldn't be anything here: Fe-Au"<<endl; */
						tmp_mat = eint1(FeAu_u1, x, y, z);
						iron_gold_up.emplace_back(tmp_mat);
						tmp_mat = eint1(FeAu_d1, x, y, z);
						iron_gold_dn.emplace_back(tmp_mat);
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

					for (int i5 = 0; i5 < Au_basis.size(); i5++){
						tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Fe_Au + Au_basis[i4] - Fe_basis[i5];
						distance = 0;
						for (int l = 0; l < 3; l++)
							distance += tmp_vec(l)*tmp_vec(l);
						distance = sqrt(distance);
						if (distance < 1e-5){
							Fe_Au_pos.emplace_back(tmp_vec);
							/* cout<<"There probably shouldn't be anything here: Fe-Au"<<endl; */
							tmp_mat = eint1(FeAu_u1, x, y, z);
							iron_gold_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d1, x, y, z);
							iron_gold_dn.emplace_back(tmp_mat);
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
							Au_Fe_pos.emplace_back(tmp_vec);
							tmp_mat = eint1(FeAu_u1, x, y, z);
							gold_iron_up.emplace_back(tmp_mat);
							tmp_mat = eint1(FeAu_d1, x, y, z);
							gold_iron_dn.emplace_back(tmp_mat);
							/* cout<<"There probably shouldn't be anything here: Au-Fe"<<endl; */
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

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Au_MgO + MgO_basis[0] - Au_basis[0];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Au_MgO_pos_11.emplace_back(tmp_vec);
					/* cout<<"There probably shouldn't be anything here: Au-MgO"<<endl; */
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_11.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nn + 1e-3){
					Au_MgO_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_11.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnn + 1e-3){
					Au_MgO_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuMg2, x, y, z);
					gold_MgO_11.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnnn + 1e-3){
					Au_MgO_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO3, x, y, z);
					gold_MgO_11.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_MgO_Fe + Fe_basis[0] - MgO_basis[0];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_Fe_pos_11.emplace_back(tmp_vec);
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_11.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_11.emplace_back(tmp_mat);
					/* cout<<"There probably shouldn't be anything here: MgO-Fe"<<endl; */
				}
				else if (distance < MgO_Fe_nn + 1e-3){
					MgO_Fe_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_11.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_11.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnn + 1e-3){
					MgO_Fe_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u2, x, y, z);
					iron_up_MgO_11.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d2, x, y, z);
					iron_dn_MgO_11.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnnn + 1e-3){
					MgO_Fe_pos_11.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u3, x, y, z);
					iron_up_MgO_11.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d3, x, y, z);
					iron_dn_MgO_11.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Au_MgO + MgO_basis[1] - Au_basis[0];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Au_MgO_pos_12.emplace_back(tmp_vec);
					/* cout<<"There probably shouldn't be anything here: Au-MgO"<<endl; */
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_12.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nn + 1e-3){
					Au_MgO_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_12.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnn + 1e-3){
					Au_MgO_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuMg2, x, y, z);
					gold_MgO_12.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnnn + 1e-3){
					Au_MgO_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO3, x, y, z);
					gold_MgO_12.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_MgO_Fe + Fe_basis[1] - MgO_basis[0];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_Fe_pos_12.emplace_back(tmp_vec);
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_12.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_12.emplace_back(tmp_mat);
					/* cout<<"There probably shouldn't be anything here: MgO-Fe"<<endl; */
				}
				else if (distance < MgO_Fe_nn + 1e-3){
					MgO_Fe_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_12.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_12.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnn + 1e-3){
					MgO_Fe_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u2, x, y, z);
					iron_up_MgO_12.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d2, x, y, z);
					iron_dn_MgO_12.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnnn + 1e-3){
					MgO_Fe_pos_12.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u3, x, y, z);
					iron_up_MgO_12.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d3, x, y, z);
					iron_dn_MgO_12.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Au_MgO + MgO_basis[0] - Au_basis[1];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Au_MgO_pos_21.emplace_back(tmp_vec);
					/* cout<<"There probably shouldn't be anything here: Au-MgO"<<endl; */
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_21.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nn + 1e-3){
					Au_MgO_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_21.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnn + 1e-3){
					Au_MgO_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuMg2, x, y, z);
					gold_MgO_21.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnnn + 1e-3){
					Au_MgO_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO3, x, y, z);
					gold_MgO_21.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_MgO_Fe + Fe_basis[0] - MgO_basis[1];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_Fe_pos_21.emplace_back(tmp_vec);
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_21.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_21.emplace_back(tmp_mat);
					/* cout<<"There probably shouldn't be anything here: MgO-Fe"<<endl; */
				}
				else if (distance < MgO_Fe_nn + 1e-3){
					MgO_Fe_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_21.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_21.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnn + 1e-3){
					MgO_Fe_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u2, x, y, z);
					iron_up_MgO_21.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d2, x, y, z);
					iron_dn_MgO_21.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnnn + 1e-3){
					MgO_Fe_pos_21.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u3, x, y, z);
					iron_up_MgO_21.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d3, x, y, z);
					iron_dn_MgO_21.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_Au_MgO + MgO_basis[1] - Au_basis[1];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					Au_MgO_pos_22.emplace_back(tmp_vec);
					/* cout<<"There probably shouldn't be anything here: Au-MgO"<<endl; */
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_22.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nn + 1e-3){
					Au_MgO_pos_22.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO1, x, y, z);
					gold_MgO_22.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnn + 1e-3){
					Au_MgO_pos_22.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuMg2, x, y, z);
					gold_MgO_22.emplace_back(tmp_mat);
				}
				else if (distance < Au_MgO_nnnn + 1e-3){
					Au_MgO_pos_22.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(AuO3, x, y, z);
					gold_MgO_22.emplace_back(tmp_mat);
				}

				tmp_vec = i1*lat_vec1 + i2*lat_vec2 + i3*lat_MgO_Fe + Fe_basis[1] - MgO_basis[1];
				distance = 0;
				for (int l = 0; l < 3; l++)
					distance += tmp_vec(l)*tmp_vec(l);
				distance = sqrt(distance);
				if (distance < 1e-5){
					MgO_Fe_pos_22.emplace_back(tmp_vec);
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_22.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_22.emplace_back(tmp_mat);
					/* cout<<"There probably shouldn't be anything here: MgO-Fe"<<endl; */
				}
				else if (distance < MgO_Fe_nn + 1e-3){
					MgO_Fe_pos_22.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u1, x, y, z);
					iron_up_MgO_22.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d1, x, y, z);
					iron_dn_MgO_22.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnn + 1e-3){
					MgO_Fe_pos_22.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u2, x, y, z);
					iron_up_MgO_22.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d2, x, y, z);
					iron_dn_MgO_22.emplace_back(tmp_mat);
				}
				else if (distance < MgO_Fe_nnnn + 1e-3){
					MgO_Fe_pos_22.emplace_back(tmp_vec);
					x = tmp_vec.dot(X)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					y = tmp_vec.dot(Y)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					z = tmp_vec.dot(Z)/sqrt(tmp_vec(0)*tmp_vec(0) + tmp_vec(1)*tmp_vec(1) + tmp_vec(2)*tmp_vec(2)); 
					tmp_mat = eint1(MgFe_u3, x, y, z);
					iron_up_MgO_22.emplace_back(tmp_mat);
					tmp_mat = eint1(MgFe_d3, x, y, z);
					iron_dn_MgO_22.emplace_back(tmp_mat);
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

	double kT = k*T;
	//set up the variables to send
	variables send;
	send.kT = kT;
	send.Ef = Ef;
	send.x = 2.532374; //for now!
	send.z = 2.532374; //for now!
	send.Au_lat_oop = &Au_lat_oop;
	send.MgO_lat_oop = &MgO_lat_oop;
	send.Fe_lat_oop = &Fe_lat_oop;
	send.lat_Au_MgO = &lat_Au_MgO;
	send.lat_MgO_Fe = &lat_MgO_Fe;
	send.lat_Fe_Au = &lat_Fe_Au;
	send.lat_Au_Fe = &lat_Au_Fe;
	send.Au_pos = &Au_pos;
	send.MgO_pos_11 = &MgO_pos_11;
	send.MgO_pos_12 = &MgO_pos_12;
	send.MgO_pos_21 = &MgO_pos_21;
	send.Fe_pos = &Fe_pos;
	send.Au_MgO_pos_11 = &Au_MgO_pos_11;
	send.Au_MgO_pos_12 = &Au_MgO_pos_12;
	send.Au_MgO_pos_21 = &Au_MgO_pos_21;
	send.Au_MgO_pos_22 = &Au_MgO_pos_22;
	send.MgO_Fe_pos_11 = &MgO_Fe_pos_11;
	send.MgO_Fe_pos_12 = &MgO_Fe_pos_12;
	send.MgO_Fe_pos_21 = &MgO_Fe_pos_21;
	send.MgO_Fe_pos_22 = &MgO_Fe_pos_22;
	send.Fe_Au_pos = &Fe_Au_pos;
	send.Au_Fe_pos = &Au_Fe_pos;
	send.Au_basis = &Au_basis;
	send.MgO_basis = &MgO_basis;
	send.Fe_basis = &Fe_basis;
	send.gold = &gold;
	send.iron_up = &iron_up;
	send.iron_dn = &iron_dn;
	send.magnesium_11 = &magnesium_11;
	send.magnesium_12 = &magnesium_12;
	send.magnesium_21 = &magnesium_21;
	send.oxide_11 = &oxide_11;
	send.oxide_21 = &oxide_21;
	send.oxide_12 = &oxide_12;
	send.gold_iron_up = &gold_iron_up;
	send.gold_iron_dn = &gold_iron_dn;
	send.iron_gold_up = &iron_gold_up;
	send.iron_gold_dn = &iron_gold_dn;
	send.iron_dn_MgO_11 = &iron_dn_MgO_11;
	send.iron_dn_MgO_12 = &iron_dn_MgO_12;
	send.iron_dn_MgO_21 = &iron_dn_MgO_21;
	send.iron_dn_MgO_22 = &iron_dn_MgO_22;
	send.iron_up_MgO_11 = &iron_up_MgO_11;
	send.iron_up_MgO_12 = &iron_up_MgO_12;
	send.iron_up_MgO_21 = &iron_up_MgO_21;
	send.iron_up_MgO_22 = &iron_up_MgO_22;
	send.gold_MgO_11 = &gold_MgO_11;
	send.gold_MgO_12 = &gold_MgO_12;
	send.gold_MgO_21 = &gold_MgO_21;
	send.gold_MgO_22 = &gold_MgO_22;
	send.V = V;
	send.lim = lim;
	send.lim2 = lim2;
	send.N = N;

	double answer;

	time_t now = time(0);
	tm *ltm = localtime(&now);
	string Mydata;
	Mydata = to_string(ltm->tm_mday);
	Mydata += "-";
	Mydata += to_string(1+ltm->tm_mon);
	Mydata += "-";
	Mydata += to_string(1900+ltm->tm_year);

	ofstream Myfile;	

	Mydata += "-Au_fermi3.txt";
	/* answer = switching(&send); */

	Myfile.open( Mydata.c_str(),ios::trunc );
	double result;
	double integrate;
	/* int n = 350;//set the number of k-points along x axis */
	int n = 1000;//set the number of k-points along x axis
	int counter2 = 0;
	double factor = 2./(n*n);
	int sumk = n*(n + 1)/2;
	int p, q;
	int product1;
	int myid, numprocs;
	time_t timer;
	int kk, l, i;
	myid = 0;
	numprocs = 1;

	double j = 0.;
	for (kk = 2*myid + 1; kk<2*sumk; kk+=2*numprocs){
		for (l = 0; l < n; l++){
			product1 = (2*n - l)*(l + 1);
			if ( kk < product1){
				p = 2*(kk/(2*n - l)) + 1;
				q = (l*(l + 1) + kk)%(2*n);
				break;
			}
		}
		send.x = M_PI*(p + q)/(2*n);//this translates the grid for fcc 1st BZ
		send.z = M_PI*(p - q)/(2*n);//this translates the grid for fcc 1st BZ
		integrate = switching(j, &send);
		if ((p==1) && (q==1))
			result = factor*0.5*integrate;
		else if (p==q)
			result = factor*0.5*integrate;
		else
			result = factor*integrate;
		counter2++;
		Myfile<<scientific<<send.x<<" "<<send.z<<" "<<result<<endl;
	}

	return 0;
}
