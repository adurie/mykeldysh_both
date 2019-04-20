#include <iostream>
#include <cmath>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include "TBdynamic.h"
/* #include "cunningham_diamond.h"//works best with error 5e-3 */
#include "cunningham_diamond_abs.h"//works best with error 1e-1
/* #include <omp.h> */
//This program calculates the realistic exchange coupling in a Co/Cu/Co(001)
//trilayer. It does so for fcc Co and fcc Cu. Interatomic spacing is considered
//identical and the interfaces are abrupt.
//use adlayer for integer thickness n, mobius method for surfaces.
//
//Use abs error as this is more efficient when summing residues

using namespace std;
using namespace Eigen;
typedef complex<double> dcomp;
typedef Matrix<dcomp, 9, 9> dmat;
typedef Matrix<dcomp, 18, 18> ddmat;
typedef Matrix<dcomp, 36, 36> dddmat;
typedef Vector3d vec;

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

VectorXcd greens(double k_x, double k_z, double a, dcomp omega, int N, dmat &u,
		dmat &t_1, dmat &t_2, dmat &t_3, dmat &t_4, dmat &t_5, dmat &t_6, dmat &t_7, 
		dmat &t_8, dmat &t_9, dmat &t_10, dmat &t_11, dmat &t_12, dmat &t_13,
	  	dmat &t_14, dmat &t_15, dmat &t_16, dmat &t_17, dmat &t_18, dmat &u_u, 
		dmat &tu_1, dmat &tu_2, dmat &tu_3, dmat &tu_4, dmat &tu_5, dmat &tu_6,
	       	dmat &tu_7, dmat &tu_8, dmat &tu_9, dmat &tu_10, dmat &tu_11, dmat &tu_12, dmat &tu_13,
		dmat &tu_14, dmat &tu_15, dmat &tu_16, dmat &tu_17, dmat &tu_18, dmat &u_d,
		dmat &td_1, dmat &td_2, dmat &td_3, dmat &td_4, dmat &td_5, dmat &td_6,
	       	dmat &td_7, dmat &td_8, dmat &td_9, dmat &td_10, dmat &td_11, dmat &td_12, dmat &td_13,
		dmat &td_14, dmat &td_15, dmat &td_16, dmat &td_17, dmat &td_18, 
		vec &d_3, vec &d_4,
	       	vec &d_9, vec &d_10, vec &d_13, vec &d_14,
	       	vec &d_17, vec &d_18){

	dcomp i;
	i = -1.;
	i = sqrt(i);
	double k_y = 0;

	Vector3d K;
	K(0) = k_x;
	K(1) = k_y;
	K(2) = k_z;

	//construct diagonalised in-plane matrices
	//Cu fcc
	Matrix<dcomp, 9, 9> u_11, u_12, u_21, T_21;
	u_11 = u + t_3*exp(i*d_3.dot(K))+ t_4*exp(i*d_4.dot(K))+ t_9*exp(i*d_9.dot(K)) + t_10*exp(i*d_10.dot(K)) + 
		t_13*exp(i*d_13.dot(K))+ t_14*exp(i*d_14.dot(K))+ t_17*exp(i*d_17.dot(K)) + t_18*exp(i*d_18.dot(K));
	u_12 = t_1 + t_5*exp(i*d_9.dot(K)) + t_7*exp(i*d_14.dot(K)) + t_12*exp(i*d_4.dot(K));
	u_21 = t_2 + t_8*exp(i*d_13.dot(K)) + t_11*exp(i*d_3.dot(K)) + t_6*exp(i*d_10.dot(K));
	Matrix<dcomp, 18, 18> U, T, OM;
	U << u_11, u_12, u_21, u_11;
	Matrix<complex<double>, 9, 9> zero = Matrix<complex<double>, 9, 9>::Zero();
	T_21 = t_7 + t_1*exp(i*d_13.dot(K)) + t_5*exp(i*d_3.dot(K)) + t_12*exp(i*d_10.dot(K));
	T << t_15, zero, T_21, t_15;
	/* cout<<U.real()<<endl<<endl; */
	/* cout<<T.real()<<endl; */

	//Co spin up fcc
	Matrix<dcomp, 9, 9> uu_11, uu_12, uu_21, Tu_21;
	uu_11 = u_u + tu_3*exp(i*d_3.dot(K))+ tu_4*exp(i*d_4.dot(K))+ tu_9*exp(i*d_9.dot(K)) + tu_10*exp(i*d_10.dot(K)) + 
		tu_13*exp(i*d_13.dot(K))+ tu_14*exp(i*d_14.dot(K))+ tu_17*exp(i*d_17.dot(K)) + tu_18*exp(i*d_18.dot(K));
	uu_12 = tu_1 + tu_5*exp(i*d_9.dot(K)) + tu_7*exp(i*d_14.dot(K)) + tu_12*exp(i*d_4.dot(K));
	uu_21 = tu_2 + tu_8*exp(i*d_13.dot(K)) + tu_11*exp(i*d_3.dot(K)) + tu_6*exp(i*d_10.dot(K));
	Matrix<dcomp, 18, 18> Uu, Tu, OMu, GLu_even, GRu, Tu_odd_init;
	Uu << uu_11, uu_12, uu_21, uu_11;
	Tu_21 = tu_7 + tu_1*exp(i*d_13.dot(K)) + tu_5*exp(i*d_3.dot(K)) + tu_12*exp(i*d_10.dot(K));
	Tu << tu_15, zero, Tu_21, tu_15;
	Tu_odd_init << tu_15, zero, Tu_21, t_15;

	//Co spin up fcc
	Matrix<dcomp, 9, 9> ud_11, ud_12, ud_21, Td_21;
	ud_11 = u_d + td_3*exp(i*d_3.dot(K))+ td_4*exp(i*d_4.dot(K))+ td_9*exp(i*d_9.dot(K)) + td_10*exp(i*d_10.dot(K)) + 
		td_13*exp(i*d_13.dot(K))+ td_14*exp(i*d_14.dot(K))+ td_17*exp(i*d_17.dot(K)) + td_18*exp(i*d_18.dot(K));
	ud_12 = td_1 + td_5*exp(i*d_9.dot(K)) + td_7*exp(i*d_14.dot(K)) + td_12*exp(i*d_4.dot(K));
	ud_21 = td_2 + td_8*exp(i*d_13.dot(K)) + td_11*exp(i*d_3.dot(K)) + td_6*exp(i*d_10.dot(K));
	Matrix<dcomp, 18, 18> Ud, Td, OMd, GLd_even, GRd, Td_odd_init;
	Ud << ud_11, ud_12, ud_21, ud_11;
	Td_21 = td_7 + td_1*exp(i*d_13.dot(K)) + td_5*exp(i*d_3.dot(K)) + td_12*exp(i*d_10.dot(K));
	Td << td_15, zero, Td_21, td_15;
	Td_odd_init << td_15, zero, Td_21, t_15;

	//Co|Cu bilayer spin up (required for odd planes)
	ddmat Ucocu_u;
	Ucocu_u << uu_11, u_12, u_21, u_11;

	//Co|Cu bilayer spin down (required for odd planes)
	ddmat Ucocu_d;
	Ucocu_d << ud_11, u_12, u_21, u_11;

      	Matrix<complex<double>, 18, 18> I = Matrix<complex<double>, 18, 18>::Identity();
	ddmat Tudagg, Tddagg, Tdagg;
	Tudagg = Tu.adjoint();
	Tddagg = Td.adjoint();
	Tdagg = T.adjoint();

	OM = omega*I-U;
	OMu = omega*I-Uu;
	OMd = omega*I-Ud;

	ddmat OMu_odd_init, OMd_odd_init;
	OMu_odd_init = omega*I-Ucocu_u;
	OMd_odd_init = omega*I-Ucocu_d;

	GLu_even = gs(OMu, Tu);
	GLd_even = gs(OMd, Td);
	GRu = gs(OMu, Tudagg);
	GRd = gs(OMd, Tddagg);

	ddmat GLu_odd, GLd_odd;
	
	ddmat Rsigma_0_u, Rsigma_0_d, Rsigma_PI_u, Rsigma_PI_d;
	dcomp Fsigma;
	VectorXcd result(N);
	result.fill(0.);

	ddmat GNu_odd, GNu_even, GNd_odd, GNd_even;

	Rsigma_0_u = (I-GRu*Tdagg*GLu_even*T);
	Rsigma_0_d = (I-GRd*Tdagg*GLd_even*T);
	Rsigma_PI_u = (I-GRd*Tdagg*GLu_even*T);
	Rsigma_PI_d = (I-GRu*Tdagg*GLd_even*T);

	Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());
	result[0] = Fsigma;

	GLu_odd = (OMu_odd_init - Tu_odd_init.adjoint()*GLu_even*Tu_odd_init).inverse();
	GLd_odd = (OMd_odd_init - Td_odd_init.adjoint()*GLd_even*Td_odd_init).inverse();

	Rsigma_0_u = (I-GRu*Tdagg*GLu_odd*T);
	Rsigma_0_d = (I-GRd*Tdagg*GLd_odd*T);
	Rsigma_PI_u = (I-GRd*Tdagg*GLu_odd*T);
	Rsigma_PI_d = (I-GRu*Tdagg*GLd_odd*T);

	Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());
	result[1] = Fsigma;

//adlayer layer 2 from layer 1 to spacer thickness, N
	for (int it=2; it < N; ++it){
		if (it%2 == 0){
			GLu_even = (OM - Tdagg*GLu_even*T).inverse();
			GLd_even = (OM - Tdagg*GLd_even*T).inverse();
			Rsigma_0_u = (I-GRu*Tdagg*GLu_even*T);
			Rsigma_0_d = (I-GRd*Tdagg*GLd_even*T);
			Rsigma_PI_u = (I-GRd*Tdagg*GLu_even*T);
			Rsigma_PI_d = (I-GRu*Tdagg*GLd_even*T);
		}
		if (it%2 == 1){
			GLu_odd = (OM -Tdagg*GLu_odd*T).inverse();
			GLd_odd = (OM -Tdagg*GLd_odd*T).inverse();
			Rsigma_0_u = (I-GRu*Tdagg*GLu_odd*T);
			Rsigma_0_d = (I-GRd*Tdagg*GLd_odd*T);
			Rsigma_PI_u = (I-GRd*Tdagg*GLu_odd*T);
			Rsigma_PI_d = (I-GRu*Tdagg*GLd_odd*T);

		}
		Fsigma = (1./M_PI)*log((Rsigma_0_d*Rsigma_0_u*Rsigma_PI_u.inverse()*Rsigma_PI_d.inverse()).determinant());
		result[it] = Fsigma;
	}

	return result;
}

int main(){

	/* cout<<"Name the data file\n"; */
	/* string Mydata; */
	/* getline(cin, Mydata); */
	/* ofstream Myfile; */	
	/* Mydata += ".txt"; */
	/* Myfile.open( Mydata.c_str(),ios::trunc ); */

	Vector3d d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8, d_9;
	Vector3d d_10, d_11, d_12, d_13, d_14, d_15, d_16;
	Vector3d d_17, d_18;
	
	double a = 1.;
	
	//position vectors of nearest neighbours in fcc
	d_1 << a/2., a/2., 0;
	d_2 << -a/2., -a/2., 0;
	d_3 << a/2., 0, a/2.;
	d_4 << -a/2., 0, -a/2.;
	d_5 << 0, a/2., a/2.;
	d_6 << 0, -a/2., -a/2.;
	d_7 << -a/2., a/2., 0;
	d_8 << a/2., -a/2., 0;
	d_9 << -a/2., 0, a/2.;
	d_10 << a/2., 0, -a/2.;
	d_11 << 0, -a/2., a/2.;
	d_12 << 0, a/2., -a/2.;

	//position vectors of next nearest neighbours
	d_13 << a, 0, 0;
	d_14 << -a, 0, 0;
	d_15 << 0, a, 0;
	d_16 << 0, -a, 0;
	d_17 << 0, 0, a;
	d_18 << 0, 0, -a;
	//initialise onsite for fcc Cu
	Matrix<dcomp, 9, 9> u;
	u = TB(2, 0, 0, 9, d_1);

	//initialise nn hopping for fcc Cu
	Matrix<dcomp, 9, 9> t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9;
	Matrix<dcomp, 9, 9> t_10, t_11, t_12;
	t_1 = TB(2, 1, 0, 9, d_1);
	t_2 = TB(2, 1, 0, 9, d_2);
	t_3 = TB(2, 1, 0, 9, d_3);
	t_4 = TB(2, 1, 0, 9, d_4);
	t_5 = TB(2, 1, 0, 9, d_5);
	t_6 = TB(2, 1, 0, 9, d_6);
	t_7 = TB(2, 1, 0, 9, d_7);
	t_8 = TB(2, 1, 0, 9, d_8);
	t_9 = TB(2, 1, 0, 9, d_9);
	t_10 = TB(2, 1, 0, 9, d_10);
	t_11 = TB(2, 1, 0, 9, d_11);
	t_12 = TB(2, 1, 0, 9, d_12);

	//initialise next nn hopping for fcc Cu
	Matrix<dcomp, 9, 9> t_13, t_14, t_15, t_16, t_17, t_18;
	t_13 = TB(2, 1, 1, 9, d_13);
	t_14 = TB(2, 1, 1, 9, d_14);
	t_15 = TB(2, 1, 1, 9, d_15);
	t_16 = TB(2, 1, 1, 9, d_16);
	t_17 = TB(2, 1, 1, 9, d_17);
	t_18 = TB(2, 1, 1, 9, d_18);

	//initialise onsite for fcc Co spin up
	Matrix<dcomp, 9, 9> u_u;
	u_u = TB(0, 0, 0, 9, d_1);

	//initialise nn hopping for fcc Co spin up
	Matrix<dcomp, 9, 9> tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8;
	Matrix<dcomp, 9, 9> tu_9, tu_10, tu_11, tu_12;
	tu_1 = TB(0, 1, 0, 9, d_1);
	tu_2 = TB(0, 1, 0, 9, d_2);
	tu_3 = TB(0, 1, 0, 9, d_3);
	tu_4 = TB(0, 1, 0, 9, d_4);
	tu_5 = TB(0, 1, 0, 9, d_5);
	tu_6 = TB(0, 1, 0, 9, d_6);
	tu_7 = TB(0, 1, 0, 9, d_7);
	tu_8 = TB(0, 1, 0, 9, d_8);
	tu_9 = TB(0, 1, 0, 9, d_9);
	tu_10 = TB(0, 1, 0, 9, d_10);
	tu_11 = TB(0, 1, 0, 9, d_11);
	tu_12 = TB(0, 1, 0, 9, d_12);

	//initialise next nn hopping for fcc Co spin up
	Matrix<dcomp, 9, 9> tu_13, tu_14, tu_15, tu_16, tu_17, tu_18;
	tu_13 = TB(0, 1, 1, 9, d_13);
	tu_14 = TB(0, 1, 1, 9, d_14);
	tu_15 = TB(0, 1, 1, 9, d_15);
	tu_16 = TB(0, 1, 1, 9, d_16);
	tu_17 = TB(0, 1, 1, 9, d_17);
	tu_18 = TB(0, 1, 1, 9, d_18);

	//initialise onsite for fcc Co spin down 
	Matrix<dcomp, 9, 9> u_d;
	u_d = TB(1, 0, 0, 9, d_1);

	//initialise nn hopping for fcc Co spin down
	Matrix<dcomp, 9, 9> td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8;
	Matrix<dcomp, 9, 9> td_9, td_10, td_11, td_12;
	td_1 = TB(1, 1, 0, 9, d_1);
	td_2 = TB(1, 1, 0, 9, d_2);
	td_3 = TB(1, 1, 0, 9, d_3);
	td_4 = TB(1, 1, 0, 9, d_4);
	td_5 = TB(1, 1, 0, 9, d_5);
	td_6 = TB(1, 1, 0, 9, d_6);
	td_7 = TB(1, 1, 0, 9, d_7);
	td_8 = TB(1, 1, 0, 9, d_8);
	td_9 = TB(1, 1, 0, 9, d_9);
	td_10 = TB(1, 1, 0, 9, d_10);
	td_11 = TB(1, 1, 0, 9, d_11);
	td_12 = TB(1, 1, 0, 9, d_12);

	//initialise next nn hopping for fcc Co spin down
	Matrix<dcomp, 9, 9> td_13, td_14, td_15, td_16, td_17, td_18;
	td_13 = TB(1, 1, 1, 9, d_13);
	td_14 = TB(1, 1, 1, 9, d_14);
	td_15 = TB(1, 1, 1, 9, d_15);
	td_16 = TB(1, 1, 1, 9, d_16);
	td_17 = TB(1, 1, 1, 9, d_17);
	td_18 = TB(1, 1, 1, 9, d_18);

	dcomp i;
	i = -1.;
	i = sqrt(i);

	//number of principle layers of spacer
	/* const int N = 50; */
	const int N = 30;
	dcomp E = 0.;
	const double Ef = 0.57553;
	const double kT = 8.617342857e-5*315.79/13.6058;
	//multiple vectors created for parallelisation
	VectorXcd result_complex(N);
	result_complex.fill(0.);

	VectorXcd result_complex_t1(N);
	result_complex_t1.fill(0.);

	VectorXcd result_complex_t2(N);
	result_complex_t2.fill(0.);

	VectorXcd result_complex_t3(N);
	result_complex_t3.fill(0.);

	double abs_error = 1e-1;
	//see integration header for what the following zero values are interpreted as
	int k_start = 0;
	int k_max = 0;
//now begin setup of parallelisation
	/* omp_set_num_threads(1); */
	/* int tid, x; */
	/* #pragma omp parallel private(E) */
	/* { */
	/* tid = omp_get_thread_num(); */
	/* if (tid == 0){ */
	/* 	x = omp_get_max_threads(); */

	/* 	if (x<5) */
	/* 		cout<<"there are "<<x<<" thread(s) working on this task"<<endl; */

	/* 	if (x>4) */
	/* 		cout<<"warning, the program can be made more efficient as num_threads > 4"<<endl; */
	/* } */
	/* for (int j=0; j!=10; j++){ */
		/* if ((j%x == omp_get_thread_num()) && (omp_get_thread_num() == 0)){ */
	int j = 0;
			E = Ef + (2.*j + 1.)*kT*M_PI*i;
			result_complex = greens(0.8, 2.2, a, E, N,
				u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9,
				t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18,
				u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9,
			       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18,
				u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9,
			       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18,
				d_3, d_4, d_9, d_10,
			       	d_13, d_14, d_17, d_18);
			/* cout<<j<<endl; */
		/* } */
		/* if ((x > 1) && (x < 5)){ */
		/* 	if ((j%x == omp_get_thread_num()) && (omp_get_thread_num() == 1)){ */
		/* 		E = Ef + (2.*j + 1.)*kT*M_PI*i; */
		/* 		result_complex_t1 = result_complex_t1 + kspace(&greens, k_start, abs_error, k_max, a, E, N, */
		/* 			u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
		/* 			t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, */
		/* 			u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9, */
		/* 		       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18, */
		/* 			u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9, */
		/* 		       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18, */
		/* 			d_3, d_4, d_9, d_10, */
		/* 		       	d_13, d_14, d_17, d_18); */
		/* 		cout<<j<<endl; */
		/* 	} */

		/* 	if ((j%x == omp_get_thread_num()) && (omp_get_thread_num() == 2)){ */
		/* 		E = Ef + (2.*j + 1.)*kT*M_PI*i; */
		/* 		result_complex_t2 = result_complex_t2 + kspace(&greens, k_start, abs_error, k_max, a, E, N, */
		/* 			u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
		/* 			t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, */
		/* 			u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9, */
		/* 		       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18, */
		/* 			u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9, */
		/* 		       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18, */
		/* 			d_3, d_4, d_9, d_10, */
		/* 		       	d_13, d_14, d_17, d_18); */
		/* 		cout<<j<<endl; */
		/* 	} */

		/* 	if ((j%x == omp_get_thread_num()) && (omp_get_thread_num() == 3)){ */
		/* 		E = Ef + (2.*j + 1.)*kT*M_PI*i; */
		/* 		result_complex_t3 = result_complex_t3 + kspace(&greens, k_start, abs_error, k_max, a, E, N, */
		/* 			u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
		/* 			t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, */
		/* 			u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9, */
		/* 		       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18, */
		/* 			u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9, */
		/* 		       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18, */
		/* 			d_3, d_4, d_9, d_10, */
		/* 		       	d_13, d_14, d_17, d_18); */
		/* 		cout<<j<<endl; */
		/* 	} */
		/* } */
	/* } */
	/* } */
	/* result_complex = result_complex + result_complex_t1 + result_complex_t2 + result_complex_t3; */

	/* VectorXd result = result_complex.real(); */
	/* result *= kT/(4.*M_PI*M_PI); */

	/* E = Ef + kT*M_PI*i; */
	/* result_complex = greens(-2.19911485751286, -0.314159265358979, a, E, N, */
	/* 		u, t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, */
	/* 		t_10, t_11, t_12, t_13, t_14, t_15, t_16, t_17, t_18, */
	/* 		u_u, tu_1, tu_2, tu_3, tu_4, tu_5, tu_6, tu_7, tu_8, tu_9, */
	/* 	       	tu_10, tu_11, tu_12, tu_13, tu_14, tu_15, tu_16, tu_17, tu_18, */
	/* 		u_d, td_1, td_2, td_3, td_4, td_5, td_6, td_7, td_8, td_9, */
	/* 	       	td_10, td_11, td_12, td_13, td_14, td_15, td_16, td_17, td_18, */
	/* 		d_3, d_4, d_9, d_10, */
	/* 	       	d_13, d_14, d_17, d_18); */
	/* VectorXd result = result_complex.real(); */

	/* Myfile<<"N Gamma"<<endl; */

	/* for (int ii=0; ii < N ; ++ii){ */
	/* 	Myfile << ii <<" ,  "<< 2.*M_PI*result[ii] << endl; */
	/* } */

	/* cout<<"finished!"<<endl; */

	/* 		Myfile.close(); */
	return 0;
}
