#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/src/Core/util/MKL_support.h>
#include <vector>
#define EIGEN_USE_MKL_ALL

using namespace Eigen;
using namespace std;
typedef complex<double> dcomp;

dcomp fermi(dcomp arg, double Ef){
	const double k = 8.617e-5/13.6058;
	const double T = 300;
	double kT = k*T;
	return 1./(1.+exp((arg-Ef)/kT));
}

vector<double> f(double x, double z, double a, dcomp E, double Ef, int N, double theta, int myswitch) {
// ...NM|ins|FM(0)|NM(n)|FM(theta)...
	dcomp i = -1;
	i = sqrt(i);
	double F = cos(x*a)+cos(z*a);

	/* const double V = 0.0; */
	const double V = 0.3;
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

	Matrix2cd OMV2=E*I-FM2-2.*T*F;

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

	Matrix2cd OM = E*I - 2.*T*F;

	Matrix2cd OMV1=E*I-NM1-2.*T*F;
	/* Matrix2cd OMV1=E*I-FM1-2.*T*F; */

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
		B = -(I-GR_dagg*T_dagg*GL.adjoint()*T).inverse();
		if (myswitch == 0){
			TOT = (GL*T*A*B*GR_dagg*T_dagg-A*B+0.5*(A+B))*Pauli;
			spincurrent = (1./(2.*M_PI))*real(TOT.trace()*(fermi(E,Ef)-fermi(E,Ef-V)));
		}
		if (myswitch == 1){
			TOT = (B-A)*Pauli;
			spincurrent = imag(TOT.trace())*(1. + V*V/4.);
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

vector<double> int_theta(double x, double z, double a, dcomp E, double Ef, int N, int myswitch){
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
		integrate = f(x, z, a, E, Ef, N, theta, myswitch);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += M_PI*(0.5/n)*integrate[i];
			else 
				result[i] += (M_PI/n)*integrate[i];
		}
	}	
	return result;
}

/* double pass(double E, void * params) */

vector<double> int_energy(double x, double z, double a, double Ef, int N){
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
	double end = 0.4;
	double start = -1.9;

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

	const int n = 20000;
	double factor = (end - start)/(n*1.);
	for (int k=0; k<n+1; k++) {
		E = start + k*(end-start)/n;
		E_send = E + 1e-6*im;
		integrate = int_theta(x, z, a, E_send, Ef, N, 0);
		for (int i = 0; i < N; i++){
			if ((k==0)||(k==n))
				result[i] += 0.5*factor*integrate[i];
			else 
				result[i] += factor*integrate[i];
		}
	}	

	return result;
}

vector<double> switching(double x, double z, double a, double Ef, int N){
	vector<double> result1, result2, integrate;
	result1.reserve(N);
	result1 = int_energy(x, z, a, Ef, N);
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
		integrate = int_theta(x, z, a, E, Ef, N, 1);
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

vector<double> int_kpoints(double a, double Ef, int N){
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
					integrate = switching(x, z, 1, Ef, N);
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
	Myfile.open( Mydata.c_str(),ios::trunc );
	const double Ef = 0.0;
	/* const double Ef = -0.3; */

	// number of spacer layers
	int N = 21;
	vector<double> answer;
	answer.reserve(N);
	/* answer = int_theta(0, 0, 1,  0.1, Ef, N); */
	answer = switching(0, 0, 1, Ef, N);
	/* answer = int_energy(0, 0, 1, Ef, N); */
	/* answer = int_kpoints(1, Ef, N); */
	/* answer = f(0, 0, 1, Ef, Ef, i, 0); */
	for (int i = 1; i < N; i++){
		Myfile<<scientific<<i<<" "<<answer[i]<<endl;
	}
	return 0;
}
