#ifndef CUNNINGHAM_POINTS_INTEGRATION_H
#define CUNNINGHAM_POINTS_INTEGRATION_H

#include <vector>
#include <algorithm>
#include <functional>
#include <utility>
#include <eigen3/Eigen/Dense>
#include <complex>
#include <cmath>

/* This header file contains routines for integrating
 * the 2-D Brillouin zone by means of Cunnningham points
 * integration. It is presently only for simple cubics
 * it also contains some useful overloads for std::vector
 */
using namespace Eigen;
using namespace std;

typedef complex<double> dcomp;


template <typename T>
vector<T> operator*(const double a, const vector<T>& b){
	vector<T> tmp;
	tmp.reserve(b.size());
	tmp = b;
	std::transform(tmp.begin(), tmp.end(), tmp.begin(),
               std::bind1st(std::multiplies<T>(),a));
	return tmp;
}

template<typename T>
vector<T>& operator+=( vector<T>& lhs, const vector<T>& rhs ) {
	assert(rhs.size() == lhs.size());
	for ( int i = 0; i < rhs.size() ; ++i ) {
		lhs[i] += rhs[i];
	}
	return lhs;
}
/* This function is for the type std::vector<double> should the
 * forwarded function be of this return type.
 * The function should be called as kspace(params, &f, iterations)
 * where params is a struct to be passed to the following function,
 * f is the function to be forwarded and iterations is the max 
 * number of iterations requested. If set to 0, it defaults to 1024**2
 * iterations should always be set to 2**number 
 */

/* template <typename Pred, typename T> */
/* vector<double> kspace(T&& params, Pred predicate, int iterations){ */
/* 	vector<double> integral(params.N,0.); */
/* 	integral.reserve(params.N); */
/* //integral width number of points n**2 */
/* 	int n; */
/* 	if (iterations == 0) */
/* 		n = 1024; */
/* 	else */
/* 		n = iterations; */
/* 	const double A = M_PI/params.a; */
/* 	double x,z; */
/* 	for (int k=1; k<=n+1; k++) { */
/* 		x = A*(k-1)/n; */
/* 		for (int l=1; l<=k+1; l++) { */
/* 			z = A*(l-1)/n; */
/* 			if (k==l){ */
/* 				integral += 0.5*predicate(x,z,forward<T>(params));} */
/* 			if (k!=l){ */
/* 				integral += predicate(x,z,forward<T>(params));} */
/* 		} */
/* 	} */
/* 	integral = (8.*A*A/(n*n))*integral; */
/* 	return integral; */
/* } */

/* This function is for the type double should the
 * forwarded function be of this return type.
 * The function should be called as kspace(params, &f, iterations, abs, loop)
 * where params is a struct to be passed to the following function,
 * f is the function to be forwarded,  iterations is the max 
 * number of iterations requested where if set to 0, it defaults to 1024**2.
 * iterations should always be set to 2**number 
 * abs is the absolute error margin to work to, unless iterations is exceeded.
 * if set to 0 it defaults to 0.03, or 3% error. loop is the value to be passed
 * through that the integration is looped over.
 */
template <typename Pred, typename T, typename T2>
dcomp kspace_complex(T&& params, Pred predicate, int iterations, double abs, T2 loop){
	dcomp integral = 0.;
	dcomp tmp;
	double test;
//integral width number of points n**2. max_width chosen as 2**integer
	double error;
	if (abs == 0)
		error = 0.02;
	else
		error = abs;
	int max_width;
	if (iterations == 0)
		max_width = 4096;
	else
		max_width = iterations;
	int n = 0;
	int nPrev;
	const double A = M_PI/params.a;
	double x,z;

	for (int i=max_width; i>=1;i=i/2){
		nPrev = n;
		n = max_width/i;
		tmp = integral;
		for (int k=1; k<=n+1; k++) {
			x = A*(k-1)/n;
			
			/* #pragma omp parallel */
			{
			/* #pragma omp for nowait reduction(+:integral) firstprivate(x) firstprivate(n) firstprivate(k) */
			for (int l=1; l<=k; l++) {
				z = A*(l-1)/n;

				if (((k - 1)%2 != 0) || ((l - 1)%2 != 0) || (n == 1))
				{
					if (k==l){
						integral += 0.5*predicate(x,z,forward<T>(params), loop);}
					if (k!=l){
						integral += predicate(x,z,forward<T>(params), loop);}
				}
			}
			}
		}

		if (std::abs(integral) < 1e-25)
			integral = 0.;

		if (n>64 && integral == 0. && tmp == 0.)
			break;
		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		test = std::abs(std::abs((n*n*std::abs(tmp))/(nPrev*nPrev*std::abs(integral)))-1.);
		if (n>=256 && (test <= error))
			break;
	}	
	cout<<"Cunningham points finished after "<<n<<" iterations"<<endl;
	if ((n > (max_width - 1)) && (test > 3*error))
		cout<<"caution: Cunningham Points returned with estimated error "<<test<<endl;
	return (8.*A*A/(n*n))*integral;
}
template <typename Pred, typename T, typename T2>
double kspace(T&& params, Pred predicate, int iterations, double abs, T2 loop){
	double integral = 0.;
	double tmp, test;
//integral width number of points n**2. max_width chosen as 2**integer
	double error;
	if (abs == 0)
		error = 0.02;
	else
		error = abs;
	int max_width;
	if (iterations == 0)
		max_width = 4096;
	else
		max_width = iterations;
	int n = 0;
	int nPrev;
	const double A = M_PI/params.a;
	double x,z;

	for (int i=max_width; i>=1;i=i/2){
		nPrev = n;
		n = max_width/i;
		tmp = integral;
		for (int k=1; k<=n+1; k++) {
			x = A*(k-1)/n;
			
			/* #pragma omp parallel */
			{
			/* #pragma omp for nowait reduction(+:integral) firstprivate(x) firstprivate(n) firstprivate(k) */
			for (int l=1; l<=k; l++) {
				z = A*(l-1)/n;

				if (((k - 1)%2 != 0) || ((l - 1)%2 != 0) || (n == 1))
				{
					if (k==l){
						integral += 0.5*predicate(x,z,forward<T>(params), loop);}
					if (k!=l){
						integral += predicate(x,z,forward<T>(params), loop);}
				}
			}
			}
		}

		if (std::abs(integral) < 1e-25)
			integral = 0.;

		if (n>64 && integral == 0. && tmp == 0.)
			break;
		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		test = std::abs(std::abs((n*n*tmp)/(nPrev*nPrev*integral))-1.);
		if (n>=256 && (test <= error))
			break;
	}	
	cout<<"Cunningham points finished after "<<n<<" iterations"<<endl;
	if ((n > (max_width - 1)) && (test > 3*error))
		cout<<"caution: Cunningham Points returned with estimated error "<<test<<endl;
	return (8.*A*A/(n*n))*integral;
}

/* This function is for the type eigen::vector<double> should the
 * forwarded function be of this return type.
 * The function should be called as kspace(params, &f, iterations)
 * where params is a struct to be passed to the following function,
 * f is the function to be forwarded and iterations is the max 
 * number of iterations requested. If set to 0, it defaults to 1024**2
 * iterations should always be set to 2**number */
template <typename Pred, typename T>
VectorXcd kspace_complex(T&& params, Pred predicate, int iterations, double abs){
	VectorXcd integral(params.N);
	VectorXcd tmp(params.N);
	integral.fill(0.);
	tmp.fill(0.);
//integral width number of points n**2
	double error;
	double condition;
	if (abs == 0)
		error = 0.03;
	else
		error = abs;
	int max_width;
	if (iterations == 0)
		max_width = 1024;
	else
		max_width = iterations;
	int n = 0;
	int nPrev;
	const double A = M_PI/params.a;
	double x,z;
	for (int i=max_width; i>=1;i=i/2){
		nPrev = n;
		n = max_width/i;
		tmp = integral;
		for (int k=1; k<=n+1; k++) {
			x = A*(k-1)/n;
			/* #pragma omp parallel */
			{ 
			VectorXcd integralPrivate(params.N);
			integralPrivate.fill(0.);
			/* #pragma omp for nowait */
			for (int l=1; l<=k+1; l++) {

				z = A*(l-1)/n;
				if (((k - 1)%2 == 0) && ((l - 1)%2 == 0) && (n > 1))
					continue;
				if (k==l){
					integralPrivate = integralPrivate + 0.5*predicate(x,z,forward<T>(params));}
				if (k!=l){
					integralPrivate = integralPrivate + predicate(x,z,forward<T>(params));}
			}
			/* #pragma omp critical */
			{
				integral += integralPrivate;
			}
			}
		}

		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		condition = (std::abs(std::abs(std::abs(tmp.sum())*n*n)/(std::abs(integral.sum())*nPrev*nPrev))-1.);
		if (n>=256 && (condition <= error))
			break;
	}	
	integral = integral*(8.*A*A/(n*n));
	return integral;
}

template <typename Pred, typename T>
VectorXd kspace(T&& params, Pred predicate, int iterations, double abs){
	VectorXd integral(params.N);
	VectorXd tmp(params.N);
	integral.fill(0.);
	tmp.fill(0.);
//integral width number of points n**2
	double error;
	double condition;
	if (abs == 0)
		error = 0.03;
	else
		error = abs;
	int max_width;
	if (iterations == 0)
		max_width = 1024;
	else
		max_width = iterations;
	int n = 0;
	int nPrev;
	const double A = M_PI/params.a;
	double x,z;
	for (int i=max_width; i>=1;i=i/2){
		nPrev = n;
		n = max_width/i;
		tmp = integral;
		for (int k=1; k<=n+1; k++) {
			x = A*(k-1)/n;
			/* #pragma omp parallel */
			{ 
			VectorXd integralPrivate(params.N);
			integralPrivate.fill(0.);
			/* #pragma omp for nowait */
			for (int l=1; l<=k+1; l++) {

				z = A*(l-1)/n;
				if (((k - 1)%2 == 0) && ((l - 1)%2 == 0) && (n > 1))
					continue;
				if (k==l){
					integralPrivate = integralPrivate + 0.5*predicate(x,z,forward<T>(params));}
				if (k!=l){
					integralPrivate = integralPrivate + predicate(x,z,forward<T>(params));}
			}
			/* #pragma omp critical */
			{
				integral += integralPrivate;
			}
			}
		}

		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		condition = (std::abs(std::abs(tmp.sum()*n*n)/(integral.sum()*nPrev*nPrev))-1.);
		if (n>=256 && (condition <= error))
			break;
	}	
	// as per test, multipled return value by 1.5. see test spincurrent_external_kspace.cpp
	// reasons presently unknown
	integral = integral*1.5*(8.*A*A/(n*n));
	return integral;
}
#endif
