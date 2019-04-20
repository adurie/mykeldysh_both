#ifndef CUNNINGHAM_H
#define CUNNINGHAM_H

#include <cmath>
#include <eigen3/Eigen/Dense>
//this is the first TRUE cunningham points algorithm I have made
//with automatic convergence checks. Parameters include the max
//number of points in one meridian (iterations), the relative 
//error value, and the number of points in one meridian to start
//with. 27-06-17
//
//Use this for abs error, the normal one for rel error

using namespace std;

double Condition(std::complex<double> result, std::complex<double> part_result){
	double condition;

	/* if ((std::abs(part_result) == 0) && (std::abs(result) == 0)) */
	/* 	condition = 0; */
	/* else */
	/* 	condition = std::abs(std::abs(9.*part_result)/std::abs(result)-1.); */

	condition = std::abs(9.*part_result - result);

	return condition;
}
double Condition(double result, double part_result){
	double condition;
	
	/* if ((std::abs(part_result) == 0) && (std::abs(result) == 0)) */
	/* 	condition = 0; */
	/* else */
	/* 	condition = std::abs(std::abs(std::abs(9.*part_result)/(result))-1.); */

	condition = std::abs(9.*part_result - result);
	
	return condition;
}
double Condition(const Eigen::VectorXd &result, const Eigen::VectorXd &part_result){
	double condition;

	/* if ((std::abs(part_result.sum()) == 0) && (std::abs(result.sum()) == 0)) */
	/* 	condition = 0; */
	/* else */
	/* 	condition = std::abs(std::abs(9.*part_result.sum())/std::abs(result.sum())-1.); */

	condition = std::abs(9.*part_result.sum() - result.sum());

	return condition;
}
double Condition(const Eigen::VectorXcd &result, const Eigen::VectorXcd &part_result){
	double condition;

	/* if ((std::abs(part_result.sum()) == 0) && (std::abs(result.sum()) == 0)) */
	/* 	condition = 0; */
	/* else */
	/* 	condition = std::abs(std::abs(9.*part_result.sum())/std::abs(result.sum())-1.); */

	condition = std::abs(9.*part_result.sum() - result.sum());

	return condition;
}

template <typename func, typename... Args>
auto kspace(func&& predicate, int iterations, double rel, int start, const double a, Args&&... params)
    -> decltype(std::forward<func>(predicate)(std::declval<double>(), std::declval<double>(), std::declval<double>(),std::forward<Args>(params)...)) {
	    typedef decltype(std::forward<func>(predicate)(std::declval<double>(), std::declval<double>(), std::declval<double>(),std::forward<Args>(params)...)) ret;

	double error;
	if (rel == 0)
		error = 0.002;
	else
		error = rel;
	int max_width;
	if (iterations == 0)
		max_width = 400;
	else
		max_width = 2*iterations; 	//2 here to cater for the doubling of N

	int N; 	//starting number of points in one meridian - editable
	if (start == 0)
		N = 1;
	else
		N = start;

	int n = 2*N;	//this caters for the special spacing but see conditions for odd numbers 
			//below. This ensures number of points in one meridian is indeed n
	
	/* string Mydata = "test.txt"; */
	/* ofstream Myfile; */	
	/* Myfile.open( Mydata.c_str(),ios::trunc ); */

	const double A = M_PI/a;
	double x,z;
	double factor;
	ret integral;
	ret tmp;
	double rel_error;

	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = A*k/n;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					z = A*l/n;
					/* Myfile<<x<<" , "<<z<<endl; */
					if ((k==1) && (l==1))
						integral = 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
					else{
						if (k==l){
							/* integral += 0.5*predicate(x,z,forward<T>(params)); */
							integral = integral + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						}
						else
							integral = integral + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
					}
				}
			}
		}
	}
	
	int Nnew;
	int condition = -2;
	double foo = 0;

	while (condition != 1){
		Nnew = 3*n;
		tmp = integral;
		for (int k = 0; k!=Nnew+1; k++){
			if (k%2!=0){
				x = A*k/Nnew;
				for (int l = 0; l!=k+1; l++){
					if ((l%2!=0) && ((k%3!=0) || (l%3!=0))){
						z = A*l/Nnew;
						/* Myfile<<x<<" , "<<z<<endl; */
						if (k==l){
							/* integral += 0.5*predicate(x,z,forward<T>(params)); */
							integral = integral + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						}
						else{
							/* integral += predicate(x,z,forward<T>(params)); */
							integral = integral + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						}
					}
				}
			}
		}
		n = Nnew;

		/* if (std::abs(integral) < 1e-25) */
		/* 	integral = 0.; */
		/* if (Nnew>64 && integral == 0. && tmp == 0.) */
		/* 	condition = 1; */

		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		rel_error = Condition(integral, tmp);//condition overloaded for all data types this header is designed for
		factor = 4*8.*A*A/(Nnew*Nnew);
		rel_error *= factor;
		cout<<rel_error<<endl;
		if (Nnew > max_width){
			cout<<"error, maximum number of iteration points reached. Estimated error = "<<rel_error<<endl;
			condition = 1;
		}
		if (Nnew>=16 && rel_error <= error)
			condition = 1;
	}
	integral = integral*factor;
	cout<<"converged for "<<Nnew/2<<" points"<<endl;
	ret integral_right;
	ret tmp_right;
	n = N*2;

	//this section of code makes the routine specific to Brillouin zones with corners (+/-2A, 0) and (0, +/-2A)
	for (int k = 0; k!=n+1; k++){
		if (k%2!=0){
			x = A*k/n + A;
			for (int l = 0; l!=k+1; l++){
				if (l%2!=0){
					z = A - A*l/n;
					/* Myfile<<x<<" , "<<z<<endl; */
					if ((k==1) && (l==1))
						integral_right = 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
					else{
						if (k==l){
							/* integral += 0.5*predicate(x,z,forward<T>(params)); */
							integral_right = integral_right + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						}
						else
							integral_right = integral_right + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
					}
				}
			}
		}
	}
	
	condition = -2;
	foo = 0;

	while (condition != 1){
		Nnew = 3*n;
		tmp_right = integral_right;
		for (int k = 0; k!=Nnew+1; k++){
			if (k%2!=0){
				x = A*k/Nnew + A;
				for (int l = 0; l!=k+1; l++){
					if ((l%2!=0) && ((k%3!=0) || (l%3!=0))){
						z = A - A*l/Nnew;
						/* Myfile<<x<<" , "<<z<<endl; */
						if (k==l){
							/* integral += 0.5*predicate(x,z,forward<T>(params)); */
							integral_right = integral_right + 0.5*std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						}
						else{
							/* integral += predicate(x,z,forward<T>(params)); */
							integral_right = integral_right + std::forward<func>(predicate)(x,z,a,std::forward<Args>(params)...);
						}
					}
				}
			}
		}
		n = Nnew;

		/* if (std::abs(integral) < 1e-25) */
		/* 	integral = 0.; */
		/* if (Nnew>64 && integral == 0. && tmp == 0.) */
		/* 	condition = 1; */

		//256 seems to be the magic number for the lower bound of sampling points
		//3% margin of error as below works well when benchmarked against fixed method 
		//of double Simpson's integral of rho-E/LDOS.cpp
		
		rel_error = Condition(integral_right, tmp_right);//condition overloaded for all data types this header is designed for
		factor = 4*8.*A*A/(Nnew*Nnew);
		rel_error *= factor;
		cout<<rel_error<<endl;
		if (Nnew > max_width){
			cout<<"error, maximum number of iteration points reached. Estimated error = "<<rel_error<<endl;
			condition = 1;
		}
		if (Nnew>=16 && rel_error <= error)
			condition = 1;
	}
	integral_right = integral_right*factor;
	cout<<"converged for "<<Nnew/2<<" points"<<endl;

	return 0.5*(integral + integral_right);
}
#endif
