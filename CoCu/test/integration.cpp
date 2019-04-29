#include <iostream>
#include <cmath>
/* #include <gsl/gsl_integration.h> */
#include "vector_integration.h"
#include <vector>

using namespace std;

double f(double x, vector<double> &vec_result, void * params){
	double alpha = *(double *) params;
	double f = log(alpha*x)/sqrt(x);
	for (int k = 0; k < vec_result.size(); k++){
		vec_result[k] = f;
		f *= x;
	}
	double result = 0;
	for (int k = 0; k < vec_result.size(); k++)
		result += vec_result[k];
	return result;
}

int main(){
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(2000);
	double result, error;
	double expected = -4.;
	double alpha = 1.;
	vector<double> vec_result;
	int n = 10; //number of dependent variables to integrate over
	vec_result.reserve(n);
	for (int k = 0; k < n; k++)
		vec_result.push_back(0.);

	my_gsl_function F;
	F.function = &f;
	F.params = &alpha;

	gsl_integration_qag(&F, 0, 1, 0, 1e-7, 1000, 1, w, vec_result, &result, &error);

	cout<<"estimated error = "<<error<<endl;
	cout<<"result from qag = "<<result<<endl;
	cout<<endl;

	gsl_integration_workspace_free (w);
	cout<<"individual results"<<endl;
	double test_result = 0;
	for (int k = 0; k < n; k++){
		cout<<k<<" "<<vec_result[k]<<endl;
		test_result += vec_result[k];
	}
	cout<<"sum of individual results = "<<test_result<<endl;

	return 0;
}
