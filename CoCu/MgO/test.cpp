#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/src/Core/util/MKL_support.h>
#include "AuMgOFe.h"

using namespace std;
typedef complex<double> dcomp;
typedef vector<Matrix<dcomp, 9, 9>, aligned_allocator<Matrix<dcomp, 9, 9>>> vM;
typedef Matrix<complex<double>, 9, 9> M9;

int main(){
	unordered_map<string, vector<vector<vM>>> umap;
	M9 mat1, mat2, mat3, mat4, zero;
	vM tmp1;
	vector<vM> tmp2;
	vector<vector<vM>> tmp3;
	mat1 = U(2,0);
	mat2 = U(1,0);
	mat3 = U(1,1);
	mat4 = U(3,0);
	zero = M9::Zero();
	tmp1.emplace_back(zero);
	for (int i = 0; i < 2; i++)
		tmp2.emplace_back(tmp1);
	for (int i = 0; i < 2; i++)
		tmp3.emplace_back(tmp2);
	umap["hey!"] = tmp3;
	umap["hi there"] = tmp3;

	umap["hey!"][0][0][0] = mat1;
	umap["hey!"][0][1][0] = mat2;
	umap["hey!"][1][0][0] = mat3;
	umap["hey!"][1][1][0] = mat4;

	umap["hi there"][0][0][0] = mat1;
	umap["hi there"][0][1][0] = mat2;
	umap["hi there"][1][0][0] = mat3;
	umap["hi there"][1][1][0] = mat4;
	umap["hi there"][0][0].emplace_back(mat1);
	umap["hi there"][0][1].emplace_back(mat2);
	umap["hi there"][1][0].emplace_back(mat3);
	umap["hi there"][1][1].emplace_back(mat4);

	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < umap["hi there"][i][j].size(); k++){
				cout<<umap["hi there"][i][j][k]<<endl<<endl;
			}
		}
	}

	return 0;
}

