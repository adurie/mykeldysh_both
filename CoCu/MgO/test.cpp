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

	vector<vector<Vector3d>> basis;
	vector<Vector3d> test_vec;
	test_vec.reserve(2);
	basis.reserve(5);
	Vector3d tmp_vec;
	tmp_vec << 0, 0, 0;
	for (int k = 0; k < 2; k++)
		test_vec.emplace_back(tmp_vec);
	for (int k = 0; k < 5; k++)
		basis.emplace_back(test_vec);
	basis[0][0] << 1, 3, 4;
	basis[1][1] << 3, 2, 5;
	basis[2][0] << 3, 2, 1;

	/* for (int i = 0; i < basis.size(); i++){ */
	/* 	for (int k = 0; k < test_vec.size(); k++) */
	/* 		cout<<basis[i][k].transpose()<<endl; */
	/* } */
	unordered_map<string, Vector3d> mappy;

	mappy["Cu"]<<0,0,1;
	mappy["Co"]<<0.5, 0.5, 0;
	mappy[species(1)]<<0.5, 0.5, 0;
	/* for (auto const& k : mappy) */
	/* 	cout<<k.first<<" "<<k.second.transpose()<<endl; */
	/* cout<<"Evalyn is a very smelly grotbag!!!"<<endl; */

	/* for (int i = 0; i < 2; i++){ */
	/* 	for (int j = 0; j < 2; j++){ */
	/* 		for (int k = 0; k < umap["hi there"][i][j].size(); k++){ */
	/* 			cout<<umap["hi there"][i][j][k]<<endl<<endl; */
	/* 		} */
	/* 	} */
	/* } */

	Matrix<double, 8, 8> blockmat;
	Matrix4d block1;
	Matrix4d block2;
	block1 = Matrix4d::Zero();
	block2 = Matrix4d::Ones();
	//topLeft
	blockmat.block<4,4>(0,0) = block1;
	//bottomLeft
	blockmat.block<4,4>(4,0) = block1;
	//topRight
	blockmat.block<4,4>(0,4) = block2;
	//bottomRight
	blockmat.block<4,4>(4,4) = block2;
	cout<<blockmat<<endl;
	
	/* unordered_map<string, M9> onsite_up, onsite_dn; */
	/* unordered_map<string, vector<vector<double>>> hop_up, hop_dn; */
	/* int numats = numatoms(); */
	/* vector<double> temporary_vector; */
	/* int numnns = numnn(); */
	/* for (int jj = 1; jj < numats + 1; jj++){ */
	/* 	onsite_up[species(jj)] = U(jj, 0); */
	/* 	onsite_dn[species(jj)] = U(jj, 1); */
	/* 	for (int kk = 1; kk < numnns + 1; kk++){ */
	/* 		temporary_vector = param(jj, kk, 0); */
	/* 		hop_up[species(jj)].emplace_back(temporary_vector); */
	/* 		temporary_vector = param(jj, kk, 1); */
	/* 		hop_dn[species(jj)].emplace_back(temporary_vector); */
	/* 	} */
	/* } */
	/* for (auto const& k : hop_dn){ */
	/* 	for (int j = 0; j < numnns; j++){ */
	/* 		cout<<k.first<<" "<<j<<endl; */
	/* 		for (int l = 0; l < 10; l++) */
	/* 			cout<<k.second[j][l]<<endl; */
	/* 		cout<<endl; */
	/* 	} */
	/* } */

	/* for (auto const& k : onsite_dn) */
	/* 	cout<<k.first<<endl<<k.second<<endl<<endl; */
	/* for (auto const& k : onsite_up) */
	/* 	cout<<k.first<<endl<<k.second<<endl<<endl; */

	return 0;
}

