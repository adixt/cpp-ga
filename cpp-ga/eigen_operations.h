#pragma once
using namespace std;
#include "vector_types.h"
#include "vector_operations.h"
#include "vector_operators.h"
#include <Eigen/Dense>

vector<float> GetXdMin(int n, int L, four_dimension_vector_float&  B, two_dimension_vector_float&  Lambda, vector<int>&  dmax) {
	// Reference stock level for full demand satisfaction
	two_dimension_vector_float summed(n, vector<float>(n));

	// table index k corresponds to the delay k
	for (int k = 0; k < L; k++)
	{
		two_dimension_vector_float multiplied = B[0][k];
		Multiply2DVectorByScalar(multiplied, static_cast<float>(k + 1));
		summed = summed + multiplied;
	}

	// Initialize the matrix as a NxN array of 0.
	two_dimension_vector_float myEye(n, vector<float>(n, 0));
	// Set the diagonal to be 1s
	for (int t = 0; t < n; t++) {
		myEye[t][t] = 1;
	}
	two_dimension_vector_float step1 = myEye + summed;

	Eigen::MatrixXf step1M(n, n);
	Eigen::MatrixXf step2M(n, n);
	Eigen::MatrixXf step3M(n, n);
	Eigen::MatrixXf step4M(n, 1);
	Eigen::MatrixXf step5M(n, 1);

	for (int i = 0; i < n; i++) {
		step4M(i, 0) = static_cast<float>(dmax[i]);
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			step1M(i, j) = step1[i][j];
			step2M(i, j) = Lambda[i][j];
		}
	}
	step3M = step1M * step2M.inverse();
	step5M = step3M * step4M;
	vector<float> xd_min(n, 0);
	for (int i = 0; i < n; i++) {
		xd_min[i] = step5M(i, 0);	// eq 23	
	}
	return xd_min;
}


vector<float> Multiple2dEigenBy1dEigen(two_dimension_vector_float&  MatrixA, vector<float>&  MatrixB) {
	int sizeA = MatrixA.size();
	int sizeB = MatrixB.size();
	Eigen::MatrixXf step1M(sizeA, sizeB);
	Eigen::MatrixXf step2M(sizeB, 1);
	Eigen::MatrixXf step3M(sizeB, 1);

	for (int i = 0; i < sizeA; i++) {
		for (int j = 0; j < sizeB; j++) {
			step1M(i, j) = MatrixA[i][j];
		}
	}

	for (int i = 0; i < sizeB; i++) {
		step2M(i, 0) = MatrixB[i];
	}

	step3M = step1M * step2M;
	vector<float> result(sizeA);
	for (int i = 0; i < sizeA; i++) {
		result[i] = step3M(i, 0);
	}

	return result;
}