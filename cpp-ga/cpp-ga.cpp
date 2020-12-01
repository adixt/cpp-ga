#include <iostream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>
#include "cpu_rnd.h"
#include "compute_time_linux.h"
using namespace std;
const unsigned int N = 3;
const unsigned int GENERATION_SIZE = 100;
const unsigned int start_xd_min[N] = { 2112,1180,417 };
unsigned int individuals[N][GENERATION_SIZE] = {};

const int n = 3;
const int ns = 2;
const int allNodes = n + ns;
int LT[allNodes][n] = {};

const int x_inf = 100000;
const int simTime = 60;
#include "vector_types.h"
#include "vector_operators.h"
#include "vector_operations.h"
#include "eigen_operations.h"
int main()
{
	ComputeTimeStart();
	for (int node = 0; node < N; node++) {
		for (int individualNo = 0; individualNo < GENERATION_SIZE; individualNo++) {
			individuals[node][individualNo] = static_cast<int>(round(rnd() * start_xd_min[node]));
		}
	}

	for (int node = 0; node < N; node++) {
		cout << "\nNODE " << node + 1 << ":\n";
		for (int individualNo = 0; individualNo < GENERATION_SIZE; individualNo++) {
			cout << individuals[node][individualNo] << " ";
		}
		cout << "\n";
	}

	LT[3][0] = 2;
	LT[4][1] = 4;
	LT[1][2] = 3;
	LT[0][2] = 3;
	LT[0][1] = 1;

	int* start = &LT[0][0];
	// max lead time
	int L = *max_element(start, start + allNodes * n);

	float LA_nom[allNodes][n] = {};
	three_dimension_vector_float LA(simTime, vector<vector<float>>(allNodes, vector<float>(n)));

	LA[0][3][0] = LA_nom[3][0] = 1;
	LA[0][4][1] = LA_nom[4][1] = 0.8f;
	LA[0][1][2] = LA_nom[1][2] = 0.6f;
	LA[0][0][2] = LA_nom[0][2] = 0.4f;
	LA[0][0][1] = LA_nom[0][1] = 0.2f;
	//Print2DVector<float>(LA[0], "delay matrix");

	// Verify if allocation correct - elements in each column should sum up to 1 or 0	
	for (int j = 0; j < n; j++) {
		float temp = 0;
		for (int i = 0; i < allNodes; i++) {
			temp = temp + LA[0][i][j];
		}

		if (temp != 0 && temp != 1) {
			throw ("Improper allocation in column: %d", j);
		}
	}

	// Initial conditions
	int time[simTime] = {};
	int u[n][simTime] = {};
	two_dimension_vector_float u_hist(n, vector<float>(simTime, 0)); // order history
	two_dimension_vector_int x(simTime, vector<int>(allNodes, 0));
	int y[allNodes][simTime] = {};
	int xd[allNodes] = { 80, 120, 98, 0, 0 };

	x[0] = { 70, 140, 88, x_inf, x_inf }; // initial stock level
	//Print2DVector<int>(x, "stock level ");

	// Demand
	vector<int> dmax{ 10, 15, 20 };
	two_dimension_vector_int d(n, vector<int>(simTime, {})); // int d[n][simTime] = {};

	for (int j = 0; j < simTime; j++) {
		double multiplier = 1;
		for (int k = 0; k < n; k++) {
			double demand = multiplier * dmax[k] * rnd();
			int randomDemand = static_cast<int>(round(demand));
			if (randomDemand > dmax[k]) d[k][j] = dmax[k];
			else d[k][j] = randomDemand;
		}
	};
	//Print2DVector<int>(d, "demand");

	// State - space description
	// System matrices
	three_dimension_vector_float B_nom(L, two_dimension_vector_float(n, vector<float>(n)));
	four_dimension_vector_float B(simTime, three_dimension_vector_float(L, two_dimension_vector_float(n, vector<float>(n))));

	// Assuming zero order processing time(eq 9)
	two_dimension_vector_float B_0(n, vector<float>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			B_0[i][j] = LA[0][i][j] * -1;
		}
	}

	// index k corresponds to delay k(eq 8)	
	for (int k = 0; k < L; k++) {
		for (int j = 0; j < n; j++) {
			float t_sum = 0;
			for (int i = 0; i < allNodes; i++) {
				if (LT[i][j] == k + 1) {
					t_sum = t_sum + LA[0][i][j];
				}
			}
			B_nom[k][j][j] = t_sum;
		}
	}
	//Print3DVector<float>(B_nom, "B matrix");
	B[0] = B_nom;

	//% Sum of delay matrices
	two_dimension_vector_float Lambda(n, vector<float>(n, {}));

	// table index k corresponds to delay k
	for (int k = 0; k < L; k++) {
		Lambda = Lambda + B[0][k];
	}

	Lambda = Lambda + B_0; // eq 11 

	vector<float> xd_min = GetXdMin(n, L, B, Lambda, dmax);
	PrintVector(xd_min, "xd min");

	double tt = ComputeTimeEnd();
	cout << "CPU time: " << tt << " ms\n";

	return 0;
}

