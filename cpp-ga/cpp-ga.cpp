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

const float x_inf = 100000;
const int simTime = 10000;
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

	two_dimension_vector_float LA_nom(allNodes, vector<float>(n));// float LA_nom[allNodes][n] = {};
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
	two_dimension_vector_float u(simTime, vector<float>(n, 0)); //order quantity
	two_dimension_vector_float ur(simTime, vector<float>(n, 0)); //orders received
	two_dimension_vector_float u_hist(simTime, vector<float>(n, 0)); // order history
	two_dimension_vector_float x(simTime + 1, vector<float>(allNodes, 0)); //stock level
	two_dimension_vector_float y(simTime, vector<float>(allNodes, 0));
	two_dimension_vector_float h(simTime, vector<float>(allNodes, 0)); //satisfied demand
	two_dimension_vector_float d(simTime, vector<float>(allNodes, 0)); //demand
	two_dimension_vector_float S_u(simTime, vector<float>(allNodes, 0));
	two_dimension_vector_float S_u_mod(simTime, vector<float>(allNodes, 0));
	vector<float> xd = { 80, 120, 98, 0, 0 };
	vector<float> xd_sub(&xd[0], &xd[n]);

	x[0] = { 70, 140, 88, x_inf, x_inf }; // initial stock level
										  //Print2DVector<int>(x, "stock level ");

										  // Demand
	vector<int> dmax{ 10, 15, 20 };

	for (int j = 0; j < simTime; j++) {
		double multiplier = 1;
		for (int k = 0; k < n; k++) {
			double demand = multiplier * dmax[k] * rnd();
			int randomDemand = static_cast<int>(round(demand));
			if (randomDemand > dmax[k]) d[j][k] = dmax[k];
			else d[j][k] = randomDemand;
			// FOR A DEBUG, YOU CAN SET CONST DEMAND = MAXDEMAND
			//d[j][k] = dmax[k];
		}
	};
	//Print2DVector<float>(d, "demand");

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


	// Main OUT policy loop
	for (int t = 0; t < simTime; t++) {
		// Calculate input proposal
		//First determine open - order quantity
		vector<float> o_q(n, 0);
		if (t < L) {
			int counter = 0;
			for (int i = 0; i < L; i++)
			{
				for (int j = i; j < L; j++)
				{
					if (t - i > 0) {
						vector<float> step1 = Multiple2dEigenBy1dEigen(B[t - i - 1][j], u_hist[t - i - 1]);
						o_q = o_q + step1;
					}
				}
			}
		}
		else {
			for (int i = 0; i < L; i++)
			{
				for (int j = i; j < L; j++)
				{
					vector<float> step1 = Multiple2dEigenBy1dEigen(B[t - i - 1][j], u_hist[t - i - 1]);
					o_q = o_q + step1;
				}
			}
		}
		// Calculate input proposal
		// Classical OUT
		vector<float> x_sub(&x[t][0], &x[t][n]);
		u[t] = xd_sub - x_sub - o_q; // (eq 14)


		for (int i = 0; i < u[t].size(); i++) {
			if (u[t][i] < 0) {
				u[t][i] = 0;
			}
		}

		u_hist[t] = u[t];

		if (t < L) {
			for (int j = 0; j < n; j++) { // j - index of recipient
				ur[t][j] = 0;
				for (int i = 0; i < allNodes; i++) // i - index of sender
				{
					if (t - LT[i][j] >= 0)
					{
						int lt = LT[i][j];
						int ts = t - lt;
						ur[t][j] = ur[t][j] + LA[ts][i][j] * u_hist[ts][j];

					}
				}

			}
		}
		else {
			for (int j = 0; j < n; j++) { // j - index of recipient
				ur[t][j] = 0;
				for (int i = 0; i < allNodes; i++) // i - index of sender
				{
					int lt = LT[i][j];
					int ts = t - lt;
					ur[t][j] = ur[t][j] + LA[ts][i][j] * u_hist[ts][j];
				}
			}
		}

		// Satisfied demand
		vector<float> d_sub(&d[t][0], &d[t][n]);
		if (d_sub < x_sub + ur[t]) {
			h[t] = d_sub;
		}
		else {
			h[t] = x_sub;
		}

		// Goods sent inside controlled network
		// Available stock for internal demand
		y[t] = x[t];

		vector<float> y_sub(&y[t][0], &y[t][n]);
		y_sub = x_sub + ur[t] - h[t];
		copy_n(y_sub.cbegin(), n, y[t].begin());

		// Total internal requests imposed at controlled nodes and external suppliers
		S_u[t] = Multiple2dEigenBy1dEigen(LA[0], u[t]);

		// Check if the available stock suffices to satisfy the internal requests, separately for each
		// node.If not - the deficiency is distributed according to the load generated by
		// request originators
		// Other methods to explore later on
		LA[t] = LA_nom;

		for (int i = 0; i < n; i++) { // i is the node index inside the controlled network
			float tmp = y[t][i] - S_u[t][i];

			if (tmp < 0) {
				// Decrease the allocation weights corresponding to requests placed at internal nodes
				// Keep input history intact
				for (int j = 0; j < n; j++) {
					LA[t][j][i] = LA_nom[j][i] * (1 + tmp / S_u[t][i]);

					if (LA[t][j][i] < 0) {
						bool flaga = true;
					}
				}
			}
		}

		// Update delay matrices B for delay 1 to L, update of B_0 is not required
		B[t] = B_nom;

		for (int k = 0; k < L; k++) {// index k corresponds to delay k
			for (int j = 0; j < n; j++) {
				float t_sum = 0;
				for (int i = 0; i < allNodes; i++) {
					if (LT[i][j] == k + 1) {
						t_sum = t_sum + LA[t][i][j];
					}
				}
				B[t][k][j][j] = t_sum;
			}
		}

		S_u_mod[t] = Multiple2dEigenBy1dEigen(LA[t], u_hist[t]);

		x[t + 1] = y[t] - S_u_mod[t];

	}
	// End of OUT Main loop

	vector<float> varianceU = GetVarianceFrom2dVector<float>(u_hist);
	PrintVector(varianceU, "Variance of U");

	vector<float> varianceD = GetVarianceFrom2dVector<float>(d);
	vector<float> varianceD_sub(&varianceD[0], &varianceD[n]);
	PrintVector(varianceD_sub, "Variance of Demand");

	vector<float> varianceX = GetVarianceFrom2dVector<float>(x);
	vector<float> varianceX_sub(&varianceX[0], &varianceX[n]);
	PrintVector(varianceX_sub, "Variance of X");


	// Calculate BI Indicators
	float maxVarU = *max_element(varianceU.begin(), varianceU.end());
	float minVarU = *min_element(varianceU.begin(), varianceU.end());
	float maxVarD = *max_element(varianceD_sub.begin(), varianceD_sub.end());
	float minVarD = *min_element(varianceD_sub.begin(), varianceD_sub.end());

	float w1 = maxVarU / minVarD;
	cout << "W1: max Var U / min Var D -> theoretical max lim BI: " << w1 << endl;

	float w2 = minVarU / maxVarD;
	cout << "W2: min Var U / max Var D -> theoretical min lim BI: " << w2 << endl;

	vector<float> biVector = varianceU / varianceD_sub;
	float w3 = *max_element(biVector.begin(), biVector.end());
	cout << "W3: max Var U / max Var D -> BI determinant, ||x||Inf, LpInf: " << w3 << endl;

	//float summedSquaresVarU = accumulate(varianceU.begin(), varianceU.end(), 0, square<float>());
	//float euclideanVarU = sqrt(summedSquaresVarU); // this gives the proper result too
	//float summedSquaresVarD = accumulate(varianceD_sub.begin(), varianceD_sub.end(), 0, square<float>());
	//float euclideanVarD = sqrt(summedSquaresVarD); // this gives the proper result too
	float euclideanVarU = vectorNorm(varianceU.begin(), varianceU.end());
	float euclideanVarD = vectorNorm(varianceD_sub.begin(), varianceD_sub.end());
	float w4 = euclideanVarU / euclideanVarD;
	cout << "W4: Euclidean Var U / Euclidean  Var D -> mean BI determinant, ||x||2, Lp 2: " << w4 << endl;

	float manhattanVarU = accumulate(varianceU.begin(), varianceU.end(), 0);
	float manhattanVarD = accumulate(varianceD_sub.begin(), varianceD_sub.end(), 0);
	float w5 = manhattanVarU / manhattanVarD;
	cout << "W5: Manhattan Var U / Manhattan  Var D -> ||x||1, Lp 1, maximum absolute column sum: " << w5 << endl;
	// End of BI calculations

	double tt = ComputeTimeEnd();
	cout << "CPU time: " << tt << " ms\n";
	system("pause");
	return 0;
}

