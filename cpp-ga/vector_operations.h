#pragma once
#include<vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <math.h> 

using namespace std;

void MultiplyVectorByScalar(vector<float>& v, float k) {
	transform(v.begin(), v.end(), v.begin(), [k](float& c) { return c * k; });
}

void Multiply2DVectorByScalar(vector<vector<float>>& v, float k) {
	for_each(v.begin(), v.end(),
		[&k](vector<float>& v) {
		MultiplyVectorByScalar(v, k);
	}
	);
}

template <typename T>
void PrintVector(vector<T>& v, string vectorName)
{
	cout << "\n" << "=============== " << vectorName << " ===============\n";
	for (int i = 0; i < v.size(); i++) {
		cout << v[i] << " ";
	}
	cout << "\n";
	cout << "==========================================\n";
}

template <typename T>
void Print2DVector(vector<vector<T>>& v, string vectorName)
{
	cout << "\n" << "=============== " << vectorName << " ===============\n";
	for (int i = 0; i < v.size(); i++) {
		vector<T> tmp = v[i];
		for (int j = 0; j < tmp.size(); j++) {
			cout << tmp[j] << " ";
		}
		cout << "\n";
	}
	cout << "==========================================\n";
}

template <typename T>
void Print3DVector(vector<vector<vector<T>>>& v, string vectorName)
{
	cout << "\n" << "=============== " << vectorName << " ===============\n";
	for (int i = 0; i < v.size(); i++) {
		vector<vector<T>> tmp2D = v[i];
		for (int j = 0; j < tmp2D.size(); j++) {
			vector<T> tmp = tmp2D[j];
			for (int k = 0; k < tmp.size(); k++) {
				cout << tmp[k] << " ";
			}
			cout << "\n";
		}
		cout << "\n\n\n";
	}
	cout << "==========================================\n";
}

template <typename T>
double avg(vector<T>& v)
{
	double return_value = 0.0;
	int n = v.size();

	for (int i = 0; i < n; i++)
	{
		return_value += v[i];
	}

	return (return_value / n);
}

template <typename T>
double variance(vector<T>& v)
{
	double mean = avg(v);
	double sum = 0.0;
	double temp = 0.0;
	double var = 0.0;
	for (int j = 0; j < v.size(); j++)
	{
		double number = v[j] - mean;
		double powed = pow(number, 2);
		sum += powed;
	}
	var = sum / (v.size() - 1);
	return var;
}

template <typename T>
vector<T> GetVarianceFrom2dVector(vector<vector<T>>& v)
{
	int columSize = v.size(); 
	int rowSize = v[0].size();
	vector<T> calculatedVariance(rowSize, 0);

	for (int j = 0; j < rowSize; j++) {
		vector<T> columns(columSize, 0);
		for (int i = 0; i < columSize; i++) {
			columns[i] = v[i][j];
		}
		double var = variance<T>(columns);
		calculatedVariance[j] = static_cast<T>(var);
	}
	return calculatedVariance;
}

template<typename T>
T varianceImproved(const vector<T> &vec) {
	const size_t sz = vec.size();
	if (sz == 1) {
		return 0.0;
	}

	// Calculate the mean
	const T mean = accumulate(vec.begin(), vec.end(), 0.0) / sz;

	// Now calculate the variance
	auto variance_func = [&mean, &sz](T accumulator, const T& val) {
		return accumulator + ((val - mean)*(val - mean) / (sz - 1));
	};

	return accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
struct square
{
	T operator()(const T& Left, const T& Right) const
	{
		return (Left + Right*Right);
	}
};

template<typename Iter_T>
long double vectorNorm(Iter_T first, Iter_T last) {
	return sqrt(inner_product(first, last, first, 0.0L));
}