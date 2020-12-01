#pragma once
#include<vector>
#include <algorithm>
#include <string>
#include <iostream>

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