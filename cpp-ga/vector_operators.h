#pragma once
#include <vector>
#include <functional>
#include <algorithm>
#include <iterator>
#include <set>
using namespace std;

template <typename T>
vector<T> operator*(const vector<T>& a, const vector<T>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), multiplies<T>());
	return result;
}

template <typename T>
vector<T> operator/(const vector<T>& a, const vector<T>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), divides<T>());
	return result;
}

template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), plus<T>());
	return result;
}

template <typename T>
vector<T> operator-(const vector<T>& a, const vector<T>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<T> result;
	result.reserve(a.size());

	transform(a.begin(), a.end(), b.begin(),
		back_inserter(result), minus<T>());
	return result;
}

template <typename T>
vector<vector<T>> operator+(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<vector<T>> result;
	result.reserve(a.size());

	for (int i = 0; i < a.size(); i++) {
		auto tmp = a[i] + b[i];
		result.push_back(tmp);
	}

	return result;
}

template <typename T>
vector<vector<T>> operator-(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<vector<T>> result;
	result.reserve(a.size());

	for (int i = 0; i < a.size(); i++) {
		auto tmp = a[i] - b[i];
		result.push_back(tmp);
	}

	return result;
}

template <typename T>
vector<vector<T>> operator*(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<vector<T>> result;
	result.reserve(a.size());

	for (int i = 0; i < a.size(); i++) {
		auto tmp = a[i] * b[i];
		result.push_back(tmp);
	}

	return result;
}

template <typename T>
vector<vector<T>> operator/(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
	if (a.size() != b.size()) {
		throw "Vectors must be equal size";
	}

	vector<vector<T>> result;
	result.reserve(a.size());

	for (int i = 0; i < a.size(); i++) {
		auto tmp = a[i] / b[i];
		result.push_back(tmp);
	}

	return result;
}
