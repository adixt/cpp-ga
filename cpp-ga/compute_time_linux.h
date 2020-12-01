#pragma once
#include <ctime>

using namespace std;

clock_t c_start;
clock_t c_end;

void ComputeTimeStart()
{
	c_start = clock();
}

double ComputeTimeEnd()
{
	c_end = clock();
	return 1000.0 * (static_cast<double>(c_end) - static_cast<double>(c_start)) / CLOCKS_PER_SEC;;
}