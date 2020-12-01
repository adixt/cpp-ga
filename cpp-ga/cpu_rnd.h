#pragma once
#include <random>
#include <functional>

using namespace std;

// uniform distribution
default_random_engine generator;
uniform_real_distribution<double> distribution(0, 1);
double dice_roll = distribution(generator);
auto rnd = bind(distribution, generator);