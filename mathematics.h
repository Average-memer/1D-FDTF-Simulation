#pragma once
#include <vector>
#include <string>

template <typename T>
double findExtremum(const std::vector<T>& data, bool maximum, bool returnIndex);

std::vector<double> calculateUpdateCoefficients(const std::vector<double>& inputArray);
std::vector<double> precomputeSource();
std::vector<double> calculateRefractiveIndexes(const std::vector<double>& permittivity, const std::vector<double>& permeability, bool approximatePermeabilityAsOne);
std::vector<double> upscaleVector(const std::vector<double>& input, size_t newSize);
double gaussianSource(double t, double tau, double a);
std::vector<double> precomputeElectricSource();
std::vector<double> precomputeMagneticSource(double epsilon_r_src = 1, double mu_r_src = 1);
