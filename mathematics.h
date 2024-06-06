#pragma once
#include <vector>
#include <string>

std::vector<double> calculateUpdateCoefficients(const std::vector<double>& inputArray);
std::vector<double> calculateRefractiveIndexes(const std::vector<double>& permittivity, const std::vector<double>& permeability, bool approximatePermeabilityAsOne);
std::vector<double> upscaleVector(const std::vector<double>& input, const int newSize);
double gaussianSource(double t, double tau, double a);
std::vector<double> precomputeElectricSource();
std::vector<double> precomputeMagneticSource(double epsilon_r_src = 1, double mu_r_src = 1);

template<typename T>
inline double findExtremum(const std::vector<T>& data, bool maximum, bool returnIndex)
{
	//if there are multiple of the same entry it will return the index of the LATEST occurance of that value
	//initialize the comparison threshold as the minimum double when looking for the max, and the maximum double when looking for the min.
	double threshold = maximum ? DBL_MIN : DBL_MAX;
	int extremeIndex = 0; //index of the most extreme value so far
	//maximum search
	if (maximum)
	{
		for (size_t i = 0; i < data.size(); i++) {
			if (data[i] > threshold)
			{
				//if the current value is the biggest one yet we record it and its index.
				threshold = data[i];
				extremeIndex = i;
			}
		}
	}
	else //we are looking for the minimum
	{
		for (size_t i = 0; i < data.size(); i++) {
			if (data[i] < threshold)
			{
				//if the current value is the biggest one yet we record it and its index.
				threshold = data[i];
				extremeIndex = i;
			}
		}
	}
	if (returnIndex)
	{
		return extremeIndex;
	}
	return threshold;
}
