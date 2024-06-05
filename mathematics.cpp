#pragma once
#include "mathematics.h"
#include "utility.h"
#include "constants.h"

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

using std::cout;
using std::endl;

//a header file for mathematical functions. The goal is to make utility.h less crowded.

//a simple gaussian source. We later extract the frequency response by calculating the FFT of the Transmitted Signal.
double gaussianSource(double t, double tau = minTau, double a = 1) {
	//t is the running time 
	//t0 is the time offset
	//tau is the width of the pulse, the higher the wider
	double t0 = 6 * minTau;
	if (tau < minTau)
	{
		tau = minTau;
	}
	return a * exp(-pow((t - t0) / tau, 2));
}

//create update coefficients to normalize the fields during the update equations
std::vector<double> calculateUpdateCoefficients(const std::vector<double>& inputArray) {
	std::vector<double> updateCoefficients(inputArray.size(), 1.);
	const double denom = c0 * dt;
	#pragma omp parallel for
	for (int i = 0; i < inputArray.size(); i++) {
		updateCoefficients[i] = denom / (inputArray[i]);
	}
	return updateCoefficients;
}

std::vector<double> precomputeSource()
{
	//precompute the source terms to speed up the main loop ever so slightly
	std::vector<double> sourceArray(steps);
	#pragma omp parallel for
	for (size_t i = 0; i < steps; i++)
	{
		sourceArray[i] = gaussianSource(i * dt, minTau, 1);
	}
	return sourceArray;
}

//calculate the index of refraction at each grid cell. 
std::vector<double> calculateRefractiveIndexes(const std::vector<double>& permittivity, const std::vector<double>& permeability, bool approximatePermeabilityAsOne)
{
	std::vector<double> indexes(permittivity.size());
	if (approximatePermeabilityAsOne)
	{
		#pragma omp parallel for
		for (size_t i = 0; i < permittivity.size(); i++)
		{
			indexes[i] = sqrt(permittivity[i]);
		}
	}
	else
	{
		#pragma omp parallel for
		for (size_t i = 0; i < permittivity.size(); i++)
		{
			indexes[i] = sqrt(permittivity[i] * permeability[i]);
		}
	}
	return indexes;
}

//upscale an array to a higher resolution while keeping intervals the same proportion
std::vector<double> upscaleVector(const std::vector<double>& input, size_t newSize)
{
	std::vector<double> upscaled(newSize);
	const size_t inputLength = input.size();
	#pragma omp parallel for
	for (size_t i = 0; i < input.size(); i++)
	{
		size_t indexInInput = round(newSize / i) * inputLength;
		upscaled[i] = input[indexInInput];
	}
}

template<typename T>
double findExtremum(const std::vector<T>& data, bool maximum, bool returnIndex)
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