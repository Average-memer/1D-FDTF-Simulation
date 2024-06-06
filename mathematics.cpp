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
std::vector<double> upscaleVector(const std::vector<double>& input, const int newSize)
{
	const int inputLength = input.size();
	if (newSize < inputLength) {
		std::cout << "New length is shorter than current length, will not scale" << std::endl;
		return input;
	}
	std::vector<double> upscaled(newSize);

	for (int i = 1; i < newSize; i++)
	{
		int indexInInput = round((i / newSize) * inputLength);
		std::cout << "i: " << i << ", indexInInput: " << indexInInput << " ,newSize: " << newSize << std::endl;
		upscaled[i] = input[indexInInput];
	}
	return upscaled;
}

std::vector<double> precomputeElectricSource()
{
	//precompute the source terms to speed up the main loop ever so slightly
	std::vector<double> sourceArray(steps);
	for (size_t i = 0; i < steps; i++)
	{
		sourceArray[i] = gaussianSource(i * dt, minTau, 1);
	}
	return sourceArray;
}

std::vector<double> precomputeMagneticSource(double epsilon_r_src, double mu_r_src)
{
	//see page 12 of learning from 1D-FDTD
	const double amplitude = -sqrt(epsilon_r_src / mu_r_src);
	//delay through one half grid cell + half a time step. The delay is needed because Ey and Hx exist at different points in space.
	const double delay = (sqrt(epsilon_r_src * mu_r_src) * ds) / (2 * c0) + dt / 2;
	std::vector<double> sourceArray(steps);
	for (size_t i = 0; i < steps; i++)
	{
		sourceArray[i] = amplitude * gaussianSource(i * dt + delay, minTau, 1);
	}
	return sourceArray;
}