#pragma once
#include "utility.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

#include "constants.h"
#include "mathematics.h"

using std::cout;
using std::endl;

//print a non-const vector of any type
template <typename T>
void printVectorNonConst(std::vector<T>& vec, std::string seperator, bool newLine){
	for (auto i : vec) {
		std::cout << i << seperator;
	}
	if (newLine) {
		std::cout << std::endl;
	}
}

//print a vector of doubles
void printVectorDouble(std::vector<double>& vec, std::string seperator, bool newLine){
for (double i : vec) {
		std::cout << i << seperator;
	}
	if (newLine) {
		std::cout << std::endl;
	}
}

//print a const vector of doubles
void printVectorDoubleConst(const std::vector<double>& vec, std::string seperator, bool newLine){
for (double i : vec) {
	std::cout << i << seperator;
}
if (newLine) {
	std::cout << std::endl;
}
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
	const double delay = (sqrt(epsilon_r_src * mu_r_src) * ds) / (2*c0) + dt/2; 
	std::vector<double> sourceArray(steps);
	for (size_t i = 0; i < steps; i++)
	{
		sourceArray[i] = amplitude * gaussianSource(i * dt + delay, minTau, 1);
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

//print information about the simulation domain, runtime, and constants
void printInformation() {
	cout << "simulation domain size: " << domainSize << "m, divided into " << cellCount << " cells." << endl;
	cout << "stepsize: " << ds << "m, timestep: " << dt << "s." << endl;
	cout << "max. frequency: " << maxF << "Hz, min Wavelength: " << lambdaMin << "m." << endl;
	cout << "Number of iterations: " << steps << ", approximate filesize: " << (saveResults ? (steps / saveResolution) * (cellCount / saveInterval) * sizeof(double) * 2 * correctionFactor / 1000 : 0) << "kilobytes." << endl;
}

void printProgress(int iteration) {
	int barWidth = 70;

	cout << "[";
	double progress = iteration / double(steps);
	int pos = floor(barWidth * progress);
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) cout << "=";
		else if (i == pos) cout << ">";
		else cout << " ";
	}
	cout << "] " << int(progress * 100.0) << " %\r";
	cout.flush();
}

//save the entered vector to a file
void saveToFile(const std::vector<double>& vecToSave, std::string filename, std::string seperator = ",") {
	std::ofstream FILE("./"+filename, std::ios::app);

	for (size_t index = 0; index < vecToSave.size(); index = index + saveResolution)
	{
		if (index <= vecToSave.size())
		{
			FILE << vecToSave.at(index);
			//don't write a comma on the final entry
			if (index + saveResolution < vecToSave.size())
			{
				FILE << seperator;
			}
		}
		else
		{
			FILE << vecToSave.back();
			break;
		}
	}
	FILE << "\n";
	FILE.close();
}

//return the amount of milliseconds that passed since the UNIX-epoch.
uint64_t timeSinceEpochMillisec()
{
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}
