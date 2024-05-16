#include "utility.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "constants.h"

using std::cout;
using std::endl;
//a simple gaussian source. We later extract the frequency response by calculating the FFT of the Transmitted Signal.
double gaussianSource(double t, double t0 = 0, double tau = minTau, double a = 1) {
	//t is the running time 
	//t0 is the time offset
	//tau is the width of the pulse, the higher the wider
	if (t0 < 3 * minTau)
	{
		double t0 = 3 * minTau;
	}
	if (tau < minTau)
	{
		tau = minTau;
	}
	return a * exp(-pow((t - t0) / tau, 2));
}

template <typename T>
void printVector(std::vector<T>& vec, std::string seperator, bool newLine){
	for (auto i : vec) {
		std::cout << i << seperator;
	}
	if (newLine) {
		std::cout << std::endl;
	}
}


void printVectorDouble(std::vector<double>& vec, std::string seperator, bool newLine){
for (auto i : vec) {
		std::cout << i << seperator;
	}
	if (newLine) {
		std::cout << std::endl;
	}
}

void printVectorDouble(const std::vector<double>& vec, std::string seperator, bool newLine){
for (auto i : vec) {
	std::cout << i << seperator;
}
if (newLine) {
	std::cout << std::endl;
}
}

//create update coefficients to normalize the fields during the update equations
std::vector<double> calculateUpdateCoefficients(const std::vector<double>& inputArray) {
	std::vector<double> updateCoefficients(inputArray.size(), 1.);
	const double denom = c0 * dt;
	std::cout << denom << std::endl;
	for (int i = 0; i < inputArray.size(); i++) {
		updateCoefficients[i] = denom / inputArray[i];
	}
	return updateCoefficients;
}

//print information about the simulation domain, runtime, and constants
void printInformation() {
	cout << "simulation domain size: " << domainSize << "m, divided into " << cellCount << " cells." << endl;
	cout << "stepsize: " << ds << "m, timestep: " << dt << "s." << endl;
	cout << "max. frequency: " << maxF << "Hz, min Wavelength: " << lambdaMin << "m." << endl;
	cout << "approximate number of iterations: " << steps * cellCount * 2 << ", approximate filesize: " << steps * cellCount * 8 *1/1000 << "kilobytes." << endl;
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