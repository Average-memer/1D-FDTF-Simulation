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


//print information about the simulation domain, runtime, and constants
void printInformation() {
	cout << "simulation domain size: " << domainSize << "m, divided into " << cellCount << " cells." << endl;
	cout << "stepsize: " << ds << "m, timestep: " << dt << "s." << endl;
	cout << "max. frequency: " << maxF << "Hz, min Wavelength: " << lambdaMin << "m." << endl;
	cout << "Number of iterations: " << steps << ", approximate filesize: " << steps * cellCount * 8 *1/1000 << "kilobytes." << endl;
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
