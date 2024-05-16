#pragma once
#include <vector>
#include <string>
#include "constants.h"

//
template <typename T>
void printVector(std::vector<T>& vec, std::string seperator, bool newLine);
void printVectorDouble(std::vector<double>& vec, std::string seperator, bool newLine);
void printVectorDouble(const std::vector<double>& vec, std::string seperator, bool newLine);

std::vector<double> calculateUpdateCoefficients(const std::vector<double>& inputArray);
double gaussianSource(double t, double t0, double tau, double a);
void printInformation();
void printProgress(int iteration);
void saveToFile(const std::vector<double>& vecToSave, std::string filename, std::string seperator);
