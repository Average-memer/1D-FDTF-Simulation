#pragma once
#include <vector>
#include <string>
#include "constants.h"

template <typename T>
void printVectorNonConst(std::vector<T>& vec, std::string seperator, bool newLine);
void printVectorDouble(std::vector<double>& vec, std::string seperator = ",", bool newLine = true);
void printVectorDoubleConst(const std::vector<double>& vec, std::string seperator = ",", bool newLine = true);

template <typename T>
double findExtremum(const std::vector<T>& data, bool maximum = true, bool returnIndex = false);

std::vector<double> calculateUpdateCoefficients(const std::vector<double>& inputArray);
std::vector<double> precomputeElectricSource();
std::vector<double> precomputeMagneticSource(double epsilon_r_src = 1, double mu_r_src = 1);
std::vector<double> calculateRefractiveIndexes(const std::vector<double>& permittivity, const std::vector<double>& permeability, bool approximatePermeabilityAsOne);
double gaussianSource(double t, double tau, double a);
void printInformation();
void printProgress(int iteration);
void saveToFile(const std::vector<double>& vecToSave, std::string filename, std::string seperator);
uint64_t timeSinceEpochMillisec();