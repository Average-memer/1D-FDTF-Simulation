#pragma once
#include <vector>
#include <string>
#include "constants.h"

template <typename T>
void printVector(std::vector<T> vec, std::string seperator, bool newLine);
std::vector<float> calculateUpdateCoefficients(const std::vector<float>& inputArray);
float gaussianSource(float t, float t0, float tau, float a);
void printInformation();
void printProgress(int iteration);
