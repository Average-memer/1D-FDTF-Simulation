#pragma once
#ifndef UTITILY_H
#define UTILITY_H

template <typename T>
void printVector(std::vector<T> vec, std::string seperator, bool newLine);
std::vector<float> calculateUpdateCoefficients(const std::vector<float>& inputArray);
float gaussianSource(float t, float t0 = 0, float tau = minTau, float a = 1);

#endif // !UTITILY_H