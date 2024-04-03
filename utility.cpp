#include "utility.h"
#include <cmath>
#include <vector>
#include <string>

#include "constants.h"

//a simple gaussian source. We later extract the frequency response by calculating the FFT of the Transmitted Signal.
float gaussianSource(float t, float t0 = 0, float tau = minTau, float a = 1) {
	//t is the running time 
	//t0 is the time offset
	//tau is the width of the pulse, the higher the wider
	if (t0 < 3 * minTau)
	{
		float t0 = 3 * minTau;
	}
	if (tau < minTau)
	{
		tau = minTau;
	}
	return a * exp(-pow((t - t0) / tau, 2));
}

template <typename T>
void printVector(std::vector<T> vec, std::string seperator, bool newLine) {
	for (auto i : vec) {
		std::cout << i << seperator;
	}
	if (newLine) {
		std::cout << std::endl;
	}
}

//create update coefficients to normalize the fields during the update equations
std::vector<float> calculateUpdateCoefficients(const std::vector<float>& inputArray) {
	std::vector<float> updateCoefficients(inputArray.size(), 1.);
	const float denom = c0 * dt;
	for (int i = 0; i < inputArray.size(); i++) {
		updateCoefficients[i] = denom / inputArray[i];
	}
	return updateCoefficients;
}