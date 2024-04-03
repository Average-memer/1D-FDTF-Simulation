#include <vector>
#include <string>
#include <stdexcept>
#include "matplotlibcpp.h"
#include <cmath>

namespace plt = matplotlibcpp;

//physical constants
const float epsilon0 = 8.854187817e-12; //electric permissivity of vacuum
const float mu0 = 1.25663706e-6; //magnetic permeability of vacuum
const float normalisation = sqrt(mu0 / epsilon0); //normalize E and H to have the same order of magnitude
const float c0 = 299792458; //speed of light in vacuum

//setup the simulation domain
const int steps = 100;
const float physize = 5; //5cm
const int resolution = 100; //100 unit cells
const float ds = physize / resolution; //size of each cell
const float dt = 1e-9;

//set epsilon and mu to be 1, aka vacuum.
const std::vector<float> permeability(resolution, 1.);
const std::vector<float> permissivity(resolution, 1.);

template <typename T>
void printVector(std::vector<T> vec, std::string seperator, bool newLine) {
	for (auto i : vec) {
		std::cout << i << seperator;
	}
	if (newLine) {
		std::cout << std::endl;
	}
}

float gaussianSource(float t, float f0, float t0, float a) {
	if (a <= 0)
	{
		throw std::invalid_argument("a has to be bigger than 0!");
	}
	return sin(f0 * t) * exp(-a * pow((t - t0), 2));
}

//create update coefficients to normalize the fields during the update equations
std::vector<float> calculateUpdateCoefficients(const std::vector<float>& inputArray) {
	std::vector<float> updateCoefficients(inputArray.size(), 1.);
	const float denom = c0 * dt;
	for (int i = 0; i < inputArray.size(); i++) {
		updateCoefficients[i] = (denom / inputArray[i]);
	}
	return updateCoefficients;
}

//modified E and H coefficients
const std::vector<float> mEy = calculateUpdateCoefficients(permeability);
const std::vector<float> mHx = calculateUpdateCoefficients(permissivity);

//E and H field vectors. These store the field strength at each point in space, given by the unit cell size
std::vector<float> Ey(resolution, 0.);
std::vector<float> Hx(resolution, 0.);

int main()
{
	//main loop
	plt::ion();
	auto fig = plt::figure();

	for (int i = 0; i < steps; i++)
	{
		//update H from E
		Ey[10] += gaussianSource(i * dt, 1e6, 0, 1);
		for (size_t j = 0; j < resolution; j++)
		{
			Hx[j] = Hx[j] + mHx[j] * ((Ey[j + 1] - Ey[j]) / ds);
		}

		//update E from H
		for (size_t k = 0; k < resolution; k++)
		{
			Ey[k] = Ey[k] + mEy[k] * ((Hx[k] - Hx[k - 1]) / ds);
		}
		
		if (i % 5 == 0)
		{
			printVector(Ey, ", ", true);
			printVector(Hx, ", ", true);
		}
		
		
	}
	plt::close();
	return 0;
}
