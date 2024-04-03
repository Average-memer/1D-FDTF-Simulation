#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>

//physical constants
const float epsilon0 = 8.854187817e-12; //electric permissivity of vacuum
const float mu0 = 1.25663706e-6; //magnetic permeability of vacuum
const float normalisation = sqrt(mu0 / epsilon0); //normalize E and H to have the same order of magnitude
const int c0 = 299792458; //speed of light in vacuum
const float PI = atan(1) * 4;

//setup the simulation domain
const int steps = 100; //amount of timesteps that will be run
const float physize = 5.0; //5cm, size of the simulation domain
const int resolution = 100; //100 unit cells
const int maxIndex = resolution - 1; //index of the final element of the grid
const int timeStepsPerPulse = 10; //amount of timesteps we want our Pulse to be resolved by?
const float ds = physize / resolution; //size of each cell
const float fmax = 1e6; //maximum frequency we care about
const float minTau = 1 / (fmax * PI); //minimum width of all gaussian sources
const float dt = minTau / timeStepsPerPulse;

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

float gaussianSource(float t, float t0 = 0, float tau = minTau) {
	//t is the running time 
	//t0 is the time offset
	//tau is the width of the pulse, the higher the wider
	if (t0 < 3*minTau)
	{
		float t0 = 3*minTau;
	}	
	if (tau < minTau)
	{
		tau = minTau;
	}
	return exp(-pow((t - t0) / tau, 2));
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
	//open files
	std::ofstream EField("./E.txt");
	std::ofstream HField("./H.txt");

	std::ostream_iterator<float> EIterator(EField, ",");
	std::ostream_iterator<float> HIterator(HField, ",");
	//main loop
	for (int i = 0; i < steps; i++)
	{
		//update H from E		
		for (size_t j = 0; j < resolution - 1; j++)
		{
			Hx[j] = Hx[j] + mHx[j] * (Ey[j + 1] - Ey[j]) / ds;
		}
		Hx[maxIndex] = Hx[maxIndex] + mHx[maxIndex] * (-Ey[maxIndex]) / ds;

		//update E from H
		Ey[0] = Ey[0] + mEy[0] * Hx[0] / ds;
		for (size_t k = 1; k < resolution; k++)
		{
			Ey[k] = Ey[k] + mEy[k] * (Hx[k] - Hx[k - 1]) / ds;
		}
		
		if (i % 5 == 0)
		{
			//printVector(Ey, ", ", true);
			//printVector(Hx, ", ", true);
			std::copy(Ey.begin(), Ey.end(), EIterator);
			EField << "\n";
			std::copy(Hx.begin(), Hx.end(), HIterator);
			HField << "\n";
		}
		Ey[10] += gaussianSource(i * dt);
	}
	EField.close();
	HField.close();
	return 0;
}
