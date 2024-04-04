#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>

#include "utility.h"
#include "constants.h"

//boundary condition variables
float H3, H2, H1 = 0;
float E3, E2, E1 = 0;

//set epsilon and mu to be 1, aka vacuum.
const std::vector<float> permeability(cellCount, 1.);
const std::vector<float> permissivity(cellCount, 1.);

//modified E and H coefficients
const std::vector<float> mEy = calculateUpdateCoefficients(permeability);
const std::vector<float> mHx = calculateUpdateCoefficients(permissivity);

//E and H field vectors. These store the field strength at each point in space, given by the unit cell size
std::vector<float> Ey(cellCount, 0.);
std::vector<float> Hx(cellCount, 0.);

int main()
{
	//open files
	std::ofstream EField("./E.txt");
	std::ofstream HField("./H.txt");

	std::ostream_iterator<float> EIterator(EField, ",");
	std::ostream_iterator<float> HIterator(HField, ",");

	printInformation();
	//main loop
	for (int i = 0; i < steps; i++)
	{
		//update H from E		
		for (size_t j = 0; j < cellCount - 1; j++)
		{
			Hx[j] = Hx[j] + mHx[j] * ((Ey[j + 1] - Ey[j]) / ds);
		}
		//boundary condition
		Hx[maxArrayIndex] = Hx[maxArrayIndex] + mHx[maxArrayIndex] * ((E3 - Ey[maxArrayIndex]) / ds);
		H3 = H2; H2 = H1; H1 = Hx[0];

		//update E from H
		Ey[0] = Ey[0] + mEy[0] * ((Hx[0] - H3)/ ds);
		for (size_t k = 1; k < cellCount; k++)
		{
			Ey[k] = Ey[k] + mEy[k] * ((Hx[k] - Hx[k - 1]) / ds);
		}
		E3 = E2; E2 = E1; E1 = Ey[maxArrayIndex];

		if (i % 2 == 0)
		{
			std::copy(Ey.begin(), Ey.end(), EIterator);
			EField << "\n";
			std::copy(Hx.begin(), Hx.end(), HIterator);
			HField << "\n";
		}
		//add the source 
		Ey[10] += gaussianSource(i * dt, 0, minTau, 1);

		if (i % 10 == 0)
		{
			printProgress(i);
		}
	}
	EField.close();
	HField.close();
	return 0;
}
