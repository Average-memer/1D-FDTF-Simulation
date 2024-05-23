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
double H2, H1 = 0;
double E2, E1 = 0;

// TODO: fix exploding values.

//set epsilon and mu to be 1, aka vacuum.
const std::vector<double> permeability(cellCount, 1.);
const std::vector<double> permissivity(cellCount, 1.);

//modified E and H coefficients
const std::vector<double> mHx = calculateUpdateCoefficients(permissivity);
const std::vector<double> mEy = calculateUpdateCoefficients(permeability);

//E and H field vectors. These store the field strength at each point in space, given by the unit cell size
std::vector<double> Ey(cellCount, 0.);
std::vector<double> Hx(cellCount, 0.);

//precompute source 
std::vector<double> source = precomputeSource();

int main()
{
	//open files
	std::ofstream EField("./E.txt");
	std::ofstream HField("./H.txt");

	std::ostream_iterator<double> EIterator(EField, ",");
	std::ostream_iterator<double> HIterator(HField, ",");

	printInformation();
	//main loop
	for (int i = 0; i < steps; i++)
	{
		//update H from E	
		H2 = H1; H1 = Hx[0];
		for (size_t j = 0; j < cellCount - 1; j++)
		{
			Hx[j] = Hx[j] + mHx[j] * ((Ey[j + 1] - Ey[j]) / ds);
		}
		//boundary condition
		Hx[maxArrayIndex] = Hx[maxArrayIndex] + mHx[maxArrayIndex] * ((E2 - Ey[maxArrayIndex]) / ds);
		
		//update E from H
		Ey[0] = Ey[0] + mEy[0] * ((Hx[0] - H2)/ ds);
		E2 = E1; E1 = Ey[maxArrayIndex];
		for (size_t k = 1; k < cellCount; k++)
		{
			Ey[k] = Ey[k] + mEy[k] * ((Hx[k] - Hx[k - 1]) / ds);
		}

		//only save every so often to conserve space
		if (i % saveInterval == 0)
		{
			saveToFile(Ey, "E.txt", ",");
			saveToFile(Hx, "H.txt", ",");
		}
		//add the source 
		Ey[sourceInjectionPoint] += source[i];

		if (i % saveInterval == 0)
		{
			printProgress(i);
		}
	}
	EField.close();
	HField.close();
	return 0;
}
