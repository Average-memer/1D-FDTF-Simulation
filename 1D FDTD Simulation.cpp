#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <chrono>

#include "utility.h"
#include "constants.h"
#include "mathematics.h"

//boundary condition variables
double H2, H1 = 0;
double E2, E1 = 0;

int main()
{
	//open files
	std::ofstream EField("./E.txt");
	std::ofstream HField("./H.txt");

	std::ostream_iterator<double> EIterator(EField, ",");
	std::ostream_iterator<double> HIterator(HField, ",");

	//initialize E and H field vectors. These store the field strength at each point in space, given by the unit cell size
	std::vector<double> Ey(cellCount, 0.);
	std::vector<double> Hx(cellCount, 0.);

	printInformation();
	//main loop
	uint64_t loopStart = timeSinceEpochMillisec();
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

		//only save or print progress every so often to conserve space and improve performance
		if (i % saveInterval == 0)
		{
			if (saveResults)
			{
				saveToFile(Ey, "E.txt", ",");
				saveToFile(Hx, "H.txt", ",");
			}
			printProgress(i);
		}
		//add the source 
		Ey[sourceInjectionPoint] += source[i];
	}
	uint64_t loopDuration = timeSinceEpochMillisec() - loopStart;
	std::cout << "\nCalculation took " << (loopDuration / 1000.0) << "s to complete." << std::endl;
	EField.close();
	HField.close();
	return 0;
}
