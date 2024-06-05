#pragma once
#include "utility.h"
#include "mathematics.h"
//physical constants
inline  double epsilon0 = 8.854187817e-12; //electric permissivity of vacuum
inline  double mu0 = 1.25663706e-6; //magnetic permeability of vacuum
inline  double normalisation = sqrt(mu0 / epsilon0); //normalize E and H to have the same order of magnitude
inline  int c0 = 299792458; //speed of light in vacuum
inline  double PI = atan(1) * 4; 

//TODO: fix this mess v

//setup the simulation properties
//set epsilon and mu to be 1, aka vacuum. This is where you would include your device
const inline int testLength = 200;
const std::vector<double> devicePermeability(testLength, 1.);
const std::vector<double> devicePermittivity(testLength, 1.);

//we can do these before scaling because the values themself don't change
const std::vector<double> indexOfRefraction = calculateRefractiveIndexes(devicePermittivity, devicePermeability, true); //IoR at every point of my device
inline const double nmax = findExtremum(indexOfRefraction, true, false); //maximum index of refraction
inline const double nmin = findExtremum(indexOfRefraction, false, false); //minimum index of refraction

//define minimum sizes
inline const double domainSize = 1; //size of the simulation domain in m. 
inline const double dmin = 2e-1; //minimum physical feature size we can resolve in m. This is the width of one element of the device arrays
inline const int Nd = 4; //amount of cells to resolve the smallest feature with
inline const double maxF = 5e9; //maximum frequency we care about
inline const int Nlambda = 10; //amount of cells to resolve the smallest wavelength

//derived values
//cell sizing and count
inline const double lambdaMin = c0 / (maxF * nmax); //smallest wavelength we care about
inline const double dsMin = fmin(lambdaMin / Nlambda, dmin / Nd); //the minimum spacing is given by the min of either the highest frequency or the smallest structure
inline const size_t cellCount = ceil(domainSize / dsMin); //number of cells needed to satisfy size requirements from frequency and spatial resolution. We will upscale the device to this size.
inline const double ds = domainSize / cellCount; //the spacing such that the domain is made up of an integer number of cells equal to the minimum size.

//time step and iteration count
inline const double minTau = 1 / (maxF * PI); //minimum width of all gaussian sources
inline const double propagationTime = (nmax * domainSize) / c0; //worst case time for a wave to propagate from one side to the other
inline const double simulationTime = 12 * minTau + 3 * propagationTime; //simulation time in seconds. Allows for five bounces of a wave plus the entire duration of the initial pulse, g(t)
inline const double dt = fmin(nmin * ds / (2 * c0), 1 / (PI * maxF)); //minimum timestep to satisfy the courant stability condition and the pulse duration
inline const size_t steps = ceil(simulationTime / dt); //how long the simulation will run, in iterations. 

inline const size_t sourceInjectionPoint = 20; //defines the TF/SF boundary. This is on the TF side.
inline const size_t maxArrayIndex = cellCount - 1; //index of the final element of the grid

//scale arrays to new cellcounts.
const std::vector<double> permeability = upscaleVector(devicePermeability, cellCount); //the actual array we want to run our simulation on.
const std::vector<double> permittivity = upscaleVector(devicePermittivity, cellCount); //the actual array we want to run our simulation on.
//modified E and H coefficients
const std::vector<double> mHx = calculateUpdateCoefficients(permittivity);
const std::vector<double> mEy = calculateUpdateCoefficients(permeability);

//precompute source 
inline const std::vector<double> source = precomputeSource();

//program settings
inline bool saveResults = false; //wether or not to save field properties during iteration
inline int saveInterval = 10; //Save fields every ith iteration
inline int saveResolution = 2; //Only save every nth cell of the fields