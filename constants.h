#pragma once
//physical inline ants
inline  double epsilon0 = 8.854187817e-12; //electric permissivity of vacuum
inline  double mu0 = 1.25663706e-6; //magnetic permeability of vacuum
inline  double normalisation = sqrt(mu0 / epsilon0); //normalize E and H to have the same order of magnitude
inline  int c0 = 299792458; //speed of light in vacuum
inline  double PI = atan(1) * 4; 

//setup the simulation properties
inline const double nmax = 3; //maximum index of refraction
inline const int steps = 200; //amount of timesteps that will be run
inline const double domainSize = 40; //size of the simulation domain
inline const double minFeatureSize = 5e-1; //5mm, minimum physical feature size we can resolve
inline const int NDres = 3; //amount of cells to resolve the smallest feature with
inline const int timeStepsPerPulse = 10; //amount of timesteps we want our Pulse to be resolved by?
inline const double maxF = 1e6; //maximum frequency we care about

//derive limits from the desired values
inline double propagationTime = domainSize / c0; //worst case time for a wave to propagate from one side to the other
inline double lambdaMin = c0 / maxF; //smallest wavelength we care about
inline double dsMin = fmin((lambdaMin * timeStepsPerPulse) / nmax, minFeatureSize / NDres); //the minimum spacing is given by the min of either the highest frequency or the smallest structure
inline double minTau = 1 / (maxF * PI); //minimum width of all gaussian sources
inline double dt = minTau / timeStepsPerPulse; //amount of time that passes between timesteps
inline size_t cellCount = ceil(domainSize / dsMin); //number of cells we divide the physical space into
inline double ds = cellCount / domainSize; //the spacing such that the domain is made up of an integer number of cells smaller or equal to the minimum size.
inline int maxArrayIndex = cellCount - 1; //index of the final element of the grid

//program settings
inline int saveInterval = 2; //Save fields every ith iteration
inline int saveResolution = 1; //Only save every nth cell of the fields