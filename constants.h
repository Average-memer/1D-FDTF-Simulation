#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

//physical constants
const float epsilon0 = 8.854187817e-12; //electric permissivity of vacuum
const float mu0 = 1.25663706e-6; //magnetic permeability of vacuum
const float normalisation = sqrt(mu0 / epsilon0); //normalize E and H to have the same order of magnitude
const int c0 = 299792458; //speed of light in vacuum
const float PI = atan(1) * 4;

//setup the simulation properties
const float nmax = 3; //maximum index of refraction
const int steps = 400; //amount of timesteps that will be run
const float domainSize = 0.05; //5cm, size of the simulation domain
const float minFeatureSize = 5e-3; //5mm, minimum physical feature size we can resolve
const int NDres = 4; //amount of cells to resolve the smallest feature with
const int timeStepsPerPulse = 10; //amount of timesteps we want our Pulse to be resolved by?
const float propagationTime = domainSize / c0; //worst case time for a wave to propagate from one side to the other
const float maxF = 1e6; //maximum frequency we care about
const float lambdaMin = 1 / maxF; //smallest wavelength we care about
const float dsMin = fmin((lambdaMin * timeStepsPerPulse) / nmax, minFeatureSize / NDres); //the minimum spacing is given by the min of either the highest frequency or the smallest structure
const float minTau = 1 / (maxF * PI); //minimum width of all gaussian sources
const float dt = minTau / (10 * timeStepsPerPulse); //amount of time that passes between timesteps
const int cellCount = ceil(domainSize / dsMin); //number of cells we divide the physical space into
const float ds = cellCount / domainSize; //the spacing such that the domain is made up of an integer number of cells smaller or equal to the minimum size.
const int maxArrayIndex = cellCount - 1; //index of the final element of the grid

#endif