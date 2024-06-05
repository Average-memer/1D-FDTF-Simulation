#pragma once
#include <vector>
#include <string>

template <typename T>
void printVectorNonConst(std::vector<T>& vec, std::string seperator, bool newLine);
void printVectorDouble(std::vector<double>& vec, std::string seperator = ",", bool newLine = true);
void printVectorDoubleConst(const std::vector<double>& vec, std::string seperator = ",", bool newLine = true);

void printInformation();
void printProgress(int iteration);
void saveToFile(const std::vector<double>& vecToSave, std::string filename, std::string seperator);
uint64_t timeSinceEpochMillisec();