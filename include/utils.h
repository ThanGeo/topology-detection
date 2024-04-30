#ifndef MAIN_UTILS_H
#define MAIN_UTILS_H

#include <sstream>
#include <vector>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>

#include "Rasterizer.h"
#include "def.h"

#define RED "\e[0;31m"
#define GREEN "\e[0;32m"
#define NC "\e[0m"

std::string getDatasetNameFromPath(std::string &datasetPath);

bool verifydatasetNickname(std::string &datasetPath);
bool verifyFileExists(std::string &filePath);
bool verifyDirectoryExists(std::string &directoryPath);

/**
 * terminal loggers
*/
void log_err(char* errorText);
void log_err(std::string errorText);
void log_err_w_text(std::string errorText, std::string txt);
void log_task(std::string text);
void log_task_w_text(std::string text, std::string txt);
void success_text(char* text);
void success_text_with_number(char* text, int number);
void success_text_with_text_and_number(char* text, std::string &extra_text, int number);
void success_text_with_double_and_unit(char *text, double val, std::string unit);
void success_text_with_percentage(char *text, double perc);
void success_text_with_time(char* text, double seconds);

/**
 * builds the APRIL interval files path
 * type: ALL or FULL
*/
std::string getAPRILFilePath(std::string &datasetName, std::string &approximationDir, std::string type, int N, int partitions, int compression);

/**
 * terminal printers
*/
void printPolygon(rasterizerlib::polygon2d &polygon, uint recID);
void printCellsOfIntervals(std::vector<uint> &intervalList, int cellsPerDim);

#endif