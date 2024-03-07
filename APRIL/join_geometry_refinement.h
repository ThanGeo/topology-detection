#pragma once

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cmath>
#include <ctime>

#include <boost/geometry/algorithms/detail/relate/de9im.hpp> 

#include "containers.h"
#include "dataset_data.h"

using namespace std;


/*
*-------------------------------------------------------
*
*     GEOMETRY RETRIEVAL
*       
*
*-------------------------------------------------------
*/

int refinementWithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS);


int refinement_DE9IM_WithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS);

int refinement_DE9IM_WithIDs_optimized(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS);