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
#include <bitset>

#include "containers.h"
#include "join_geometry_refinement.h"


#include "../libvbyte-master/vbyte.h"
#include "../libvbyte-master/varintdecode.h"

using namespace std;



/*
*-------------------------------------------------------
*
*     COMMON
*       
*
*-------------------------------------------------------
*/


int compareIntervals(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2);
