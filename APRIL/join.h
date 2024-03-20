#ifndef JOIN_H
#define JOIN_H

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




enum{
    R_INTERVALS_CONTAINED_IN_S_INTERVALS,
    S_INTERVALS_CONTAINED_IN_R_INTERVALS,
    NO_CONTAINMENT,
    SYMMETRICAL_CONTAINMENT,
    MATCH,
};

/**
 * compares intervals of R and S 
 * checks for containment (returns 1 if intervals of R completely contained in intervals of S) 
 * while also flagging for intersection (bool intersect)
*/
int compareIntervalsEfficient(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2, bool &intersects);

/**
 * compares intervals 
 * checks for containment (returns 1 if intervals1 are completely contained in intervals2) 
 * while also flagging for intersection (bool intersect)
*/
int compareIntervalsHybrid(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2, bool &intersects);

/**
 * compares intervals of R and S and returns 1 if at least one intersection is found
 * otherwise returns 0
*/
int compareIntervalsForIntersection(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2);

/**
 * 
 * 
*/
int compareIntervalsSymmetricalContainmentAndIntersection(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2, bool &intersects);

/**
 * compares intervals of R with intervals of S and returns 1 if:
 *  ALL of intervals of R are contained completely in the intervals of S
 *  otherwise returns 0
*/
int compareIntervalsForContainment(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2);

/**
 * compares intervals of R with intervals of S and returns 1 if they match 100%
 * otherwise returns 0
*/
int compareIntervalsForMatch(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2);



#endif