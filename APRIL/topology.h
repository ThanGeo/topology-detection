#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <stdio.h>
#include <stdlib.h>
#include "omp.h"


#include "containers.h"
#include "join.h"

using namespace std;

extern int (*checkPredicateFunctionPtr)(int AA, int AF, int FA, int FF);

extern int checkAllRelations(Polygon *polA, Polygon *polB);
extern int checkSpecificRelation(Polygon *polA, Polygon *polB);

// topology check functions
int checkTopologyDisjoint(int AA, int AF, int FA, int FF);
int checkTopologyEqual(int AA, int AF, int FA, int FF);
int checkTopologyRContainedInS(int AA, int AF, int FA, int FF);
int checkTopologyRContainsS(int AA, int AF, int FA, int FF);
int checkTopologyRCoveredByS(int AA, int AF, int FA, int FF);
int checkTopologyRCoversS(int AA, int AF, int FA, int FF);
int checkTopologyCrosses(int AA, int AF, int FA, int FF);
int checkTopologyMeet(int AA, int AF, int FA, int FF);

int checkTopologyIntersects(Polygon *polA, Polygon *polB);   //todo: use normal april here instead of this?


#endif