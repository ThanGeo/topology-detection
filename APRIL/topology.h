#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <stdio.h>
#include <stdlib.h>
#include "omp.h"


#include "containers.h"
#include "join.h"

using namespace std;


enum {
	REF_MARKED,
	REF_UNMARKED,
	REF_BLOCKED
};

extern int findRelationUsingAPRIL(Polygon *polA, Polygon *polB, bool &markedForEqual);
extern int checkAllRelations_old(Polygon *polA, Polygon *polB);

// topology check functions
int checkTopologyDisjoint(Polygon *polA, Polygon *polB);
int checkTopologyEqual(Polygon *polA, Polygon *polB);
int checkTopologyRwithinS(Polygon *polA, Polygon *polB);
int checkTopologyRContainsS(Polygon *polA, Polygon *polB);
int checkTopologyRCoveredByS(Polygon *polA, Polygon *polB);
int checkTopologyRCoversS(Polygon *polA, Polygon *polB);
int checkTopologyCrosses(Polygon *polA, Polygon *polB);
int checkTopologyMeet(Polygon *polA, Polygon *polB);

int checkTopologyIntersects(Polygon *polA, Polygon *polB);   //todo: use normal april here instead of this?


#endif