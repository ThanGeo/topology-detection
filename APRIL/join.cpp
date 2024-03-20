#include "join.h"


/*
*-------------------------------------------------------
*
*     COMMON
*       
*
*-------------------------------------------------------
*/

void nextinterval(CONTAINER *ar, uint *offset, ID *st, ID *end)
{
	
	(*st) = *(ID *)(ar+(*offset));
    (*offset)+=1;
    (*end) = *(ID *)(ar+(*offset));
    (*offset)+=1;

}

int compareIntervalsForIntersection(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2)
{
    //they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	
	// ID st1,st2,end1,end2;
	uint cur1=0;
	uint cur2=0;
	
	auto st1 = ar1.begin();
	auto end1 = ar1.begin() + 1;
    cur1++;
    
	auto st2 = ar2.begin();
	auto end2 = ar2.begin() + 1;
    cur2++;

    do {
		if (*st1<=*st2)
		{
			if (*end1>*st2) // overlap, end1>=st2 if intervals are [s,e] and end1>st2 if intervals are [s,e)
			{
				//they overlap, return 1
				return 1;
			}	
			else
			{
				st1 += 2;
				end1 += 2;
				cur1++;
			}
		}
		else // st2<st1
		{
			if (*end2>*st1) // overlap, end2>=st1 if intervals are [s,e] and end2>st1 if intervals are [s,e)
			{

				//they overlap, return 1
				return 1;
			}
			else
			{
				st2 += 2;
				end2 += 2;
				cur2++;
			}
		}
	} while(cur1<=numintervals1 && cur2<=numintervals2);
	
	//no overlap, return 0
	return 0;
}
static int checkIntervalsForContainmentAndIntersection(uint st1, uint end1, uint st2, uint end2, 
							bool &currentRcontainedInS, bool &currentScontainedInR, bool &intersects) {
	if(st1 >= st2) {
		if (end1 <= end2) {
			// if the current interval of R is contained completely in interval of S
			currentRcontainedInS = true;
			intersects = true;
		}
		if (st1 < end2){
			// there is an intersection
			intersects = true;
		}
	}
	if (st2 >= st1) {
		if (end2 <= end1) {
			// if the current interval of S is contained completely in interval of R
			currentScontainedInR = true;
			intersects = true;
		}
		if (end1 > st2) {
			// there is an intersection
			intersects = true;
		}
	}
	if (currentRcontainedInS && currentScontainedInR) {
		return MATCH;
	} else if(currentRcontainedInS) {
		return R_INTERVALS_CONTAINED_IN_S_INTERVALS;
	} else if (currentScontainedInR) {
		return S_INTERVALS_CONTAINED_IN_R_INTERVALS;
	}
	return NO_CONTAINMENT;
}
int compareIntervalsSymmetricalContainmentAndIntersection(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2, bool &intersects) {
	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return NO_CONTAINMENT;
	}

	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	
	// ID st1,st2,end1,end2;
	uint cur1=0;
	uint cur2=0;
	
	auto st1 = ar1.begin();
	auto end1 = ar1.begin() + 1;
    cur1++;
    
	auto st2 = ar2.begin();
	auto end2 = ar2.begin() + 1;
    cur2++;

	bool RinS = true;
	bool SinR = true;

	bool currentRintervalIsContained = false;
	bool currentSintervalIsContained = false;

    do {
		// printf("Intervals:\n");
		// printf("R[%d]: [%d,%d): containing interval %d\n", cur1, *st1, *end1, currentRintervalIsContained);
		// printf("S[%d]: [%d,%d): containing interval %d\n", cur2, *st2, *end2, currentSintervalIsContained);

		if (*st1 >= *st2) {
			if (*end1 <= *end2) {
				// R is contained
				currentRintervalIsContained = true;
				intersects = true;
			} else {
				if (currentRintervalIsContained){
					// R not in S
					RinS = false;
				}
			}
			if (*end2 > *st1) {
				// intersection
				intersects = true;
			}
			
		} else {
			if (!currentRintervalIsContained && cur2 >= numintervals2){
				// R not in S
				RinS = false;
			}
		}

		if (*st2 >= *st1) {
			if (*end2 <= *end1) {
				// S is contained
				currentSintervalIsContained = true;
				intersects = true;
			} else {
				if (currentSintervalIsContained){
					// S not in R
					SinR = false;
				}
			}
			if (*end1 > *st2) {
				// intersection
				intersects = true;
			}
		} else {
			if (!currentSintervalIsContained && cur1 >= numintervals1){
				// S not in R
				SinR = false;
			}
		}

		if (!RinS && !SinR && intersects) {
			// early break
			break;
		}

		if (*end1 < *end2) {
			st1 += 2;
			end1 += 2;
			cur1++;
			if (cur1 <= numintervals1) {
				currentRintervalIsContained = false;
			} else {
				break;
			}
		} else {
			st2 += 2;
			end2 += 2;
			cur2++;
			if (cur2 <= numintervals2) {
				currentSintervalIsContained = false;
			} else {
				break;
			}
		}

	} while(cur1<=numintervals1 && cur2<=numintervals2);

	// symmetrical containment
	if (RinS && SinR) {
		if(cur1 >= numintervals1 && cur2 >= numintervals2 && currentRintervalIsContained && currentSintervalIsContained) {
			return SYMMETRICAL_CONTAINMENT;
		}
	}
	// R in S
	if (cur1 >= numintervals1) {
		if (RinS && currentRintervalIsContained) {
			return R_INTERVALS_CONTAINED_IN_S_INTERVALS;
		}
	}
	// S in R
	if (cur2 >= numintervals2) {
		if (SinR && currentSintervalIsContained) {
			return S_INTERVALS_CONTAINED_IN_R_INTERVALS;
		}
	}
	// no containment
	return NO_CONTAINMENT;
}


int compareIntervalsForContainment(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2){
	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	
	// ID st1,st2,end1,end2;
	uint cur1=0;
	uint cur2=0;
	
	auto st1 = ar1.begin();
	auto end1 = ar1.begin() + 1;
    cur1++;
    
	auto st2 = ar2.begin();
	auto end2 = ar2.begin() + 1;
    cur2++;

    bool intervalRcontained = false;

    do {
    	//check if the current interval of R is contained completely in interval of S
    	if(*st1 >= *st2 && *end1 <= *end2){
    		intervalRcontained = true;
    	}

		if (*end1<=*end2)
		{
			if(!intervalRcontained){
				//we are skipping this interval because it was never contained, so return false (not within)
				return 0;
			}
			st1 += 2;
			end1 += 2;
			cur1++;
			intervalRcontained = false;
		}
		else 
		{
			st2 += 2;
			end2 += 2;
			cur2++;
		}
	} while(cur1<=numintervals1 && cur2<=numintervals2);
			
	//if we didnt check all of the R intervals
	if(cur1 <= numintervals1){	
		return 0;
	}
	//all intervals R were contained
	return 1;
}

int compareIntervalsForMatch(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2){
	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	
	// ID st1,st2,end1,end2;
	uint cur1=0;
	uint cur2=0;
	
	auto st1 = ar1.begin();
	auto end1 = ar1.begin() + 1;
    cur1++;
    
	auto st2 = ar2.begin();
	auto end2 = ar2.begin() + 1;
    cur2++;

    do {
    	//check if the current interval of R is contained completely in interval of S
    	if(*st1 != *st2 || *end1 != *end2){
    		return 0;
    	}
		// move both lists' indexes
		st1 += 2;
		end1 += 2;
		cur1++;
		
		st2 += 2;
		end2 += 2;
		cur2++;
	} while(cur1<=numintervals1 && cur2<=numintervals2);
			
	//if we didnt check all of the intervals, then the lists do not match
	if(cur1 <= numintervals1 || cur2 <= numintervals2){	
		return 0;
	}
	//all intervals match
	return 1;
}

int compareIntervalsEfficient(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2, bool &intersects){
	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	
	// ID st1,st2,end1,end2;
	uint cur1=0;
	uint cur2=0;
	
	auto st1 = ar1.begin();
	auto end1 = ar1.begin() + 1;
    cur1++;
    
	auto st2 = ar2.begin();
	auto end2 = ar2.begin() + 1;
    cur2++;

    bool intervalRcontained = false;

    do {
    	//check if the current interval of R is contained completely in interval of S
    	if(*st1 >= *st2 && *end1 <= *end2){
    		intervalRcontained = true;
			// containment == intersection
			intersects = true;
    	}else if(*st1<=*st2 && *end1>*st2){
			// there is an intersection
			intersects = true;
		}else if(*st1>*st2 && *st1<*end2){
			// there is an intersection
			intersects = true;
		}

		if (*end1<=*end2)
		{
			if(!intervalRcontained && intersects){
				//we are skipping this interval because it was never contained, so return false (not within)
				// while also we have detected at least one intersection, so there's no point in continuing
				return 0;
			}
			st1 += 2;
			end1 += 2;
			cur1++;
			intervalRcontained = false;
		}
		else 
		{
			st2 += 2;
			end2 += 2;
			cur2++;
		}
	} while(cur1<=numintervals1 && cur2<=numintervals2);
			
	//if we didnt check all of the R intervals
	if(cur1 <= numintervals1 || !intervalRcontained){	
		return 0;
	}
	//all intervals R were contained
	return 1;
}

int compareIntervalsHybrid(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2, bool &intersects){
	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	
	// ID st1,st2,end1,end2;
	uint cur1=0;
	uint cur2=0;
	
	auto st1 = ar1.begin();
	auto end1 = ar1.begin() + 1;
    cur1++;
    
	auto st2 = ar2.begin();
	auto end2 = ar2.begin() + 1;
    cur2++;

    bool intervalRcontained = false;

    do {
    	//check if the current interval of R is contained completely in interval of S
    	if(*st1 >= *st2 && *end1 <= *end2){
    		intervalRcontained = true;
			intersects = true;
			// printf("set intervalRcontained and intersects to true\n");
    	}else if(*st1<=*st2 && *end1>*st2){
			// there is an intersection
			// printf("set intersects to true\n");
			intersects = true;
		}else if(*st1>*st2 && *st1<*end2){
			// there is an intersection
			// printf("set intersects to true\n");
			intersects = true;
		}

		if (*end1<=*end2)
		{
			if(!intervalRcontained){
				// found an interval that is not contained, continue looking for intersections
				// if one was not found yet
				if (!intersects) {
					// get next interval
					st1 += 2;
					end1 += 2;
					cur1++;
					// and jump to continue looking for intersection
					// printf("jumped to look for overlap\n");
					goto LOOK_FOR_OVERLAP;
				}
				// we have an intersection and non-containment, so return non-containment
				return 0;
			}
			st1 += 2;
			end1 += 2;
			cur1++;
			intervalRcontained = false;
		}
		else 
		{
			st2 += 2;
			end2 += 2;
			cur2++;
		}
	} while(cur1<=numintervals1 && cur2<=numintervals2);

	//if we didnt check all of the R intervals
	if(cur1 <= numintervals1){	
		// printf("returning non-containment\n");
		return 0;
	}
	//all intervals R were contained
	// printf("returning containment\n");
	return 1;

LOOK_FOR_OVERLAP:
	do {
		if (*st1<=*st2)
		{
			if (*end1>*st2) // overlap, end1>=st2 if intervals are [s,e] and end1>st2 if intervals are [s,e)
			{
				//they intersect
				intersects = true;
				// printf("found intersection\n");
				break;
			}	
			else
			{
				st1 += 2;
				end1 += 2;
				cur1++;
			}
		}
		else // st2<st1
		{
			if (*end2>*st1) // overlap, end2>=st1 if intervals are [s,e] and end2>st1 if intervals are [s,e)
			{

				//they intersect
				intersects = true;
				// printf("found intersection\n");
				break;
			}
			else
			{
				st2 += 2;
				end2 += 2;
				cur2++;
			}
		}
	} while(cur1<=numintervals1 && cur2<=numintervals2);

	// guaranteed non-containment at this point
	// printf("returning non-containment\n");
	// exit(0);
	return 0;
}

