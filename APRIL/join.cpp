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

int compareIntervalsWithin(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2){
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
    	//check if it is contained completely
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


/*
*	
*	INTERVALS CAN EITHER:
*	
*	no hits			= 0
*	at least 1 hit 	= 1
*	match 100%		= 2
*	all R, any S	= 3
*	all S, any R	= 4
*
*/

int compareIntervals(vector<ID> &ar1, uint &numintervals1, vector<ID> &ar2, uint &numintervals2){
	//they may not have any intervals of this type
	if(numintervals1 == 0 || numintervals2 == 0){
		return 0;
	}
	

	bool zero_overlaps = true;
	bool at_least_one_hit = false;
	bool equal = true;
	bool all_of_R = true;
	bool all_of_S = true;

	uint R_last_overlap = 0;
	uint S_last_overlap = 0;

	//different amount of intervals = not equal polygons
	if(numintervals1 != numintervals2){
		equal = false;
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

    // cout << "total intervals: " << numintervals1 << "," << numintervals2 << endl;

    do {

    	// cout << *st1 << "," << *end1 << " with " << *st2 << "," << *end2 << "(" << cur1 << " with " << cur2 << ")" << endl;

    	if(cur1 == cur2){
    		if(*st1 != *st2 || *end1 != *end2){
    			equal = false;
    		}
    	}

		if (*st1<=*st2)
		{
			if (*end1>*st2) // overlap, end1>=st2 if intervals are [s,e] and end1>st2 if intervals are [s,e)
			{
				//they overlap
				zero_overlaps = false;
				at_least_one_hit = true;
				// cout << "1. overlap!" << endl;

				if((R_last_overlap != cur1 - 1 && R_last_overlap != cur1) || (*st1 < *st2 || *end1 > *end2)){
					all_of_R = false;
				}
				if((S_last_overlap != cur2 - 1 && S_last_overlap != cur2) || (*st2 < *st1 || *end2 > *end1)){
					all_of_S = false;
				}

				R_last_overlap = cur1;
				S_last_overlap = cur2;

				if(*end1 < *end2){
					st1 += 2;
					end1 += 2;
					cur1++;
				}else{
					st2 += 2;
					end2 += 2;
					cur2++;
				}


			}else{
				st1 += 2;
				end1 += 2;
				cur1++;
			}
		}
		else // st2<st1
		{
			if (*end2>*st1) // overlap, end2>=st1 if intervals are [s,e] and end2>st1 if intervals are [s,e)
			{

				//they overlap
				zero_overlaps = false;
				at_least_one_hit = true;	
				// cout << "2. overlap!" << endl;
			
				if((R_last_overlap != cur1 - 1 && R_last_overlap != cur1) || (*st1 < *st2 || *end1 > *end2)){
					all_of_R = false;
				}
				if((S_last_overlap != cur2 - 1 && S_last_overlap != cur2) || (*st2 < *st1 || *end2 > *end1)){
					all_of_S = false;
				}

				R_last_overlap = cur1;
				S_last_overlap = cur2;

				if(*end1 < *end2){
					st1 += 2;
					end1 += 2;
					cur1++;
				}else{
					st2 += 2;
					end2 += 2;
					cur2++;
				}
			}else{
				st2 += 2;
				end2 += 2;
				cur2++;
			}
		}
	} while(cur1<=numintervals1 && cur2<=numintervals2);

	
	// cout << R_last_overlap << "," << cur1 << endl;
	// cout << S_last_overlap << "," << cur2 << endl;

	//end checks
	
	if(cur1 == numintervals1 && R_last_overlap != cur1){
		all_of_R = false;
	}
	if(cur2 == numintervals2 && S_last_overlap != cur2){
		all_of_S = false;
	}
	

	if(R_last_overlap != cur1 - 1 && R_last_overlap != cur1){
		all_of_R = false;
	}
	if(S_last_overlap != cur2 - 1 && S_last_overlap != cur2){
		all_of_S = false;
	}


	if(all_of_R && cur1 < numintervals1){
		//if we didnt check all of R's intervals, then it cannot be contained in S
		all_of_R = false;
	}
	if(all_of_S && cur2 < numintervals2){
		//if we didnt check all of S's intervals, then it cannot be contained in R
		all_of_S = false;
	}

	// cout << "equal: " << equal << endl;
	// cout << "zero_overlaps: " << zero_overlaps << endl;
	// cout << "at_least_one_hit: " << at_least_one_hit << endl;
	// cout << "all_of_R: " << all_of_R << endl;
	// cout << "all_of_S: " << all_of_S << endl;


	//no overlaps have been detected, disjoint
	if(zero_overlaps){
		return 0;
	}

	//if they are the same
	if(equal){
		return 2;
	}

	//if all of R overlap with any of S
	if(all_of_R){
		return 3;
	}

	//if all of S overlap with any of R
	if(all_of_S){
		return 4;
	}

	//else they intersect at least once
	return 1;
}
