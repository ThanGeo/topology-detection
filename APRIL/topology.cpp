#include "topology.h"

typedef enum {
	NO_OVERLAPS = 0,
	AT_LEAST_ONE_OVERLAP,
	LISTS_MATCH,
	ALL_R_INTERVALS_IN_S_INTERVALS,
	ALL_S_INTERVALS_IN_R_INTERVALS,
} A4IMvalueE;

int (*checkPredicateFunctionPtr)(int AA, int AF, int FA, int FF);

//same granularity
static void createTopologyTableUncompressed(Polygon *polA, Polygon *polB, int &AA, int &AF, int &FA, int &FF){
	// cout << "entered createTopologyTableUncompressed..." << endl;
	//parallel
	// #pragma omp parallel sections
	// {	
	// 	//ALL - ALL
	// 	#pragma omp section
	// 	{
	// 		// cout << "A-A" << endl;
	// 		AA = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	// 	}
	// 	//ALL - F
	// 	#pragma omp section
	// 	{
	// 		// cout << "hi from thread " << omp_get_thread_num() << endl;
	// 		if(polB->F){
	// 			// cout << "A-F" << endl;
	// 			AF = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
	// 		}else{
	// 			AF = -1;
	// 		}
	// 	}
	// 	//F - ALL
	// 	#pragma omp section
	// 	{
	// 		// cout << "hi from thread " << omp_get_thread_num() << endl;
	// 		if(polA->F){
	// 			// cout << "F-A" << endl;
	// 			FA = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL);
	// 		}else{
	// 			FA = -1;
	// 		}
	// 	}
	// 	//F - F
	// 	#pragma omp section
	// 	{
	// 		// cout << "hi from thread " << omp_get_thread_num() << endl;
	// 		if(polA->F && polB->F){
	// 			// cout << "F-F" << endl;
	// 			FF = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedF, polB->numIntervalsF);
	// 		}else{
	// 			FF = -1;
	// 		}
	// 	}
	// }
	//serial
	// cout << "polygons " << polA->recID << " and " << polB->recID << endl;
	// cout << "A-A" << endl;
	AA = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	if(polB->F){
		// cout << "A-F" << endl;
		AF = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
	}else{
		AF = -1;
	}
	if(polA->F){
		// cout << "F-A" << endl;
		FA = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL);
	}else{
		FA = -1;
	}
	if(polA->F && polB->F){
		// cout << "F-F" << endl;
		FF = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedF, polB->numIntervalsF);
	}else{
		FF = -1;
	}
	// cout << "---- Table ----" << endl;
	// cout << AA << "   |   " << AF << endl;
	// cout << FA << "   |   " << FF << endl;
	// cout << "---------------" << endl;
}

int checkAllRelations(Polygon *polA, Polygon *polB){

	// cout << "CHECKING POLYGONS " << polA->recID << " and " << polB->recID << endl;
	int AA,AF,FA,FF;
	
	//parallel
	// #pragma omp parallel sections
	// {	
	// 	//ALL - ALL
	// 	#pragma omp section
	// 	{
	// 		// cout << "A-A" << endl;
	// 		AA = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	// 	}
	// 	//ALL - F
	// 	#pragma omp section
	// 	{
	// 		// cout << "hi from thread " << omp_get_thread_num() << endl;
	// 		if(polB->F){
	// 			// cout << "A-F" << endl;
	// 			AF = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
	// 		}else{
	// 			AF = -1;
	// 		}
	// 	}
	// 	//F - ALL
	// 	#pragma omp section
	// 	{
	// 		// cout << "hi from thread " << omp_get_thread_num() << endl;
	// 		if(polA->F){
	// 			// cout << "F-A" << endl;
	// 			FA = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL);
	// 		}else{
	// 			FA = -1;
	// 		}
	// 	}
	// 	//F - F
	// 	#pragma omp section
	// 	{
	// 		// cout << "hi from thread " << omp_get_thread_num() << endl;
	// 		if(polA->F && polB->F){
	// 			// cout << "F-F" << endl;
	// 			FF = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedF, polB->numIntervalsF);
	// 		}else{
	// 			FF = -1;
	// 		}
	// 	}
	// }

	//serial
	// cout << "A-A" << endl;
	AA = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	if(polB->F){
		// cout << "A-F" << endl;
		AF = compareIntervals(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
	}else{
		AF = -1;
	}
	if(polA->F){
		// cout << "F-A" << endl;
		FA = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL);
	}else{
		FA = -1;
	}
	if(polA->F && polB->F){
		// cout << "F-F" << endl;
		FF = compareIntervals(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedF, polB->numIntervalsF);
	}else{
		FF = -1;
	}

	// cout << "---- Table ----" << endl;
	// cout << AA << "   |   " << AF << endl;
	// cout << FA << "   |   " << FF << endl;
	// cout << "---------------" << endl;


	if(AA == NO_OVERLAPS){
		//disjoint
		return DISJOINT;
	}

	if((AA == LISTS_MATCH && FF == LISTS_MATCH) || (AA == LISTS_MATCH && (!polA->F && !polB->F))){
		//equal - needs refinement
		// return NEEDS_REFINEMENT;
		return EQUAL;
	}

	if(AA == AT_LEAST_ONE_OVERLAP){
		if(AF >= AT_LEAST_ONE_OVERLAP || FA >= AT_LEAST_ONE_OVERLAP){
			//overlap
			return OVERLAP;
		}else if(AF + FA <= NO_OVERLAPS){			
			//meet - needs refinement
			// return NEEDS_REFINEMENT;
			return MEET;
		}
	}

	if(AA == ALL_R_INTERVALS_IN_S_INTERVALS){
		if(AF == ALL_R_INTERVALS_IN_S_INTERVALS){
			//R inside S
			return R_CONTAINED_IN_S;
		}else{
			//R covered_by S - needs refinement
			// return NEEDS_REFINEMENT;
			return R_COVERED_BY_S;
		}
	}

	if(AA == ALL_S_INTERVALS_IN_R_INTERVALS){
		if(FA == ALL_S_INTERVALS_IN_R_INTERVALS){
			//S inside R
			return R_CONTAINS_S;
		}else{
			//S covered_by R - needs refinement
			// return NEEDS_REFINEMENT;
			return R_COVERS_S;
		}
	}

	//error: could not identify topological relationship
	return -1;

}

int checkTopologyDisjoint(int AA, int AF, int FA, int FF) {
	if(AA == NO_OVERLAPS){
		//disjoint
		return 1;
	}
	// APRIL has some false negatives in the disjoint case, so send to refinement
	return 2;
}

int checkTopologyEqual(int AA, int AF, int FA, int FF) {
	if((AA == LISTS_MATCH && FF == LISTS_MATCH) || (AA == LISTS_MATCH && (FA == -1 && AF == -1))){
		//candidate for equal - needs refinement though
		return 2;
	}
	// guaranteed not equal
	return 0;
}

int checkTopologyRContainedInS(int AA, int AF, int FA, int FF) {
	if(AA == ALL_R_INTERVALS_IN_S_INTERVALS){
		if(AF == ALL_R_INTERVALS_IN_S_INTERVALS){
			// R inside S guaranteed
			return 1;
		}else{
			// candidate, needs refinement
			return 2;
		}
	}
	// guaranteed non-containment
	return 0;
}

int checkTopologyRContainsS(int AA, int AF, int FA, int FF) {
	if(AA == ALL_S_INTERVALS_IN_R_INTERVALS){
		if(FA == ALL_S_INTERVALS_IN_R_INTERVALS){
			// S inside R guaranteed
			return 1;
		}else{
			// candidate - needs refinement
			return 2;
		}
	}
	// guaranteed non-containment
	return 0;
}

int checkTopologyRCoveredByS(int AA, int AF, int FA, int FF) {
	if (AA <= NO_OVERLAPS) {
		// disjoint, guaranteed non-coverage
		return 0;
	} else if(AA == ALL_R_INTERVALS_IN_S_INTERVALS){
		if(AF == ALL_R_INTERVALS_IN_S_INTERVALS){
			// R contained in S (subset of covered by)
			return 1;
		}
	}
	// candidate, needs refinement
	return 2;
}

int checkTopologyRCoversS(int AA, int AF, int FA, int FF) {
	if (AA <= NO_OVERLAPS) {
		// disjoint, guaranteed non-coverage
		return 0;
	} else if(AA == ALL_S_INTERVALS_IN_R_INTERVALS){
		if(FA == ALL_S_INTERVALS_IN_R_INTERVALS){
			// R contains S (subset of covers)
			return 1;
		}
	}
	// candidate, needs refinement
	return 2;
}

int checkTopologyCrosses(int AA, int AF, int FA, int FF) {
	// TODO
}

int checkTopologyMeet(int AA, int AF, int FA, int FF) {
	if(AA >= AT_LEAST_ONE_OVERLAP){
		if(AF >= AT_LEAST_ONE_OVERLAP || FA >= AT_LEAST_ONE_OVERLAP){
			//overlap, guaranteed non-meet
			return 0;
		}else if(FF > NO_OVERLAPS) {
			// interiors intersect, guaranteed non-meet
			return 0;
		} else {		
			// case: AF + FA <= NO_OVERLAPS	
			// candidate for meet - needs refinement
			return 2;
		} 
	}
	// disjoint, guaranteed non-meet
	return 0;
}

//join two uncompressed APRIL approximations
int checkTopologyIntersects(Polygon *polA, Polygon *polB){
	//check ALL - ALL
	if(compareIntervalsForIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL) == 0){
		//guaranteed not hit
		return 0;
	}

	//check ALL - F (if any F in B)
	if(polB->F){
		if(compareIntervalsForIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF)){
			//hit
			return 1;
		}
	}

	//check F - ALL (if any F in A)
	if(polA->F){
		if(compareIntervalsForIntersection(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL)){
			//hit
			return 1;
		}
	}

	//The weak have not been checked
	return 2; //send to refinement
}

int checkSpecificRelation(Polygon *polA, Polygon *polB) {
	int AA = -1;
	int AF = -1;
	int FA = -1;
	int FF = -1;
	// cout << "creating topology table..." << endl;
	// create the topology table using APRIL
	createTopologyTableUncompressed(polA, polB, AA, AF, FA, FF);

	// check the specified predicate relation
	return (*checkPredicateFunctionPtr)(AA, AF, FA, FF);
}