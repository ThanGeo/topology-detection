#include "topology.h"




vector<int> relationRefinementBitMap = {REF_UNMARKED, REF_UNMARKED, REF_UNMARKED, REF_UNMARKED, REF_UNMARKED, REF_UNMARKED, REF_UNMARKED, REF_UNMARKED, REF_UNMARKED};

void printTable(int AA, int AF, int FA, int FF) {
	cout << "---- Table ----" << endl;
	cout << AA << "   |   " << AF << endl;
	cout << FA << "   |   " << FF << endl;
	cout << "---------------" << endl;
}

typedef struct testCase
{
	vector<uint> ar1;
	uint numintervals1;
	vector<uint> ar2;
	uint numintervals2;
}testCaseT;

void testAlgorithmContainmentIntersection() {
	vector<testCaseT> testCases;
	vector<uint> containmentExpectedResults;
	vector<bool> intersectionExpectedResults;
	// INIT TEST CASES
	// test 1
	testCaseT test = {{0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7, {10,11,12,13,14,15,201,202}, 4};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 2
	test = {{10,11,12,13,14,15,201,202}, 4, {0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(R_INTERVALS_CONTAINED_IN_S_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 3
	test = {{0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7, {201,202}, 1};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 4
	test = {{201,202}, 1, {0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(R_INTERVALS_CONTAINED_IN_S_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 5
	test = {{0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7, {0,1}, 1};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 6
	test = {{0,1}, 1, {0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(R_INTERVALS_CONTAINED_IN_S_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 7
	test = {{0,10,20,50,60,100,800,900}, 4, {0,5,800,899}, 2};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 8
	test = {{0,5,800,899}, 2, {0,10,20,50,60,100,800,900}, 4};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(R_INTERVALS_CONTAINED_IN_S_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 9
	test = {{0,10,20,50,60,100,800,900}, 4, {800,900}, 1};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 10
	test = {{800,900}, 1, {0,10,20,50,60,100,800,900}, 4};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(R_INTERVALS_CONTAINED_IN_S_INTERVALS);
	intersectionExpectedResults.emplace_back(true);
	// test 11
	test = {{0,10,20,50,60,100,800,900}, 4, {800,901}, 1};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(NO_CONTAINMENT);
	intersectionExpectedResults.emplace_back(true);
	// test 12
	test = {{800,901}, 1, {0,10,20,50,60,100,800,900}, 4};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(NO_CONTAINMENT);
	intersectionExpectedResults.emplace_back(true);
	// test 13
	test = {{0,10,20,50,60,100,800,900}, 4, {900,1000}, 1};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(NO_CONTAINMENT);
	intersectionExpectedResults.emplace_back(false);
	// test 14
	test = {{900,1000}, 1, {0,10,20,50,60,100,800,900}, 4};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(NO_CONTAINMENT);
	intersectionExpectedResults.emplace_back(false);
	// test 15 - equal
	test = {{1,2,10,20,40,100}, 3, {1,2,10,20,40,100}, 3};
	testCases.emplace_back(test);
	containmentExpectedResults.emplace_back(SYMMETRICAL_CONTAINMENT);
	intersectionExpectedResults.emplace_back(true);

	// TEST ALL CASES
	for (int i=0; i<testCases.size(); i++) {
		// if (i == 4) {
		// }
		bool intersection = false;
		int containment = compareIntervalsSymmetricalContainmentAndIntersection(testCases.at(i).ar1, testCases.at(i).numintervals1, testCases.at(i).ar2, testCases.at(i).numintervals2, intersection);
		if (containment == containmentExpectedResults.at(i) && intersection == intersectionExpectedResults.at(i)) {
			printf("Test %d passed \n", i+1);
		} else {
			printf("Test %d failed \n", i+1);
			printf("	Containment: %d, intersection: %d\n", containment, intersection);
		}
	}

}

void testAlgorithmMatch() {
	vector<testCaseT> testCases;
	vector<int> matchExpectedResults;
	vector<int> containmentExpectedResults;
	// INIT TEST CASES
	// test 1
	testCaseT test = {{0,1,2,3,4,5,6,7,8,9,10,100,200,900}, 7, {10,11,12,13,14,15,201,202}, 4};
	testCases.emplace_back(test);
	matchExpectedResults.emplace_back(0);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);
	// test 2
	test = {{0,1,2,3,4,5,6,7}, 4, {0,1,2,3,4,5,6,7}, 4};
	testCases.emplace_back(test);
	matchExpectedResults.emplace_back(1);
	containmentExpectedResults.emplace_back(SYMMETRICAL_CONTAINMENT);
	// test 3
	test = {{0,1}, 1, {0,1}, 1};
	testCases.emplace_back(test);
	matchExpectedResults.emplace_back(1);
	containmentExpectedResults.emplace_back(SYMMETRICAL_CONTAINMENT);
	// test 4
	test = {{0,1,2,3,4,5,6,7}, 4, {0,1,2,3,4,5,6,7,8,9}, 5};
	testCases.emplace_back(test);
	matchExpectedResults.emplace_back(0);
	containmentExpectedResults.emplace_back(R_INTERVALS_CONTAINED_IN_S_INTERVALS);
	// test 5
	test = {{0,1,2,3,4,5,6,7,8,9}, 5, {0,1,2,3,4,5,6,7}, 4};
	testCases.emplace_back(test);
	matchExpectedResults.emplace_back(0);
	containmentExpectedResults.emplace_back(S_INTERVALS_CONTAINED_IN_R_INTERVALS);

	// TEST ALL CASES
	for (int i=0; i<testCases.size(); i++) {
		// if (i == 4) {
		// }
		bool intersects;
		int containment = compareIntervalsSymmetricalContainmentAndIntersection(testCases.at(i).ar1, testCases.at(i).numintervals1, testCases.at(i).ar2, testCases.at(i).numintervals2, intersects);
		int matchRes = compareIntervalsForMatch(testCases.at(i).ar1, testCases.at(i).numintervals1, testCases.at(i).ar2, testCases.at(i).numintervals2);
		if (matchRes == matchExpectedResults.at(i) && containment == containmentExpectedResults.at(i)) {
			printf("Test %d passed \n", i+1);
		} else {
			printf("Test %d failed \n", i+1);
			printf("	Match: %d\n", matchRes);
			printf("	Containment: %d\n", containment);
		}
	}

}

int findRelationUsingAPRIL(Polygon *polA, Polygon *polB, bool &markedForEqual) {
	bool AAintersect = false;
	bool AFintersect = false;
	bool FAintersect = false;

	int AAcontainment = compareIntervalsSymmetricalContainmentAndIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL, AAintersect);
	if (AAcontainment == SYMMETRICAL_CONTAINMENT) {
		// AA symmetrical containment, happens in EQUAL
		// check if candidate for EQUAL
		if (polA->F && polB->F) {
			int FF = compareIntervalsForMatch(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedF, polB->numIntervalsF);
			if (FF == 1) {
				// if full lists match, then the polygons might be EQUAL
				markedForEqual = true;
			}
		} else if (!polA->F && !polB->F) {
			// if either has full, also possibly EQUAL
			markedForEqual = true;
		}
	}
	if (AAcontainment == R_INTERVALS_CONTAINED_IN_S_INTERVALS || AAcontainment == SYMMETRICAL_CONTAINMENT) {
		// R_ALL contained in S_ALL guaranteed
		if(polB->F) {
			int AFcontainment = compareIntervalsHybrid(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF, AFintersect);
			if (AFcontainment) {
				// true hit within
				if (markedForEqual) {
					// if it is an equal candidate, return this
					return REFINE_EQUAL_POST_WITHIN;
				}
				// can't be equal, its within
				return WITHIN;
			}

			if (AFintersect) {
				// intersection between ALL-FULL
				// no true hit for containment doesnt mean NO containment
				// refine
				return REFINE_CONTAINMENT;
			} else {
				// if no true hits, need refinement for all
				// printf("Containment + No intersection\n");
				return REFINE_WITHIN_PLUS;
				// return REFINE_ALL;
			}
		}
		// no full cells in S, TODO: do what?	
	}
	
	if(AAcontainment == S_INTERVALS_CONTAINED_IN_R_INTERVALS || AAcontainment == SYMMETRICAL_CONTAINMENT) {
		// S_ALL contained in R_ALL guaranteed
		if (polA->F) {
			int FAcontainment = compareIntervalsHybrid(polB->uncompressedALL, polB->numIntervalsALL, polA->uncompressedF, polA->numIntervalsF, FAintersect);
			if (FAcontainment) {
				// true hit contains
				if (markedForEqual) {
					// if it is an equal candidate, return this
					return REFINE_EQUAL_POST_CONTAIN;
				}

				// true hit contains
				return CONTAINS;
			} 

			if (FAintersect) {
				// intersection between FULL-ALL
				// no true hit for containment doesnt mean NO containment
				// refine
				return REFINE_CONTAINMENT;
			} else {
				// if no true hits, need refinement for all
				// printf("Containment + No intersection\n");
				return REFINE_CONTAINS_PLUS;
				// return REFINE_ALL;
			}
		}
		// no full cells in R, TODO: do what?
	} 
	
	if (AAcontainment == NO_CONTAINMENT) {
		// NO-CONTAINMENT IN AA
		if (!AAintersect) {
			// true hit for disjoint
			return DISJOINT;
		} else{
			if(polB->F) {
				AFintersect = compareIntervalsForIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
				if (AFintersect) {
					// true hit intersect, but candidate for some relations
					return REFINE_EQUAL_INTERSECT;
				}
			}
			if(polA->F) {
				FAintersect = compareIntervalsForIntersection(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL);
				if (FAintersect) {
					// true hit intersect, but candidate for some relations
					return REFINE_EQUAL_INTERSECT;
				}
			}
			// no containment
			return REFINE_NO_CONTAINMENT;
		}
	}

	// probably error if we reach here
	// printf("Reached here: \n");
	// printf("	AAcontainment: %d\n", AAcontainment);
	// printf("	AAintersects: %d\n", AAintersect);
	return REFINE_ALL;
}

int checkTopologyDisjoint(Polygon *polA, Polygon *polB) {
	// compute AA
	int res = checkTopologyIntersects(polA, polB);
	switch (res) {
		case TRUE_NEGATIVE:
			return TRUE_HIT;
			break;
		case TRUE_HIT:
			return TRUE_NEGATIVE;
			break;
	}
	return INCONCLUSIVE;
}

int checkTopologyEqual(Polygon *polA, Polygon *polB) {
	// optimized
	// compute AA for match
	int res = compareIntervalsForMatch(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	if (res == 0) {
		return TRUE_NEGATIVE;
	}
	if(polA->F && polB->F) {
		int res = compareIntervalsForMatch(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedF, polB->numIntervalsF);
		if (res == 0) {
			return TRUE_NEGATIVE;
		}
	}
	return INCONCLUSIVE;
}

int checkTopologyRwithinS(Polygon *polA, Polygon *polB) {
	// optimized
	// compute AA for containment
	int res = compareIntervalsForContainment(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	if (res == 0) {
		return TRUE_NEGATIVE;
	}
	// compute AF for containment (if any F in polB)
	if(polB->F){
		res = compareIntervalsForContainment(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
		if (res == 1) {
			return TRUE_HIT;
		}
	}
	return INCONCLUSIVE;
}

int checkTopologyRContainsS(Polygon *polA, Polygon *polB) {
	// optimized
	// compute AA for containment
	int res = compareIntervalsForContainment(polB->uncompressedALL, polB->numIntervalsALL, polA->uncompressedALL, polA->numIntervalsALL);
	if (res == 0) {
		return TRUE_NEGATIVE;
	}
	// compute FA for containment (if any F in polA)
	if(polA->F){
		res = compareIntervalsForContainment(polB->uncompressedALL, polB->numIntervalsALL, polA->uncompressedF, polA->numIntervalsF);
		if (res == 1) {
			// if all polB->ALL in polA->F
			return TRUE_HIT;
		}
	}
	return INCONCLUSIVE;
}

int checkTopologyRCoveredByS(Polygon *polA, Polygon *polB) {
	// optimized
	// compute AA for containment
	int res = compareIntervalsForContainment(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL);
	if (res == 0) {
		return TRUE_NEGATIVE;
	}
	// compute AF for containment (if any F in polB)
	if(polB->F){
		res = compareIntervalsForContainment(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
		if (res == 1) {
			return TRUE_HIT;
		}
	}
	return INCONCLUSIVE;
}

int checkTopologyRCoversS(Polygon *polA, Polygon *polB) {
	// optimized
	// compute AA for containment
	int res = compareIntervalsForContainment(polB->uncompressedALL, polB->numIntervalsALL, polA->uncompressedALL, polA->numIntervalsALL);
	if (res == 0) {
		return TRUE_NEGATIVE;
	}
	// compute FA for containment (if any F in polA)
	if(polA->F){
		res = compareIntervalsForContainment(polB->uncompressedALL, polB->numIntervalsALL, polA->uncompressedF, polA->numIntervalsF);
		if (res == 1) {
			// if all polB->ALL in polA->F
			return TRUE_HIT;
		}
	}
	return INCONCLUSIVE;
}

int checkTopologyCrosses(Polygon *polA, Polygon *polB) {
	// TODO
}

int checkTopologyMeet(Polygon *polA, Polygon *polB) {
	// optimized
	// compute AA for intersection
	int res = compareIntervalsForIntersection(polB->uncompressedALL, polB->numIntervalsALL, polA->uncompressedALL, polA->numIntervalsALL);
	if (res == 0) {
		return TRUE_NEGATIVE;
	}
	// compute AF for intersection (if any F in polB)
	if(polB->F){
		res = compareIntervalsForIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF);
		if (res == 1) {
			return TRUE_NEGATIVE;
		}
	}
	// compute FA for intersection (if any F in polA)
	if(polA->F){
		res = compareIntervalsForIntersection(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL);
		if (res == 1) {
			return TRUE_NEGATIVE;
		}
	}
	return INCONCLUSIVE;
}

//join two uncompressed APRIL approximations
int checkTopologyIntersects(Polygon *polA, Polygon *polB){
	//check ALL - ALL
	if(compareIntervalsForIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedALL, polB->numIntervalsALL) == 0){
		//guaranteed not hit
		return TRUE_NEGATIVE;
	}

	//check ALL - F (if any F in B)
	if(polB->F){
		if(compareIntervalsForIntersection(polA->uncompressedALL, polA->numIntervalsALL, polB->uncompressedF, polB->numIntervalsF)){
			//hit
			return TRUE_HIT;
		}
	}

	//check F - ALL (if any F in A)
	if(polA->F){
		if(compareIntervalsForIntersection(polA->uncompressedF, polA->numIntervalsF, polB->uncompressedALL, polB->numIntervalsALL)){
			//hit
			return TRUE_HIT;
		}
	}

	//The weak have not been checked
	return INCONCLUSIVE; //send to refinement
}

