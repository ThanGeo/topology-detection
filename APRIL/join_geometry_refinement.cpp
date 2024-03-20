#include "join_geometry_refinement.h"

char withinCode[] = "T*F**F***";
char coveredbyCode1[] = "T*F**F***";
char coveredbyCode2[] = "*TF**F***";
char coveredbyCode3[] = "**FT*F***";
char coveredbyCode4[] = "**F*TF***";
char containsCode[] = "T*****FF*";
char coversCode1[] = "T*****FF*";
char coversCode2[] = "*T****FF*";
char coversCode3[] = "***T**FF*";
char coversCode4[] = "****T*FF*";
char meetCode1[] = "FT*******"; 
char meetCode2[] = "F**T*****"; 
char meetCode3[] = "F***T****"; 
char equalCode[] = "T*F**FFF*"; 
char disjointCode[] = "FF*FF****";
char overlapCode1[] = "T********";
char overlapCode2[] = "*T*******";
char overlapCode3[] = "***T*****";
char overlapCode4[] = "****T****";

//define topological masks for refinement
// a within b
boost::geometry::de9im::mask withinMask(withinCode); 
// a contains b
boost::geometry::de9im::mask containsMask(containsCode); 
// a covered by b
vector<boost::geometry::de9im::mask> coveredByMaskList = {
				boost::geometry::de9im::mask(coveredbyCode1),
				boost::geometry::de9im::mask(coveredbyCode2),
				boost::geometry::de9im::mask(coveredbyCode3),
				boost::geometry::de9im::mask(coveredbyCode4)};
// a covers b
vector<boost::geometry::de9im::mask> coversMaskList = {
				boost::geometry::de9im::mask(coversCode1),
				boost::geometry::de9im::mask(coversCode2),
				boost::geometry::de9im::mask(coversCode3),
				boost::geometry::de9im::mask(coversCode4)};
// a and b meet
boost::geometry::de9im::mask meetMask1(meetCode1); 
boost::geometry::de9im::mask meetMask2(meetCode2); 
boost::geometry::de9im::mask meetMask3(meetCode3); 
// a and b are equal
boost::geometry::de9im::mask equalMask(equalCode); 
// a and b are disjoint
boost::geometry::de9im::mask disjointMask(disjointCode); 
// a overlaps b
vector<boost::geometry::de9im::mask> overlapMaskList = {
				boost::geometry::de9im::mask(overlapCode1),
				boost::geometry::de9im::mask(overlapCode2),
				boost::geometry::de9im::mask(overlapCode3),
				boost::geometry::de9im::mask(overlapCode4)};

static polygon loadPolygonGeometryBOOST(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	polygon pol;
	int readID;
	int vertexCount, polygonCount;
	double x,y;

	//search the map for the specific polygon offset
	unordered_map<uint,unsigned long>::const_iterator got = offsetMap.find(recID);
	if(got != offsetMap.end()){ 
		//set read offset
		fin.seekg(got->second-fin.tellg(), fin.cur);		
		//read rec ID
		fin.read((char*) &readID, sizeof(int));
		//read vertex count
		fin.read((char*) &vertexCount, sizeof(int));
		for(int i=0; i<vertexCount; i++){
			fin.read((char*) &x, sizeof(double));
			fin.read((char*) &y, sizeof(double));

			pol.outer().push_back(point_xy(x,y));
		}
	}

	boost::geometry::correct(pol);
	return pol;
}

static void printPolygonFromDisk(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	int readID;
	int vertexCount, polygonCount;
	double x,y;

	//search the map for the specific polygon offset
	unordered_map<uint,unsigned long>::const_iterator got = offsetMap.find(recID);
	if(got != offsetMap.end()){ 
		//set read offset
		fin.seekg(got->second-fin.tellg(), fin.cur);		
		//read rec ID
		fin.read((char*) &readID, sizeof(int));
		//read vertex count
		fin.read((char*) &vertexCount, sizeof(int));
		for(int i=0; i<vertexCount; i++){
			fin.read((char*) &x, sizeof(double));
			fin.read((char*) &y, sizeof(double));
			printf("(%f,%f),",x,y);
		}
	}
	printf("\n");
}

int refinement_DE9IM_WithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);

    
    //disjoint
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, disjointMask)){
    	return DISJOINT;
    }

    //equal
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, equalMask)){
    	return EQUAL;
    }

    //R within S
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, withinMask)){
    	return WITHIN;
    }
	// R contains S
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, containsMask)){
    	return CONTAINS;
    }

	// R covered by S
	for(auto &it: coveredByMaskList){
		if(boost::geometry::relate(boostPolygonR, boostPolygonS, it)){
    		return R_COVERED_BY_S;
		}
	}
	// R covers S
	for(auto &it: coversMaskList){
		if(boost::geometry::relate(boostPolygonR, boostPolygonS, it)){
    		return R_COVERS_S;
		}
	}

    //meet
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask1) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask2) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask3)){
    	return MEET;
    }

    //else return overlap/intersects
    return INTERSECT;

}

// 'mask_char' character is always going to be T or F
// 'character' character can be 0,1,2 which count as T
static bool compareDe9imChars(char character, char char_mask) {
	if (character != 'F' && char_mask == 'T') {
		// character is 0,1,2 and char_mask is T
		return true;
	} else if (character == 'F' && char_mask == 'F'){
		// both are F
		return true;
	} else {
		// no match
		return false;
	}
}

static bool compareMasks(std::string &de9imCode, char* maskCode) {
	for(int i=0; i<9; i++) {
		if (de9imCode[i] == '*' || maskCode[i] == '*' || compareDe9imChars(de9imCode[i], maskCode[i])){
			continue;
		} else {
			return false;
		}
	}
	return true;
}

static void createMaskCodes(polygon &boostPolygonR, polygon &boostPolygonS, std::string &codeRS, std::string &codeSR) {
	boost::geometry::de9im::matrix matrixRS = boost::geometry::relation(boostPolygonR, boostPolygonS);
	boost::geometry::de9im::matrix matrixSR = boost::geometry::relation(boostPolygonS, boostPolygonR);
   	codeRS = matrixRS.str();
    codeSR = matrixSR.str();
}

static void createMaskCode(polygon &boostPolygonR, polygon &boostPolygonS, std::string &code) {
	boost::geometry::de9im::matrix matrix = boost::geometry::relation(boostPolygonR, boostPolygonS);
   	code = matrix.str();
}

bool refinementDisjoint(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
	
	// // get the mask codes
    // std::string codeRS,codeSR;
	// createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
	// // compare masks
	// if(compareMasks(codeRS, disjointCode)){
    // 	return TRUE_HIT;
    // }
	// return TRUE_NEGATIVE;

	return boost::geometry::disjoint(boostPolygonR, boostPolygonS);
}
bool refinementEqual(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    
	// // get the mask codes
	// std::string codeRS,codeSR;
	// createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
	// // compare masks
	// if(compareMasks(codeRS, equalCode)){
    // 	return TRUE_HIT;
    // }
	// return TRUE_NEGATIVE;
	
	return boost::geometry::equals(boostPolygonR, boostPolygonS);
}
bool refinementIntersect(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    // std::string codeRS,codeSR;
	// // get the mask codes
	// createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
	// // compare masks
	// if(compareMasks(codeSR, overlapCode1) || 
	// 			compareMasks(codeSR, overlapCode2) || 
	// 			compareMasks(codeSR, overlapCode3) || 
	// 			compareMasks(codeSR, overlapCode4)){
	// 	return TRUE_HIT;
	// }
	// return TRUE_NEGATIVE;

	return boost::geometry::intersects(boostPolygonR, boostPolygonS);
}
bool refinementMeet(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    // std::string codeRS,codeSR;
	// // get the mask codes
	// createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
	// // compare masks
	// if(compareMasks(codeRS, meetCode1) || 
	// 			compareMasks(codeRS, meetCode2) || 
	// 			compareMasks(codeRS, meetCode3)){
	// 	return TRUE_HIT;
	// }
	// return TRUE_NEGATIVE;

	return boost::geometry::touches(boostPolygonR, boostPolygonS);
}
bool refinementCrosses(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// todo
}
bool refinementRcoversS(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
	
	// // get the mask codes
    // std::string codeRS,codeSR;
	// createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
	// // compare masks
	// if(compareMasks(codeRS, coversCode1) || 
	// 			compareMasks(codeRS, coversCode2) || 
	// 			compareMasks(codeRS, coversCode3) || 
	// 			compareMasks(codeRS, coversCode4)){
	// 	return TRUE_HIT;
	// }
	// return TRUE_NEGATIVE;

	// R covers S == S covered by R
	return boost::geometry::covered_by(boostPolygonS, boostPolygonR);
}
bool refinementRcoveredByS(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    
	// get the mask codes
	std::string code;
	createMaskCode(boostPolygonR, boostPolygonS, code);
	// compare masks
	if(compareMasks(code, coveredbyCode1) || 
				compareMasks(code, coveredbyCode2) || 
				compareMasks(code, coveredbyCode3) || 
				compareMasks(code, coveredbyCode4)){
		return TRUE_HIT;
	}
	return TRUE_NEGATIVE;

	// return boost::geometry::covered_by(boostPolygonR, boostPolygonS) && !(boost::geometry::within(boostPolygonR, boostPolygonS));
}
bool refinementRcontainsS(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
   
	// // get the mask codes
    // std::string codeRS,codeSR;
	// createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
	// // compare masks
	// if(compareMasks(codeRS, containsCode)){
	// 	return TRUE_HIT;
	// }
	// return TRUE_NEGATIVE;

	// R contains S == S within R, and NOT equal
	return (boost::geometry::within(boostPolygonS, boostPolygonR) && !boost::geometry::equals(boostPolygonS, boostPolygonR));
}
bool refinementRwithinS(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS) {
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    
	// // using DE9IM masks
	// std::string code;
	// createMaskCode(boostPolygonR, boostPolygonS, code);
	// // compare masks
	// if(compareMasks(code, withinCode)){
	// 	return TRUE_HIT;
	// }
	// return TRUE_NEGATIVE;

	// within and NOT equal
	return (boost::geometry::within(boostPolygonR, boostPolygonS) && !boost::geometry::equals(boostPolygonR, boostPolygonS));
}


int refinementAllTopologyRelations(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS, bool markedForEqual){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    std::string code;
	// get the mask codes
	createMaskCode(boostPolygonR, boostPolygonS, code);
    
	//disjoint
	if(compareMasks(code, disjointCode)){
    	return DISJOINT;
    }

	/**
	 * it definitely intersects at this point
	*/

	if (markedForEqual) {
		// check equality first because it is a subset of covers and covered by
		if(compareMasks(code, equalCode)){
			return EQUAL;
		}
	}

	// covers
	if(compareMasks(code, coversCode1) || compareMasks(code, coversCode2) || compareMasks(code, coversCode3)|| compareMasks(code, coversCode4)){
		// first check contains because it is a subset of covers
		if(compareMasks(code, containsCode)){
			return CONTAINS;
		}
		return R_COVERS_S;
    }

	// covered by
	if(compareMasks(code, coveredbyCode1) || compareMasks(code, coveredbyCode2) || compareMasks(code, coveredbyCode3)|| compareMasks(code, coveredbyCode4)){
		// first check within because it is a subset of covered by
		if(compareMasks(code, withinCode)){
			return WITHIN;
		}
		return R_COVERED_BY_S;
    }

	// meet
	if(compareMasks(code, meetCode1) || 
				compareMasks(code, meetCode2) || 
				compareMasks(code, meetCode3)){
		return MEET;
	}

	// else return intersects
	return INTERSECT;
}

int refinementContainmentPostAPRIL(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS, bool markedForEqual){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    std::string code;
	// get the mask codes
	createMaskCode(boostPolygonR, boostPolygonS, code);
   
	// equal
	if (markedForEqual) {
		if (compareMasks(code, equalCode)) {
			return EQUAL;
		}
	}

	// covers
	if(compareMasks(code, coversCode1) || compareMasks(code, coversCode2) || compareMasks(code, coversCode3)|| compareMasks(code, coversCode4)){
		// it DOES cover
		
		// check contains
		if(compareMasks(code, containsCode)){
			return CONTAINS;
		}

		// it covers
		return R_COVERS_S;

    }
	// covered by
	if(compareMasks(code, coveredbyCode1) || compareMasks(code, coveredbyCode2) || compareMasks(code, coveredbyCode3)|| compareMasks(code, coveredbyCode4)){
		// it IS covered by

		// check within
		if(compareMasks(code, withinCode)){
			return WITHIN;
		}

		// it covers
		return R_COVERED_BY_S;

    }

	// else return intersect
	return INTERSECT;
}

int refinementContainsPlusPostAPRIL(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS, bool markedForEqual){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    std::string code;
	// get the mask codes
	createMaskCode(boostPolygonR, boostPolygonS, code);
    
	//disjoint
	if(compareMasks(code, disjointCode)){
    	return DISJOINT;
    }

	if (markedForEqual) {
		// equal
		if(compareMasks(code, equalCode)){
			return EQUAL;
		}
	}

	// covers
	if(compareMasks(code, coversCode1) || compareMasks(code, coversCode2) || compareMasks(code, coversCode3)|| compareMasks(code, coversCode4)){
		// it DOES cover
		
		// check contains
		if(compareMasks(code, containsCode)){
			return CONTAINS;
		}

		// it covers
		return R_COVERS_S;

    }

	// meet
	if(compareMasks(code, meetCode1) || 
				compareMasks(code, meetCode2) || 
				compareMasks(code, meetCode3)){
		return MEET;
	}

	// else return intersects
	return INTERSECT;
}

int refinementWithinPlusPostAPRIL(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS, bool markedForEqual){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    std::string code;
	// get the mask codes
	createMaskCode(boostPolygonR, boostPolygonS, code);
    
	//disjoint
	if(compareMasks(code, disjointCode)){
    	return DISJOINT;
    }

	if (markedForEqual) {
		// equal
		if(compareMasks(code, equalCode)){
			return EQUAL;
		}
	}

	// covered by
	if(compareMasks(code, coveredbyCode1) || compareMasks(code, coveredbyCode2) || compareMasks(code, coveredbyCode3)|| compareMasks(code, coveredbyCode4)){
		// first check within because it is a subset of covered by
		if(compareMasks(code, withinCode)){
			return WITHIN;
		}
		return R_COVERED_BY_S;
    }

	// meet
	if(compareMasks(code, meetCode1) || 
				compareMasks(code, meetCode2) || 
				compareMasks(code, meetCode3)){
		return MEET;
	}

	// else return intersects
	return INTERSECT;
}

int refinementNoContainmentPostAPRIL(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS, bool markedForEqual){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    std::string code;
	// get the mask codes
	createMaskCode(boostPolygonR, boostPolygonS, code);
    
	//disjoint
	if(compareMasks(code, disjointCode)){
    	return DISJOINT;
    }

	if (markedForEqual) {
		// equal
		if(compareMasks(code, equalCode)){
			return EQUAL;
		}
	}

	// meet
	if(compareMasks(code, meetCode1) || 
				compareMasks(code, meetCode2) || 
				compareMasks(code, meetCode3)){
		return MEET;
	}

	// else return intersects
	return INTERSECT;
}

int refinementAllTopologyRelations_old(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// load boost polygons
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
    std::string codeRS,codeSR;
	// get the mask codes
	createMaskCodes(boostPolygonR, boostPolygonS, codeRS, codeSR);
    
    //disjoint
	if(compareMasks(codeRS, disjointCode)){
    	return DISJOINT;
    }

    //equal
	if(compareMasks(codeRS, equalCode)){
    	return EQUAL;
    }

    //inside
	if(compareMasks(codeRS, withinCode)){
    	return WITHIN;
    }
	if(compareMasks(codeSR, withinCode)){
    	return CONTAINS;
    }

    //covered by
	if(compareMasks(codeRS, coveredbyCode1) || 
				compareMasks(codeRS, coveredbyCode2) || 
				compareMasks(codeRS, coveredbyCode3) || 
				compareMasks(codeRS, coveredbyCode4)){
		return R_COVERED_BY_S;
	}
	if(compareMasks(codeSR, coveredbyCode1) || 
				compareMasks(codeSR, coveredbyCode2) || 
				compareMasks(codeSR, coveredbyCode3) || 
				compareMasks(codeSR, coveredbyCode4)){
		return R_COVERS_S;
	}

    //meet
	if(compareMasks(codeRS, meetCode1) || 
				compareMasks(codeRS, meetCode2) || 
				compareMasks(codeRS, meetCode3)){
		return MEET;
	}

    //else return overlap/intersects
    return INTERSECT;

}

