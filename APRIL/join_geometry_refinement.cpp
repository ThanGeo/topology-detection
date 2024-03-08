#include "join_geometry_refinement.h"

char withinCode[] = "T*F*FF***";
char coveredbyCode1[] = "**F*TF***";
char coveredbyCode2[] = "*TF**F***";
char coveredbyCode3[] = "**FT*F***";
char coveredbyCode4[] = "**F*TF***";
char meetCode1[] = "FT*******"; 
char meetCode2[] = "F**T*****"; 
char meetCode3[] = "F***T****"; 
char equalCode[] = "T*F**FFF*"; 
char disjointCode[] = "FF*FF****";

//define topological masks for refinement
// a within b
boost::geometry::de9im::mask withinMask(withinCode); 
// a covered by b
vector<boost::geometry::de9im::mask> coveredByMaskList = {
				boost::geometry::de9im::mask(coveredbyCode1),
				boost::geometry::de9im::mask(coveredbyCode2),
				boost::geometry::de9im::mask(coveredbyCode3),
				boost::geometry::de9im::mask(coveredbyCode4)};

// a and b meet
boost::geometry::de9im::mask meetMask1(meetCode1); 
boost::geometry::de9im::mask meetMask2(meetCode2); 
boost::geometry::de9im::mask meetMask3(meetCode3); 
// a and b are equal
boost::geometry::de9im::mask equalMask(equalCode); 
// a and b are disjoint
boost::geometry::de9im::mask disjointMask(disjointCode); 

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

    //inside
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, withinMask)){
    	return R_INSIDE_S;
    }
    if(boost::geometry::relate(boostPolygonS, boostPolygonR, withinMask)){
    	return S_INSIDE_R;
    }

    //covered by
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, coveredbyMask)){
    // 	return R_COVERED_BY_S;
    // }
    // if(boost::geometry::relate(boostPolygonS, boostPolygonR, coveredbyMask)){
    // 	return S_COVERED_BY_R;
    // }

	for(auto &it: coveredByMaskList){
		if(boost::geometry::relate(boostPolygonR, boostPolygonS, it)){
    		return R_COVERED_BY_S;
		}
		if(boost::geometry::relate(boostPolygonS, boostPolygonR, it)){
			return S_COVERED_BY_R;
		}
	}

    //meet
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask1) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask2) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask3)){
    	return MEET;
    }

    //else return overlap/intersects
    return OVERLAP;

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

int refinement_DE9IM_WithIDs_optimized(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);
	
	boost::geometry::de9im::matrix matrixRS = boost::geometry::relation(boostPolygonR, boostPolygonS);
	boost::geometry::de9im::matrix matrixSR = boost::geometry::relation(boostPolygonS, boostPolygonR);
    std::string codeRS = matrixRS.str();
    std::string codeSR = matrixSR.str();
    
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
    	return R_INSIDE_S;
    }
	if(compareMasks(codeSR, withinCode)){
    	return S_INSIDE_R;
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
		return S_COVERED_BY_R;
	}

    //meet
	if(compareMasks(codeRS, meetCode1) || 
				compareMasks(codeRS, meetCode2) || 
				compareMasks(codeRS, meetCode3)){
		return MEET;
	}

    //else return overlap/intersects
    return OVERLAP;

}