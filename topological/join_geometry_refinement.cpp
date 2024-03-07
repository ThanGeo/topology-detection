#include "join_geometry_refinement.h"

//define topological masks for refinement
char equalMask[] = "T*F**FFF*";
char disjointMask[] = "FF*FF****";
char meetMask1[] = "FT*******";
char meetMask2[] = "F**T*****";
char meetMask3[] = "F***T****";
char withinMask[] = "T*F*FF***";
char coveredbyMask1[] = "**F*TF***";
char coveredbyMask2[] = "**F*TF***";
char coveredbyMask3[] = "**F*TF***";
char coveredbyMask4[] = "**F*TF***";
vector<char*> coveredByMaskList = {
				coveredbyMask1,
				coveredbyMask2,
				coveredbyMask3,
				coveredbyMask4};
char intersectsMask[] = "T********";

GEOSGeometry* MakePoly(std::vector<double> const& xCoords, std::vector<double> const& yCoords)
{
	
	size_t seqSize = xCoords.size();
	GEOSCoordSequence* seq = GEOSCoordSeq_copyFromArrays(
			xCoords.data(),
			yCoords.data(),
			NULL,  /* Zs */
			NULL,  /* Ms */
			seqSize);
	
	GEOSGeometry *ring = GEOSGeom_createLinearRing(seq);
	GEOSGeometry *polygon = GEOSGeom_createPolygon(ring, NULL, 0);

	return polygon;
}

GEOSGeometry* loadPolygonGeometryGEOS(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
	// polygon pol;
	int readID;
	int vertexCount, polygonCount;
	double x,y;

	std::vector<double> x_coords;
	std::vector<double> y_coords;

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
			x_coords.emplace_back(x);
			y_coords.emplace_back(y);
		}
	}

	GEOSGeometry* polygon = MakePoly(x_coords, y_coords);
	return polygon;
}

int refinement_DE9IM_WithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT TO DETECT TOPOLOGICAL RELATIONSHIP */
	// polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	// polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);

    
    //disjoint
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, disjointMask)){
    // 	return DISJOINT;
    // }

    // //equal
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, equalMask)){
    // 	return EQUAL;
    // }

    // //inside
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, withinMask)){
    // 	return R_INSIDE_S;
    // }
    // if(boost::geometry::relate(boostPolygonS, boostPolygonR, withinMask)){
    // 	return S_INSIDE_R;
    // }

    // //covered by
    // // if(boost::geometry::relate(boostPolygonR, boostPolygonS, coveredbyMask)){
    // // 	return R_COVERED_BY_S;
    // // }
    // // if(boost::geometry::relate(boostPolygonS, boostPolygonR, coveredbyMask)){
    // // 	return S_COVERED_BY_R;
    // // }

	// for(auto &it: coveredByMaskList){
	// 	if(boost::geometry::relate(boostPolygonR, boostPolygonS, it)){
    // 		return R_COVERED_BY_S;
	// 	}
	// 	if(boost::geometry::relate(boostPolygonS, boostPolygonR, it)){
	// 		return S_COVERED_BY_R;
	// 	}
	// }

    // //meet
    // if(boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask1) || 
    // 	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask2) || 
    // 	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask3)){
    // 	return MEET;
    // }

    // //else return overlap/intersects
    // return OVERLAP;

}



int refinement_GEOS_WithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	// printf("Checking objects %d and %d \n", idA, idB);
	// load geos objects
	GEOSGeometry* geosPolygonR = loadPolygonGeometryGEOS(idA, offsetMapR, finR);
	GEOSGeometry* geosPolygonS = loadPolygonGeometryGEOS(idB, offsetMapS, finS);
	int returnRelation;

	// printf("loaded GEOSGeometries\n");

	// use the relate function to get the DE-9IM
    char* resRS = GEOSRelate(geosPolygonR, geosPolygonS);
    char* resSR = GEOSRelate(geosPolygonS, geosPolygonR);

    // disjoint
    if (GEOSRelatePatternMatch(resRS, disjointMask) == 1) {
		returnRelation = DISJOINT;
		goto FREE_AND_EXIT;
    }
	// equal
    if (GEOSRelatePatternMatch(resRS, equalMask) == 1) {
		returnRelation = EQUAL;
		goto FREE_AND_EXIT;
    }
    // meet
    if (GEOSRelatePatternMatch(resRS, meetMask1) == 1 || GEOSRelatePatternMatch(resRS, meetMask2) == 1 || GEOSRelatePatternMatch(resRS, meetMask3) == 1) {
		returnRelation = MEET;
		goto FREE_AND_EXIT;
    }
    // within
    if (GEOSRelatePatternMatch(resRS, withinMask) == 1) {
		returnRelation = R_INSIDE_S;
		goto FREE_AND_EXIT;
    }
	if (GEOSRelatePatternMatch(resSR, withinMask) == 1) {
		returnRelation = S_INSIDE_R;
		goto FREE_AND_EXIT;
    }
    // covers
	for(auto &it: coveredByMaskList) {
		if (GEOSRelatePatternMatch(resRS, it) == 1) {
			returnRelation = R_COVERED_BY_S;
			goto FREE_AND_EXIT;
		}
		if (GEOSRelatePatternMatch(resSR, it) == 1) {
			returnRelation = S_COVERED_BY_R;
			goto FREE_AND_EXIT;
		}
	}
    // intersects (TODO: make just returning intersects instead of checking)
    if (GEOSRelatePatternMatch(resRS, intersectsMask) == 1) {
		returnRelation = OVERLAP;
		goto FREE_AND_EXIT;
    }
	
FREE_AND_EXIT:
	GEOSGeom_destroy(geosPolygonR);
	GEOSGeom_destroy(geosPolygonS);
	return returnRelation;
}

