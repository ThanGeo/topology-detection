#include "join_geometry_refinement.h"

//define topological masks for refinement
// a within b
boost::geometry::de9im::mask withinMask("T*F*FF***"); 
// a covered by b
boost::geometry::de9im::mask coveredbyMask("**F*TF***"); 
// a and b meet
boost::geometry::de9im::mask meetMask1("FT*******"); 
boost::geometry::de9im::mask meetMask2("F**T*****"); 
boost::geometry::de9im::mask meetMask3("F***T****"); 
// a and b are equal
boost::geometry::de9im::mask equalMask("T*F**FFF*"); 
// a and b are disjoint
boost::geometry::de9im::mask disjointMask("FF*FF****"); 




polygon loadPolygonGeometryBOOST(uint &recID, unordered_map<uint,unsigned long> &offsetMap, ifstream &fin){
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

int refinementWithIDs(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS){
	/*  BOOST GEOMETRY REFINEMENT TO CLEAR IF TOPOLOGICAL RELATIONSHIP IS INDEED MEET OR COVERED BY*/
	polygon boostPolygonR = loadPolygonGeometryBOOST(idA, offsetMapR, finR);
	polygon boostPolygonS = loadPolygonGeometryBOOST(idB, offsetMapS, finS);

	//equal
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, equalMask)){
    	return EQUAL;
    }

    //disjoint
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, disjointMask)){
    	return DISJOINT;
    }

	//inside
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, withinMask)){
    	return R_INSIDE_S;
    }
    if(boost::geometry::relate(boostPolygonS, boostPolygonR, withinMask)){
    	return S_INSIDE_R;
    }

	//covered by
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, coveredbyMask)){
    	return R_COVERED_BY_S;
    }
    if(boost::geometry::relate(boostPolygonS, boostPolygonR, coveredbyMask)){
    	return S_COVERED_BY_R;
    }

	//meet
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask1) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask2) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask3)){
    	return MEET;
    }

    

    //else return overlap
    return OVERLAP;   
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
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, coveredbyMask)){
    	return R_COVERED_BY_S;
    }
    if(boost::geometry::relate(boostPolygonS, boostPolygonR, coveredbyMask)){
    	return S_COVERED_BY_R;
    }

    //meet
    if(boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask1) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask2) || 
    	boost::geometry::relate(boostPolygonR, boostPolygonS, meetMask3)){
    	return MEET;
    }

    //else return overlap
    return OVERLAP;

}