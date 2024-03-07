#include "topology.h"

//same granularity
int create_topology_table_uncompressed(Polygon *polA, Polygon *polB){

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


	if(AA == 0){
		//disjoint
		return DISJOINT;
	}

	if((AA == 2 && FF == 2) || (AA == 2 && (!polA->F && !polB->F))){
		//equal - needs refinement
		// return NEEDS_REFINEMENT;
		return EQUAL;
	}

	if(AA == 1){
		if(AF >= 1 || FA >= 1){
			//overlap
			return OVERLAP;
		}else if(AF + FA <= 0){			
			//meet - needs refinement
			// return NEEDS_REFINEMENT;
			return MEET;
		}
	}

	if(AA == 3){
		if(AF == 3){
			//R inside S
			return R_INSIDE_S;
		}else{
			//R covered_by S - needs refinement
			// return NEEDS_REFINEMENT;
			return R_COVERED_BY_S;
		}
	}

	if(AA == 4){
		if(FA == 4){
			//S inside R
			return S_INSIDE_R;
		}else{
			//S covered_by R - needs refinement
			// return NEEDS_REFINEMENT;
			return S_COVERED_BY_R;
		}
	}

	//error: could not identify topological relationship
	return -1;

}