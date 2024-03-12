#ifndef PIPE_H
#define PIPE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <fstream>
#include <future>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include "omp.h"
#include <functional>
#include <array>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>

#include "./APRIL/april-main.h"

using namespace std;


Dataset rasterIntervalsR;
Dataset rasterIntervalsS;

//geometry filenames
string geometryFileNameR;
string geometryFileNameS;
//offset maps for binary geometries
unordered_map<uint,unsigned long> offsetMapR;
unordered_map<uint,unsigned long> offsetMapS;
//binary geometry files
ifstream finR;
ifstream finS;

string argument1, argument2;

string result_filename = "results_pairs.csv";

double intermediateFilterTime = 0;
double refinementTime = 0;

int postMBRCandidates = 0;
int refinementCandidates = 0;

clock_t timer;

int (*find_relation_function_ptr)(Polygon *polA, Polygon *polB);
int (*refinement_function_ptr)(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS);
void (*main_topology_function_ptr)(uint &idA, uint &idB);

//-----------------------------
//
//          UTILITY
//
//-----------------------------

void initTopologyResultMap() {
        result_count_map.clear();

        result_count_map.insert(make_pair<int,int>(DISJOINT, 0));
        result_count_map.insert(make_pair<int,int>(OVERLAP, 0));
        result_count_map.insert(make_pair<int,int>(EQUAL, 0));
        result_count_map.insert(make_pair<int,int>(R_COVERED_BY_S, 0));
        result_count_map.insert(make_pair<int,int>(R_COVERS_S, 0));
        result_count_map.insert(make_pair<int,int>(R_CONTAINED_IN_S, 0));
        result_count_map.insert(make_pair<int,int>(R_CONTAINS_S, 0));
        result_count_map.insert(make_pair<int,int>(MEET, 0));
        result_count_map.insert(make_pair<int,int>(CROSSES, 0));
}

void resetMetricParameters(){
        intermediateFilterTime = 0;
        refinementTime = 0;
        postMBRCandidates = 0;
        refinementCandidates = 0;

        initTopologyResultMap();

}

void printSections(){
        cout << "DATA SPACE (" << DATA_SPACE.xMin << "," << DATA_SPACE.yMin << "),(" << DATA_SPACE.xMax << "," << DATA_SPACE.yMax << ") " << DATA_SPACE.sections.size() << " sections: " << endl;

        for(auto &it : DATA_SPACE.sections){
                cout << "       Section " << it.sectionID << " (" << it.x << "," << it.y << "): " << endl;
                cout << "               interest area: (" << it.interestxMin << "," << it.interestyMin << "),(" << it.interestxMax << "," << it.interestyMax << ")" << endl;
                cout << "               rasterin area: (" << it.rasterxMin << "," << it.rasteryMin << "),(" << it.rasterxMax << "," << it.rasteryMax << ")" << endl;
        }
}

void saveResultPair(uint &idA, uint &idB){
        ofstream fout(result_filename, fstream::out | ios_base::binary | fstream::app);

        fout << idA << " " << idB << endl; 

        fout.close();
}


//-----------------------------
//
//          MAIN FUNCTIONS
//
//-----------------------------

void gatewayAllRelations(uint &idA, uint &idB) {
        double timing;
        bool refFlag = false;
        int result=-1;
        postMBRCandidates++;
        if(INTERMEDIATE_FILTER){
                timing = omp_get_wtime();
                //find common sections
                vector<uint> commonSections = DATA_SPACE.getCommonSectionsIDOfObjects(idA, idB);
                //check in each common section
                for(auto &secID : commonSections){
                        result = (*find_relation_function_ptr)(rasterIntervalsR.getPolygonByIDAndSection(secID, idA), rasterIntervalsS.getPolygonByIDAndSection(secID, idB));
                        // if(result == EQUAL || result == R_COVERED_BY_S || result == R_COVERS_S || result == MEET){
                        if(result >= EQUAL){
                                refFlag = true;
                        }
                        if(result == -1){
                                //error
                                cout << "Error for " << idA << "," << idB << ": could not identify topological relationship." << endl;
                                exit(0);
                        }

                }
                //if it isn't marked for refinement, return
                if(!refFlag){
                        //save result
                        result_count_map.at(result) += 1;
                        intermediateFilterTime += (omp_get_wtime() - timing);
                        return;
                }                           
                intermediateFilterTime += (omp_get_wtime() - timing);
        }
        refinementCandidates++;  
        if(REFINEMENT){
                //needs refinement
                timing = omp_get_wtime();
                result = (*refinement_function_ptr)(idA, idB, offsetMapR, offsetMapS, finR, finS);

                //save result
                result_count_map.at(result) += 1;

                refinementTime += (omp_get_wtime() - timing);
        }
}

void gatewaySpecificRelation(uint &idA, uint &idB) {
        double timing;
        int result = -1;
        postMBRCandidates++;
        // cout << idA << "," << idB << endl;
        if(INTERMEDIATE_FILTER){
                timing = omp_get_wtime();
                //find common sections
                vector<uint> commonSections = DATA_SPACE.getCommonSectionsIDOfObjects(idA, idB);
                // cout << commonSections.size() << " common sections." << endl;
                //check in each common section
                for(auto &secID : commonSections){
                        /**
                         * -1 : error
                         * 0 : true negative
                         * 1 : true positive
                         * 2 : inconclusive (needs refinement)
                        */        
                        result = (*find_relation_function_ptr)(rasterIntervalsR.getPolygonByIDAndSection(secID, idA), rasterIntervalsS.getPolygonByIDAndSection(secID, idB));
                        if(result == -1){
                                //error
                                cout << "Error for " << idA << "," << idB << ": could not identify topological relationship." << endl;
                                exit(0);
                        }
                }
                //if it isn't marked for refinement, return
                if(result == 1){
                        //equal, save result
                        result_count_map.at(TOPOLOGY_PREDICATE) += 1;
                        intermediateFilterTime += (omp_get_wtime() - timing);
                        return;
                } else if (result == 0) {
                        // not equal
                        intermediateFilterTime += (omp_get_wtime() - timing);
                        return;
                }                 
                intermediateFilterTime += (omp_get_wtime() - timing);
        }
        refinementCandidates++;  
        if(REFINEMENT){
                //needs refinement
                timing = omp_get_wtime();
                result = (*refinement_function_ptr)(idA, idB, offsetMapR, offsetMapS, finR, finS);
                //save result
                if (result == 1) {
                        result_count_map.at(TOPOLOGY_PREDICATE) += 1;
                }
                refinementTime += (omp_get_wtime() - timing);
        }
}

//-----------------------------
//
//          INITIALIZE
//
//-----------------------------

void setTopologyFunctions() {
        if (TOPOLOGY_PREDICATE == NONE) {
                /**
                 * ALL RELATIONS
                */
                find_relation_function_ptr = &checkAllRelations;
                main_topology_function_ptr = &gatewayAllRelations;
                refinement_function_ptr = &refinementAllTopologyRelations;
                return;
        }
        /**
         * A SINGLE SPECIFIC RELATION
        */

        find_relation_function_ptr = &checkSpecificRelation;
        main_topology_function_ptr = &gatewaySpecificRelation;
        switch (TOPOLOGY_PREDICATE) {
                case DISJOINT:
                        checkPredicateFunctionPtr = &checkTopologyDisjoint;
                        refinement_function_ptr = &refinementDisjoint;
                        break;
                case EQUAL:
                        checkPredicateFunctionPtr = &checkTopologyEqual;
                        refinement_function_ptr = &refinementEqual;
                        break;
                case OVERLAP:
                        // Use regular APRIL for this instead of the matrix, override find_relation_function
                        find_relation_function_ptr = &checkTopologyIntersects;
                        refinement_function_ptr = &refinementOverlap;
                        break;
                case MEET:
                        checkPredicateFunctionPtr = &checkTopologyMeet;
                        refinement_function_ptr = &refinementMeet;
                        break;
                case CROSSES:
                        checkPredicateFunctionPtr = &checkTopologyCrosses;
                        refinement_function_ptr = &refinementCrosses;
                        break;
                case R_CONTAINED_IN_S:
                        checkPredicateFunctionPtr = &checkTopologyRContainedInS;
                        refinement_function_ptr = &refinementRcontainedInS;
                        break;
                case R_CONTAINS_S:
                        checkPredicateFunctionPtr = &checkTopologyRContainsS;
                        refinement_function_ptr = &refinementRcontainsS;
                        break;
                case R_COVERED_BY_S:
                        checkPredicateFunctionPtr = &checkTopologyRCoveredByS;
                        refinement_function_ptr = &refinementRcoveredByS;
                        break;
                case R_COVERS_S:
                        checkPredicateFunctionPtr = &checkTopologyRCoversS;
                        refinement_function_ptr = &refinementRcoversS;
                        break;
        }
}

void initialize(string &arg1, string &arg2){
        string info_message = "";
        argument1 = arg1;
        argument2 = arg2;

        //build file paths
        buildFilePaths(arg1, arg2);

        geometryFileNameR = getBinaryGeometryFilename(0);
        geometryFileNameS = getBinaryGeometryFilename(1);
        offsetMapR = loadOffsetMap(0);
        offsetMapS = loadOffsetMap(1);
        finR.open(geometryFileNameR, fstream::in | ios_base::binary);
        finS.open(geometryFileNameS, fstream::in | ios_base::binary);

        if(!finR){
                cout << "Error opening: " << geometryFileNameR << endl;
                exit(-1);
        }
        if(!finS){
                cout << "Error opening: " << geometryFileNameS << endl;
                exit(-1);
        }

        //set universal coordinates
        // cout << "Setting universal min/max..." << endl;

        //JOIN 
        if((argument1.at(0) == 'T' && argument2.at(0) == 'T')){
                getUniversalCoordinates(0);
        }else{
                string continent = argument1.substr(argument1.find("_") + 1);
                if(continent == "Oceania"){
                        getUniversalCoordinates(1);
                }else if(continent == "Asia"){
                        getUniversalCoordinates(2);
                }else if(continent == "Europe"){
                        getUniversalCoordinates(3);
                }else if(continent == "NorthAmerica"){
                        getUniversalCoordinates(4);                        
                }else if(continent == "Africa"){
                        getUniversalCoordinates(5);                        
                }else if(continent == "SouthAmerica"){
                        getUniversalCoordinates(6);                        
                }
        }


        //resize the data space and create its sections
        DATA_SPACE.resize();
        //pass it to the data space object
        DATA_SPACE.setUniversalCoordinates();

        // cout << fixed << setprecision(6) << "Done: " << universalMinX << " " << universalMinY << "," << universalMaxX << " " << universalMaxY << endl;
        //return to begining
        finR.seekg(0, ios::beg);
        finS.seekg(0, ios::beg);

        // ofstream fout(result_filename, fstream::out | ios_base::binary);
        // fout.close();

        //create the data space sections
        for(int i=0; i<H; i++){
                for(int j=0; j<H; j++){
                        Section sec(i,j);
                        //cout << "created section " << sec.sectionID << " with interest area: " << sec.interestxMin << "," << sec.interestyMin << " and " << sec.interestxMax << "," << sec.interestyMax << endl; 
                        DATA_SPACE.sections.at(sec.sectionID) = sec;
                }
        }

        info_message += "-Datasets (R & S): \t\t" + argument1 + " & " + argument2 + "\n";
        info_message += "-Partitions: \t\t\t" + to_string(H) + "\n";
        info_message += "-Data: \t\t\t\tpolygon-polygon\n";

        //set the join function (based on whether it has to join compressed APRIL intervals)        
        if(COMPRESSION == 0){
                info_message += "-Compression: \t\t\tno\n";
                //UNCOMPRESSED APRIL
                if(DIFF_GRANULARITY_FIXED){
                        //DIFFERENT GRANULARITY
                        info_message += "-Granularity: \t\t\t16-" + to_string(DESIGNATED_ORDER) + "\n";
                }else{
                        //SAME GRANULARITY
                        setTopologyFunctions();
                        info_message += "-Granularity: \t\t\t16-16\n";
                }
        }else{
                info_message += "-Compression: \t\t\tyes\n";
                //COMPRESSED APRIL
                if(DIFF_GRANULARITY_FIXED){
                        //DIFFERENT GRANULARITY
                        info_message += "-Granularity: \t\t\t16-" + to_string(DESIGNATED_ORDER) + "\n";
                }else{
                        //SAME GRANULARITY
                        info_message += "-Granularity: \t\t\t16-16\n";
                }
        }
        
        //set ID type
        setIDtype();

        //reset result file
        ofstream fout(result_filename, fstream::out | ios_base::binary);
        fout.close();

        //check if APRIL dirs exists
        DIR* dirUncompressed = opendir("APRIL/interval_data/uncompressed/");
        if(dirUncompressed) {
                /* Directory exists. */
                closedir(dirUncompressed);
        }else if(ENOENT == errno) {
                /* Directory does not exist. */
                mkdir("APRIL/interval_data/uncompressed/", 0700);
        }else{
                /* opendir() failed for some other reason. */
                cout << "Init error: Cannot open directory 'APRIL/interval_data/uncompressed/'" << endl;
                exit(-1);
        }
        DIR* dirCompressed = opendir("APRIL/interval_data/compressed/");
        if(dirCompressed) {
                /* Directory exists. */
                closedir(dirCompressed);
        }else if(ENOENT == errno) {
                /* Directory does not exist. */
                mkdir("APRIL/interval_data/compressed/", 0700);
        }else{
                /* opendir() failed for some other reason. */
                cout << "Init error: Cannot open directory 'APRIL/interval_data/compressed/'" << endl;
                exit(-1);
        }

        // init result map
        initTopologyResultMap();
        
        cout << "***************************************************" << endl;
        cout << "*************** Pipeline configuration: **************" << endl;
        cout << info_message << "***************************************************" << endl;
        cout << "MBR-join";
        if(INTERMEDIATE_FILTER){
                cout << " -> Intermediate Filter";
        }
        if(REFINEMENT){
                cout << " -> Refinement Stage";
        }
        cout << " -> Result" << endl;
        cout << "***************************************************" << endl << endl;

}

//-----------------------------
//
//          ENABLERS
//
//-----------------------------

void enableIntermediateFilter(string &argument1, string &argument2){        
        rasterIntervalsR.argument = argument1;
        rasterIntervalsR.letterID = "A";
        rasterIntervalsS.argument = argument2;
        rasterIntervalsS.letterID = "B";

        cout << "Loading datasets' APRIL... " << endl;
        loadApproximations(rasterIntervalsR, argument1, 0);
        loadApproximations(rasterIntervalsS, argument2, 1);
        cout << "Finished." << endl;
}


void initiateRasterIntervalsCreation(string &argument1, string &argument2){
        // initializeRasterizationOnOpenGL();
        createApproximations(argument1, 0);
        createApproximations(argument2, 1);
}


//-----------------------------
//
//          CONNECTORS
//
//-----------------------------

void forwardCandidatePair(uint idA, uint idB){
        
        (*main_topology_function_ptr)(idA, idB);
}

#endif
