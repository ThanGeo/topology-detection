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

#include <stdarg.h> /* for va_list - GEOS */

#include "topological/dataset_data.h"
#include "topological/binary_files.h"
#include "topological/join_geometry_refinement.h"

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

int (*create_topology_table_function_pointer)(Polygon *polA, Polygon *polB);
int (*refinement_function)(uint &idA, uint &idB, unordered_map<uint,unsigned long> &offsetMapR, unordered_map<uint,unsigned long> &offsetMapS, ifstream &finR, ifstream &finS);


//-----------------------------
//
//          UTILITY
//
//-----------------------------

void resetMetricParameters(){
        intermediateFilterTime = 0;
        refinementTime = 0;
        postMBRCandidates = 0;
        refinementCandidates = 0;

        result_count_map.insert(make_pair<int,int>(DISJOINT, 0));
        result_count_map.insert(make_pair<int,int>(OVERLAP, 0));
        result_count_map.insert(make_pair<int,int>(EQUAL, 0));
        result_count_map.insert(make_pair<int,int>(R_COVERED_BY_S, 0));
        result_count_map.insert(make_pair<int,int>(S_COVERED_BY_R, 0));
        result_count_map.insert(make_pair<int,int>(R_INSIDE_S, 0));
        result_count_map.insert(make_pair<int,int>(S_INSIDE_R, 0));
        result_count_map.insert(make_pair<int,int>(MEET, 0));

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
//          INITIALIZE
//
//-----------------------------

static void geos_msg_handler(const char* fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf (fmt, ap);
    va_end(ap);
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

        //no intermediate filter, so complete refinement
        refinement_function = &refinement_GEOS_WithIDs;
        
        //set ID type
        setIDtype();

        // initialize GEOS
        /* Send notice and error messages to the terminal */
        initGEOS(geos_msg_handler, geos_msg_handler);

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


        result_count_map.insert(make_pair<int,int>(DISJOINT, 0));
        result_count_map.insert(make_pair<int,int>(OVERLAP, 0));
        result_count_map.insert(make_pair<int,int>(EQUAL, 0));
        result_count_map.insert(make_pair<int,int>(R_COVERED_BY_S, 0));
        result_count_map.insert(make_pair<int,int>(S_COVERED_BY_R, 0));
        result_count_map.insert(make_pair<int,int>(R_INSIDE_S, 0));
        result_count_map.insert(make_pair<int,int>(S_INSIDE_R, 0));
        result_count_map.insert(make_pair<int,int>(MEET, 0));


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
//          CONNECTORS
//
//-----------------------------

void forwardCandidatePair(uint idA, uint idB){
        double timing;

        bool refFlag = false;
        int result=-1;
        postMBRCandidates++;

        refinementCandidates++;  
        if(REFINEMENT){
                //needs refinement
                timing = omp_get_wtime();
                result = (*refinement_function)(idA, idB, offsetMapR, offsetMapS, finR, finS);
        
                //save result
                result_count_map.at(result) += 1;

                refinementTime += (omp_get_wtime() - timing);
        }
}

#endif
