#ifndef SPATIAL_LIB_DEF_H
#define SPATIAL_LIB_DEF_H

#include <vector>
#include <sys/types.h>
#include <string>
#include <unordered_map>
#include <fstream>

#include "data.h"

namespace spatial_lib
{
    typedef enum QueryResult {
        TRUE_NEGATIVE,
        TRUE_HIT,
        INCONCLUSIVE,
    } QueryResultT;

    typedef struct PolygonStats {
        std::unordered_map<uint,uint> THappearancesMapR, THappearancesMapS; // recID -> appearances in true hit pairs
        std::unordered_map<uint,uint> RFappearancesMapR, RFappearancesMapS; // recID -> appearances in refinement pairs
        std::unordered_map<uint,uint> polVerticesMapR, polVerticesMapS; // recID -> number of vertices
    } PolygonStatsT;

    typedef struct QueryOutput {
        // for regular query rsesults
        int queryResults;
        // for topology relations results
        std::unordered_map<int,uint> topologyRelationsResultMap;
        // statistics
        int postMBRFilterCandidates;
        int refinementCandidates;
        int trueHits;
        int trueNegatives;
        // times
        double totalTime;
        double mbrFilterTime;
        double iFilterTime;
        double refinementTime;
        // on the fly april
        uint rasterizationsDone;
        // development only
        PolygonStatsT polygonStats;
    } QueryOutputT;

    // global query output variable
    extern QueryOutputT g_queryOutput;

    // just for testing, experiments about scalability
    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1,T2> &p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);

            // Mainly for demonstration purposes, i.e. works but is overly simple
            // In the real world, use sth. like boost.hash_combine
            return h1 ^ h2;  
        }
    };
    typedef struct ScalabilityContainer {
        std::string bucketContentsPath;
        std::string bucketDataPath;
        uint numberOfBuckets;
        std::unordered_map<uint,uint> bucketPolygonCount;
        std::unordered_map<uint,std::pair<uint,uint>> bucketIDTobucketRangeMap;
        std::unordered_map<std::pair<uint,uint>,uint,pair_hash> pairToBucketIDMap;
        std::unordered_map<uint,double> bucketIfilterTime;
        std::unordered_map<uint,double> bucketRefinementTime;
        std::unordered_map<uint,uint> bucketInconclusiveCount;
    } ScalabilityContainerT;

    extern ScalabilityContainerT g_scalContainer;
    uint getBucketOfPair(uint idR, uint idS);
    void addIFilterTimeToBucket(uint bucketID, double time);
    void addRefinementTimeToBucket(uint bucketID, double time);

    void countInconclusiveToBucket(uint bucketID);
    typedef enum QueryType{
        RANGE,
        Q_INTERSECT,
        Q_INSIDE,
        Q_DISJOINT,
        Q_EQUAL,
        Q_MEET,
        Q_CONTAINS,
        Q_COVERS,
        Q_COVERED_BY,
        Q_FIND_RELATION,    // find what type of topological relation is there
    }QueryTypeE;


    typedef struct Section {
        uint sectionID;
        // axis position indexes
        uint i,j;
        //objects that intersect this MBR will be assigned to this area
        double interestxMin, interestyMin, interestxMax, interestyMax;
        // double normInterestxMin, normInterestyMin, normInterestxMax, normInterestyMax;
        //this MBR defines the rasterization area (widened interest area to include intersecting polygons completely)
        double rasterxMin, rasteryMin, rasterxMax, rasteryMax;
        // double normRasterxMin, normRasteryMin, normRasterxMax, normRasteryMax;
        // APRIL data
        uint objectCount = 0;
        std::unordered_map<uint, spatial_lib::AprilDataT> aprilData;
    } SectionT;

    typedef struct DataspaceInfo {
        double xMinGlobal, yMinGlobal, xMaxGlobal, yMaxGlobal;  // global bounds based on dataset bounds
        double xExtent, yExtent;
    } DataspaceInfoT;

    typedef struct Dataset{
        spatial_lib::DataTypeE dataType;
        std::string path;
        std::string offsetMapPath;
        // derived from the path
        std::string datasetName;
        // as given by arguments and specified by datasets.ini config file
        std::string nickname;
        // map: recID -> vector data (polygon, linestring etc.)
        std::unordered_map<uint, spatial_lib::VectorDataT> vectorData;
        // double xMinGlobal, yMinGlobal, xMaxGlobal, yMaxGlobal;  // global bounds based on dataset bounds
        DataspaceInfoT dataspaceInfo;
        /**
         * Approximations
        */
        ApproximationTypeE approxType;
        // APRIL
        spatial_lib::AprilConfigT aprilConfig;
        std::unordered_map<uint, SectionT> sectionMap;           // map: k,v = sectionID,(unordered map of k,v = recID,aprilData)
        std::unordered_map<uint,std::vector<uint>> recToSectionIdMap;         // map: k,v = recID,vector<sectionID>: maps recs to sections          
    }DatasetT;

    typedef struct Query{
        spatial_lib::QueryTypeE type;
        int numberOfDatasets;
        spatial_lib::DatasetT R;         // R: left dataset
        spatial_lib::DatasetT S;         // S: right dataset
        bool boundsSet = false;
        // double xMinGlobal, yMinGlobal, xMaxGlobal, yMaxGlobal;  // global bounds based on dataset bounds
        DataspaceInfoT dataspaceInfo;
    }QueryT;


    void resetQueryOutput();
    void setupScalabilityTesting();
    void countAPRILResult(int result);
    void countResult();
    void countTopologyRelationResult(int relation);

    void addAprilDataToApproximationDataMap(DatasetT &dataset, uint sectionID, uint recID, AprilDataT* aprilData);

    void addObjectToSectionMap(DatasetT &dataset, uint sectionID, uint recID);
    /**
     * @brief returns the APRIL data of an object based on section and rec IDs from a given dataset
     * 
     * @param dataset 
     * @param sectionID 
     * @param recID 
     * @return AprilDataT* 
     */
    AprilDataT* getAprilDataBySectionAndObjectIDs(Dataset &dataset, uint sectionID, uint recID);

    std::unordered_map<uint,unsigned long> loadOffsetMap(std::string &offsetMapPath);
    bg_polygon loadPolygonFromDiskBoostGeometry(uint recID, std::ifstream &fin, std::unordered_map<uint,unsigned long> &offsetMap);

    /**
     * @brief Returns a pointer to the section that corresponds to the argument ID from the dataset
     * 
     * @param dataset 
     * @param sectionID 
     * @return SectionT* 
     */
    SectionT* getSectionByID(DatasetT &dataset, uint sectionID);

    /**
     * @brief Returns a list of pointers to the sections that intersect with the given MBR
     * 
     * @return std::vector<uint> 
     */
    std::vector<SectionT*> getSectionsOfMBR(spatial_lib::DatasetT &dataset, double xMin, double yMin, double xMax, double yMax);

    inline uint getSectionIDFromIdxs(uint i, uint j, uint partitionsNum) {
        return (j * partitionsNum) + i;
    }

    /**
     * @brief returns a list of the common section ids between two objects of two datasets
     * 
     * @param datasetR 
     * @param datasetS 
     * @param idR 
     * @param idS 
     * @return std::vector<uint> 
     * 
     */
    std::vector<uint> getCommonSectionIDsOfObjects(Dataset &datasetR, Dataset &datasetS, uint idR, uint idS);

    /**
     * @brief Get the Section IDs Of an object by ID
     * 
     * @param dataset 
     * @param id 
     * @return std::vector<uint>* 
     */
    std::vector<uint>* getSectionIDsOfObject(DatasetT &dataset, uint id);


    void printBoostPolygon(bg_polygon &polygon);
}

#endif