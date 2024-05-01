#include "config.h"

/**
 * Global configuration variable
*/
ConfigT g_config;

static std::unordered_map<std::string, spatial_lib::QueryTypeE> queryStringToIntMap = {
    {"range",spatial_lib::RANGE},
    {"intersect",spatial_lib::Q_INTERSECT},
    {"inside",spatial_lib::Q_INSIDE},
    {"disjoint",spatial_lib::Q_DISJOINT},
    {"equal",spatial_lib::Q_EQUAL},
    {"meet",spatial_lib::Q_MEET},
    {"contains",spatial_lib::Q_CONTAINS},
    {"covered_by",spatial_lib::Q_COVERED_BY},
    {"covers",spatial_lib::Q_COVERS},
    {"find_relation",spatial_lib::Q_FIND_RELATION},
    };

static std::unordered_map<std::string, spatial_lib::PipelineSettingE> pipelineSettingsStringToIntMap = {
    {"ST2", spatial_lib::P_ST2},
    {"ST3", spatial_lib::P_ST3},
    {"OP2", spatial_lib::P_OP2},
    {"OP3", spatial_lib::P_OP3},
    {"SCALABILITY_OP2", spatial_lib::P_SCALABILITY_OP2},
    {"SCALABILITY_OP3", spatial_lib::P_SCALABILITY_OP3},
    {"OTF", spatial_lib::P_OTF},
};

static std::string iFilterTypeIntToText(int val){
    switch(val) {
        case spatial_lib::IF_NONE: return "NONE";
        case spatial_lib::IF_APRIL_STANDARD: return "APRIL STANDARD";
        case spatial_lib::IF_APRIL_OPTIMIZED: return "APRIL FIND RELATION";
        case spatial_lib::IF_APRIL_OTF: return "APRIL ON THE FLY";
        case spatial_lib::IF_APRIL_SCALABILITY: return "APRIL FIND RELATION SCALABILITY TEST";
    }
}

static std::string mbrFilterTypeIntToText(int val){
    switch(val) {
        case spatial_lib::MF_STANDARD: return "STANDARD MBR FILTER";
        case spatial_lib::MF_OPTIMIZED: return "FIND RELATION MBR FILTER";
        case spatial_lib::MF_SCALABILITY: return "FIND RELATION SCALABILITY MBR FILTER";
    }
}

static std::string pipelineSettingIntToText(spatial_lib::PipelineSettingE val) {
    switch(val) {
        case spatial_lib::P_ST2: return "ST2";
        case spatial_lib::P_OP2: return "OP2";
        case spatial_lib::P_ST3: return "ST3";
        case spatial_lib::P_OP3: return "OP3";
        case spatial_lib::P_SCALABILITY_OP2: return "SCALABILITY_OP2";
        case spatial_lib::P_SCALABILITY_OP3: return "SCALABILITY_OP3";
        case spatial_lib::P_OTF: return "OTF";
    }
}

static std::string dataComboIntToText(int val){
    switch(val) {
        case spatial_lib::POLYGON_POLYGON: return "POLYGON_POLYGON";
        case spatial_lib::POLYGON_LINESTRING: return "POLYGON_LINESTRING";
    }
}

void printConfig() {
    std::cout << "---------- Configuration ----------" << std::endl;
    std::cout << "Configuration file path: \n\t" << g_config.dirPaths.configFilePath << std::endl;
    std::cout << "Query: " << std::endl;
    std::cout << "\ttype: " << queryTypeIntToText(g_config.queryData.type) << std::endl;
    std::cout << "\tR: " << g_config.queryData.R.path << " (" << dataTypeIntToText(g_config.queryData.R.dataType) <<  ")" << std::endl;
    if (g_config.queryData.numberOfDatasets == 2) {
        std::cout << "\tS: " << g_config.queryData.S.path << " (" << dataTypeIntToText(g_config.queryData.S.dataType) <<  ")" << std::endl;
    }
    std::cout << "Pipeline: " << std::endl;
    std::cout << "\tMBR Filter: " << mbrFilterTypeIntToText(g_config.pipeline.MBRFilterType) << std::endl;
    std::cout << "\tIntermediate Filter: " << iFilterTypeIntToText(g_config.pipeline.iFilterType) << std::endl;
    std::cout << "\tRefinement: " << g_config.pipeline.RefinementEnabled << std::endl;

    std::cout << "------------------------------------" << std::endl;
}

std::string getConfigSettingsStr() {
    std::string str = "";
    str += "Query: " + queryTypeIntToText(g_config.queryData.type) + "|";
    str += "Datasets: " + g_config.queryData.R.nickname + "," + g_config.queryData.S.nickname + "|";
    // str += "Pipeline: " + mbrFilterTypeIntToText(g_config.pipeline.MBRFilterType) + "->" + iFilterTypeIntToText(g_config.pipeline.iFilterType) + "->Refinement\n";
    str += "Pipeline: " + pipelineSettingIntToText(g_config.pipeline.setting) + "\n";
    return str;
}

static bool verifyIntermediateFilterType(std::string queryType) {
    if (auto it = queryStringToIntMap.find(queryType) != queryStringToIntMap.end()) {
        return true;
    }
    return false;
}


static bool verifyQuery(QueryStatementT *queryStmt) {
    if (queryStmt->queryType == "") {
        log_err("No query type specified.");
        return false;
    }
    
    auto itqt = queryStringToIntMap.find(queryStmt->queryType);
    // verify type
    if (itqt == queryStringToIntMap.end()) {
        log_err("Failed when verifying query type.");
        // todo: print available query types strings
        return false;
    }

    // check if its a scalability experiment and verify
    if ((g_config.pipeline.setting == spatial_lib::P_SCALABILITY_OP2 || g_config.pipeline.setting == spatial_lib::P_SCALABILITY_OP3) && 
        (itqt->second != spatial_lib::Q_FIND_RELATION || 
        queryStmt->datasetNicknameR != "O5EU" || queryStmt->datasetNicknameS != "O6EU")) {

        log_err("Scalability pipeline is available only for find relation queries and datasets O5EU,O6EU");
        return false;
    }
    
    // verify OTF setting
    if (g_config.pipeline.setting == spatial_lib::P_OTF && itqt->second != spatial_lib::Q_FIND_RELATION) {
        log_err("On-the-fly rasterization during query evaluation has only been implemented for 'find relation' queries.");
        return false;
    }

    return true;
}

bool verifyAndbuildQuery(QueryStatementT *queryStmt) {
    // verify
    if (!verifyQuery(queryStmt)) {
        log_err("Failed when verifying query config.");
        return false;
    }

    // build
    auto it = queryStringToIntMap.find(queryStmt->queryType);
    spatial_lib::QueryT query;
    query.type = (spatial_lib::QueryTypeE) it->second;
    query.numberOfDatasets = queryStmt->datasetCount;
    query.boundsSet = queryStmt->boundsSet;
    spatial_lib::DataspaceInfoT dataspaceInfo;
    dataspaceInfo.xMinGlobal = queryStmt->xMinGlobal;
    dataspaceInfo.yMinGlobal = queryStmt->yMinGlobal;
    dataspaceInfo.xMaxGlobal = queryStmt->xMaxGlobal;
    dataspaceInfo.yMaxGlobal = queryStmt->yMaxGlobal;
    dataspaceInfo.xExtent = dataspaceInfo.xMaxGlobal - dataspaceInfo.xMinGlobal;
    dataspaceInfo.yExtent = dataspaceInfo.yMaxGlobal - dataspaceInfo.yMinGlobal;
    query.dataspaceInfo = dataspaceInfo;
    // printf("Global bounds: (%f,%f),(%f,%f)\n", query.xMinGlobal, query.yMinGlobal, query.xMaxGlobal, query.yMaxGlobal);
    
    // verify dataset paths
    if (!verifydatasetNickname(queryStmt->datasetPathR) ){
        log_err("Failed when verifying dataset R path.");
        return false;
    }
    if (!verifydatasetNickname(queryStmt->datasetPathS) ){
        log_err("Failed when verifying dataset S path.");
        return false;
    }
    // verify offset map paths
    if (!verifydatasetNickname(queryStmt->offsetMapPathR) ){
        log_err("Failed when verifying offset map path of R.");
        return false;
    }
    if (!verifydatasetNickname(queryStmt->offsetMapPathS) ){
        log_err("Failed when verifying offset map path of S.");
        return false;
    }

    // build dataset objects
    spatial_lib::DatasetT R,S;
    R.dataspaceInfo = query.dataspaceInfo;
    S.dataspaceInfo = query.dataspaceInfo;
    R.dataType = queryStmt->datatypeR;
    S.dataType = queryStmt->datatypeS;
    R.path = queryStmt->datasetPathR;
    R.nickname = queryStmt->datasetNicknameR;
    S.path = queryStmt->datasetPathS;
    S.nickname = queryStmt->datasetNicknameS;
    R.offsetMapPath = queryStmt->offsetMapPathR;
    S.offsetMapPath = queryStmt->offsetMapPathS;
    R.datasetName = getDatasetNameFromPath(R.path);
    S.datasetName = getDatasetNameFromPath(S.path);

    // add datasets to query
    query.R = R;
    query.S = S;

    // add query to global config
    g_config.queryData = query;

    // all is well
    return true;
}


bool verifyPipelineSettingsAndBuild(pipelineStatementT &pipelineStmt) {
    // verify setting
    auto it = pipelineSettingsStringToIntMap.find(pipelineStmt.settingStr);
    if (it == pipelineSettingsStringToIntMap.end()) {
        log_err_w_text("Invalid pipeline settings", pipelineStmt.settingStr);
        // todo: print availables
        return false;
    }
    // add setting
    g_config.pipeline.setting = it->second;
    // based on setting, add the corresponding pipeline parts
    switch(g_config.pipeline.setting) {
        case spatial_lib::P_ST2:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_STANDARD;
            g_config.pipeline.iFilterType = spatial_lib::IF_NONE;
            g_config.pipeline.RefinementEnabled = true;
            break;
        case spatial_lib::P_ST3:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_STANDARD;
            g_config.pipeline.iFilterType = spatial_lib::IF_APRIL_STANDARD;
            g_config.pipeline.RefinementEnabled = true;
            break;
        case spatial_lib::P_OP2:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_OPTIMIZED;
            g_config.pipeline.iFilterType = spatial_lib::IF_NONE;
            g_config.pipeline.RefinementEnabled = true;
            break;
        case spatial_lib::P_OP3:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_OPTIMIZED;
            g_config.pipeline.iFilterType = spatial_lib::IF_APRIL_OPTIMIZED;
            g_config.pipeline.RefinementEnabled = true;
            break;
        case spatial_lib::P_SCALABILITY_OP2:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_SCALABILITY;
            g_config.pipeline.iFilterType = spatial_lib::IF_NONE;
            g_config.pipeline.RefinementEnabled = true;
            break;
        case spatial_lib::P_SCALABILITY_OP3:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_SCALABILITY;
            g_config.pipeline.iFilterType = spatial_lib::IF_APRIL_SCALABILITY;
            g_config.pipeline.RefinementEnabled = true;
            break;
        case spatial_lib::P_OTF:
            g_config.pipeline.MBRFilterType = spatial_lib::MF_OPTIMIZED;
            g_config.pipeline.iFilterType = spatial_lib::IF_APRIL_OTF;
            g_config.pipeline.RefinementEnabled = true;
            break;
    }

    return true;
}


/**
 * creates APRIL for a dataset using the rasterizer lib
 * returns true if completed successsfully
*/
static bool createAPRIL(spatial_lib::DatasetT &dataset) {
    uint polygonCount;
    int vertexCount;
    double x,y;
    uint recID;
    // check dataset file path
    ifstream datasetFile(dataset.path, fstream::in | ios_base::binary);
    if(!datasetFile) {
        log_err("Cannot open dataset file.");
        return false;
    }
    // check approximation file paths
    ofstream foutALL(dataset.aprilConfig.ALL_intervals_path, fstream::out | ios_base::binary);
    if(!foutALL) {
        log_err("Cannot open output approximation file.");
        return false;
    }
    ofstream foutFULL(dataset.aprilConfig.FULL_intervals_path, fstream::out | ios_base::binary);
    if(!foutFULL) {
        log_err("Cannot open output approximation file.");
        return false;
    }
    
    //first read the total polygon count
    datasetFile.read((char*) &polygonCount, sizeof(int));
    // write it to the out files
    foutALL.write((char*)&polygonCount, sizeof(uint));
    foutFULL.write((char*)&polygonCount, sizeof(uint));
    //read polygons
    for(int j=0; j<polygonCount; j++){
        //read/write the polygon id
        datasetFile.read((char*) &recID, sizeof(uint)); 
        //read the vertex count
        datasetFile.read((char*) &vertexCount, sizeof(int));
        std::vector<spatial_lib::PointT> vertices;
        vertices.reserve(vertexCount);
        // for MBR
        double xMin = std::numeric_limits<int>::max();
		double yMin = std::numeric_limits<int>::max();
		double xMax = -std::numeric_limits<int>::max();
		double yMax = -std::numeric_limits<int>::max();
        // read the points
        for(int i=0; i<vertexCount; i++){
            datasetFile.read((char*) &x, sizeof(double));
            datasetFile.read((char*) &y, sizeof(double));
            spatial_lib::PointT p;
            p.x = x;
            p.y = y;
            vertices.emplace_back(p);
            // store MBR
            xMin = min(xMin, x);
            yMin = min(yMin, y);
            xMax = max(xMax, x);
            yMax = max(yMax, y);
        }
        // create a polygon
        rasterizerlib::polygon2d polygon = rasterizerlib::createPolygon(vertices);

        // get sections for the polygon based on how many have been set (stored in dataset object) 
        std::vector<spatial_lib::SectionT*> sections = spatial_lib::getSectionsOfMBR(dataset, xMin, yMin, xMax, yMax);

        // for each section
        for (auto &section : sections) {      
            // printf("rec %u at section %u\n", recID, section->sectionID);
            // printf("Setion interest area: (%f,%f),(%f,%f) \n", section->interestxMin, section->interestyMin, section->interestxMax, section->interestyMax);  
            // printf("Setion raster area: (%f,%f),(%f,%f) \n", section->rasterxMin, section->rasteryMin, section->rasterxMax, section->rasteryMax);   
            // set rasterizer lib dataspace to the section's data space
            rasterizerlib::setDataspace(section->rasterxMin, section->rasteryMin, section->rasterxMax, section->rasteryMax);

            // generate the raster based on the configuration
            spatial_lib::AprilDataT aprilData = rasterizerlib::generate(polygon, rasterizerlib::GT_APRIL);
            
            if (dataset.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                // save on disk
                APRIL::saveAPRILonDisk(foutALL, foutFULL, recID, section->sectionID, &aprilData);
            } else{
                // compress and save
                APRIL::compressAndSaveAPRILonDisk(foutALL, foutFULL, recID, section->sectionID, &aprilData);
            }

        }
    }

    datasetFile.close();
    foutALL.close();
    foutFULL.close();

    return true;
}


static void initializeSections(spatial_lib::QueryT &queryData, spatial_lib::DatasetT &dataset) {
    // partitions/sections per dimension
    uint partitionsNum = dataset.aprilConfig.partitions;
    for (int i=0; i<partitionsNum; i++) {
        for (int j=0; j<partitionsNum; j++) {
            spatial_lib::SectionT section;
            section.i = i;
            section.j = j;
            section.sectionID = spatial_lib::getSectionIDFromIdxs(i, j, partitionsNum);
            section.interestxMin = queryData.dataspaceInfo.xMinGlobal + ((i) * (queryData.dataspaceInfo.xExtent / (double) partitionsNum));
            section.interestxMax = queryData.dataspaceInfo.xMinGlobal + ((i+1) * (queryData.dataspaceInfo.xExtent / (double) partitionsNum));
            section.interestyMin = queryData.dataspaceInfo.yMinGlobal + ((j) * (queryData.dataspaceInfo.yExtent / (double) partitionsNum));
            section.interestyMax = queryData.dataspaceInfo.yMinGlobal + ((j+1) * (queryData.dataspaceInfo.yExtent / (double) partitionsNum));
            section.rasterxMin = section.interestxMin;
            section.rasteryMin = section.interestyMin;
            section.rasterxMax = section.interestxMax;
            section.rasteryMax = section.interestyMax;
            section.objectCount = 0;
            // store to dataset
            dataset.sectionMap.insert(std::make_pair(section.sectionID, section));
            // printf("Initialized section interest area: (%f,%f),(%f,%f) \n", section.interestxMin, section.interestyMin, section.interestxMax, section.interestyMax);  
            
        }
    }
}

static void adjustSections(spatial_lib::DatasetT &dataset, uint recID, double xMin, double yMin, double xMax, double yMax) {
    std::vector<spatial_lib::SectionT*> sections = spatial_lib::getSectionsOfMBR(dataset, xMin, yMin,xMax, yMax);
    for (auto &section : sections) {
        if(xMin < section->interestxMin){
            section->rasterxMin = std::min(section->rasterxMin, xMin);
        }
        if(yMin < section->interestyMin){
            section->rasteryMin = std::min(section->rasteryMin, yMin);
        }
        if(xMax > section->interestxMax){
            section->rasterxMax = std::max(section->rasterxMax, xMax);
        }
        if(yMax > section->interestyMax){
            section->rasteryMax = std::max(section->rasteryMax, yMax);
        }

        // printf("Adjusted section raster area from: (%f,%f),(%f,%f) \n", section->interestxMin, section->interestyMin, section->interestxMax, section->interestyMax);  
        // printf("To: (%f,%f),(%f,%f) \n", section->rasterxMin, section->rasteryMin, section->rasterxMax, section->rasteryMax);  

        // if on the fly
        if (g_config.pipeline.iFilterType == spatial_lib::IF_APRIL_OTF) {
            // store the rec to section map
            spatial_lib::addObjectToSectionMap(dataset, section->sectionID, recID);
        }

    }
    
}

/**
 * @brief sets the raster areas of the sections of both datasets to match
 * 
 * @param datasetR 
 * @param datasetS 
 */
static void matchSections(spatial_lib::DatasetT &datasetR, spatial_lib::DatasetT &datasetS) {    
    for(int secID = 0; secID < datasetR.sectionMap.size(); secID++) {
        auto itR = datasetR.sectionMap.find(secID);
        auto itS = datasetS.sectionMap.find(secID);

        double xMin = std::min(itR->second.rasterxMin,itS->second.rasterxMin);
        itR->second.rasterxMin = xMin;
        itS->second.rasterxMin = xMin;
        double yMin = std::min(itR->second.rasteryMin,itS->second.rasteryMin);
        itR->second.rasteryMin = yMin;
        itS->second.rasteryMin = yMin;
        double xMax = std::max(itR->second.rasterxMax,itS->second.rasterxMax);
        itR->second.rasterxMax = xMax;
        itS->second.rasterxMax = xMax;
        double yMax = std::max(itR->second.rasteryMax,itS->second.rasteryMax);
        itR->second.rasteryMax = yMax;
        itS->second.rasteryMax = yMax;

        // printf("Finalized Section %d: (%f,%f),(%f,%f)\n", secID, itR->second.rasterxMin, itR->second.rasteryMin, itR->second.rasterxMax, itR->second.rasteryMax);
    }
}

/**
 * @brief works only for custom testing EU. Hardcoded file paths etc.
 * 
 */
static void initScalabilityTesting() {
    spatial_lib::setupScalabilityTesting();
}

static void initAPRIL() {

    // create if requested
    if (g_config.actions.createApproximations) {
        // init rasterizer lib
        rasterizerlib::init(g_config.queryData.dataspaceInfo.xMinGlobal, g_config.queryData.dataspaceInfo.yMinGlobal, g_config.queryData.dataspaceInfo.xMaxGlobal, g_config.queryData.dataspaceInfo.yMaxGlobal);
        // printf("Initialized rasterizer to (%f,%f),(%f,%f)\n",g_config.queryData.dataspaceInfo.xMinGlobal, g_config.queryData.dataspaceInfo.yMinGlobal, g_config.queryData.dataspaceInfo.xMaxGlobal, g_config.queryData.dataspaceInfo.yMaxGlobal);
        // create APRIL for R
        clock_t timer = clock();
        log_task("Creating APRIL for R...");
        bool ret = createAPRIL(g_config.queryData.R);
        if (!ret) {
            log_err("Create APRIL failed for dataset R.");
            exit(-1);
        }
        success_text_with_time("Finished in ", spatial_lib::time::getElapsedTime(timer));

        // create APRIL for S
        timer = clock();
        log_task("Creating APRIL for S...");
        ret = createAPRIL(g_config.queryData.S);
        if (!ret) {
            log_err("Create APRIL failed for dataset S.");
            exit(-1);
        }

        success_text_with_time("Finished in ", spatial_lib::time::getElapsedTime(timer));

    }

    // verify paths
    if (!verifyFileExists(g_config.queryData.R.aprilConfig.ALL_intervals_path)) {
        log_err_w_text("Error verifying APRIL filepath",g_config.queryData.R.aprilConfig.ALL_intervals_path);
        exit(-1);
    }
    if (!verifyFileExists(g_config.queryData.R.aprilConfig.FULL_intervals_path)) {
        log_err_w_text("Error verifying APRIL filepath",g_config.queryData.R.aprilConfig.FULL_intervals_path);
        exit(-1);
    }
    if (!verifyFileExists(g_config.queryData.S.aprilConfig.ALL_intervals_path)) {
        log_err_w_text("Error verifying APRIL filepath",g_config.queryData.S.aprilConfig.ALL_intervals_path);
        exit(-1);
    }
    if (!verifyFileExists(g_config.queryData.S.aprilConfig.FULL_intervals_path)) {
        log_err_w_text("Error verifying APRIL filepath",g_config.queryData.S.aprilConfig.FULL_intervals_path);
        exit(-1);
    }

    // load APRIL from disk
    clock_t timer = clock();
    if (g_config.queryData.R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
        // only check R april data because they have to match. Todo: join compressed with uncompressed
        APRIL::loadAPRILfromDisk(g_config.queryData.R);
        APRIL::loadAPRILfromDisk(g_config.queryData.S);
    } else {
        // load compressed APRIL
        APRIL::loadCompressedAPRILfromDisk(g_config.queryData.R);
        APRIL::loadCompressedAPRILfromDisk(g_config.queryData.S);
    }
    // success_text_with_time("Loaded APRIL", spatial_lib::time::getElapsedTime(timer));

    // setup filter (same april config for both datasets)
    // todo: extend to allow different april configs for R and S
    APRIL::setupAPRILIntermediateFilter(&g_config.queryData);

    // success_text("APRIL intermediate filter set.");
}

static void initAPRILOTF() {

    // setup filter (same april config for both datasets)
    // todo: extend to allow different april configs for R and S
    APRIL::setupAPRILIntermediateFilter(&g_config.queryData);

    success_text("APRIL ON THE FLY intermediate filter set.");
}

static void loadDataset(std::string &filepath, bool left) {
    int polygonCount;
    int vertexCount;
    int recID;
    two_layer::Coord x,y;

    ifstream file( filepath, fstream::in | ios_base::binary);
    if(!file)
    {
        cerr << "Cannot open the File: " << filepath << endl;
        exit(1);
    }

    //first read the total polygon count
    file.read((char*) &polygonCount, sizeof(int));
    //read polygons
    for(int j=0; j<polygonCount; j++){

        two_layer::Coord minXmbr, minYmbr, maxXmbr, maxYmbr;
        minXmbr = std::numeric_limits<two_layer::Coord>::max();
        maxXmbr = -std::numeric_limits<two_layer::Coord>::max();
        minYmbr = std::numeric_limits<two_layer::Coord>::max();
        maxYmbr = -std::numeric_limits<two_layer::Coord>::max();

        //read/write the polygon id
        file.read((char*) &recID, sizeof(int)); 

        //read the vertex count
        file.read((char*) &vertexCount, sizeof(int));

        for(int i=0; i<vertexCount; i++){
            file.read((char*) &x, sizeof(double));
            file.read((char*) &y, sizeof(double));

            minXmbr = std::min(minXmbr, x);
            maxXmbr = std::max(maxXmbr, x);
            minYmbr = std::min(minYmbr, y);
            maxYmbr = std::max(maxYmbr, y);

            // if (recID == 206434) {
            //     printf("(%f,%f),",x,y);
            // }
            // if (recID == 206434 || recID == 9173759 || recID == 220594 || recID == 357429) {
            //     printf("(%f,%f),",x,y);
            // }
        }
        
        // if (recID == 206434) {
        //     printf("\n^ polygon %u\n\n", recID);    
        // }
        // if (recID == 206434 || recID == 9173759 || recID == 220594 || recID == 357429) {
        //     printf("\n^ polygon %u\n\n", recID);    
        // }

        // if (recID == 1252559 || recID == 313074) {
        //     printf("Pol %u has vertices: %u\n\n", recID, vertexCount);    
        // }
        
        // MBR
        // if (recID == 95697 || recID == 113621) {
        //     printf("(%f,%f),(%f,%f),(%f,%f),(%f,%f)",minXmbr,minYmbr,maxXmbr,minYmbr,maxXmbr,maxYmbr,minXmbr,maxYmbr);
        //     printf("\n^ polygon %u\n\n", recID); 
        // }
        // add to relation for the MBR filter
        two_layer::addObjectToDataset(left, recID, minXmbr, minYmbr, maxXmbr, maxYmbr);

        // if APRIL is set
        if(g_config.pipeline.iFilterType >= spatial_lib::IF_MARK_APRIL_BEGIN && g_config.pipeline.iFilterType < spatial_lib::IF_MARK_APRIL_END) {
            // adjust sections
            if (left) {
                adjustSections(g_config.queryData.R, recID, minXmbr, minYmbr, maxXmbr, maxYmbr);
            } else {
                adjustSections(g_config.queryData.S, recID, minXmbr, minYmbr, maxXmbr, maxYmbr);
            }
        }
    }
    // success_text_with_text_and_number("Loaded dataset", filepath, polygonCount);
    file.close();
}

// TODO: optimize the section stuff because its too slow
static void setupDataspace(std::string &Rfilepath, std::string &Sfilepath) {
    // if APRIL is set
    if(g_config.pipeline.iFilterType >= spatial_lib::IF_MARK_APRIL_BEGIN && g_config.pipeline.iFilterType < spatial_lib::IF_MARK_APRIL_END) {
        // create sections
        initializeSections(g_config.queryData, g_config.queryData.R);
        initializeSections(g_config.queryData, g_config.queryData.S);
    }
    // load datasets
    loadDataset(Rfilepath, true);
    loadDataset(Sfilepath, false);

    // if APRIL is set
    if(g_config.pipeline.iFilterType >= spatial_lib::IF_MARK_APRIL_BEGIN && g_config.pipeline.iFilterType < spatial_lib::IF_MARK_APRIL_END) {
        // match sections built by R and S
        matchSections(g_config.queryData.R, g_config.queryData.S);
    }
    
    // success_text("Setup data space.");
}

void initConfig() {
    // log_task("Initializing..."); 
    // set up the dataspace
    setupDataspace(g_config.queryData.R.path, g_config.queryData.S.path);
    
    // initialize MBR filter
    switch(g_config.pipeline.MBRFilterType) {
        // check bounds
        if (!g_config.queryData.boundsSet) {
            two_layer::getDatasetGlobalBounds(g_config.queryData.dataspaceInfo.xMinGlobal, g_config.queryData.dataspaceInfo.yMinGlobal, g_config.queryData.dataspaceInfo.xMaxGlobal, g_config.queryData.dataspaceInfo.yMaxGlobal);
            g_config.queryData.boundsSet = true;
            success_text("Calculated global bounds based on dataset bounds.");
        }
        case spatial_lib::MF_STANDARD:
            two_layer::initTwoLayer(1000, spatial_lib::MF_STANDARD, g_config.queryData.type);
            break;
        case spatial_lib::MF_OPTIMIZED:
            two_layer::initTwoLayer(1000, spatial_lib::MF_OPTIMIZED, g_config.queryData.type);
            break;
        case spatial_lib::MF_SCALABILITY:
            two_layer::initTwoLayer(1000, spatial_lib::MF_SCALABILITY, g_config.queryData.type);
            // init scalability stuff
            initScalabilityTesting();
            break;
        default:
            log_err_w_text("Invalid selection of MBR filter.", std::to_string(g_config.pipeline.MBRFilterType));
            exit(-1);
            break;
    }

    // initialize intermediate filter, if any
    switch (g_config.pipeline.iFilterType) {
        case spatial_lib::IF_APRIL_OPTIMIZED:
            // init APRIL
            initAPRIL();
            // set as next stage after two-layer MBR filter
            two_layer::setNextStage(spatial_lib::IF_APRIL_OPTIMIZED);
            break;
        case spatial_lib::IF_APRIL_STANDARD:
            // init APRIL
            initAPRIL();
            // set as next stage after two-layer MBR filter
            two_layer::setNextStage(spatial_lib::IF_APRIL_STANDARD);
            break;
        case spatial_lib::IF_APRIL_OTF:
            // init ON ThE FLY APRIL
            initAPRILOTF();
            // set as next stage after two-layer MBR filter
            two_layer::setNextStage(spatial_lib::IF_APRIL_OTF);
            break;
        case spatial_lib::IF_APRIL_SCALABILITY:
            // init APRIL
            initAPRIL();
            // set as next stage after two-layer MBR filter
            two_layer::setNextStage(spatial_lib::IF_APRIL_SCALABILITY);
            break;
        case spatial_lib::IF_NONE:
            // set refinement as next stage
            two_layer::setNextStage(spatial_lib::IF_NONE);
            break;
        default:
            log_err("Unrecognized intermediate filter.");
            return;
    }

    // initialize refinement, if enabled
    if (g_config.pipeline.RefinementEnabled) {
        // setup refinement
        spatial_lib::setupRefinement(g_config.queryData);
    }
    // reset query output
    spatial_lib::resetQueryOutput();

    // success_text("Init Done!");
}