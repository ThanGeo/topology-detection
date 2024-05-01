#include "include/def.h"
#include "include/config.h"
#include "include/parsing.h"

static void exportCSV() {
    // check if file exists
    ofstream foutCSV;
    if (verifyFileExists(g_config.actions.csvFilepath)) {
        // exists, so append contents
        foutCSV.open(g_config.actions.csvFilepath, fstream::out | ios_base::binary | ios_base::app);
    } else {
        // doesnt exist, create new file
        foutCSV.open(g_config.actions.csvFilepath, fstream::out | ios_base::binary);
    }
    
    // get settings
    std::string settingsStr = getConfigSettingsStr();
    foutCSV << settingsStr;
    std::string header,contents = "";

    switch (g_config.queryData.type) {
        case spatial_lib::Q_INTERSECT:
        case spatial_lib::Q_INSIDE:
        case spatial_lib::Q_EQUAL:
        case spatial_lib::Q_MEET:
        case spatial_lib::Q_CONTAINS:
        case spatial_lib::Q_COVERED_BY:
        case spatial_lib::Q_COVERS:
        case spatial_lib::Q_DISJOINT:
            header = "Accepted%,Rejected%,Inconclusive%,MBR Filter time,Intermediate Filter time,Refinement time\n";
            contents += std::to_string(spatial_lib::g_queryOutput.trueHits / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100) + "%,";
            contents += std::to_string(spatial_lib::g_queryOutput.trueNegatives / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100) + "%,";
            contents += std::to_string(spatial_lib::g_queryOutput.refinementCandidates / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100) + "%,";
            contents += std::to_string(spatial_lib::g_queryOutput.totalTime - (spatial_lib::g_queryOutput.iFilterTime + spatial_lib::g_queryOutput.refinementTime)) + ",";
            contents += std::to_string(spatial_lib::g_queryOutput.iFilterTime) + ",";
            contents += std::to_string(spatial_lib::g_queryOutput.refinementTime);
            break;
        case spatial_lib::Q_FIND_RELATION:
            header = "Inconclusive%,MBR Filter time,Intermediate Filter time,Refinement time\n";
            contents += std::to_string(spatial_lib::g_queryOutput.refinementCandidates / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100) + "%,";
            contents += std::to_string(spatial_lib::g_queryOutput.totalTime - (spatial_lib::g_queryOutput.iFilterTime + spatial_lib::g_queryOutput.refinementTime)) + ",";
            contents += std::to_string(spatial_lib::g_queryOutput.iFilterTime) + ",";
            contents += std::to_string(spatial_lib::g_queryOutput.refinementTime);
            break;
    }

    // write
    foutCSV << header;
    foutCSV << contents;
    foutCSV << std::endl << std::endl;

    // close
    foutCSV.close();
}

static void printPolStats() {
    // true hit pairs
    unsigned long long vertexSumTH = 0;
    unsigned long long vertexSumRF = 0;
    unsigned long long averageAppearanceOfObjectInTH = 0;
    unsigned long long averageAppearanceOfObjectInRF = 0;
    for (auto &kv : spatial_lib::g_queryOutput.polygonStats.polVerticesMapR) {
        
        auto f = spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.find(kv.first);
        if (f != spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.end()) {
            averageAppearanceOfObjectInTH += f->second;
            
            vertexSumTH += kv.second;
            // printf("%u,%u,%u\n",f->first, kv.second, f->second);
            // vertexSumTH += f->second * kv.second;
        }
        
        f = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.find(kv.first);
        if (f != spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.end()) {
            averageAppearanceOfObjectInRF += f->second;

            vertexSumRF += kv.second;
            // printf("%u,%u\n",f->first, f->second);
            // vertexSumRF += f->second * kv.second;
        }
    }
    double averageTHVertexCountR = (double) vertexSumTH / spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.size();
    double averageRFVertexCountR = (double) vertexSumRF / spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.size();
    printf("%d unique ids in true hits from dataset %s.\n", spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.size(), g_config.queryData.R.nickname.c_str());
    printf("Average vertex count in true hit unique polygons of %s: %0.1f vertices\n", averageTHVertexCountR, g_config.queryData.R.nickname.c_str());
    printf("Average number of appearances in TH per unique polygon of %s: %0.1f times\n", (double) averageAppearanceOfObjectInTH  / spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.size(), g_config.queryData.R.nickname.c_str());
    printf("%d unique ids in refinement from dataset %s.\n", spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.size(), g_config.queryData.R.nickname.c_str());
    printf("Average vertex count in refinement unique polygons of %s: %0.1f veftices\n", averageRFVertexCountR, g_config.queryData.R.nickname.c_str());
    printf("Average number of appearances in RF per unique polygon of %s: %0.1f times\n", (double) averageAppearanceOfObjectInRF  / spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.size(), g_config.queryData.R.nickname.c_str());

    vertexSumTH = 0;
    vertexSumRF = 0;
    averageAppearanceOfObjectInTH = 0;
    averageAppearanceOfObjectInRF = 0;
    for (auto &kv : spatial_lib::g_queryOutput.polygonStats.polVerticesMapS) {
        auto f = spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.find(kv.first);
        if (f != spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.end()) {
            averageAppearanceOfObjectInTH += f->second;


            vertexSumTH += kv.second;
            // vertexSumTH += f->second * kv.second;
        }
        f = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.find(kv.first);
        if (f != spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.end()) {
            averageAppearanceOfObjectInRF += f->second;


            vertexSumRF += kv.second;
            // vertexSumRF += f->second * kv.second;
        }
    }
    double averageTHVertexCountS = (double) vertexSumTH / spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.size();
    double averageRFVertexCountS = (double) vertexSumRF / spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.size();
    printf("%d unique ids in true hits from dataset %s.\n", spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.size(), g_config.queryData.S.nickname.c_str());
    printf("Average vertex count in true hit unique polygons of %s: %0.1f vertices\n", averageTHVertexCountS, g_config.queryData.S.nickname.c_str());
    printf("Average number of appearances in TH per unique polygon of %s: %0.1f times\n", (double) averageAppearanceOfObjectInTH  / spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.size(), g_config.queryData.S.nickname.c_str());
    printf("%d unique ids in refinement from dataset %s.\n", spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.size(), g_config.queryData.S.nickname.c_str());
    printf("Average vertex count in refinement unique polygons of %s: %0.1f veftices\n", averageRFVertexCountS, g_config.queryData.S.nickname.c_str());
    printf("Average number of appearances in RF per unique polygon of %s: %0.1f times\n", (double) averageAppearanceOfObjectInRF  / spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.size(), g_config.queryData.S.nickname.c_str());

}

static void printAnalysis() {
    if (g_config.pipeline.setting != spatial_lib::P_SCALABILITY_OP2 && g_config.pipeline.setting != spatial_lib::P_SCALABILITY_OP3) {
        // relate or find relation
        success_text_with_double_and_unit("Throughput", (double) spatial_lib::g_queryOutput.postMBRFilterCandidates / spatial_lib::g_queryOutput.totalTime, "pairs/sec");
        success_text_with_percentage("Inconclusive pairs (% of MBR candidates)", spatial_lib::g_queryOutput.refinementCandidates / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100);
    }
    else if (g_config.pipeline.setting == spatial_lib::P_SCALABILITY_OP2) {
        for(int i=0; i<spatial_lib::g_scalContainer.numberOfBuckets; i++) {
            success_text_with_number("Complexity level", i+1);
            success_text_with_number("\tRefinement pairs", spatial_lib::g_scalContainer.bucketPolygonCount[i]);
            success_text_with_time("\tRefinement time", spatial_lib::g_scalContainer.bucketRefinementTime[i]);
        }
    }
    else if (g_config.pipeline.setting == spatial_lib::P_SCALABILITY_OP3) {
        for(int i=0; i<spatial_lib::g_scalContainer.numberOfBuckets; i++) {
            success_text_with_number("Complexity level", i+1);
            success_text_with_percentage("\tInconclusive pairs", (spatial_lib::g_scalContainer.bucketInconclusiveCount[i] / (double) spatial_lib::g_scalContainer.bucketPolygonCount[i]) * 100);
            // printf("%u\n", spatial_lib::g_scalContainer.bucketInconclusiveCount[i]);
            success_text_with_time("\tIntermediate filter time", spatial_lib::g_scalContainer.bucketIfilterTime[i]);
            success_text_with_time("\tRefinement time", spatial_lib::g_scalContainer.bucketRefinementTime[i]);
        }
    }

}

static void printResults(int runTimes) {
    spatial_lib::g_queryOutput.queryResults += spatial_lib::g_queryOutput.trueHits;
    printf("Query '%s' on datasets %s and %s:\n",spatial_lib::queryTypeIntToText(g_config.queryData.type).c_str(), g_config.queryData.R.nickname.c_str(),g_config.queryData.S.nickname.c_str());
    printf("Repeats: %d\n", runTimes);
    if (g_config.pipeline.iFilterType >= spatial_lib::IF_MARK_APRIL_BEGIN && g_config.pipeline.iFilterType < spatial_lib::IF_MARK_APRIL_END) {
        printf("Partitions: %d\n", g_config.queryData.R.aprilConfig.partitions);
        if(g_config.pipeline.iFilterType == spatial_lib::IF_APRIL_OTF) {
            printf("Rasterizations done: %u\n", spatial_lib::g_queryOutput.rasterizationsDone);
        }
    }
    printf("Total Time:\t\t %0.4f sec.\n", spatial_lib::g_queryOutput.totalTime);
    printf("- MBR filter:\t\t %0.4f sec.\n", spatial_lib::g_queryOutput.totalTime - (spatial_lib::g_queryOutput.iFilterTime + spatial_lib::g_queryOutput.refinementTime));
    printf("- Intermediate filter:\t %0.4f sec.\n", spatial_lib::g_queryOutput.iFilterTime);
    printf("- Refinement:\t\t %0.4f sec.\n", spatial_lib::g_queryOutput.refinementTime);
    printf("Post MBR filter:\t %d pairs.\n", spatial_lib::g_queryOutput.postMBRFilterCandidates);
    printf("Refinement Candidates:\t %d pairs (%0.2f%).\n", spatial_lib::g_queryOutput.refinementCandidates, spatial_lib::g_queryOutput.refinementCandidates / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100);
    switch (g_config.queryData.type) {
        case spatial_lib::Q_INTERSECT:
        case spatial_lib::Q_INSIDE:
        case spatial_lib::Q_EQUAL:
        case spatial_lib::Q_MEET:
        case spatial_lib::Q_CONTAINS:
        case spatial_lib::Q_COVERED_BY:
        case spatial_lib::Q_COVERS:
        case spatial_lib::Q_DISJOINT:
            printf("Filter Breakdown:\n");
            printf("- True Hits:\t\t %0.2f%\n", spatial_lib::g_queryOutput.trueHits / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100);
            printf("- True Negatives:\t %0.2f%\n", spatial_lib::g_queryOutput.trueNegatives / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100);
            printf("- Inconclusive:\t\t %0.2f%\n", spatial_lib::g_queryOutput.refinementCandidates / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100);
            printf("\n");
            printf("Query Results  :\t %d pairs.\n", spatial_lib::g_queryOutput.queryResults);
            break;
        case spatial_lib::Q_FIND_RELATION:
            printf("Query Results\n");
            printf("- Disjoint  :\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_DISJOINT]);
            printf("- Intersects:\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_INTERSECT]);
            printf("- Meet      :\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_MEET]);
            printf("- Contains  :\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_CONTAINS]);
            printf("- Inside    :\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_INSIDE]);
            printf("- Covers    :\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_COVERS]);
            printf("- Covered by:\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_COVERED_BY]);
            printf("- Equal     :\t\t %d\n", spatial_lib::g_queryOutput.topologyRelationsResultMap[spatial_lib::TR_EQUAL]);

            break;
    }
    if (g_config.pipeline.iFilterType == spatial_lib::IF_APRIL_SCALABILITY || g_config.pipeline.MBRFilterType == spatial_lib::MF_SCALABILITY) {
        for(int i=0; i<spatial_lib::g_scalContainer.numberOfBuckets; i++) {
            printf("Bucket %d:\n",i);
            printf("\tInconclusive pairs: %u\n", spatial_lib::g_scalContainer.bucketInconclusiveCount[i]);
            printf("\tIntermediate filter time: %f\n", spatial_lib::g_scalContainer.bucketIfilterTime[i]);
            printf("\tRefinement time: %f\n", spatial_lib::g_scalContainer.bucketRefinementTime[i]);
        }
    }
}

static void freeMemory() {
    
}

int main(int argc, char *argv[]) {
    // args and configuration file
    parseArgumentsAndConfigurationFile(argc, argv);
    // print config
    // printConfig();

    // init
    initConfig();

    // if experiments option is selected, run 10 times and count average time
    int runTimes = 1;
    if (g_config.actions.runExperiments) {
        runTimes = 10;
    }

    // begin evaluation
    clock_t timer;
    double totalTime = 0;
    double mbrTime = 0;
    double iFilterTime = 0;
    double refinementTime = 0;
    for (int i=0; i<runTimes; i++) {
        log_task_w_text("Running iteration", std::to_string(i+1));
        // reset
        spatial_lib::resetQueryOutput();
        
        // evaluate
        timer = spatial_lib::time::getNewTimer();
        long long totalMBRFilterResults = two_layer::evaluateTwoLayer();
        totalTime += spatial_lib::time::getElapsedTime(timer);
        mbrTime += spatial_lib::g_queryOutput.mbrFilterTime;
        iFilterTime += spatial_lib::g_queryOutput.iFilterTime;
        refinementTime += spatial_lib::g_queryOutput.refinementTime;
    }
    
    // set average time
    spatial_lib::g_queryOutput.totalTime = totalTime / (double) runTimes;
    spatial_lib::g_queryOutput.mbrFilterTime = mbrTime / (double) runTimes;
    spatial_lib::g_queryOutput.iFilterTime = iFilterTime / (double) runTimes;
    spatial_lib::g_queryOutput.refinementTime = refinementTime / (double) runTimes;

    // print results
    // printResults(runTimes);

    printAnalysis();

    // dev only
    // printPolStats();

    // output to CSV
    if (g_config.actions.exportCSV) {
        exportCSV();
    }

    // free any memory
    freeMemory();

    // success_text("Finished successfuly!");
    return 0;
}