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

static void printAnalysis() {
    success_text_with_double_and_unit("Throughput", (double) spatial_lib::g_queryOutput.postMBRFilterCandidates / spatial_lib::g_queryOutput.totalTime, "pairs/sec");
    success_text_with_percentage("Inconclusive pairs (% of MBR candidates)", spatial_lib::g_queryOutput.refinementCandidates / (double) spatial_lib::g_queryOutput.postMBRFilterCandidates * 100);
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
        // log_task_w_text("Running iteration ", std::to_string(i+1));
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

    // output to CSV
    if (g_config.actions.exportCSV) {
        exportCSV();
    }

    // free any memory
    freeMemory();

    // success_text("Finished successfuly!");
    return 0;
}