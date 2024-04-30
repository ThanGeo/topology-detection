#ifndef MAIN_DEF_H
#define MAIN_DEF_H

#include "APRIL.h"
#include "SpatialLib.h"
#include "Rasterizer.h"
#include <vector>
#include <unordered_map>

typedef struct Section {
    uint sectionID;
    // x,y indexes
    uint x,y;
} SectionT;

typedef struct DirectoryPath {
    std::string configFilePath = "../config.ini";
    std::string datasetsConfigPath = "../datasets.ini";
    std::string approximationPath = "../data/";
} DirectoryPathT;

typedef struct Pipeline{
    spatial_lib::PipelineSettingE setting;
    spatial_lib::MBRFilterTypeE MBRFilterType;
    spatial_lib::IntermediateFilterTypeE iFilterType;
    bool RefinementEnabled;
} PipelineT;










#endif