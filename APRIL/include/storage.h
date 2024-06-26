#ifndef APRIL_STORAGE_H
#define APRIL_STORAGE_H

#include "SpatialLib.h"
#include "utils.h"

#include "../libvbyte-master/vbyte.h"

#include <fstream>
#include <string.h>

namespace APRIL
{

    /**
     * save function for APRIL binary file
    */
    void saveAPRILonDisk(std::ofstream &foutALL, std::ofstream &foutFULL, uint recID, uint sectionID, spatial_lib::AprilDataT* aprilData);
    /**
     * compress and save function for APRIL binary file
    */
    void compressAndSaveAPRILonDisk(std::ofstream &foutALL, std::ofstream &foutFULL, uint recID, uint sectionID, spatial_lib::AprilDataT* aprilData);
    /**
     * load function for APRIL binary file
    */
    void loadAPRILfromDisk(spatial_lib::DatasetT &dataset);
    /**
     * load function for compressed APRIL binary file
    */
    void loadCompressedAPRILfromDisk(spatial_lib::DatasetT &dataset);
}


#endif