#include "storage.h"

namespace APRIL
{
    uint8_t loadSpace8[10000000];
    uint32_t loadSpace32[10000000];

    void saveAPRILonDisk(std::ofstream &foutALL, std::ofstream &foutFULL, uint recID, uint sectionID, spatial_lib::AprilDataT* aprilData) {
        // std::cout << "Polygon " << recID << " ALL intervals: " << polygon.rasterData.data.listA.size() << std::endl; 
        //write ID
        foutALL.write((char*)(&recID), sizeof(uint));
        //write section
        foutALL.write((char*)(&sectionID), sizeof(uint));
        //write number of intervals
        foutALL.write((char*)(&aprilData->numIntervalsALL), sizeof(uint));
        //write the interval data (numintervals * 2 uint values)
        foutALL.write((char*)(&aprilData->intervalsALL.data()[0]), aprilData->numIntervalsALL * 2 * sizeof(uint));

        // same for FULL (if any)
        if(aprilData->numIntervalsFULL > 0){
            foutFULL.write((char*)(&recID), sizeof(uint));
            foutFULL.write((char*)(&sectionID), sizeof(uint));
            foutFULL.write((char*)(&aprilData->numIntervalsFULL), sizeof(uint));
            foutFULL.write((char*)(&aprilData->intervalsFULL.data()[0]), aprilData->numIntervalsFULL * 2 * sizeof(uint));	
        }
    }

    void compressAndSaveAPRILonDisk(std::ofstream &foutALL, std::ofstream &foutFULL, uint recID, uint sectionID, spatial_lib::AprilDataT* aprilData) {
        // std::cout << "Polygon " << recID << " ALL intervals: " << polygon.rasterData.data.listA.size() << std::endl; 
        //write ID
        foutALL.write((char*)(&recID), sizeof(uint));
        //write section
        foutALL.write((char*)(&sectionID), sizeof(uint));
        // write number of intervals
        foutALL.write((char*)(&aprilData->numIntervalsALL), sizeof(uint));
        
        // compress interval data
        std::copy(aprilData->intervalsALL.begin(), aprilData->intervalsALL.end(), loadSpace32);
        uint bytesALL = vbyte_compress_sorted32((uint32_t*)&loadSpace32[0], &loadSpace8[0], 0, aprilData->numIntervalsALL*2);
    
        //write number of bytes of the compressed interval data
        foutALL.write((char*)(&bytesALL), sizeof(uint));
        //write the compressed data
        foutALL.write((char*)(&loadSpace8[0]), bytesALL*sizeof(uint8_t));	

        // printf("Saved Polygon %u, Section %u, intervals: %u, bytes %u\n", recID, sectionID, aprilData->numIntervalsALL, bytesALL);
        // for(int i=0; i<aprilData->numIntervalsALL*2; i++) {
        //     printf("%u\n",loadSpace32[i]);
        // }
        // exit(0);
        // do the same same for FULL intervals (if any)
        if(aprilData->numIntervalsFULL > 0){
            foutFULL.write((char*)(&recID), sizeof(uint));
            foutFULL.write((char*)(&sectionID), sizeof(uint));
            foutFULL.write((char*)(&aprilData->numIntervalsFULL), sizeof(uint));

            // compress interval data
            std::copy(aprilData->intervalsFULL.begin(), aprilData->intervalsFULL.end(), loadSpace32);
            uint bytesFULL = vbyte_compress_sorted32((uint32_t*)&loadSpace32[0], &loadSpace8[0], 0, aprilData->numIntervalsFULL*2);
        
            //write number of bytes of the compressed interval data
            foutFULL.write((char*)(&bytesFULL), sizeof(uint));
            //write the compressed data
            foutFULL.write((char*)(&loadSpace8[0]), bytesFULL * sizeof(uint8_t));	
        }
    }

    void loadAPRILfromDisk(spatial_lib::DatasetT &dataset) {
        long fileSize;
        uint8_t* buffer;
        size_t result;
        int sectionID;
        uint totalPolygons;
        uint recID;
        /* ---------- ALL ---------- */
        FILE* pFile = fopen(dataset.aprilConfig.ALL_intervals_path.c_str(), "rb");
        // obtain file size:
        fseek (pFile , 0 , SEEK_END);
        fileSize = ftell (pFile);
        rewind (pFile);
        // allocate memory to contain the whole file:
        buffer = (uint8_t*) malloc (sizeof(uint8_t)*fileSize);
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
        // copy the file into the buffer:
        result = fread (buffer, 1, fileSize, pFile);
        if (result != fileSize) {fputs ("Reading error",stderr); exit (3);}
        //init read indices
        unsigned long bufferIndex = 0;
        unsigned long loadIndex;

        //read total polygons
        memcpy(&totalPolygons, &buffer[bufferIndex], sizeof(uint));
        bufferIndex += sizeof(uint);

        //read loop
        uint polCounter = 0;
        while(bufferIndex < fileSize){
            loadIndex = 0;
            //id
            memcpy(&recID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            // make approximation object
            spatial_lib::AprilDataT aprilData = spatial_lib::createEmptyAprilDataObject();
            
            //section
            memcpy(&sectionID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);
            
            //total intervals in polygon
            memcpy(&aprilData.numIntervalsALL, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            //copy uncompressed data
            uint totalValues = aprilData.numIntervalsALL * 2;
            memcpy(&loadSpace32[loadIndex], &buffer[bufferIndex], totalValues * sizeof(uint));
            loadIndex += totalValues;
            bufferIndex += totalValues * sizeof(uint);

            //add intervals to april data
            aprilData.intervalsALL.insert(aprilData.intervalsALL.begin(), &loadSpace32[0], loadSpace32 + loadIndex);

            // add to dataset map
            spatial_lib::addAprilDataToApproximationDataMap(dataset, sectionID, recID, &aprilData);
            
            polCounter++;
        }
        // terminate this file
        fclose (pFile);
        
        //free memory
        free (buffer);

        /* ---------- FULL ---------- */
        // obtain file size:
        pFile = fopen(dataset.aprilConfig.FULL_intervals_path.c_str(), "rb");
        fseek (pFile , 0 , SEEK_END);
        fileSize = ftell (pFile);
        rewind (pFile);
        // allocate memory to contain the whole file:
        buffer = (uint8_t*) malloc (sizeof(uint8_t)*fileSize);
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
        // copy the file into the buffer:
        result = fread (buffer, 1, fileSize, pFile);
        if (result != fileSize) {fputs ("Reading error",stderr); exit (3);}
        //init read indices
        bufferIndex = 0;
        loadIndex;

        //read total polygons
        memcpy(&totalPolygons, &buffer[bufferIndex], sizeof(uint));
        bufferIndex += sizeof(uint);

        //read loop
        polCounter = 0;
        while(bufferIndex < fileSize){
            loadIndex = 0;
            //id
            memcpy(&recID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            //section
            memcpy(&sectionID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);
            
            // fetch approximation object
            spatial_lib::AprilDataT* aprilData = spatial_lib::getAprilDataBySectionAndObjectIDs(dataset, sectionID, recID);
            
            //total intervals in polygon
            memcpy(&aprilData->numIntervalsFULL, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            //copy uncompressed data
            uint totalValues = aprilData->numIntervalsFULL * 2;
            memcpy(&loadSpace32[loadIndex], &buffer[bufferIndex], totalValues * sizeof(uint));
            loadIndex += totalValues;
            bufferIndex += totalValues * sizeof(uint);		

            //add intervals to the april data reference
            aprilData->intervalsFULL.insert(aprilData->intervalsFULL.begin(), &loadSpace32[0], loadSpace32 + loadIndex);


            // if (recID == 206434) {
            //     printAPRIL(*aprilData);
            //     printf("^Dataset: %s polygon %u\n", dataset.nickname.c_str(), recID);
            // }

            polCounter++;
        }
        // terminate this file
        fclose (pFile);
        
        //free memory
        free (buffer);
    }

    void loadCompressedAPRILfromDisk(spatial_lib::DatasetT &dataset) {
        long fileSize;
        uint8_t* buffer;
        size_t result;
        int sectionID;
        uint totalPolygons;
        uint recID;
        /* ---------- ALL ---------- */
        FILE* pFile = fopen(dataset.aprilConfig.ALL_intervals_path.c_str(), "rb");
        // obtain file size:
        fseek (pFile , 0 , SEEK_END);
        fileSize = ftell (pFile);
        rewind (pFile);
        // allocate memory to contain the whole file:
        buffer = (uint8_t*) malloc (sizeof(uint8_t)*fileSize);
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
        // copy the file into the buffer:
        result = fread (buffer, 1, fileSize, pFile);
        if (result != fileSize) {fputs ("Reading error",stderr); exit (3);}
        //init read indices
        unsigned long bufferIndex = 0;
        unsigned long loadIndex;

        //read total polygons
        memcpy(&totalPolygons, &buffer[bufferIndex], sizeof(uint));
        bufferIndex += sizeof(uint);

        //read loop
        uint polCounter = 0;
        while(bufferIndex < fileSize){
            loadIndex = 0;
            //id
            memcpy(&recID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            // make approximation object
            spatial_lib::AprilDataT aprilData = spatial_lib::createEmptyAprilDataObject();
            
            //section
            memcpy(&sectionID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);
            
            //total intervals in polygon
            memcpy(&aprilData.numIntervalsALL, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            //total bytes for the compressed data
            
            memcpy(&aprilData.compressedLengthALL, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            // printf("Polygon %u, Section %u, intervals: %u, bytes %u\n", recID, sectionID, aprilData.numIntervalsALL, bytesALL);
            
            //copy compressed data
            memcpy(&loadSpace8[loadIndex], &buffer[bufferIndex], aprilData.compressedLengthALL * sizeof(uint8_t));
            loadIndex += aprilData.compressedLengthALL;
            bufferIndex += aprilData.compressedLengthALL * sizeof(uint8_t);
            
            //add to polygon
            aprilData.intervalsALLcompressed.insert(aprilData.intervalsALLcompressed.end(), loadSpace8, loadSpace8 + loadIndex);
            // aprilData.intervalsALLcompressed.insert(aprilData.intervalsALLcompressed.end(), buffer, buffer + aprilData.compressedLengthALL);
            // bufferIndex += aprilData.compressedLengthALL * sizeof(uint8_t);

            // add to dataset map
            spatial_lib::addAprilDataToApproximationDataMap(dataset, sectionID, recID, &aprilData);
            
            // printf("Loaded polygon %u \n", recID);

            polCounter++;
        }
        // terminate this file
        fclose (pFile);
        
        //free memory
        free (buffer);

        /* ---------- FULL ---------- */
        // obtain file size:
        pFile = fopen(dataset.aprilConfig.FULL_intervals_path.c_str(), "rb");
        fseek (pFile , 0 , SEEK_END);
        fileSize = ftell (pFile);
        rewind (pFile);
        // allocate memory to contain the whole file:
        buffer = (uint8_t*) malloc (sizeof(uint8_t)*fileSize);
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}
        // copy the file into the buffer:
        result = fread (buffer, 1, fileSize, pFile);
        if (result != fileSize) {fputs ("Reading error",stderr); exit (3);}
        //init read indices
        bufferIndex = 0;
        loadIndex;

        //read total polygons
        memcpy(&totalPolygons, &buffer[bufferIndex], sizeof(uint));
        bufferIndex += sizeof(uint);

        //read loop
        polCounter = 0;
        while(bufferIndex < fileSize){
            loadIndex = 0;
            //id
            memcpy(&recID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            //section
            memcpy(&sectionID, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);
            
            // fetch approximation object
            spatial_lib::AprilDataT* aprilData = spatial_lib::getAprilDataBySectionAndObjectIDs(dataset, sectionID, recID);
            
            //total intervals in polygon
            memcpy(&aprilData->numIntervalsFULL, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);


            //total bytes for the compressed data
            memcpy(&aprilData->compressedLengthFULL, &buffer[bufferIndex], sizeof(uint));
            bufferIndex += sizeof(uint);

            //copy compressed data
            memcpy(&loadSpace8[loadIndex], &buffer[bufferIndex], aprilData->compressedLengthFULL * sizeof(uint8_t));
            loadIndex += aprilData->compressedLengthFULL;
            bufferIndex += aprilData->compressedLengthFULL;
            
            //add to polygon
            aprilData->intervalsFULLcompressed.insert(aprilData->intervalsFULLcompressed.end(), &loadSpace8[0], loadSpace8 + loadIndex);
            // aprilData->intervalsFULLcompressed.insert(aprilData->intervalsFULLcompressed.end(), buffer, buffer + aprilData->compressedLengthFULL);
            // bufferIndex += aprilData->compressedLengthFULL * sizeof(uint8_t);

            polCounter++;
        }
        // terminate this file
        fclose (pFile);
        
        //free memory
        free (buffer);
    }
}