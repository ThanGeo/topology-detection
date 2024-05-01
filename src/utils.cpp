#include "utils.h"

void log_err(char* errorText) {
    std::string errTextStr(errorText);
    std::string msg = RED "[ERROR]" NC ": " + errTextStr;
    fprintf(stderr, "%s\n", msg.c_str());
}

void log_err(std::string errorText) {
    std::string msg = RED "[ERROR]" NC ": " + errorText;
    fprintf(stderr, "%s\n", msg.c_str());
}

void log_err_w_text(std::string errorText, std::string txt) {
    std::string msg = RED "[ERROR]" NC ": " + errorText + ": " + txt;
    fprintf(stderr, "%s\n", msg.c_str());
}

void success_text(char* text) {
    std::string textStr(text);
    std::string msg = GREEN "[SUCCESS]" NC ": " + textStr;
    fprintf(stderr, "%s\n", msg.c_str());
}

void success_text_with_number(char* text, int number) {
    std::string textStr(text);
    std::string msg = GREEN "[SUCCESS]" NC ": " + textStr + ": " + std::to_string(number);
    fprintf(stderr, "%s\n", msg.c_str());
}

void success_text_with_text_and_number(char* text, std::string &extra_text, int number) {
    std::string textStr(text);
    std::string msg = GREEN "[SUCCESS]" NC ": " + textStr + " '" + extra_text + "': " + std::to_string(number);
    fprintf(stderr, "%s\n", msg.c_str());
}

void success_text_with_time(char* text, double seconds) {
    std::string textStr(text);
    std::string msg = GREEN "[SUCCESS]" NC ": " + textStr;
    fprintf(stderr, "%s: %0.2f seconds.\n", msg.c_str(), seconds);
}

void success_text_with_double_and_unit(char *text, double val, std::string unit) {
    std::string textStr(text);
    std::string msg = GREEN "[SUCCESS]" NC ": " + textStr;
    fprintf(stderr, "%s: %0.1f %s.\n", msg.c_str(), val, unit.c_str());
}

void success_text_with_percentage(char *text, double perc) {
    std::string textStr(text);
    std::string msg = GREEN "[SUCCESS]" NC ": " + textStr;
    fprintf(stderr, "%s: %0.2f%%.\n", msg.c_str(), perc);
}

void log_task(std::string text) {
    std::string msg = "[  ...  ]: " + text;
    fprintf(stderr, "%s\n", msg.c_str());
}

void log_task_w_text(std::string text, std::string txt) {
    std::string msg = "[  ...  ]: " + text + ": " + txt;
    fprintf(stderr, "%s\n", msg.c_str());
}

std::string getDatasetNameFromPath(std::string &datasetPath) {
    std::stringstream ss(datasetPath);
    std::string token;

    // set delimiter of directories
    std::size_t found = datasetPath.find('/');
    char delimiter;
    if (found!=std::string::npos){
        delimiter = '/';
    } else if (datasetPath.find('\\')) {
        delimiter = '\\';
    } else {
        // no delimiter
        delimiter = '\0';
    }

    if (delimiter != '\0') {
        while(std::getline(ss, token, delimiter)) {
            // do nothing
        }
    }
    // last item is stored in variable token
    std::string datasetName = token.substr(0, token.length() - 4);
    return datasetName;
}

static bool checkDatafileTypeExtension(std::string &fileExtension) {
    if (fileExtension == "dat") {
        // dat (custom binary)
    } else if (fileExtension == "csv") {
        // csv
        log_err("CSV not supported yet.");
        return false;
    } else {
        // error
        log_err("Unsupported data type:");
        printf("\t%s\n", fileExtension.c_str());
        return false;
    }
    return true;
}

bool verifydatasetNickname(std::string &datasetPath){
    std::string fileExtension = datasetPath.substr(datasetPath.length() - 3);
    // check file extension
    if (!checkDatafileTypeExtension(fileExtension)) {
        return false;
    }

    // check if file exists
    if (!verifyFileExists(datasetPath)) {
        log_err_w_text("Error verifying dataset filepath", datasetPath);
        return false;
    }

    return true;
}

bool verifyFileExists(std::string &filePath) {
    if (FILE *file = fopen(filePath.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}

bool verifyDirectoryExists(std::string &directoryPath) {
    //check if APRIL dirs exists
    DIR* dir = opendir(directoryPath.c_str());
    if(dir) {
        /* Directory exists. */
        closedir(dir);
        return true;
    }else if(ENOENT == errno) {
        /* Directory does not exist. */
        // mkdir(directoryPath.c_str(), 0700);
        return false;
    }else{
        /* opendir() failed for some other reason. */
        return false;
    }
}


std::string getAPRILFilePath(std::string &datasetName, std::string &approximationDir, std::string type, int N, int partitions, int compression) {
    if (type != "ALL" && type != "FULL") {
        log_err("Error specifying april interval type");
        return "";
    }
    std::string approximationFilePath = approximationDir + datasetName;

    approximationFilePath += "_" + std::to_string(N) + "_" + std::to_string(partitions);
    if (compression) {
        approximationFilePath += "_compressed";
    }

    approximationFilePath += "_" + type + ".dat";
    return approximationFilePath;
}

void printPolygon(rasterizerlib::polygon2d &polygon, uint recID) {
    printf("Polygon id %u: \n", recID);
    auto it = polygon.vertices.begin();
    printf("(%f,%f)", it->x, it->y);
    it++;
    for (it; it != polygon.vertices.end(); it++) {
        printf(",(%f,%f)", it->x, it->y);
    }
    printf("\n");
    printf("MBR: \n");
    printf("(%f,%f),(%f,%f)\n", polygon.mbr.minPoint.x, polygon.mbr.minPoint.y, polygon.mbr.maxPoint.x, polygon.mbr.maxPoint.y);
}

void printCellsOfIntervals(std::vector<uint> &intervalList, int cellsPerDim){
    uint x,y;
    for(auto it = intervalList.begin(); it != intervalList.end(); it +=2) {
        uint start = *it;
        uint end = *(it + 1);
        for(uint i=start; i<end; i++) {
            
            rasterizerlib::d2xy(cellsPerDim, i, x, y);
            printf("(%u,%u),",x,y);
        }
    }
    printf("\n");
}
