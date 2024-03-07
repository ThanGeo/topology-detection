#include "binary_files.h"

/*
*-------------------------------------------------------
*
*     GEOMETRY BINARY DATA
*       
*
*-------------------------------------------------------
*/

//offset map
unordered_map<uint,unsigned long> loadOffsetMap(int flag){
	unsigned long offset;
	uint lineCounter = 0;

	uint recID;

	ifstream fin(getOffsetMap(flag), fstream::out | ios_base::binary);

	unordered_map<uint,unsigned long> offset_map;

	int totalLines;

	//read total lines
	fin.read((char*) &totalLines, sizeof(int));

	while(lineCounter < totalLines){
		//read rec id
		fin.read((char*) &recID, sizeof(int));
		//read byte offset
		fin.read((char*) &offset, sizeof(unsigned long));

		offset_map.insert(make_pair(recID, offset));		
		lineCounter++;
	}
	

	fin.close();

	return offset_map;
}