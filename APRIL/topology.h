#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "omp.h"


#include "containers.h"
#include "join.h"

using namespace std;

int create_topology_table_uncompressed(Polygon *polA, Polygon *polB);