cmake_minimum_required(VERSION 3.5.0)
project(SpatialLib VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(ProjectType SHARED) # make STATIC or SHARED

set(SrcFiles
  src/def.cpp
  src/refinement.cpp
  src/mapppings.cpp
  src/data.cpp
  src/timer.cpp
)

# supress warnings
add_definitions(-w)

# find boost
find_package(Boost 1.50 REQUIRED COMPONENTS)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(include)
add_library(${PROJECT_NAME} ${ProjectType} ${SrcFiles})