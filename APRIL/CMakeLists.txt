cmake_minimum_required(VERSION 3.5.0)
project(APRIL VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(ProjectType SHARED) # make STATIC or SHARED

include_directories(libvbyte-master)
if(NOT TARGET libvbyte-master)
    add_subdirectory(libvbyte-master)
endif()

include_directories(SpatialLib)
if(NOT TARGET SpatialLib)
    add_subdirectory(SpatialLib)
endif()

set(SrcFiles
  src/controller.cpp
  src/storage.cpp
  src/join.cpp
  src/utils.cpp
  src/topology.cpp
)

# find boost
find_package(Boost 1.50 REQUIRED COMPONENTS)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

include_directories(include)
add_library(${PROJECT_NAME} ${ProjectType} ${SrcFiles})