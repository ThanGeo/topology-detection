cmake_minimum_required(VERSION 3.5.0)
project(libvbyte-master VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(ProjectType STATIC) # make STATIC or SHARED


set(SrcFiles
  vbyte.cc
)

include_directories(.)
add_library(${PROJECT_NAME} ${ProjectType} ${SrcFiles})