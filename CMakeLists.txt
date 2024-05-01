cmake_minimum_required(VERSION 3.0.0)
project(Main VERSION 1.0.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(MSVC)
    add_compile_options(/std::c++latest)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
    set(BUILD_SHARED_LIBS TRUE)
endif()

# supress warnings
add_definitions(-w)

include_directories(libvbyte-master)
if(NOT TARGET libvbyte-master)
    add_subdirectory(libvbyte-master)
endif()

include_directories(SpatialLib)
if(NOT TARGET SpatialLib)
    add_subdirectory(SpatialLib)
endif()

include_directories(Rasterizer)
if(NOT TARGET Rasterizer)
    add_subdirectory(Rasterizer)
endif()

include_directories(APRIL)
if(NOT TARGET APRIL)
    add_subdirectory(APRIL)
endif()

include_directories(TwoLayerFilter)
if(NOT TARGET TwoLayerFilter)
    add_subdirectory(TwoLayerFilter)
endif()

set(SrcFiles
    src/def.cpp
    src/config.cpp
    src/parsing.cpp
    src/utils.cpp
)

# supress warnings
add_definitions(-w)
include_directories(include)
add_library(${PROJECT_NAME} ${ProjectType} ${SrcFiles})

# find boost
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.50 COMPONENTS REQUIRED)

#main
add_executable(main main.cpp)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

# link all
target_link_libraries(main ${PROJECT_NAME})
target_link_libraries(main ${Boost_LIBRARIES})
target_link_libraries(main APRIL)
target_link_libraries(main libvbyte-master)
target_link_libraries(main TwoLayerFilter)
target_link_libraries(main SpatialLib)
target_link_libraries(main Rasterizer)