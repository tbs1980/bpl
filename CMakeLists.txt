project(bpl CXX)
cmake_minimum_required(VERSION 2.8)

# guard against in-source builds (got this from Eigen)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

#add the customised package searches to the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(NOT cmake_build_type_tolower STREQUAL "debug" AND NOT cmake_build_type_tolower STREQUAL "release")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug and Release.")
endif()

include(set_compiler_flags)

find_package(Boost 1.53.0 REQUIRED COMPONENTS
    unit_test_framework
    program_options
)
find_package(SHARP REQUIRED)
find_package(HEALPix REQUIRED)
# TODO make sure that we are using the same fitsio for the healpix

include_directories(include)
include_directories(
    ${Boost_INCLUDE_DIRS}
    ${SHARP_INCLUDE_DIRS}
)

enable_testing()
add_subdirectory(test)
add_subdirectory(perf)
add_subdirectory(src)

# cmake ../ -DSHARP_ROOT=/home/sbalan/github/libsharp/auto_libs -DHEALPix_ROOT=/home/sbalan/github/healpix-code/src/cxx/optimized_gcc
