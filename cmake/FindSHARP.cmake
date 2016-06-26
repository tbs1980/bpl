#.rst:
# FindSHARP
# -----------
#
# Try to find the SHARP library using the search path SHARP_ROOT
#
# Once done this will define
#
# ::
#
#   SHARP_FOUND - System has SHARP
#   SHARP_INCLUDE_DIRS - The SHARP include directory
#   SHARP_LIBRARIES - The libraries needed to use SHARP
#

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_path(SHARP_INCLUDE_DIRS NAMES sharp.h
    HINTS ${SHARP_ROOT} /usr /usr/local
    PATH_SUFFIXES "include"
)

find_library(
    SHARP_LIBRARY_sharp
    NAMES sharp
    HINTS ${SHARP_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

find_library(
    SHARP_LIBRARY_c_utils
    NAMES c_utils
    HINTS ${SHARP_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

find_library(
    SHARP_LIBRARY_fftpack
    NAMES fftpack
    HINTS ${SHARP_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

set(
    SHARP_LIBRARIES
    ${SHARP_LIBRARY_sharp}
    ${SHARP_LIBRARY_c_utils}
    ${SHARP_LIBRARY_fftpack}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    SHARP
    DEFAULT_MSG
    SHARP_INCLUDE_DIRS
    SHARP_LIBRARIES
)
mark_as_advanced(
    SHARP_INCLUDE_DIRS
    SHARP_LIBRARIES
)
