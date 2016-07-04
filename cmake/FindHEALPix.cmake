#.rst:
# FindHEALPix
# -----------
#
# Try to find the HEALPix library using the search path HEALPix_ROOT
#
# Once done this will define
#
# ::
#
#   HEALPix_FOUND - System has HEALPix
#   HEALPix_INCLUDE_DIRS - The HEALPix include directory
#   HEALPix_LIBRARIES - The libraries needed to use HEALPix
#   HEALPix_LIBRARY_healpix_cxx - HEALPix cxx library
#   HEALPix_LIBRARY_cxxsupport - HEALPix cxxsupport library
#   HEALPix_LIBRARY_sharp - HEALPix sharp library
#   HEALPix_LIBRARY_c_utils - HEALPix c_utils library
#   HEALPix_LIBRARY_fftpack - HEALPix fftpack library
#

find_package(CFITSIO REQUIRED)

find_path(HEALPix_INCLUDE_DIRS NAMES healpix_base.h
    HINTS ${HEALPix_ROOT} /usr /usr/local
    PATH_SUFFIXES "include"
)

find_library(
    HEALPix_LIBRARY_healpix_cxx
    NAMES healpix_cxx
    HINTS ${HEALPix_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

find_library(
    HEALPix_LIBRARY_cxxsupport
    NAMES cxxsupport
    HINTS ${HEALPix_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

find_library(
    HEALPix_LIBRARY_sharp
    NAMES sharp
    HINTS ${HEALPix_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

find_library(
    HEALPix_LIBRARY_c_utils
    NAMES c_utils
    HINTS ${HEALPix_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

find_library(
    HEALPix_LIBRARY_fftpack
    NAMES fftpack
    HINTS ${HEALPix_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

set(
    HEALPix_LIBRARIES
    ${HEALPix_LIBRARY_healpix_cxx}
    ${HEALPix_LIBRARY_cxxsupport}
    ${HEALPix_LIBRARY_sharp}
    ${HEALPix_LIBRARY_c_utils}
    ${HEALPix_LIBRARY_fftpack}
    ${CFITSIO_LIBRARIES}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    HEALPix
    DEFAULT_MSG
    HEALPix_INCLUDE_DIRS
    HEALPix_LIBRARIES
    HEALPix_LIBRARY_healpix_cxx
    HEALPix_LIBRARY_cxxsupport
    HEALPix_LIBRARY_sharp
    HEALPix_LIBRARY_c_utils
    HEALPix_LIBRARY_fftpack
)
mark_as_advanced(
    HEALPix_INCLUDE_DIRS
    HEALPix_LIBRARIES
    HEALPix_LIBRARY_healpix_cxx
    HEALPix_LIBRARY_cxxsupport
    HEALPix_LIBRARY_sharp
    HEALPix_LIBRARY_c_utils
    HEALPix_LIBRARY_fftpack
)
