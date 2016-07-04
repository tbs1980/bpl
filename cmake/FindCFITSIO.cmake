#.rst:
# FindCFITSIO
# -----------
#
# Try to find the CFITSIO library using the search path CFITSIO_ROOT
#
# Once done this will define
#
# ::
#
#   CFITSIO_FOUND - System has CFITSIO
#   CFITSIO_INCLUDE_DIRS - The CFITSIO include directory
#   CFITSIO_LIBRARIES - The libraries needed to use CFITSIO
#

find_path(CFITSIO_INCLUDE_DIRS NAMES fitsio.h
    HINTS ${CFITSIO_ROOT} /usr /usr/local
    PATH_SUFFIXES "include" "include/cfitsio"
)

find_library(
    CFITSIO_LIBRARIES
    NAMES cfitsio
    HINTS ${CFITSIO_ROOT} /usr /usr/local
    PATH_SUFFIXES lib lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    CFITSIO
    DEFAULT_MSG
    CFITSIO_INCLUDE_DIRS
    CFITSIO_LIBRARIES
)
mark_as_advanced(
    CFITSIO_INCLUDE_DIRS
    CFITSIO_LIBRARIES
)
