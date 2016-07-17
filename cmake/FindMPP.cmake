#.rst:
# FindMPP
# -----------
#
# Try to find the MPP library using the search path MPP_ROOT
#
# Once done this will define
#
# ::
#
#   MPP_FOUND - System has MPP
#   MPP_INCLUDE_DIRS - The MPP include directory
#   MPP_LIBRARIES - The libraries needed to use MPP
#   MPP_ROOT - The search path
#

find_path(MPP_INCLUDE_DIRS NAMES mpp.hpp
    HINTS ${MPP_ROOT} /usr /usr/local
    PATH_SUFFIXES "include"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    MPP
    DEFAULT_MSG
    MPP_INCLUDE_DIRS
)
mark_as_advanced(
    MPP_INCLUDE_DIRS
)
