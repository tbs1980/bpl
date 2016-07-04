include_directories(../include)
include_directories(${HEALPix_INCLUDE_DIRS})
add_executable(bpl bpl.cpp)
target_link_libraries(bpl
    ${Boost_LIBRARIES}
    ${SHARP_LIBRARIES}
    ${HEALPix_LIBRARY_healpix_cxx}
    ${HEALPix_LIBRARY_cxxsupport}
    ${CFITSIO_LIBRARIES}
)
