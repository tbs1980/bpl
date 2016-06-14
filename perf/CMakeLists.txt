include_directories(../include)

function(add_blackpearl_perf_test TEST_NAME TEST_SOURCE)
    get_filename_component(TEST_TARGET ${TEST_SOURCE} NAME_WE)
    add_executable(perf_${TEST_TARGET} ${TEST_SOURCE})
    target_link_libraries(
        perf_${TEST_TARGET}
        ${Boost_LIBRARIES}
        ${SHARP_LIBRARIES}
    )
endfunction()

add_blackpearl_perf_test("num_sph_hrm_coeffs" num_sph_hrm_coeffs.cpp)
add_blackpearl_perf_test("sht" sht.cpp)
