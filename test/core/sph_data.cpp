#define BOOST_TEST_MODULE shperical data
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <iostream>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_data.hpp>

template<typename size_type>
size_type num_sph_hrm_coeffs(size_type const l_max,size_type const m_max) {
    return (m_max+size_type(1))*(l_max+size_type(2))/size_type(2)+(l_max+size_type(1))*(l_max-m_max);
}

// template<typename real_scalar_type>
// void test_shp_data(){
//     using namespace blackpearl::core;
//     std::size_t num_pixels = 12*4096*4096;
//     std::vector<size_t> spins = {0,0,2,2,2,2};
//     shp_data<real_scalar_type> data(num_pixels,spins);

// }

// BOOST_AUTO_TEST_CASE(shp_data){
//     test_shp_data<float>();
//     test_shp_data<double>();
// }

BOOST_AUTO_TEST_CASE(harm_coeffs)
{
    typedef size_t size_type;
    size_type const l_max(4096);
    size_type const m_max(l_max);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for(size_t i=0;i<100000;++i){
        size_type num_alms = num_sph_hrm_coeffs(l_max,m_max);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout<<"Elapsed time = "<<elapsed_seconds.count()<<std::endl;
}
