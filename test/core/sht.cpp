#define BOOST_TEST_MODULE shperical harmonic transforms
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sht.hpp>

template<typename real_scalar_type>
void test_sht(){
    using namespace blackpearl::core;
    size_t l_max(6143);
    size_t m_max(l_max);
    size_t n_side(2048);
    size_t num_pixels = 12*n_side*n_side;
    size_t num_fields(6);
    sht<real_scalar_type> sht_test(l_max,m_max,num_pixels,num_fields);
}

BOOST_AUTO_TEST_CASE(shp_data){
    test_sht<float>();
    test_sht<double>();
}

