#define BOOST_TEST_MODULE shperical harmonic coefficients
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>

template<typename real_scalar_type>
void test_shp_hrm_coeffs(){
    using namespace blackpearl::core;
    size_t l_max(4096);
    size_t m_max(l_max);
    size_t num_fields(6);
    sph_hrm_coeffs<real_scalar_type> sh(l_max,m_max,num_fields);
}

BOOST_AUTO_TEST_CASE(shp_hrm_coeffs){
    test_shp_hrm_coeffs<float>();
    test_shp_hrm_coeffs<double>();
}

