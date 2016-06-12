#define BOOST_TEST_MODULE power spectrum
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/pow_spec.hpp>

template<typename real_scalar_type>
void test_pow_spec(){
    using namespace blackpearl::core;
    std::size_t l_max = 2048;
    std::size_t num_fields = 6;
    pow_spec<real_scalar_type> ps(l_max,num_fields);
}

BOOST_AUTO_TEST_CASE(shp_data){
    test_pow_spec<float>();
    test_pow_spec<double>();
}
