#define BOOST_TEST_MODULE shperical data
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_data.hpp>

template<typename real_scalar_type>
void test_shp_data(){
    using namespace blackpearl::core;
    std::size_t n_side = 2048;
    std::size_t num_pixels = 12*n_side*n_side;
    std::vector<size_t> spins = {0,0,2,2,2,2};
    shp_data<real_scalar_type> data(num_pixels,spins);
}

BOOST_AUTO_TEST_CASE(shp_data){
    test_shp_data<float>();
    test_shp_data<double>();
}

