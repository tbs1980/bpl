#define BOOST_TEST_MODULE shperical harmonic transforms
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_data.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <blackpearl/core/sht.hpp>

template<typename real_scalar_type>
void test_sht(){
    using namespace blackpearl::core;
    size_t n_side(512);
    size_t l_max = 2*n_side;
    size_t m_max(l_max);
    size_t num_pixels = 12*n_side*n_side;
    std::vector<size_t> spins = {0};//,0};//,2,2,2,2};
    size_t num_fields = spins.size();
    sph_hrm_coeffs<real_scalar_type> sh(l_max,m_max,num_fields);
    shp_data<real_scalar_type> data(num_pixels,spins);
    std::mt19937 rng(1234);
    std::uniform_real_distribution<real_scalar_type> uni_real_dist;
    for(size_t i=0;i<num_pixels;++i)
    {
        for(size_t j=0;j<num_fields;++j)
        {
            data(i,j) = uni_real_dist(rng);
        }
    }
    sht<real_scalar_type> sht_test(l_max,m_max,num_pixels,num_fields);
    sht_test.analyse(data,sh);
    shp_data<real_scalar_type> data_test(num_pixels,spins);
    sht_test.synthesise(sh,data_test);
    real_scalar_type exp_rms
        = std::numeric_limits<real_scalar_type>::epsilon()*1e12;
    real_scalar_type rms = 0.;
    for(size_t i=0;i<num_pixels;++i)
    {
        for(size_t j=0;j<num_fields;++j)
        {
            rms += (data(i,j) - data_test(i,j))*(data(i,j) - data_test(i,j));
        }
    }
    rms = std::sqrt(rms)/(real_scalar_type)num_pixels;
    BOOST_CHECK(rms < exp_rms);
}

BOOST_AUTO_TEST_CASE(shp_data){
    test_sht<float>();
    test_sht<double>();
}
