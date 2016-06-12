#define BOOST_TEST_MODULE shperical harmonic transforms
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <limits>
#include <complex>
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
    for(size_t i=0;i<num_pixels;++i) {
        for(size_t j=0;j<num_fields;++j) {
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
    for(size_t i=0;i<num_pixels;++i) {
        for(size_t j=0;j<num_fields;++j) {
            rms += (data(i,j) - data_test(i,j))*(data(i,j) - data_test(i,j));
        }
    }
    rms = std::sqrt(rms)/(real_scalar_type)num_pixels;
    BOOST_CHECK(rms < exp_rms);
}

template<typename real_scalar_type>
void test_sht_unit_alms(){
    typedef std::complex<real_scalar_type> complex_scalar_type;
    using namespace blackpearl::core;
    size_t n_side(8);
    size_t l_max = 2*n_side;
    size_t m_max(l_max);
    size_t num_pixels = 12*n_side*n_side;
    std::vector<size_t> spins = {0};//,0};//,2,2,2,2};
    size_t num_fields = spins.size();
    sph_hrm_coeffs<real_scalar_type> shcfs(l_max,m_max,num_fields);
    for(size_t m=0;m<=m_max;++m) {
        for(size_t l=m;l<=l_max;++l) {
            shcfs(l,m,0) = complex_scalar_type (real_scalar_type(1),real_scalar_type(0) );
        }
    }
    shp_data<real_scalar_type> data(num_pixels,spins);
    sht<real_scalar_type> sht_test(l_max,m_max,num_pixels,num_fields);
    sht_test.synthesise(shcfs,data);

    sph_hrm_coeffs<real_scalar_type> shcfs_test(l_max,m_max,num_fields);
    sht_test.analyse(data,shcfs_test);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    for(size_t m=0;m<=m_max;++m) {
        for(size_t l=m;l<=l_max;++l) {
            // NOTE we can conly achieve limited amount of accuracy here
            // Y^-1 is will not give you the value back with great accuracy :(
            BOOST_CHECK(
                std::abs( shcfs_test(l,m,0).real() - real_scalar_type(1.) )
                <= eps*1e15
            );
            BOOST_CHECK(
                std::abs( shcfs_test(l,m,0).imag() - real_scalar_type(0.) )
                <= eps*1e15
            );
        }
    }
}



BOOST_AUTO_TEST_CASE(shp_data){
    test_sht<float>();
    test_sht<double>();
}

BOOST_AUTO_TEST_CASE(shp_data_unit_alms){
    test_sht_unit_alms<float>();
    test_sht_unit_alms<double>();
}
