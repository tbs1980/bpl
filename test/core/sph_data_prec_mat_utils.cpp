#define BOOST_TEST_MODULE shperical data prec mat utils
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_data.hpp>
#include <blackpearl/core/sph_data_prec_mat_utils.hpp>

template<typename real_scalar_type>
void test_gauss_noise_data(){
    using namespace blackpearl::core;
    size_t const n_side = 512;
    size_t const num_pixels = 12*n_side*n_side;
    std::vector<size_t> spins = {0,0,2,2};
    size_t const  num_fields = spins.size();
    std::mt19937 rng(1234);
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1;
        }
    }
    sph_data<real_scalar_type> data
        = gen_gauss_noise_data(p_mat,spins,rng);
}

BOOST_AUTO_TEST_CASE(gauss_noise_data){
    test_gauss_noise_data<float>();
    test_gauss_noise_data<double>();
}
