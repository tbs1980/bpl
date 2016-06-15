#define BOOST_TEST_MODULE shperical precision matrix
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_prec_mat.hpp>

template<typename real_scalar_type>
void test_sph_diag_prec_mat(){
    using namespace blackpearl::core;
    std::size_t n_side = 2048;
    std::size_t num_pixels = 12*n_side*n_side;
    size_t num_fields = 6;
    sph_diag_prec_mat<real_scalar_type> prec_mat(num_pixels,num_fields);
}

BOOST_AUTO_TEST_CASE(sph_diag_prec_mat){
    test_sph_diag_prec_mat<float>();
    test_sph_diag_prec_mat<double>();
}
