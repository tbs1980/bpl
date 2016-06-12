#define BOOST_TEST_MODULE power spectrum to spherical harmonics
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <complex>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <blackpearl/core/shc_ps_utils.hpp>

template<typename real_scalar_type>
void test_create_shc(){
    using namespace blackpearl::core;
    using namespace boost::numeric::ublas;

    std::size_t l_max = 32;
    std::size_t num_fields = 2;
    pow_spec<real_scalar_type> ps(l_max,num_fields);
    for(size_t mtpl_l=0;mtpl_l<=l_max;++mtpl_l){
        ps(mtpl_l) = identity_matrix<double>(num_fields)*4.;
    }

    sph_hrm_coeffs<real_scalar_type> shc = create_gauss_sph_hrm_coeffs(ps);
}

template<typename real_scalar_type>
void test_extract_pow_spec() {
    using namespace blackpearl::core;
    using namespace boost::numeric::ublas;
    typedef std::complex<real_scalar_type> complex_scalar_type;

    size_t const l_max = 32;
    size_t const m_max = l_max;
    size_t const num_fields = 1;
    pow_spec<real_scalar_type> ps(l_max,num_fields);
    for(size_t mtpl_l=0;mtpl_l<=l_max;++mtpl_l){
        ps(mtpl_l) = identity_matrix<double>(num_fields)*4.;
    }

    sph_hrm_coeffs<real_scalar_type> shc(l_max,m_max,num_fields);
    for(size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m){
        for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
            shc(mtpl_l,mtpl_m,0) = complex_scalar_type(1.,0.);
        }
    }

    pow_spec<real_scalar_type> ps_test =
        extract_pow_spec<real_scalar_type>(shc,shc);

    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
        std::cout<<mtpl_l<<"\t"<<ps_test(mtpl_l,0,0)<<std::endl;
    }
}

// BOOST_AUTO_TEST_CASE(create_shc){
//     test_create_shc<float>();
//     test_create_shc<double>();
// }

BOOST_AUTO_TEST_CASE(extract_pow_spec){
    //test_extract_pow_spec<float>();
    test_extract_pow_spec<double>();
}
