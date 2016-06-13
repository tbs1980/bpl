#define BOOST_TEST_MODULE power spectrum to spherical harmonics
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <complex>
#include <limits>
#include <random>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <blackpearl/core/shc_ps_utils.hpp>

template<typename real_scalar_type>
void test_create_shc(){
    using namespace blackpearl::core;
    using namespace boost::numeric::ublas;

    std::mt19937 rng(1234);
    std::size_t l_max = 32;
    std::size_t num_fields = 2;
    pow_spec<real_scalar_type> ps(l_max,num_fields);
    pow_spec<real_scalar_type> ps_test_cum(l_max,num_fields);
    for(size_t mtpl_l=0;mtpl_l<=l_max;++mtpl_l){
        ps(mtpl_l) = identity_matrix<real_scalar_type>(num_fields);
        ps_test_cum(mtpl_l) = zero_matrix<real_scalar_type>(num_fields);
    }
    size_t num_realzns = 1000;
    for(size_t rlzn_i = 0; rlzn_i < num_realzns; ++rlzn_i ){
        sph_hrm_coeffs<real_scalar_type> shc =
            create_gauss_sph_hrm_coeffs(ps,rng);
        pow_spec<real_scalar_type> ps_test =
            extract_pow_spec<real_scalar_type>(shc,shc);
        for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
            ps_test_cum(mtpl_l) += ps_test(mtpl_l);
        }
    }
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        ps_test_cum(mtpl_l) /= real_scalar_type(num_realzns);
    }
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        real_scalar_type err_mc(0.1);
        for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
            for(size_t fld_j = 0; fld_j < num_fields; ++fld_j){
                if( fld_i == fld_j) {
                    BOOST_REQUIRE(
                        std::abs(
                            ps_test_cum(mtpl_l,fld_i,fld_j)
                            - real_scalar_type(1)
                        )
                        <= err_mc
                    );
                }
                else {
                    BOOST_REQUIRE(
                        std::abs(
                            ps_test_cum(mtpl_l,fld_i,fld_j)
                            - real_scalar_type(0)
                        )
                        <= err_mc
                    );
                }
            }
        }
    }
}

template<typename real_scalar_type>
void test_extract_pow_spec() {
    using namespace blackpearl::core;
    using namespace boost::numeric::ublas;
    typedef std::complex<real_scalar_type> complex_scalar_type;

    size_t const l_max = 32;
    size_t const m_max = l_max;
    size_t const num_fields = 2;
    pow_spec<real_scalar_type> ps(l_max,num_fields);
    for(size_t mtpl_l=0;mtpl_l<=l_max;++mtpl_l){
        ps(mtpl_l) = identity_matrix<real_scalar_type>(num_fields)*4.;
    }

    sph_hrm_coeffs<real_scalar_type> shc(l_max,m_max,num_fields);
    for(size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m){
        for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
            shc(mtpl_l,mtpl_m,0) = complex_scalar_type(1.,0.);
            shc(mtpl_l,mtpl_m,1) = complex_scalar_type(1.,0.);
        }
    }

    pow_spec<real_scalar_type> ps_test =
        extract_pow_spec<real_scalar_type>(shc,shc);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
        for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
            for(size_t fld_j = 0; fld_j < num_fields; ++fld_j){
                BOOST_REQUIRE(
                    std::abs(
                        ps_test(mtpl_l,fld_i,fld_j)
                        - real_scalar_type(1)
                    )
                    <= eps
                );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(create_shc){
    test_create_shc<float>();
    test_create_shc<double>();
}

BOOST_AUTO_TEST_CASE(extract_pow_spec){
    test_extract_pow_spec<float>();
    test_extract_pow_spec<double>();
}
