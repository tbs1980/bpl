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
    matrix<real_scalar_type> ps_mat_l(num_fields,num_fields);
    ps_mat_l(0,0) = 1;
    ps_mat_l(0,1) = 0;
    ps_mat_l(1,1) = 2;
    for(size_t mtpl_l=0;mtpl_l<=l_max;++mtpl_l){
        ps.set_mtpl(mtpl_l,ps_mat_l);
        ps_test_cum.set_mtpl(mtpl_l,zero_matrix<real_scalar_type>(num_fields));
    }

    size_t num_realzns = 1000;
    for(size_t rlzn_i = 0; rlzn_i < num_realzns; ++rlzn_i ){
        sph_hrm_coeffs<real_scalar_type> shc =
            create_gauss_sph_hrm_coeffs(ps,rng);
        pow_spec<real_scalar_type> ps_test =
            extract_pow_spec<real_scalar_type>(shc,shc);
        for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
            for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
                for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                    ps_test_cum(fld_i,fld_j,mtpl_l) +=
                        ps_test(fld_i,fld_j,mtpl_l);
                }
            }
        }
    }

    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                ps_test_cum(fld_i,fld_j,mtpl_l) /=
                    real_scalar_type(num_realzns);
            }
        }
    }

    real_scalar_type err_mc = 0.1;
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
        BOOST_REQUIRE(
            std::abs( ps_test_cum(0,0,mtpl_l) - real_scalar_type(1) )
            <= err_mc
        );
        BOOST_REQUIRE(
            std::abs( ps_test_cum(0,1,mtpl_l) - real_scalar_type(0) )
            <= err_mc
        );
        BOOST_REQUIRE(
            std::abs( ps_test_cum(1,1,mtpl_l) - real_scalar_type(2) )
            <= err_mc
        );
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

    sph_hrm_coeffs<real_scalar_type> shc(l_max,m_max,num_fields);
    for(size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m){
        for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
            shc(0,mtpl_l,mtpl_m) = complex_scalar_type(1.,0.);
            shc(1,mtpl_l,mtpl_m) = complex_scalar_type(1.,0.);
        }
    }

    pow_spec<real_scalar_type> ps_test =
        extract_pow_spec<real_scalar_type>(shc,shc);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                BOOST_REQUIRE(
                    std::abs(
                        ps_test(fld_i,fld_j,mtpl_l)
                        - real_scalar_type(1)
                    )
                    <= eps
                );
            }
        }
    }
}

template<typename real_scalar_type>
void test_convert(){
    using namespace blackpearl::core;
    using namespace boost::numeric::ublas;
    typedef std::complex<real_scalar_type> complex_scalar_type;

    size_t const l_max = 1024;
    size_t const m_max = l_max;
    size_t const num_fields = 5;

    sph_hrm_coeffs<real_scalar_type> shc(l_max,m_max,num_fields);
    pow_spec<real_scalar_type> pspec(l_max,num_fields);

    size_t num_rac = shc.num_real_indep_coeffs(num_fields,l_max,m_max);
    size_t num_rpc = pspec.num_real_indep_coeffs(num_fields,l_max);
    vector<real_scalar_type> pos_q(num_rac + num_rpc);

    size_t index = 0;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i)  {
        for(size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                real_scalar_type re = (real_scalar_type) index;
                real_scalar_type im = (real_scalar_type) index+1;
                shc(fld_i,mtpl_l,mtpl_m) = complex_scalar_type(re,im);
                index += 2;
            }
        }
    }

    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                pspec(fld_i,fld_j,mtpl_l) = (real_scalar_type) index;
                ++index;
            }
        }
    }

    convert_to_real_vector<real_scalar_type>(shc,pspec,pos_q);
    for(size_t ind_i = 0; ind_i < pos_q.size(); ++ind_i){
        BOOST_REQUIRE(ind_i == (size_t) pos_q(ind_i));
    }

    sph_hrm_coeffs<real_scalar_type> shc_test(l_max,m_max,num_fields);
    pow_spec<real_scalar_type> pspec_test(l_max,num_fields);
    convert_to_coeffs<real_scalar_type>(pos_q,shc_test,pspec_test);

    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i)  {
        for(size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                BOOST_REQUIRE(
                    shc(fld_i,mtpl_l,mtpl_m).real()
                    == shc_test(fld_i,mtpl_l,mtpl_m).real()
                );
                BOOST_REQUIRE(
                    shc(fld_i,mtpl_l,mtpl_m).imag()
                    == shc_test(fld_i,mtpl_l,mtpl_m).imag()
                );
            }
        }
    }

    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                BOOST_REQUIRE(
                    pspec(fld_i,fld_j,mtpl_l)
                    == pspec_test(fld_i,fld_j,mtpl_l)
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

BOOST_AUTO_TEST_CASE(conversion){
    test_convert<float>();
    test_convert<double>();
}
