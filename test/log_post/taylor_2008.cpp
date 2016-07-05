#define BOOST_TEST_MODULE taylor 2008
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <complex>
#include <limits>
#include <random>
#include <cstddef>
#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <blackpearl/core/shc_ps_utils.hpp>
#include <blackpearl/core/pow_spec.hpp>
#include <blackpearl/core/win_func.hpp>
#include <blackpearl/core/sph_prec_mat.hpp>
#include <blackpearl/core/sph_data.hpp>
#include <blackpearl/core/sht.hpp>
#include <blackpearl/core/sph_data_prec_mat_utils.hpp>
#include <blackpearl/log_post/taylor_2008.hpp>
#include <blackpearl/utils/lin_alg_utils.hpp>

template<typename real_scalar_type>
void test_taylor_2008_init(){
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;

    std::vector<size_t> spins = {0,0,2,2};
    size_t const num_fields = spins.size();
    size_t const l_max = 8;
    size_t const m_max = l_max;
    size_t const num_sides = l_max/2;
    size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_type> cls(num_fields,l_max);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                cls(fld_i,fld_j,mtpl_l) = 0;
                if(fld_i == fld_j ){
                    cls(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_type> alms =
        create_gauss_sph_hrm_coeffs<real_scalar_type>(cls,rng);
    sht<real_scalar_type> sh_trans(num_fields,l_max,m_max,num_pixels);
    sph_data<real_scalar_type> maps(spins,num_pixels);
    sh_trans.synthesise(alms,maps);
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    real_scalar_type const sigma_pixel = 1e-3;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1./(sigma_pixel*sigma_pixel);
        }
    }
    sph_data<real_scalar_type> noise_maps
        = gen_gauss_noise_data(p_mat,spins,rng);
    maps.data() += noise_maps.data();
    win_func<real_scalar_type> w_func(num_fields,l_max);
    taylor_2008<real_scalar_type> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,cls,pos_q);
    lp_t08.log_post(pos_q);
}

template<typename real_scalar_type>
void test_taylor_2008_log_post_case_1(){
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;
    typedef std::complex<real_scalar_type> complex_scalar_type;
    typedef matrix<real_scalar_type> real_matrix_type;

    std::vector<size_t> spins = {0,0,2,2};
    size_t const num_fields = spins.size();
    size_t const l_max = 8;
    size_t const m_max = l_max;
    size_t const num_sides = l_max/2;
    size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_type> cls(num_fields,l_max);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                cls(fld_i,fld_j,mtpl_l) = 0;
                if(fld_i == fld_j ){
                    cls(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
    pow_spec<real_scalar_type> gls(num_fields,l_max);
    for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        real_matrix_type c_l(num_fields,num_fields);
        cls.get_mtpl(mtpl_l,c_l);
        real_matrix_type g_l = compute_matrix_log<real_scalar_type>(c_l);
        gls.set_mtpl(mtpl_l,g_l);
    }
    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_type> alms(num_fields,l_max,m_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            for(std::size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                alms(fld_i,mtpl_l,mtpl_m) = complex_scalar_type(0,0);
            }
        }
    }
    sph_data<real_scalar_type> maps(spins,num_pixels);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t pix_i = 0; pix_i < num_pixels; ++pix_i) {
            maps(fld_i,pix_i) = 0;
        }
    }
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(std::size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1;
        }
    }
    win_func<real_scalar_type> w_func(num_fields,l_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
            w_func(fld_i,mtpl_l) = 1;
        }
    }
    taylor_2008<real_scalar_type> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,gls,pos_q);
    real_scalar_type log_post_val = lp_t08.log_post(pos_q);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    BOOST_REQUIRE(std::abs(log_post_val) <= eps);
}

template<typename real_scalar_type>
void test_taylor_2008_grad_log_post_case_1(){
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;
    typedef std::complex<real_scalar_type> complex_scalar_type;
    typedef matrix<real_scalar_type> real_matrix_type;

    std::vector<size_t> spins = {0,0,2,2};
    size_t const num_fields = spins.size();
    size_t const l_max = 512;
    size_t const m_max = l_max;
    size_t const num_sides = l_max/2;
    size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_type> cls(num_fields,l_max);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                cls(fld_i,fld_j,mtpl_l) = 0;
                if(fld_i == fld_j ){
                    cls(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
    pow_spec<real_scalar_type> gls(num_fields,l_max);
    for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        real_matrix_type c_l(num_fields,num_fields);
        cls.get_mtpl(mtpl_l,c_l);
        real_matrix_type g_l = compute_matrix_log<real_scalar_type>(c_l);
        gls.set_mtpl(mtpl_l,g_l);
    }
    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_type> alms(num_fields,l_max,m_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            for(std::size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                alms(fld_i,mtpl_l,mtpl_m) = complex_scalar_type(0,0);
            }
        }
    }
    sph_data<real_scalar_type> maps(spins,num_pixels);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t pix_i = 0; pix_i < num_pixels; ++pix_i) {
            maps(fld_i,pix_i) = 0;
        }
    }
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(std::size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1;
        }
    }
    win_func<real_scalar_type> w_func(num_fields,l_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
            w_func(fld_i,mtpl_l) = 1;
        }
    }
    taylor_2008<real_scalar_type> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,gls,pos_q);
    vector<real_scalar_type> grad_pos_q = lp_t08.grad_log_post(pos_q);

    pow_spec<real_scalar_type> grad_gls(num_fields,l_max);
    sph_hrm_coeffs<real_scalar_type> grad_alms(num_fields,l_max,m_max);
    convert_to_coeffs(grad_pos_q,grad_alms,grad_gls);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){

                if(fld_i == fld_j ){
                    real_scalar_type fact_l = -0.5*(2.*mtpl_l - 1.);
                    BOOST_REQUIRE(
                        std::abs( grad_gls(fld_i,fld_j,mtpl_l) - fact_l) <= eps
                    );
                }else {
                    BOOST_REQUIRE(
                        std::abs( grad_gls(fld_i,fld_j,mtpl_l) ) <= eps
                    );
                }
            }
        }
    }

    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            for(std::size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                BOOST_REQUIRE(
                    std::abs( grad_alms(fld_i,mtpl_l,mtpl_m).real() ) <= eps
                );
                BOOST_REQUIRE(
                    std::abs( grad_alms(fld_i,mtpl_l,mtpl_m).imag() ) <= eps
                );
            }
        }
    }
}

template<typename real_scalar_type>
void test_taylor_2008_log_post_case_2(){
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;
    typedef matrix<real_scalar_type> real_matrix_type;

    std::vector<size_t> spins = {0,2,2};
    std::size_t const num_fields = spins.size();
    std::size_t const l_max = 512;
    std::size_t const m_max = l_max;
    std::size_t const num_sides = l_max/2;
    std::size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_type> cls(num_fields,l_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(std::size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                cls(fld_i,fld_j,mtpl_l) = 0;
                if(fld_i == fld_j ){
                    cls(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
    pow_spec<real_scalar_type> gls(num_fields,l_max);
    for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        real_matrix_type c_l(num_fields,num_fields);
        cls.get_mtpl(mtpl_l,c_l);
        real_matrix_type g_l = compute_matrix_log<real_scalar_type>(c_l);
        gls.set_mtpl(mtpl_l,g_l);
    }

    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_type> alms =
        create_gauss_sph_hrm_coeffs<real_scalar_type>(cls,rng);
    sht<real_scalar_type> sh_trans(num_fields,l_max,m_max,num_pixels);
    sph_data<real_scalar_type> maps(spins,num_pixels);
    sh_trans.synthesise(alms,maps);
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1;
        }
    }

    win_func<real_scalar_type> w_func(num_fields,l_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
            w_func(fld_i,mtpl_l) = 1;
        }
    }
    taylor_2008<real_scalar_type> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,gls,pos_q);
    real_scalar_type log_post_val = lp_t08.log_post(pos_q);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    pow_spec<real_scalar_type> sigma =  extract_pow_spec(alms);
    real_scalar_type log_prior(0);
    for(size_t mtpl_l =0; mtpl_l <= l_max; ++mtpl_l)
    {
        real_scalar_type fact_l(2*mtpl_l+1);
        real_matrix_type sigma_l(num_fields,num_fields);
        sigma.get_mtpl(mtpl_l,sigma_l);
        real_scalar_type trace_sig_l = compute_trace<real_scalar_type>(sigma_l);
        log_prior += -0.5*fact_l*trace_sig_l;
    }
    BOOST_REQUIRE( std::abs(log_post_val - log_prior) <= eps );
}

template<typename real_scalar_type>
struct grad_log_post_case_2_scale;

template<>
struct grad_log_post_case_2_scale<double>{
    static constexpr double const val_for_gls = 1e16;
    static constexpr double const pd_alms = 1e-2;
};

template<typename real_scalar_type>
void test_taylor_2008_grad_log_post_case_2(){
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;
    typedef matrix<real_scalar_type> real_matrix_type;

    std::vector<size_t> spins = {0};//,0,2,2};
    size_t const num_fields = spins.size();
    size_t const l_max = 8;
    size_t const m_max = l_max;
    size_t const num_sides = 128;//l_max/2;
    size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_type> cls(num_fields,l_max);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                cls(fld_i,fld_j,mtpl_l) = 0;
                if(fld_i == fld_j ){
                    cls(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
    pow_spec<real_scalar_type> gls(num_fields,l_max);
    for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        real_matrix_type c_l(num_fields,num_fields);
        cls.get_mtpl(mtpl_l,c_l);
        real_matrix_type g_l = compute_matrix_log<real_scalar_type>(c_l);
        gls.set_mtpl(mtpl_l,g_l);
    }
    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_type> alms =
        create_gauss_sph_hrm_coeffs<real_scalar_type>(cls,rng);
    sph_data<real_scalar_type> maps(spins,num_pixels);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t pix_i = 0; pix_i < num_pixels; ++pix_i) {
            maps(fld_i,pix_i) = 0;
        }
    }
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(std::size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1;
        }
    }
    win_func<real_scalar_type> w_func(num_fields,l_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
            w_func(fld_i,mtpl_l) = 1;
        }
    }
    taylor_2008<real_scalar_type> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,gls,pos_q);
    vector<real_scalar_type> grad_pos_q = lp_t08.grad_log_post(pos_q);

    pow_spec<real_scalar_type> grad_gls(num_fields,l_max);
    sph_hrm_coeffs<real_scalar_type> grad_alms(num_fields,l_max,m_max);
    convert_to_coeffs(grad_pos_q,grad_alms,grad_gls);

    pow_spec<real_scalar_type> ps_sigma =  extract_pow_spec(alms);
    real_scalar_type eps_gls = std::numeric_limits<real_scalar_type>::epsilon();
    eps_gls *= grad_log_post_case_2_scale<real_scalar_type>::val_for_gls;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                real_scalar_type fact_1 = 0.5*(2.*mtpl_l + 1.);
                real_scalar_type fact_2 = 0.5*(2.*mtpl_l - 1.);
                if(fld_i == fld_j ){
                    real_scalar_type exp_val
                        = fact_1*ps_sigma(fld_i,fld_j,mtpl_l) - fact_2;
                    real_scalar_type comp_val = grad_gls(fld_i,fld_j,mtpl_l);
                    BOOST_REQUIRE( std::abs( exp_val - comp_val ) <= eps_gls );
                }else {
                    real_scalar_type exp_val
                        = fact_1*ps_sigma(fld_i,fld_j,mtpl_l);
                    real_scalar_type comp_val = grad_gls(fld_i,fld_j,mtpl_l);
                    BOOST_REQUIRE( std::abs( exp_val - comp_val ) <= eps_gls );
                }
            }
        }
    }

    real_scalar_type eps_alms
        = std::numeric_limits<real_scalar_type>::epsilon();
    eps_alms *= grad_log_post_case_2_scale<real_scalar_type>::val_for_gls;
    real_scalar_type omega_pix = 4.*M_PI/(real_scalar_type) num_pixels;
    real_scalar_type grad_fact = -(1. + omega_pix)/omega_pix;
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
            real_scalar_type exp_val_re
                = grad_fact*alms(fld_i,mtpl_l,0).real();
            real_scalar_type comp_val_re
                = grad_alms(fld_i,mtpl_l,0).real();
            BOOST_REQUIRE(
                std::abs( (exp_val_re - comp_val_re )/exp_val_re )
                <= grad_log_post_case_2_scale<real_scalar_type>::pd_alms
            );
        }
        for(std::size_t mtpl_m = 1; mtpl_m <= m_max; ++mtpl_m) {
            for(std::size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                real_scalar_type exp_val_re
                    = 2*grad_fact*alms(fld_i,mtpl_l,mtpl_m).real();
                real_scalar_type comp_val_re
                    = grad_alms(fld_i,mtpl_l,mtpl_m).real();
                BOOST_REQUIRE(
                    std::abs((exp_val_re -comp_val_re)/exp_val_re)
                    <= grad_log_post_case_2_scale<real_scalar_type>::pd_alms
                );
                real_scalar_type exp_val_im
                    = 2*grad_fact*alms(fld_i,mtpl_l,mtpl_m).imag();
                real_scalar_type comp_val_im
                    = grad_alms(fld_i,mtpl_l,mtpl_m).imag();
                BOOST_REQUIRE(
                    std::abs((exp_val_im -comp_val_im)/exp_val_im)
                    <= grad_log_post_case_2_scale<real_scalar_type>::pd_alms
                );
            }
        }
    }
}

template<typename real_scalar_type>
struct grad_log_post_case_3_scale;

template<>
struct grad_log_post_case_3_scale<double>{
    static constexpr double const val_num_grads = 1e15;
};

template<class real_scalar_type>
void test_taylor_2008_grad_log_post_case_3(){
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;

    std::vector<std::size_t> spins = {0,0,2,2};
    std::size_t const num_fields = spins.size();
    std::size_t const l_max = 4;
    std::size_t const m_max = l_max;
    std::size_t const num_sides = l_max/2;
    std::size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_type> cls(num_fields,l_max);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                cls(fld_i,fld_j,mtpl_l) = 0;
                if(fld_i == fld_j ){
                    cls(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
    // std::random_device rng_dev;
    // std::size_t rng_seed = 1234;//rng_dev();
    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_type> alms =
        create_gauss_sph_hrm_coeffs<real_scalar_type>(cls,rng);
    pow_spec<real_scalar_type> sigma =  extract_pow_spec(alms);
    sht<real_scalar_type> sh_trans(num_fields,l_max,m_max,num_pixels);
    sph_data<real_scalar_type> maps(spins,num_pixels);
    sh_trans.synthesise(alms,maps);
    sph_diag_prec_mat<real_scalar_type> p_mat(num_fields,num_pixels);
    real_scalar_type const omega_pix = 4.*M_PI/(real_scalar_type) num_pixels;
    real_scalar_type const sigma_pixel = std::sqrt(1e-2/omega_pix);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1./(sigma_pixel*sigma_pixel);
        }
    }
    sph_data<real_scalar_type> noise_maps
        = gen_gauss_noise_data(p_mat,spins,rng);
    maps.data() += noise_maps.data();
    win_func<real_scalar_type> w_func(num_fields,l_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
            w_func(fld_i,mtpl_l) = 1;
        }
    }

    taylor_2008<real_scalar_type> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,sigma,pos_q);
    vector<real_scalar_type> grad_pos_q = lp_t08.grad_log_post(pos_q);

    vector<real_scalar_type> grad_pos_q_num(grad_pos_q.size());
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    vector<real_scalar_type> p_plus_h(pos_q);
    vector<real_scalar_type> p_minus_h(pos_q);
    for(std::size_t ind_i = 0; ind_i<pos_q.size(); ++ind_i){
        if(pos_q(ind_i) != 0){
            real_scalar_type h = pos_q(ind_i)*std::sqrt(eps)*10;
            p_plus_h(ind_i) += h;
            p_minus_h(ind_i) -= h;
            real_scalar_type f_plus_h = lp_t08.log_post(p_plus_h);
            real_scalar_type f_minus_h = lp_t08.log_post(p_minus_h);
            grad_pos_q_num(ind_i) = 0.5*(f_plus_h - f_minus_h)/h;
            p_plus_h(ind_i) = pos_q(ind_i);
            p_minus_h(ind_i) = pos_q(ind_i);
        }
        else {
            grad_pos_q_num(ind_i) = 0.;
        }
    }

    for(std::size_t ind_i = num_cls; ind_i<num_alms+num_cls; ++ind_i) {
        if(pos_q(ind_i) > 0) {
            BOOST_REQUIRE(
                std::abs( grad_pos_q(ind_i) - grad_pos_q_num(ind_i) ) <=
                eps*grad_log_post_case_3_scale<real_scalar_type>::val_num_grads
            );
        }
    }
}

BOOST_AUTO_TEST_CASE(taylor_2008_init){
    // test_taylor_2008_init<float>();
    test_taylor_2008_init<double>();
}

BOOST_AUTO_TEST_CASE(taylor_2008_log_post_case_1){
    // test_taylor_2008_log_post_case_1<float>();
    test_taylor_2008_log_post_case_1<double>();
}

BOOST_AUTO_TEST_CASE(taylor_2008_grad_log_post_case_1){
    // test_taylor_2008_grad_log_post_case_1<float>();
    test_taylor_2008_grad_log_post_case_1<double>();
}

BOOST_AUTO_TEST_CASE(taylor_2008_log_post_case_2){
    // test_taylor_2008_log_post_case_2<float>();
    test_taylor_2008_log_post_case_2<double>();
}

BOOST_AUTO_TEST_CASE(taylor_2008_grad_log_post_case_2){
    // test_taylor_2008_grad_log_post_case_2<float>();
    test_taylor_2008_grad_log_post_case_2<double>();
}

BOOST_AUTO_TEST_CASE(taylor_2008_grad_log_post_case_3){
    test_taylor_2008_grad_log_post_case_3<double>();
}
