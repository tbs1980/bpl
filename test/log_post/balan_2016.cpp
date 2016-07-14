#define BOOST_TEST_MODULE balan 2016
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
#include <blackpearl/log_post/balan_2016.hpp>
#include <blackpearl/utils/lin_alg_utils.hpp>

template<typename real_scalar_t>
void test_balan_2016_init() {
    using namespace blackpearl::core;
    using namespace blackpearl::utils;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;
    typedef matrix<real_scalar_t> real_matrix_t;

    std::vector<size_t> spins = {0,0,2,2};
    size_t const num_fields = spins.size();
    size_t const l_max = 8;
    size_t const m_max = l_max;
    size_t const num_sides = l_max/2;
    size_t const num_pixels = 12*num_sides*num_sides;

    pow_spec<real_scalar_t> cls(num_fields,l_max);
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

    pow_spec<real_scalar_t> fls(num_fields,l_max);
    pow_spec<real_scalar_t> cinvls(num_fields,l_max);
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l) {
        real_matrix_t c_l(num_fields, num_fields);
        cls.get_mtpl(mtpl_l,c_l);
        real_matrix_t f_l = compute_matrix_exp(c_l);
        fls.set_mtpl(mtpl_l,f_l);
        real_matrix_t c_inv_l(num_fields,num_fields);
        compute_inverse<real_scalar_t>(c_l,c_inv_l);
        cinvls.set_mtpl(mtpl_l,c_inv_l);
    }

    std::mt19937 rng(1234);
    sph_hrm_coeffs<real_scalar_t> alms =
        create_gauss_sph_hrm_coeffs<real_scalar_t>(cls,rng);
    sph_hrm_coeffs<real_scalar_t> xlms(num_fields,l_max,m_max);
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_l=0; mtpl_l <= l_max; ++mtpl_l) {
            xlms(fld_i,mtpl_l,0)
                = cinvls(fld_i,0,mtpl_l)*alms(0,mtpl_l,0);
        }
        for(std::size_t mtpl_m = 1; mtpl_m <= m_max; ++mtpl_m) {
            for(std::size_t mtpl_l=mtpl_m; mtpl_l<=l_max; ++mtpl_l) {
                xlms(fld_i,mtpl_l,mtpl_m)
                    = real_scalar_t(0.5)*cinvls(fld_i,0,mtpl_l)
                        *alms(0,mtpl_l,mtpl_m);
            }
        }

        for(std::size_t fld_j = 1; fld_j < num_fields; ++fld_j) {
            for(std::size_t mtpl_l=0; mtpl_l <= l_max; ++mtpl_l) {
                xlms(fld_i,mtpl_l,0)
                    = cinvls(fld_i,fld_j,mtpl_l)*alms(fld_j,mtpl_l,0);
            }
            for(std::size_t mtpl_m = 1; mtpl_m <= m_max; ++mtpl_m) {
                for(std::size_t mtpl_l=mtpl_m; mtpl_l<=l_max; ++mtpl_l) {
                    xlms(fld_i,mtpl_l,mtpl_m)
                        = real_scalar_t(0.5)*cinvls(fld_i,0,mtpl_l)
                            *alms(fld_j,mtpl_l,mtpl_m);
                }
            }
        }
    }

    sht<real_scalar_t> sh_trans(num_fields,l_max,m_max,num_pixels);
    sph_data<real_scalar_t> maps(spins,num_pixels);
    sh_trans.synthesise(alms,maps);
    sph_diag_prec_mat<real_scalar_t> p_mat(num_fields,num_pixels);
    real_scalar_t const sigma_pixel = 1e-3;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            p_mat(fld_i,pix_i) = 1./(sigma_pixel*sigma_pixel);
        }
    }
    sph_data<real_scalar_t> noise_maps
        = gen_gauss_noise_data(p_mat,spins,rng);
    maps.data() += noise_maps.data();
    win_func<real_scalar_t> w_func(num_fields,l_max);
    balan_2016<real_scalar_t> lp_t08(maps,p_mat,w_func);
    std::size_t num_cls = fls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = xlms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_t> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_t>(xlms,fls,pos_q);
    lp_t08.log_post(pos_q);
}

BOOST_AUTO_TEST_CASE(balan_2016_init){
    test_balan_2016_init<double>();
}
