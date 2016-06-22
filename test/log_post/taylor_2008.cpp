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

template<typename real_scalar_type>
void test_taylor_2008_init(){
    using namespace blackpearl::core;
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
    convert_to_coeffs<real_scalar_type>(alms,cls,pos_q);
    lp_t08.log_post(pos_q);
}

BOOST_AUTO_TEST_CASE(taylor_2008_init){
    test_taylor_2008_init<float>();
    test_taylor_2008_init<double>();
}
