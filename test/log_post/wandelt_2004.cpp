#define BOOST_TEST_MODULE wandelt 2004
#define BOOST_TEST_DYN_LINK

#include <cstddef>
#include <vector>
#include <random>

#include <boost/test/unit_test.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/assert.hpp>

#include <blackpearl/core/pow_spec.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <blackpearl/core/sht.hpp>
#include <blackpearl/core/sph_data.hpp>
#include <blackpearl/core/shc_ps_utils.hpp>
#include <blackpearl/core/sph_prec_mat.hpp>
#include <blackpearl/core/sph_data_prec_mat_utils.hpp>
#include <blackpearl/core/win_func.hpp>
#include <blackpearl/log_post/wandelt_2004.hpp>

template<typename real_scalar_type>
void test_wandelt_2004_init(){
    using namespace blackpearl::core;
    using namespace blackpearl::log_post;
    using namespace boost::numeric::ublas;
    typedef vector<real_scalar_type> real_vector_type;

    std::vector<std::size_t> spins = {0}; // FIXME more than one filed is necessary
    std::size_t const num_fields = spins.size();
    BOOST_ASSERT(num_fields == 1);
    std::size_t const l_max = 8;
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

    std::size_t num_cls = cls.num_real_indep_coeffs(num_fields,l_max);
    std::size_t num_alms = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
    vector<real_scalar_type> pos_q(num_cls+num_alms);
    convert_to_real_vector<real_scalar_type>(alms,cls,pos_q);

    wandelt_2004<real_scalar_type> lp_w04(maps,p_mat,w_func);
    lp_w04.log_post(pos_q);

    convert_to_real_vector<real_scalar_type>(alms,cls,pos_q);
    real_vector_type grad_pos_q = lp_w04.grad_log_post(pos_q);
}

BOOST_AUTO_TEST_CASE(wandelt_2004_init){
    test_wandelt_2004_init<double>();
}
