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
#include <blackpearl/core/pow_spec.hpp>
#include <blackpearl/core/shc_ps_utils.hpp>

template<typename real_scalar_type>
void test_sht(){
    using namespace blackpearl::core;
    size_t n_side(512);
    size_t l_max = 2*n_side;
    size_t m_max(l_max);
    size_t num_pixels = 12*n_side*n_side;
    std::vector<size_t> spins = {0};//,0};//,2,2,2,2};
    size_t num_fields = spins.size();
    sph_hrm_coeffs<real_scalar_type> sh(num_fields,l_max,m_max);
    sph_data<real_scalar_type> data(spins,num_pixels);
    std::mt19937 rng(1234);
    std::uniform_real_distribution<real_scalar_type> uni_real_dist;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(size_t fld_j = 0; fld_j < num_pixels; ++fld_j) {
            data(fld_i,fld_j) = uni_real_dist(rng);
        }
    }
    sht<real_scalar_type> sht_test(num_fields,l_max,m_max,num_pixels);
    sht_test.analyse(data,sh);
    sph_data<real_scalar_type> data_test(spins,num_pixels);
    sht_test.synthesise(sh,data_test);
    real_scalar_type exp_rms
        = std::numeric_limits<real_scalar_type>::epsilon()*1e12;
    real_scalar_type rms = 0.;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(size_t fld_j = 0; fld_j < num_pixels; ++fld_j) {
            rms += (data(fld_i,fld_j) - data_test(fld_i,fld_j))*
                (data(fld_i,fld_j) - data_test(fld_i,fld_j));
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
    sph_hrm_coeffs<real_scalar_type> shcfs(num_fields,l_max,m_max);
    for(size_t mtpl_m = 0; mtpl_m <= m_max;++mtpl_m) {
        for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max;++mtpl_l) {
            shcfs(0,mtpl_l,mtpl_m) = complex_scalar_type (1,0);
        }
    }
    sph_data<real_scalar_type> data(spins,num_pixels);
    sht<real_scalar_type> sht_test(num_fields,l_max,m_max,num_pixels);
    sht_test.synthesise(shcfs,data);

    sph_hrm_coeffs<real_scalar_type> shcfs_test(num_fields,l_max,m_max);
    sht_test.analyse(data,shcfs_test);
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    for(size_t mtpl_m = 0; mtpl_m <= m_max;++mtpl_m) {
        for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max;++mtpl_l) {
            // NOTE we can conly achieve limited amount of accuracy here
            // Y^-1 is will not give you the value back with great accuracy :(
            BOOST_CHECK(
                std::abs(
                    shcfs_test(0,mtpl_l,mtpl_m).real() - real_scalar_type(1.)
                )
                <= eps*1e15
            );
            BOOST_CHECK(
                std::abs(
                    shcfs_test(0,mtpl_l,mtpl_m).imag() - real_scalar_type(0.)
                )
                <= eps*1e15
            );
        }
    }
}

template<typename real_scalar_type>
struct sht_unit_pow_spec_scale;

// template<>
// struct sht_unit_pow_spec_scale<float>{
//     static constexpr float const val = 1e7;
// };
//
template<>
struct sht_unit_pow_spec_scale<double>{
    static constexpr double const val = 1e16;
};

template<typename real_scalar_type>
void test_sht_unit_pow_spec() {
    using namespace blackpearl::core;
    using namespace boost::numeric::ublas;

    std::vector<size_t> spins = {0,0,2,2};
    size_t const num_fields = spins.size();
    size_t const l_max = 32;
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
    sph_hrm_coeffs<real_scalar_type> alms_test(num_fields,l_max,m_max);
    sh_trans.analyse(maps,alms_test);

    real_scalar_type omega_pix = 1.;//4*M_PI/(real_scalar_type) num_pixels;
    real_scalar_type eps = std::numeric_limits<real_scalar_type>::epsilon();
    eps *= sht_unit_pow_spec_scale<real_scalar_type>::val;
    for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i) {
        for(std::size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            std::size_t l_min
                = spins[fld_i] == 0 ? mtpl_m :  std::max(size_t(2),mtpl_m);
            for(std::size_t mtpl_l = l_min ; mtpl_l <= l_max; ++mtpl_l) {
                // std::cout<< alms(fld_i,mtpl_l,mtpl_m) << "\t"
                // << alms_test(fld_i,mtpl_l,mtpl_m)/omega_pix
                // << "\t" << std::abs(
                //     alms(fld_i,mtpl_l,mtpl_m).real()
                //     - alms_test(fld_i,mtpl_l,mtpl_m).real() /omega_pix
                // ) << std::endl;
                BOOST_CHECK(
                    std::abs(
                        alms(fld_i,mtpl_l,mtpl_m).real()
                        - alms_test(fld_i,mtpl_l,mtpl_m).real()/omega_pix
                    )
                    <= eps
                );
                BOOST_CHECK(
                    std::abs(
                        alms(fld_i,mtpl_l,mtpl_m).imag()
                        - alms_test(fld_i,mtpl_l,mtpl_m).imag()/omega_pix
                    )
                    <= eps
                );
            }
        }
    }
}


// BOOST_AUTO_TEST_CASE(sph_data){
//     test_sht<float>();
//     test_sht<double>();
// }
//
// BOOST_AUTO_TEST_CASE(sph_data_unit_alms){
//     test_sht_unit_alms<float>();
//     test_sht_unit_alms<double>();
// }

BOOST_AUTO_TEST_CASE(sht_unit_pow_spec){
    // test_sht_unit_pow_spec<float>();
    test_sht_unit_pow_spec<double>();
}
