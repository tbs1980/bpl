#ifndef BLACKPEARL_SPH_DATA_PREC_MAT_UTILS_HPP
#define BLACKPEARL_SPH_DATA_PREC_MAT_UTILS_HPP

#include <cstddef>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <random>
#include <cmath>
#include <boost/assert.hpp>
#include "sph_data.hpp"
#include "sph_prec_mat.hpp"

namespace blackpearl{ namespace core {

template<typename real_scalar_type, class rng_type>
blackpearl::core::sph_data<real_scalar_type> gen_gauss_noise_data(
    blackpearl::core::sph_diag_prec_mat<real_scalar_type> const & p_mat,
    std::vector<size_t> const & spins,
    rng_type & rng
) {
    BOOST_ASSERT( spins.size() == p_mat.num_fields() );
    size_t const num_fields =  p_mat.num_fields();
    size_t const num_pixels = p_mat.num_pixels();
    blackpearl::core::sph_data<real_scalar_type> maps(spins,num_pixels);
    std::normal_distribution<real_scalar_type> norm_dist;
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
            maps(fld_i,pix_i) = norm_dist(rng)/std::sqrt( p_mat(fld_i,pix_i) );
        }
    }
    return maps;
}

}}

#endif //BLACKPEARL_SPH_DATA_PREC_MAT_UTILS_HPP
