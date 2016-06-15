#ifndef BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
#define BLACKPEARL_LOG_POST_TAYLOR_2008_HPP

#include <memory>
#include <cmath>
#include <complex>
#include <exception>
#include <sstream>
#include <string>

#include "sph_data.hpp"
#include "sph_hrm_coeffs.hpp"
#include "sph_prec_mat.hpp"

namespace blackpearl{ namespace log_post {

template<typename real_scalar_type>
class taylor_2008{
public:
    static_assert(std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type");

    static const indexType num_mc_samples = 1000;

    taylor_2008(
        blackpearl::core::shp_data<real_scalar_type> const & data,
        sph_diag_prec_mat<real_scalar_type> const & prec_mat,
        windowFuncType const & R
    ){

    }
};

}}

#endif //BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
