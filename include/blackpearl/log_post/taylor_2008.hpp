#ifndef BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
#define BLACKPEARL_LOG_POST_TAYLOR_2008_HPP

#include <memory>
#include <cmath>
#include <complex>
#include <exception>
#include <sstream>
#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include "sph_data.hpp"
#include "sph_hrm_coeffs.hpp"
#include "sph_prec_mat.hpp"

namespace blackpearl{ namespace log_post {

template<typename real_scalar_type>
class taylor_2008{
public:
    static_assert(std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type");

    typedef boost::numeric::ublas::vector<real_scalar_type> real_vector_type;

    static const indexType num_mc_samples = 1000;

    taylor_2008(
        blackpearl::core::shp_data<real_scalar_type> const & data,
        blackpearl::core::sph_diag_prec_mat<real_scalar_type> const & p_mat,
        blackpearl::core::win_func<real_scalar_type> const & w_func
    ){

    }

    real_scalar_type log_post(real_vector_type const & pos_q){

    }

    real_vector_type grad_log_post(real_vector_type const & pos_q) {

    }

private:
    blackpearl::core::shp_data<real_scalar_type> m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_type> m_prec_mat;
    blackpearl::core::win_func<real_scalar_type> m_win_func;
};

}}

#endif //BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
