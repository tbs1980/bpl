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
    ) throw() {
        if(data.num_fields() != p_mat.num_fields()){
            std::stringstream msg;
            msg << "The number of fields"
                << " in the data "
                << data.num_fields()
                << " does not match the ones from precision matrix "
                << p_mat.num_fields()
                << std::endl;
            throw std::length_error(msg.str());
        }
        if(data.num_fields() != w_func.num_fields()){
            std::stringstream msg;
            msg << "The number of fields"
                << " in the data "
                << data.num_fields()
                << " does not match the ones from window function "
                << w_func.num_fields()
                << std::endl;
            throw std::length_error(msg.str());
        }
        m_num_fields = data.num_fields();
        m_l_max = w_func.l_max();
        m_m_max = w_func.l_max();
        size_t const num_real_alms
            = blackpearl::core::sph_hrm_coeffs<real_scalar_type>(
                m_num_fields,
                m_l_max,
                m_m_max
            );
        size_t const num_real_cls
            = blackpearl::core::pow_spec<real_scalar_type>(
                m_num_fields,
                m_l_max
            );
        m_num_real_coeffs = num_real_alms + num_real_cls;
    }

    real_scalar_type log_post(real_vector_type const & pos_q){
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);

    }

    real_vector_type grad_log_post(real_vector_type const & pos_q) {
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        blackpearl::core::sph_hrm_coeffs<real_scalar_type> alms()
    }

private:
    blackpearl::core::shp_data<real_scalar_type> m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_type> m_prec_mat;
    blackpearl::core::win_func<real_scalar_type> m_win_func;
    size_t m_num_fields;
    size_t m_l_max;
    size_t m_m_max;
    size_t m_num_real_coeffs;
};

}}

#endif //BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
