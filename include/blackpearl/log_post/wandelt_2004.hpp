#ifndef BLACKPEARL_LOG_POST_WANDELT_2004_HPP
#define BLACKPEARL_LOG_POST_WANDELT_2004_HPP

#include <sstream>
#include <string>
#include <exception>

#include <boost/assert.hpp>

#include "../core/sph_data.hpp"
#include "../core/sph_prec_mat.hpp"
#include "../core/win_func.hpp"
#include "../core/sph_hrm_coeffs.hpp"
#include "../core/pow_spec.hpp"
#include "../core/shc_ps_utils.hpp"

namespace blackpearl{ namespace log_post {

template<typename real_scalar_type>
class wandelt_2004{
public:
    static_assert(std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type");
    typedef boost::numeric::ublas::vector<real_scalar_type> real_vector_type;

    wandelt_2004(
        blackpearl::core::sph_data<real_scalar_type> const & data,
        blackpearl::core::sph_diag_prec_mat<real_scalar_type> const & p_mat,
        blackpearl::core::win_func<real_scalar_type> const & w_func
    )
    : m_data(data)
    , m_prec_mat(p_mat)
    , m_win_func(w_func) {
        using namespace blackpearl::core;
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
        std::size_t const num_real_alms
            = sph_hrm_coeffs<real_scalar_type>::num_real_indep_coeffs(
                m_num_fields,
                m_l_max,
                m_m_max
            );
        std::size_t const num_real_cls
            = pow_spec<real_scalar_type>::num_real_indep_coeffs(
                m_num_fields,
                m_l_max
            );
        m_num_real_coeffs = num_real_alms + num_real_cls;
    }

    ~wandelt_2004(){

    }

    real_scalar_type log_post(real_vector_type const & pos_q){
        using namespace blackpearl::core;
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_type> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_type> ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_type>(pos_q, shc_a, ps_c);
        pow_spec<real_scalar_type> ps_sigma =  extract_pow_spec(shc_a);

        real_scalar_type log_prior = 0;
        for(std::size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l) {
            
        }
    }

private:
    blackpearl::core::sph_data<real_scalar_type> const & m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_type> const & m_prec_mat;
    blackpearl::core::win_func<real_scalar_type> const & m_win_func;
    std::size_t m_num_fields;
    std::size_t m_l_max;
    std::size_t m_m_max;
    std::size_t m_num_real_coeffs;
};

}}

#endif // BLACKPEARL_LOG_POST_WANDELT_2004_HPP
