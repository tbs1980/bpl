#ifndef BLACKPEARL_LOG_POST_BALAN_2016_HPP
#define BLACKPEARL_LOG_POST_BALAN_2016_HPP

namespace blackpearl{ namespace log_post {

template<typename real_scalar_type>
class balan_2016{
public:
    balan_2016(
        blackpearl::core::sph_data<real_scalar_type> const & data,
        blackpearl::core::sph_diag_prec_mat<real_scalar_type> const & p_mat,
        blackpearl::core::win_func<real_scalar_type> const & w_func
    ) throw()
    :m_data(data)
    ,m_prec_mat(p_mat)
    ,m_win_func(w_func){
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
        size_t const num_real_alms
            = sph_hrm_coeffs<real_scalar_type>::num_real_indep_coeffs(
                m_num_fields,
                m_l_max,
                m_m_max
            );
        size_t const num_real_cls
            = pow_spec<real_scalar_type>::num_real_indep_coeffs(
                m_num_fields,
                m_l_max
            );
        m_num_real_coeffs = num_real_alms + num_real_cls;
    }
    
    real_scalar_type log_post(real_vector_type const & pos_q){
        using namespace blackpearl::core;
        using namespace boost::numeric::ublas;
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_type> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_type> ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_type>(pos_q, shc_a, ps_c);
        pow_spec<real_scalar_type> ps_sigma =  extract_pow_spec(shc_a);
        real_scalar_type log_post = 0;
        for(size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l){
            matrix<real_scalar_type> c_l(m_num_fields,m_num_fields);
            matrix<real_scalar_type> sigma_l(m_num_fields,m_num_fields);
            matrix<real_scalar_type> c_inv_l(m_num_fields,m_num_fields);
            matrix<real_scalar_type> c_inv_sigma_l(m_num_fields,m_num_fields);
            ps_c.get_mtpl(mtpl_l,c_l);
            ps_sigma.get_mtpl(mtpl_l,sigma_l);
            real_scalar_type fact_l = real_scalar_type(2*mtpl_l+1);
            real_scalar_type log_det_c_l
                = compute_inverse<real_scalar_type>(c_l);
            real_scalar_type tr_c_l
                = compute_trace<real_scalar_type>(c_l);
            bool has_inv = compute_inverse<real_scalar_type>(c_l,c_inv_l);
            BOOST_ASSERT(has_inv == true);
            c_inv_sigma_l = prod(c_inv_l,sigma_l);
            real_scalar_type tr_cl_inv_sigma_l
                = compute_trace<real_scalar_type>(c_inv_sigma_l);
            log_post += (
                fact_l*log_det_c_l
                + 0.5*fact_l*tr_cl_inv_sigma_l
                + fact_l*tr_c_l
            );
        }

    }

    real_vector_type grad_log_post(real_vector_type const & pos_q) {
        using namespace blackpearl::core;
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);

        sph_hrm_coeffs<real_scalar_type> alms(
            m_num_fields,
            m_l_max,
            m_m_max
        );
        pow_spec<real_scalar_type> cls(m_num_fields,m_l_max);
        convert_to_coeffs<real_scalar_type>(pos_q, alms, cls);
    }

private:
    blackpearl::core::sph_data<real_scalar_type> m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_type> m_prec_mat;
    blackpearl::core::win_func<real_scalar_type> m_win_func;
    size_t m_num_fields;
    size_t m_l_max;
    size_t m_m_max;
    size_t m_num_real_coeffs;
};

}}

#endif //BLACKPEARL_LOG_POST_BALAN_2016_HPP
