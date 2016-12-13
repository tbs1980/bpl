#ifndef BLACKPEARL_LOG_POST_WANDELT_2004_HPP
#define BLACKPEARL_LOG_POST_WANDELT_2004_HPP

#include <cstddef>
#include <sstream>
#include <string>
#include <exception>
#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/assert.hpp>

#include "../core/sph_data.hpp"
#include "../core/sph_prec_mat.hpp"
#include "../core/win_func.hpp"
#include "../core/sph_hrm_coeffs.hpp"
#include "../core/pow_spec.hpp"
#include "../core/shc_ps_utils.hpp"
#include "../utils/lin_alg_utils.hpp"
#include "../core/win_func_shc_utils.hpp"
#include "../core/sht.hpp"

namespace blackpearl{ namespace log_post {

template<typename real_scalar_t>
class wandelt_2004{
public:
    static_assert(std::is_floating_point<real_scalar_t>::value,
        "The real_scalar_t should be a floating point type");
    typedef boost::numeric::ublas::vector<real_scalar_t> real_vector_t;
    typedef boost::numeric::ublas::matrix<real_scalar_t> real_matrix_t;
    typedef boost::numeric::ublas::unbounded_array<real_matrix_t> real_matrix_array_t;

    wandelt_2004(
        blackpearl::core::sph_data<real_scalar_t> const & data,
        blackpearl::core::sph_diag_prec_mat<real_scalar_t> const & p_mat,
        blackpearl::core::win_func<real_scalar_t> const & w_func
    )
    : m_data(data)
    , m_prec_mat(p_mat)
    , m_win_func(w_func)
    , m_sh_trans(
        data.num_fields(),
        w_func.l_max(),
        w_func.l_max(),
        data.num_pixels()
    ) {
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
            = sph_hrm_coeffs<real_scalar_t>::num_real_indep_coeffs(
                m_num_fields,
                m_l_max,
                m_m_max
            );
        std::size_t const num_real_cls
            = pow_spec<real_scalar_t>::num_real_indep_coeffs(
                m_num_fields,
                m_l_max
            );
        m_num_real_coeffs = num_real_alms + num_real_cls;
        m_num_pixels = data.num_pixels();
    }

    ~wandelt_2004(){

    }

    real_scalar_t log_post(real_vector_t const & pos_q){
        using namespace blackpearl::core;
        using namespace blackpearl::utils;

        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_t> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_t> ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(pos_q, shc_a, ps_c);

        pow_spec<real_scalar_t> ps_sigma =  extract_pow_spec(shc_a);
        real_scalar_t log_prior = 0;
        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
            real_matrix_t c_l(m_num_fields, m_num_fields);
            ps_c.get_mtpl(mtpl_l,c_l);
            real_matrix_t sigma_l(m_num_fields,m_num_fields);
            ps_sigma.get_mtpl(mtpl_l,sigma_l);
            real_scalar_t const fact_l = 0.5*(2*mtpl_l+1);
            real_scalar_t const det_cl
                = compute_determinant<real_scalar_t>(c_l);
            BOOST_ASSERT(det_cl > 0);
            real_matrix_t c_inv_l(m_num_fields,m_num_fields);
            compute_inverse<real_scalar_t>(c_l,c_inv_l);
            real_matrix_t c_inv_sigma_l = prod(c_inv_l,sigma_l);
            real_scalar_t const trace_cl_inv_sig_l
                = compute_trace<real_scalar_t>(c_inv_sigma_l);
            log_prior += -fact_l*( std::log(det_cl) + trace_cl_inv_sig_l );
        }

        sph_hrm_coeffs<real_scalar_t> shc_a_fwd(shc_a);
        apply_win_func<real_scalar_t>(m_win_func,shc_a_fwd);
        sph_data<real_scalar_t> data_fwd(m_data.spins(),m_data.num_pixels());
        real_scalar_t log_lik = 0;
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
            for(std::size_t pix_i = 0; pix_i < m_num_pixels; ++pix_i) {
                real_scalar_t const diff
                    = m_data(fld_i,pix_i) - data_fwd(fld_i,pix_i);
                log_lik -= diff*diff*m_prec_mat(fld_i,pix_i);
            }
        }
        log_lik *= 0.5;

        return log_prior + log_lik;
    }

    real_vector_t grad_log_post(real_vector_t const & pos_q) {
        using namespace blackpearl::core;
        using namespace blackpearl::utils;
        using namespace boost::math::constants;

        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_t> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_t> ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(pos_q, shc_a, ps_c);

        pow_spec<real_scalar_t> ps_sigma =  extract_pow_spec(shc_a);
        pow_spec<real_scalar_t> ps_dc(m_num_fields, m_l_max);
        pow_spec<real_scalar_t> ps_c_inv(m_num_fields, m_l_max);
        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
            real_matrix_t c_l(m_num_fields, m_num_fields);
            ps_c.get_mtpl(mtpl_l,c_l);
            real_matrix_t sigma_l(m_num_fields,m_num_fields);
            ps_sigma.get_mtpl(mtpl_l,sigma_l);
            real_scalar_t const fact_l = 0.5*(2*mtpl_l+1);
            real_matrix_t c_inv_l(m_num_fields,m_num_fields);
            compute_inverse<real_scalar_t>(c_l,c_inv_l);
            ps_c_inv.set_mtpl(mtpl_l,c_inv_l);
            real_matrix_t const c_inv_sigma_l = prod(c_inv_l,sigma_l);
            real_matrix_t const id_I = real_identity_matrix_t(m_num_fields);
            real_matrix_t const phi_l = prod((c_inv_sigma_l - id_I), c_inv_l);
            real_matrix_t const dc_l
                = fact_l*( 2.*phi_l - element_prod(phi_l,id_I) );
            ps_dc.set_mtpl(mtpl_l,dc_l);
        }

        sph_hrm_coeffs<real_scalar_t> shc_a_fwd(shc_a);
        apply_win_func<real_scalar_t>(m_win_func,shc_a_fwd);
        sph_data<real_scalar_t> data_fwd(m_data.spins(),m_data.num_pixels());
        m_sh_trans.synthesise(shc_a_fwd, data_fwd);
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
            for(std::size_t pix_i = 0; pix_i < m_num_pixels; ++pix_i) {
                data_fwd(fld_i,pix_i)
                    = ( m_data(fld_i,pix_i) - data_fwd(fld_i,pix_i) )
                        *m_prec_mat(fld_i,pix_i);
            }
        }
        m_sh_trans.analyse(data_fwd,shc_a_fwd);
        // FIXME need to see if this is really reaquired
        real_scalar_t const omega_pix
            = 4.*pi<real_scalar_t>()/(real_scalar_t) m_num_pixels;
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
            for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
                shc_a_fwd(fld_i,mtpl_l,0)
                    *=  m_win_func(fld_i,mtpl_l)/omega_pix;
            }
            for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                for(std::size_t mtpl_l = mtpl_m; mtpl_l <= m_l_max; ++mtpl_l) {
                    shc_a_fwd(fld_i,mtpl_l,mtpl_m)
                        *= 2.*m_win_func(fld_i,mtpl_l)/omega_pix;
                }
            }

            for(std::size_t fld_j = 0; fld_j < m_num_fields; ++fld_j) {
                for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
                    shc_a_fwd(fld_i,mtpl_l,0)
                        -=  ps_c_inv(fld_i,fld_j,mtpl_l)*shc_a(fld_j,mtpl_l,0);
                }
                for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                    for(std::size_t mtpl_l = mtpl_m; mtpl_l<=m_l_max; ++mtpl_l) {
                        shc_a_fwd(fld_i,mtpl_l,mtpl_m)
                            -= real_scalar_t(2)*ps_c_inv(fld_i,fld_j,mtpl_l)
                                *shc_a(fld_j,mtpl_l,mtpl_m);
                    }
                }
            }
        }

        real_vector_t d_pos_q( pos_q.size() );
        convert_to_real_vector<real_scalar_t>(shc_a_fwd,ps_dc,d_pos_q);
        return d_pos_q;
    }

    real_matrix_t metric_tensor_log_posterior_dense (
        real_vector_t const & pos_q
    ) const {
        using namespace blackpearl::core;
        using namespace blackpearl::utils;
        using namespace boost::numeric::ublas;
        typedef unit_vector<real_scalar_t> real_unit_vector_t;
        typedef sph_hrm_coeffs<real_scalar_t> sph_hrm_coeffs_t;
        typedef pow_spec<real_scalar_t> pow_spec_t;
        typedef sph_data<real_scalar_t> sph_data_t;

        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_t> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_t> ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(pos_q, shc_a, ps_c);

        pow_spec<real_scalar_t> ps_c_inv(m_num_fields, m_l_max);
        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
            real_matrix_t c_l(m_num_fields, m_num_fields);
            ps_c.get_mtpl(mtpl_l, c_l);
            real_matrix_t c_inv_l(m_num_fields, m_num_fields);
            compute_inverse<real_scalar_t>(c_l, c_inv_l);
            ps_c_inv.set_mtpl(mtpl_l, c_inv_l);
        }

        std::size_t const num_real_alms
            = sph_hrm_coeffs<real_scalar_t>::num_real_indep_coeffs(
            m_num_fields,
            m_l_max,
            m_m_max
        );
        std::size_t const num_real_cls
            = pow_spec<real_scalar_t>::num_real_indep_coeffs(
            m_num_fields,
            m_l_max
        );

        BOOST_ASSERT(num_real_alms + num_real_cls == m_num_real_coeffs);

        real_matrix_t mtrc_tnsr_g(m_num_real_coeffs,m_num_real_coeffs);

        for(std::size_t dim_i = 0; dim_i < num_real_alms; ++dim_i) {
            real_unit_vector_t unit_vect(m_num_real_coeffs, dim_i);
            sph_hrm_coeffs_t unit_shc_a(m_num_fields, m_l_max, m_m_max);
            pow_spec_t unit_ps_c(m_num_fields, m_l_max);
            convert_to_coeffs<real_scalar_t>(unit_vect, unit_shc_a, unit_ps_c);
            sph_hrm_coeffs_t unit_shc_a_temp(unit_shc_a);
            apply_win_func<real_scalar_t>(m_win_func, unit_shc_a);
            sph_data_t unit_data(m_data.spins(), m_data.num_pixels());
            m_sh_trans.synthesise(unit_shc_a, unit_data);
            for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
                for(std::size_t pix_i = 0; pix_i < m_num_pixels; ++pix_i) {
                    unit_data(fld_i, pix_i)
                        = unit_data(fld_i, pix_i)*m_prec_mat(fld_i, pix_i);
                }
            }
            m_sh_trans.analyse(unit_data, unit_shc_a);
            // FIXME do we need omega_pix here?
            apply_win_func<real_scalar_t>(m_win_func, unit_shc_a);
            for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
                for(std::size_t fld_j = 0; fld_j < m_num_fields; ++fld_j) {
                    for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
                        unit_shc_a(fld_i, mtpl_l, 0)
                            += ps_c_inv(fld_i, fld_j, mtpl_l)
                                *unit_shc_a_temp(fld_j, mtpl_l, 0);
                    }
                    for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                        for(std::size_t mtpl_l = mtpl_m; mtpl_l <= m_l_max; ++mtpl_l) {
                            unit_shc_a(fld_i, mtpl_l, mtpl_m)
                                += real_scalar_t(2)*ps_c_inv(fld_i,fld_j,mtpl_l)
                                    *unit_shc_a_temp(fld_j, mtpl_l, mtpl_m);
                        }
                    }
                }
            }
            real_vector_t mtrc_tnsr_g_col(m_num_real_coeffs);
            convert_to_real_vector<real_scalar_t>(unit_shc_a, unit_ps_c, mtrc_tnsr_g_col);
            for(std::size_t dim_j = 0; dim_j < num_real_alms; ++dim_j){
                mtrc_tnsr_g(dim_i, dim_j) = mtrc_tnsr_g_col(dim_j);
            }
        }

        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
            std::size_t const stride = num_real_alms + mtpl_l;
            for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
                for(std::size_t fld_j = 0; fld_j < m_num_fields; ++fld_j) {
                    mtrc_tnsr_g(stride + fld_i, stride + fld_j)
                        = 2.*ps_c_inv(fld_i, fld_j, mtpl_l);
                }
            }
        }

        return mtrc_tnsr_g;
    }

    real_matrix_array_t deriv_metric_tensor_log_posterior_dense(real_vector_t const & pos_q ) const {
        using namespace blackpearl::core;
        using namespace blackpearl::utils;
        using namespace boost::numeric::ublas;
        typedef unit_vector<real_scalar_t> real_unit_vector_t;
        typedef zero_vector<real_scalar_t> real_zero_vector_t;
        typedef sph_hrm_coeffs<real_scalar_t> sph_hrm_coeffs_t;
        typedef pow_spec<real_scalar_t> pow_spec_t;
        typedef sph_data<real_scalar_t> sph_data_t;

        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_t> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_t> ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(pos_q, shc_a, ps_c);

        pow_spec<real_scalar_t> ps_c_inv_2(m_num_fields, m_l_max);
        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
            real_matrix_t c_l(m_num_fields, m_num_fields);
            ps_c.get_mtpl(mtpl_l, c_l);
            real_matrix_t c_inv_l(m_num_fields, m_num_fields);
            compute_inverse<real_scalar_t>(c_l,c_inv_l);
            real_matrix_t c_inv_2_l = prod(c_inv_l, c_inv_l);
            ps_c_inv_2.set_mtpl(mtpl_l, c_inv_2_l);
        }

        real_matrix_array_t d_mtrc_tnsr_g(
            pos_q.size(),
            zero_matrix<real_scalar_t>( pos_q.size(), pos_q.size() )
        );

        std::size_t const num_real_alms
            = sph_hrm_coeffs<real_scalar_t>::num_real_indep_coeffs(
            m_num_fields,
            m_l_max,
            m_m_max
        );
        std::size_t const num_real_cls
            = pow_spec<real_scalar_t>::num_real_indep_coeffs(
            m_num_fields,
            m_l_max
        );

        BOOST_ASSERT(num_real_alms + num_real_cls == m_num_real_coeffs);

        std::size_t blk_pos = 0;
        for(std::size_t mtpl_blk_l = 0; mtpl_blk_l <= m_l_max; ++mtpl_blk_l) {
            std::size_t const num_real_alms_l
                = sph_hrm_coeffs<real_scalar_t>::num_real_indep_coeffs(
                m_num_fields,
                (mtpl_blk_l+1),
                (mtpl_blk_l+1)
            ) - sph_hrm_coeffs<real_scalar_t>::num_real_indep_coeffs(
                m_num_fields,
                mtpl_blk_l,
                mtpl_blk_l
            );
            for(std::size_t dim_i = 0; dim_i < num_real_alms_l; ++dim_i) {
                real_unit_vector_t unit_vect(m_num_real_coeffs, blk_pos + dim_i);
                sph_hrm_coeffs_t unit_shc_a(m_num_fields, m_l_max, m_m_max);
                pow_spec_t unit_ps_c(m_num_fields, m_l_max);
                convert_to_coeffs<real_scalar_t>(unit_vect, unit_shc_a, unit_ps_c);

                for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
                    for(std::size_t fld_j = 0; fld_j < m_num_fields; ++fld_j) {
                        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
                            unit_shc_a(fld_i, mtpl_l, 0)
                                = -ps_c_inv_2(fld_i, fld_j, mtpl_l)
                                    *unit_shc_a(fld_j, mtpl_l, 0);
                        }
                        for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                            for(std::size_t mtpl_l = mtpl_m; mtpl_l <= m_l_max; ++mtpl_l) {
                                unit_shc_a(fld_i, mtpl_l, mtpl_m)
                                    = -2.*ps_c_inv_2(fld_i,fld_j,mtpl_l)
                                        *unit_shc_a(fld_j, mtpl_l, mtpl_m);
                            }
                        }
                    }
                }
                real_vector_t d_mtrc_tnsr_g_col(m_num_real_coeffs);
                convert_to_real_vector<real_scalar_t>(unit_shc_a, unit_ps_c, d_mtrc_tnsr_g_col);
                for(std::size_t dim_j = 0; dim_j < num_real_alms; ++dim_j){
                    d_mtrc_tnsr_g[num_real_alms + mtpl_blk_l](blk_pos + dim_i, dim_j)
                        = mtrc_tnsr_g_col(dim_j);
                }
            }
            blk_pos += num_real_alms_l;
        }


        real_zero_vector_t zero_vect(m_num_real_coeffs);
        sph_hrm_coeffs_t zero_shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec_t zero_ps_c(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(zero_vect, zero_shc_a, zero_ps_c);
        pow_spec<real_scalar_t> ps_neg_2_over_c_inv_2(m_num_fields, m_l_max);
        for(std::size_t mtpl_l = 0; mtpl_l <= m_l_max; ++mtpl_l) {
            real_matrix_t c_inv_2_l(m_num_fields, m_num_fields);
            ps_c_inv_2.get_mtpl(mtpl_l, c_inv_2_l);
            c_inv_2_l *= -2.;
            ps_neg_2_over_c_inv_2.set_mtpl(mtpl_l, c_inv_2_l);
        }

        std::size_t blk_pos = num_real_alms;
        for(std::size_t mtpl_blk_l = 0; mtpl_blk_l <= m_l_max; ++mtpl_blk_l) {
            real_vector_t d_mtrc_tnsr_g_col(m_num_real_coeffs);
            convert_to_real_vector<real_scalar_t>(zero_shc_a, ps_neg_2_over_c_inv_2, d_mtrc_tnsr_g_col);
        }
    }

private:
    blackpearl::core::sph_data<real_scalar_t> const & m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_t> const & m_prec_mat;
    blackpearl::core::win_func<real_scalar_t> const & m_win_func;
    std::size_t m_num_fields;
    std::size_t m_l_max;
    std::size_t m_m_max;
    std::size_t m_num_real_coeffs;
    blackpearl::core::sht<real_scalar_t> m_sh_trans;
    std::size_t m_num_pixels;
};

}}

#endif // BLACKPEARL_LOG_POST_WANDELT_2004_HPP
