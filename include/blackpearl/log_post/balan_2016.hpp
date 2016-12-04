#ifndef BLACKPEARL_LOG_POST_BALAN_2016_HPP
#define BLACKPEARL_LOG_POST_BALAN_2016_HPP

#include <memory>
#include <cmath>
#include <complex>
#include <exception>
#include <sstream>
#include <string>
#include <cstddef>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "../core/pow_spec.hpp"
#include "../core/shc_ps_utils.hpp"
#include "../core/sht.hpp"
#include "../core/sph_data.hpp"
#include "../core/sph_hrm_coeffs.hpp"
#include "../core/sph_prec_mat.hpp"
#include "../core/win_func.hpp"
#include "../core/win_func_shc_utils.hpp"
#include "../utils/lin_alg_utils.hpp"

namespace blackpearl{ namespace log_post {

template<typename real_scalar_t>
class balan_2016{
public:
    static_assert(std::is_floating_point<real_scalar_t>::value,
        "The real_scalar_t should be a floating point type");
    typedef std::complex<real_scalar_t> complex_scalar_type;
    typedef boost::numeric::ublas::vector<real_scalar_t> real_vector_t;
    typedef boost::numeric::ublas::vector<complex_scalar_type>
        complex_vector_type;

    balan_2016(
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
    ){
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

    real_scalar_t log_post(real_vector_t const & pos_q){
        using namespace blackpearl::core;
        using namespace blackpearl::utils;
        using namespace boost::numeric::ublas;
        typedef matrix<real_scalar_t> real_matrix_t;
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_t> shc_x(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_t> ps_f(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(pos_q, shc_x, ps_f);

        pow_spec<real_scalar_t> ps_c(m_num_fields, m_l_max);
        for(std::size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l){
            real_matrix_t f_l(m_num_fields, m_num_fields);
            ps_f.get_mtpl(mtpl_l,f_l);
            real_matrix_t c_l = compute_matrix_log(f_l);
            ps_c.set_mtpl(mtpl_l,c_l);
        }
        sph_hrm_coeffs<real_scalar_t> shc_a(m_num_fields, m_l_max, m_m_max);
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
            for(std::size_t mtpl_l=0; mtpl_l <= m_l_max; ++mtpl_l) {
                shc_a(fld_i,mtpl_l,0)
                    = ps_c(fld_i,0,mtpl_l)*shc_x(0,mtpl_l,0);
            }
            for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                for(std::size_t mtpl_l=mtpl_m; mtpl_l<=m_l_max; ++mtpl_l) {
                    shc_a(fld_i,mtpl_l,mtpl_m)
                        = real_scalar_t(0.5)*ps_c(fld_i,0,mtpl_l)
                            *shc_x(0,mtpl_l,mtpl_m);
                }
            }

            for(std::size_t fld_j = 1; fld_j < m_num_fields; ++fld_j) {
                for(std::size_t mtpl_l=0; mtpl_l <= m_l_max; ++mtpl_l) {
                    shc_a(fld_i,mtpl_l,0)
                        = ps_c(fld_i,fld_j,mtpl_l)*shc_x(fld_j,mtpl_l,0);
                }
                for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                    for(std::size_t mtpl_l=mtpl_m; mtpl_l<=m_l_max; ++mtpl_l) {
                        shc_a(fld_i,mtpl_l,mtpl_m)
                            = real_scalar_t(0.5)*ps_c(fld_i,0,mtpl_l)
                                *shc_x(fld_j,mtpl_l,mtpl_m);
                    }
                }
            }
        }

        pow_spec<real_scalar_t> ps_sigma =  extract_pow_spec(shc_a);

        real_scalar_t log_prior = 0;
        for(std::size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l){
            real_matrix_t c_l(m_num_fields,m_num_fields);
            ps_c.get_mtpl(mtpl_l,c_l);
            real_matrix_t sigma_l(m_num_fields,m_num_fields);
            ps_sigma.get_mtpl(mtpl_l,sigma_l);
            real_scalar_t const fact_l(2*mtpl_l+1);
            real_scalar_t const det_cl
                = compute_determinant<real_scalar_t>(c_l);
            real_matrix_t c_inv_l(m_num_fields,m_num_fields);
            compute_inverse<real_scalar_t>(c_l,c_inv_l);
            real_matrix_t c_inv_sigma_l = prod(c_inv_l,sigma_l);
            real_scalar_t const trace_cl_inv_sig_l =
                compute_trace<real_scalar_t>(c_inv_sigma_l);
            real_scalar_t const trace_cl =
                compute_trace<real_scalar_t>(c_l);
            log_prior += 0.5*fact_l*std::log(det_cl)
                - 0.5*fact_l*trace_cl_inv_sig_l
                - fact_l*trace_cl;
        }

        sph_hrm_coeffs<real_scalar_t> shc_a_fwd(shc_a);
        apply_win_func<real_scalar_t>(m_win_func,shc_a_fwd);
        sph_data<real_scalar_t> data_fwd(m_data.spins(),m_data.num_pixels());
        m_sh_trans.synthesise(shc_a_fwd,data_fwd);
        real_scalar_t log_lik = 0;
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
            for(std::size_t pix_i = 0; pix_i < m_num_pixels; ++pix_i){
                real_scalar_t diff
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
        using namespace boost::numeric::ublas;
        typedef matrix<real_scalar_t> real_matrix_t;
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_t> shc_x(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_t> ps_f(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_t>(pos_q, shc_x, ps_f);

        pow_spec<real_scalar_t> ps_c(m_num_fields, m_l_max);
        for(std::size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l){
            real_matrix_t f_l(m_num_fields, m_num_fields);
            ps_f.get_mtpl(mtpl_l,f_l);
            real_matrix_t c_l = compute_matrix_log(f_l);
            ps_c.set_mtpl(mtpl_l,c_l);
        }
        sph_hrm_coeffs<real_scalar_t> shc_a(m_num_fields, m_l_max, m_m_max);
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i) {
            for(std::size_t mtpl_l=0; mtpl_l <= m_l_max; ++mtpl_l) {
                shc_a(fld_i,mtpl_l,0)
                    = ps_c(fld_i,0,mtpl_l)*shc_x(0,mtpl_l,0);
            }
            for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                for(std::size_t mtpl_l=mtpl_m; mtpl_l<=m_l_max; ++mtpl_l) {
                    shc_a(fld_i,mtpl_l,mtpl_m)
                        = real_scalar_t(0.5)*ps_c(fld_i,0,mtpl_l)
                            *shc_x(0,mtpl_l,mtpl_m);
                }
            }

            for(std::size_t fld_j = 1; fld_j < m_num_fields; ++fld_j) {
                for(std::size_t mtpl_l=0; mtpl_l <= m_l_max; ++mtpl_l) {
                    shc_a(fld_i,mtpl_l,0)
                        = ps_c(fld_i,fld_j,mtpl_l)*shc_x(fld_j,mtpl_l,0);
                }
                for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m) {
                    for(std::size_t mtpl_l=mtpl_m; mtpl_l<=m_l_max; ++mtpl_l) {
                        shc_a(fld_i,mtpl_l,mtpl_m)
                            = real_scalar_t(0.5)*ps_c(fld_i,0,mtpl_l)
                                *shc_x(fld_j,mtpl_l,mtpl_m);
                    }
                }
            }
        }

        pow_spec<real_scalar_t> ps_sigma =  extract_pow_spec(shc_a);

        pow_spec<real_scalar_t> ps_dg(m_num_fields, m_l_max);
        pow_spec<real_scalar_t> ps_c_inv(m_num_fields, m_l_max);
        for(std::size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l){
            real_matrix_t g_l(m_num_fields, m_num_fields);
            // ps_g.get_mtpl(mtpl_l,g_l); //FIXME
            real_matrix_t c_l = compute_matrix_exp(g_l);
            real_matrix_t sigma_l(m_num_fields,m_num_fields);
            ps_sigma.get_mtpl(mtpl_l,sigma_l);
            real_matrix_t c_inv_l(m_num_fields,m_num_fields);
            compute_inverse<real_scalar_t>(c_l,c_inv_l);
            real_matrix_t c_inv_sigma_l = prod(c_inv_l,sigma_l);
            real_matrix_t id_I
                = identity_matrix<real_scalar_t>(m_num_fields);
            real_matrix_t phi_l = 0.5*(2.*mtpl_l+1.)*(c_inv_sigma_l - id_I);
            real_matrix_t dg_l
                = 2.*phi_l - element_prod(phi_l,id_I) + id_I;
            ps_dg.set_mtpl(mtpl_l,dg_l);
            ps_c_inv.set_mtpl(mtpl_l,c_inv_l);
        }

        // sph_hrm_coeffs<real_scalar_t> shc_a_fwd(shc_a);
        // apply_win_func<real_scalar_t>(m_win_func,shc_a_fwd);
        // sph_data<real_scalar_t> data_fwd(m_data.spins(),m_data.num_pixels());
        // m_sh_trans.synthesise(shc_a_fwd, data_fwd);
        // for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
        //     for(std::size_t pix_i = 0; pix_i < m_num_pixels; ++pix_i){
        //         data_fwd(fld_i,pix_i)
        //             = ( m_data(fld_i,pix_i) - data_fwd(fld_i,pix_i) )
        //                 *m_prec_mat(fld_i,pix_i);
        //     }
        // }
        // m_sh_trans.analyse(data_fwd,shc_a_fwd);
        // real_scalar_t omega_pix = 4.*M_PI/(real_scalar_t) m_num_pixels;
        // for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
        //     for(std::size_t mtpl_l=0; mtpl_l <= m_l_max; ++mtpl_l){
        //         shc_a_fwd(fld_i,mtpl_l,0)
        //             *=  m_win_func(fld_i,mtpl_l)/omega_pix;
        //     }
        //     for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m){
        //         for(std::size_t mtpl_l= mtpl_m; mtpl_l <= m_l_max; ++mtpl_l){
        //             shc_a_fwd(fld_i,mtpl_l,mtpl_m)
        //                 *= 2.*m_win_func(fld_i,mtpl_l)/omega_pix;
        //         }
        //     }
        //
        //     for(std::size_t fld_j = 0; fld_j < m_num_fields; ++fld_j){
        //         for(std::size_t mtpl_l=0; mtpl_l <= m_l_max; ++mtpl_l){
        //             shc_a_fwd(fld_i,mtpl_l,0)
        //                 -=  ps_c_inv(fld_i,fld_j,mtpl_l)*shc_a(fld_j,mtpl_l,0);
        //         }
        //         for(std::size_t mtpl_m = 1; mtpl_m <= m_m_max; ++mtpl_m){
        //             for(std::size_t mtpl_l=mtpl_m; mtpl_l<=m_l_max; ++mtpl_l){
        //                 shc_a_fwd(fld_i,mtpl_l,mtpl_m)
        //                     -= real_scalar_t(2)*ps_c_inv(fld_i,fld_j,mtpl_l)
        //                         *shc_a(fld_j,mtpl_l,mtpl_m);
        //             }
        //         }
        //     }
        // }
        //
        // real_vector_t d_pos_q( pos_q.size() );
        // convert_to_real_vector<real_scalar_t>(shc_a_fwd,ps_dg,d_pos_q);
        // return d_pos_q;
    }

private:
    blackpearl::core::sph_data<real_scalar_t> m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_t> m_prec_mat;
    blackpearl::core::win_func<real_scalar_t> m_win_func;
    blackpearl::core::sht<real_scalar_t> m_sh_trans;
    std::size_t m_num_fields;
    std::size_t m_l_max;
    std::size_t m_m_max;
    std::size_t m_num_real_coeffs;
    std::size_t m_num_pixels;
};

}}

#endif //BLACKPEARL_LOG_POST_BALAN_2016_HPP
