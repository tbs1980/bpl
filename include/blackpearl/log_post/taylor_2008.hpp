#ifndef BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
#define BLACKPEARL_LOG_POST_TAYLOR_2008_HPP

#include <memory>
#include <cmath>
#include <complex>
#include <exception>
#include <sstream>
#include <string>
#include <cstddef>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "../core/pow_spec.hpp"
#include "../core/shc_ps_utils.hpp"
#include "../core/sht.hpp"
#include "../core/sph_data.hpp"
#include "../core/sph_hrm_coeffs.hpp"
#include "../core/sph_prec_mat.hpp"
#include "../core/win_func.hpp"
#include "../core/win_func_shc_utils.hpp"
#include "../utils/EigenvalueDecomposition.hpp"

namespace blackpearl{ namespace log_post {

template<typename real_scalar_type>
real_scalar_type compute_determinant(
    boost::numeric::ublas::matrix<real_scalar_type>  m
) {
    // http://www.richelbilderbeek.nl/CppUblasMatrixExample7.htm
    using namespace boost::numeric::ublas;
    BOOST_ASSERT( m.size1() == m.size2() );
    permutation_matrix<std::size_t> pivots( m.size1() );
    int const is_singular = lu_factorize(m, pivots);
    if (is_singular){
        return real_scalar_type(0);
    }
    real_scalar_type d(1);
    std::size_t const sz = pivots.size();
    for (std::size_t i=0; i != sz; ++i){
        if ( pivots(i) != i ){
            d *= -1.;
        }
        d *= m(i,i);
    }
    return d;
}

template<typename real_scalar_type>
real_scalar_type compute_trace(
    boost::numeric::ublas::matrix<real_scalar_type> const & m
) {
    BOOST_ASSERT( m.size1() == m.size2() );
    real_scalar_type trace(0);
    for(std::size_t i=0; i<m.size1(); ++i){
        trace += m(i,i);
    }
    return trace;
}

template<typename real_scalar_type>
bool compute_inverse(
    boost::numeric::ublas::matrix<real_scalar_type> const & input,
    boost::numeric::ublas::matrix<real_scalar_type> & inverse
) {
    // https://gist.github.com/lilac/2464434
    BOOST_ASSERT( input.size1() == input.size2() );
    BOOST_ASSERT( input.size1() == inverse.size1() );
    BOOST_ASSERT( inverse.size1() == inverse.size2() );
    using namespace boost::numeric::ublas;
    matrix<real_scalar_type> A(input);
    permutation_matrix<std::size_t> pm(A.size1());
    int const res = lu_factorize(A,pm);
    if( res != 0 ) {
        return false;
    }
    inverse.assign( identity_matrix<real_scalar_type>( A.size1() ) );
    lu_substitute(A, pm, inverse);
    return true;
}

template<typename real_scalar_type>
boost::numeric::ublas::matrix<real_scalar_type> compute_matrix_exp (
    boost::numeric::ublas::matrix<real_scalar_type> const & mat_G
) {
    using namespace boost::numeric::ublas;
    BOOST_ASSERT( mat_G.size1() == mat_G.size2() );
    EigenvalueDecomposition<real_scalar_type> eig_decomp(mat_G);
    matrix<real_scalar_type> mat_D = eig_decomp.getD();
    matrix<real_scalar_type> mat_U = eig_decomp.getV();
    for(std::size_t ind_i = 0; ind_i < mat_G.size1(); ++ind_i) {
        mat_D(ind_i, ind_i) = std::exp( mat_D(ind_i, ind_i) );
    }
    matrix<real_scalar_type> mat_DU = prod(mat_D,trans (mat_U));
    matrix<real_scalar_type> mat_exp_G =  prod( mat_U, mat_DU);
    return mat_exp_G;
}

template<typename real_scalar_type>
boost::numeric::ublas::matrix<real_scalar_type> compute_matrix_log (
    boost::numeric::ublas::matrix<real_scalar_type> const & mat_C
) {
    using namespace boost::numeric::ublas;
    BOOST_ASSERT( mat_C.size1() == mat_C.size2() );
    EigenvalueDecomposition<real_scalar_type> eig_decomp(mat_C);
    matrix<real_scalar_type> mat_D = eig_decomp.getD();
    matrix<real_scalar_type> mat_U = eig_decomp.getV();
    for(std::size_t ind_i = 0; ind_i < mat_C.size1(); ++ind_i) {
        BOOST_ASSERT(mat_D(ind_i, ind_i) > 0);
        mat_D(ind_i, ind_i) = std::log( mat_D(ind_i, ind_i) );
    }
    matrix<real_scalar_type> mat_DU = prod(mat_D,trans (mat_U));
    matrix<real_scalar_type> mat_log_C =  prod( mat_U, mat_DU);
    return mat_log_C;
}

template<typename real_scalar_type>
class taylor_2008{
public:
    static_assert(std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type");

    typedef boost::numeric::ublas::vector<real_scalar_type> real_vector_type;

    static const size_t num_mc_samples = 1000;

    taylor_2008 (
        blackpearl::core::sph_data<real_scalar_type> const & data,
        blackpearl::core::sph_diag_prec_mat<real_scalar_type> const & p_mat,
        blackpearl::core::win_func<real_scalar_type> const & w_func
    ) throw()
    :m_data(data)
    ,m_prec_mat(p_mat)
    ,m_win_func(w_func)
    ,m_sh_trans(
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
        m_num_pixels = data.num_pixels();
    }

    real_scalar_type log_post(real_vector_type const & pos_q){
        using namespace blackpearl::core;
        using namespace boost::numeric::ublas;
        typedef matrix<real_scalar_type> real_matrix_type;
        BOOST_ASSERT(pos_q.size() == m_num_real_coeffs);
        sph_hrm_coeffs<real_scalar_type> shc_a(m_num_fields, m_l_max, m_m_max);
        pow_spec<real_scalar_type> ps_g(m_num_fields, m_l_max);
        convert_to_coeffs<real_scalar_type>(pos_q, shc_a, ps_g);
        pow_spec<real_scalar_type> ps_sigma =  extract_pow_spec(shc_a);

        real_scalar_type log_prior = 0;
        for(size_t mtpl_l =0; mtpl_l <= m_l_max; ++mtpl_l){
            real_matrix_type g_l(m_num_fields, m_num_fields);
            ps_g.get_mtpl(mtpl_l,g_l);
            real_matrix_type c_l = compute_matrix_exp(g_l);
            real_matrix_type sigma_l(m_num_fields,m_num_fields);
            ps_sigma.get_mtpl(mtpl_l,sigma_l);
            real_scalar_type fact_l(2*mtpl_l+1);
            real_scalar_type det_cl
                = compute_determinant<real_scalar_type>(c_l);
            real_matrix_type c_inv_l(m_num_fields,m_num_fields);
            compute_inverse<real_scalar_type>(c_l,c_inv_l);
            real_matrix_type c_inv_sigma_l = prod(c_inv_l,sigma_l);
            real_scalar_type trace_cl_inv_sig_l =
                compute_trace<real_scalar_type>(c_inv_sigma_l);
            log_prior += -0.5*fact_l*( std::log(det_cl) + trace_cl_inv_sig_l)
                + std::log(det_cl);
        }

        sph_hrm_coeffs<real_scalar_type> shc_a_fwd(shc_a);
        apply_win_func<real_scalar_type>(m_win_func,shc_a_fwd);
        sph_data<real_scalar_type> data_fwd(m_data.spins(),m_data.num_pixels());
        m_sh_trans.synthesise(shc_a_fwd,data_fwd);
        real_scalar_type log_lik = 0;
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
            for(std::size_t pix_i = 0; pix_i <= m_num_pixels; ++pix_i){
                real_scalar_type diff
                    = m_data(fld_i,pix_i) - data_fwd(fld_i,pix_i);
                log_lik -= diff*diff*m_prec_mat(fld_i,pix_i);
            }
        }
        log_lik *= 0.5;

        return log_prior + log_lik;
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
    blackpearl::core::sht<real_scalar_type> m_sh_trans;
    size_t m_num_fields;
    size_t m_l_max;
    size_t m_m_max;
    size_t m_num_real_coeffs;
    std::size_t m_num_pixels;
};

}}

#endif //BLACKPEARL_LOG_POST_TAYLOR_2008_HPP
