#ifndef BLACKPEARL_CORE_SHC_PS_UTILS_HPP
#define BLACKPEARL_CORE_SHC_PS_UTILS_HPP

#include <cmath>
#include <iostream>
#include <complex>
#include <random>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "pow_spec.hpp"
#include "sph_hrm_coeffs.hpp"

namespace blackpearl { namespace core {

template <typename real_scalar_type>
size_t cholesky_decompose(
    boost::numeric::ublas::matrix<real_scalar_type> const & mat_A,
    boost::numeric::ublas::matrix<real_scalar_type> &  mat_L
) {
    using namespace boost::numeric::ublas;
    typedef matrix<real_scalar_type> real_matrix_type;
    // http://www.guwi17.de/ublas/examples/
    BOOST_ASSERT( mat_A.size1() == mat_A.size2() );
    BOOST_ASSERT( mat_A.size1() ==  mat_L.size1() );
    BOOST_ASSERT( mat_A.size2() ==  mat_L.size2() );

    size_t const num_rows = mat_A.size1();
    for (size_t k=0 ; k < num_rows; k++) {
        real_scalar_type qL_kk = mat_A(k,k) - inner_prod(
            project( row( mat_L, k), range(0, k) ),
            project( row( mat_L, k), range(0, k) )
        );

        if (qL_kk <= 0) {
            return 1 + k;
        }
        else {
            real_scalar_type L_kk = std::sqrt( qL_kk );
             mat_L(k,k) =  L_kk;
            matrix_column<real_matrix_type> cLk(mat_L, k);
            project( cLk, range(k+1, num_rows) ) = (
                project( column( mat_A, k), range(k+1, num_rows) ) - prod(
                    project( mat_L, range(k+1, num_rows), range(0, k) ) ,
                    project( row(mat_L, k), range(0, k) )
                )
            ) / L_kk;
        }
    }
    return 0;
}

template<typename real_scalar_type>
pow_spec<real_scalar_type> extract_pow_spec(
    sph_hrm_coeffs<real_scalar_type> const & shc
) {
    size_t const l_max = shc.l_max();
    size_t const m_max = shc.m_max();
    size_t const num_fields = shc.num_fields();
    pow_spec<real_scalar_type> ps(num_fields,l_max);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                ps(fld_i,fld_j,mtpl_l) =
                    shc(fld_i,mtpl_l,0).real()*shc(fld_j,mtpl_l,0).real();
            }
            for(size_t mtpl_m = 1; mtpl_m <= m_max; ++mtpl_m){
                for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l){
                    ps(fld_i,fld_j,mtpl_l) +=
                        2*(
                            shc(fld_i,mtpl_l,mtpl_m).real()*
                            shc(fld_j,mtpl_l,mtpl_m).real()
                            + shc(fld_i,mtpl_l,mtpl_m).imag()*
                            shc(fld_j,mtpl_l,mtpl_m).imag()
                        );
                }
            }
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                ps(fld_i,fld_j,mtpl_l) /= real_scalar_type(2*mtpl_l+1);
            }
        }
    }
    return ps;
}

template<typename real_scalar_type,class rng_type>
sph_hrm_coeffs<real_scalar_type> create_gauss_sph_hrm_coeffs(
    pow_spec<real_scalar_type> const & ps,
    rng_type & rng
){
    using namespace boost::numeric::ublas;
    typedef std::complex<real_scalar_type>  complex_scalar_type;
    typedef matrix<real_scalar_type> real_matrix_type;
    typedef vector<complex_scalar_type> complex_vector_type;
    size_t const l_max = ps.l_max();
    size_t const m_max = l_max;
    size_t const num_fields = ps.num_fields();
    sph_hrm_coeffs<real_scalar_type> shc(num_fields,l_max,m_max);
    const real_scalar_type sqrt_half = std::sqrt(0.5);
    std::normal_distribution<real_scalar_type> norm_dist;
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
        real_matrix_type ps_mat_l(num_fields,num_fields);
        ps.get_mtpl(mtpl_l,ps_mat_l);
        real_matrix_type ps_chol_l
            = zero_matrix<real_scalar_type>(num_fields,num_fields);
        size_t info = cholesky_decompose( ps_mat_l, ps_chol_l );
        BOOST_ASSERT_MSG(
            info == 0,
            "The power spectrum is not positive definite"
        );
        if (info == 0) {
            complex_vector_type shc_lm(num_fields);
            for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
                real_scalar_type re = norm_dist(rng);
                real_scalar_type im = 0;
                shc_lm(fld_i) = complex_scalar_type(re,im);
            }
            shc_lm = prod(ps_chol_l,shc_lm);
            shc.set_mtpl(mtpl_l,0,shc_lm);
            for(size_t mtpl_m = 1; mtpl_m <= mtpl_l; ++mtpl_m) {
                for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
                    real_scalar_type re = norm_dist(rng);
                    real_scalar_type im = norm_dist(rng);
                    shc_lm(fld_i) = complex_scalar_type(re,im);
                }
                shc_lm = sqrt_half*prod(ps_chol_l,shc_lm);
                shc.set_mtpl(mtpl_l,mtpl_m,shc_lm);
            }
        }
    }
    return shc;
}

template<typename real_scalar_type>
void convert_to_real_vector(
    blackpearl::core::sph_hrm_coeffs<real_scalar_type> const & shc,
    blackpearl::core::pow_spec<real_scalar_type> const & pspec,
    boost::numeric::ublas::vector<real_scalar_type> & pos_q
){
    using namespace boost::numeric::ublas;
    using namespace blackpearl::core;
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );

    BOOST_ASSERT(shc.l_max() == pspec.l_max());
    BOOST_ASSERT(shc.num_fields() == pspec.num_fields());

    size_t num_real_alms = shc.num_real_indep_coeffs(
        shc.num_fields(),
        shc.l_max(),
        shc.m_max()
    );

    size_t num_real_cls = pspec.num_real_indep_coeffs(
        pspec.num_fields(),
        pspec.l_max()
    );

    BOOST_ASSERT(pos_q.size() == num_real_alms + num_real_cls);

    std::complex<real_scalar_type> const * p_shc_data = & shc(0,0,0);
    real_scalar_type const * p_real_shc_data =
        reinterpret_cast<real_scalar_type const *>(p_shc_data);
    real_scalar_type const * p_pspec_data = & pspec(0,0,0);
    real_scalar_type * p_pos_q = & pos_q(0);
    std::copy(
        p_real_shc_data,
        p_real_shc_data + num_real_alms,
        p_pos_q
    );
    std::copy(
        p_pspec_data,
        p_pspec_data + num_real_cls,
        p_pos_q + num_real_alms
    );
}

template<typename real_scalar_type>
void convert_to_coeffs(
    boost::numeric::ublas::vector<real_scalar_type> const & pos_q,
    blackpearl::core::sph_hrm_coeffs<real_scalar_type> & shc,
    blackpearl::core::pow_spec<real_scalar_type> & pspec
){
    using namespace boost::numeric::ublas;
    using namespace blackpearl::core;
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );

    BOOST_ASSERT(shc.l_max() == pspec.l_max());
    BOOST_ASSERT(shc.num_fields() == pspec.num_fields());

    size_t num_real_alms = shc.num_real_indep_coeffs(
        shc.num_fields(),
        shc.l_max(),
        shc.m_max()
    );

    size_t num_real_cls = pspec.num_real_indep_coeffs(
        pspec.num_fields(),
        pspec.l_max()
    );

    BOOST_ASSERT(pos_q.size() == num_real_alms + num_real_cls);

    real_scalar_type const * p_pos_q = & pos_q(0);
    std::complex<real_scalar_type> const * p_cplx_q =
        reinterpret_cast< std::complex<real_scalar_type> const *>(p_pos_q);
    std::complex<real_scalar_type> * p_shc_data = & shc(0,0,0);
    real_scalar_type * p_pspec_data = & pspec(0,0,0);

    std::copy(
        p_cplx_q,
        p_cplx_q + num_real_alms/2,
        p_shc_data
    );
    std::copy(
        p_pos_q + num_real_alms,
        p_pos_q + num_real_alms + num_real_cls,
        p_pspec_data
    );
}

}}

#endif //BLACKPEARL_CORE_SHC_PS_UTILS_HPP
