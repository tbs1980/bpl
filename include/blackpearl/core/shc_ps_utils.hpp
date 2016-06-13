#ifndef BLACKPEARL_CORE_SHC_PS_UTILS_HPP
#define BLACKPEARL_CORE_SHC_PS_UTILS_HPP

#include <cmath>
#include <iostream>
#include <complex>
#include <random>
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
    sph_hrm_coeffs<real_scalar_type> const & shc_1,
    sph_hrm_coeffs<real_scalar_type> const & shc_2
) {
    size_t const l_max = shc_1.l_max();
    size_t const m_max = shc_1.m_max();
    size_t const num_fields = shc_1.num_fields();
    pow_spec<real_scalar_type> ps(l_max,num_fields);
    for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                ps(mtpl_l,fld_i,fld_j) =
                    shc_1(mtpl_l,0,fld_i).real()*shc_2(mtpl_l,0,fld_j).real();
            }
            for(size_t mtpl_m = 1; mtpl_m <= m_max; ++mtpl_m){
                for(size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l){
                    ps(mtpl_l,fld_i,fld_j) +=
                        2*(
                            shc_1(mtpl_l,mtpl_m,fld_i).real()*
                            shc_2(mtpl_l,mtpl_m,fld_j).real()
                            + shc_1(mtpl_l,mtpl_m,fld_i).imag()*
                            shc_2(mtpl_l,mtpl_m,fld_j).imag()
                        );
                }
            }
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                ps(mtpl_l,fld_i,fld_j) /= real_scalar_type(2*mtpl_l+1);
                if(fld_i != fld_j){
                    ps(mtpl_l,fld_j,fld_i) = ps(mtpl_l,fld_i,fld_j);
                }
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
    sph_hrm_coeffs<real_scalar_type> shc(l_max,m_max,num_fields);
    const real_scalar_type sqrt_half = std::sqrt(0.5);
    std::normal_distribution<real_scalar_type> norm_dist;
    for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
        real_matrix_type ps_mat_l = ps(mtpl_l);
        real_matrix_type ps_chol_l
            = zero_matrix<real_scalar_type>(num_fields,num_fields);
        size_t info = cholesky_decompose( ps_mat_l, ps_chol_l );
        BOOST_ASSERT_MSG(
            info == 0,
            "The power spectrum is not positive definite"
        );
        complex_vector_type shc_lm(num_fields);
        for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
            real_scalar_type re = norm_dist(rng);
            real_scalar_type im = 0;
            shc_lm(fld_i) = complex_scalar_type(re,im);
        }
        shc_lm = prod(ps_chol_l,shc_lm);
        shc.set_row(mtpl_l,0,shc_lm);
        for(size_t mtpl_m = 1; mtpl_m <= mtpl_l; ++mtpl_m) {
            for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
                real_scalar_type re = norm_dist(rng);
                real_scalar_type im = norm_dist(rng);
                shc_lm(fld_i) = complex_scalar_type(re,im);
            }
            shc_lm = sqrt_half*prod(ps_chol_l,shc_lm);
            shc.set_row(mtpl_l,mtpl_m,shc_lm);
        }
    }
    return shc;
}

}}

#endif //BLACKPEARL_CORE_SHC_PS_UTILS_HPP
