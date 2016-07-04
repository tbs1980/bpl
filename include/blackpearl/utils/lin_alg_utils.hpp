#ifndef BLACKPEARL_UTILS_LIN_ALG_UTILS_HPP
#define BLACKPEARL_UTILS_LIN_ALG_UTILS_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "../utils/EigenvalueDecomposition.hpp"

namespace blackpearl{ namespace utils{

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

}}

#endif //BLACKPEARL_UTILS_LIN_ALG_UTILS_HPP
