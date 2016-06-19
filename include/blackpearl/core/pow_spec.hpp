#ifndef BLACKPEARL_CORE_POW_SPEC_HPP
#define BLACKPEARL_CORE_POW_SPEC_HPP

#include <cstddef>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include "../config.hpp"

namespace blackpearl { namespace core {

template<typename real_scalar_type>
class pow_spec {
public:
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );
    typedef boost::numeric::ublas::matrix<real_scalar_type> real_matrix_type;

    static size_t num_power_coeffs(size_t const l_max){
        return size_t(l_max+1);
    }

    static size_t num_real_indep_coeffs(
        size_t const num_fields,
        size_t const l_max
    ) {
        return num_power_coeffs(l_max)*num_fields*(num_fields+1)/2;
    }

    static size_t get_tri_index(
        size_t const ind_i,
        size_t const ind_j,
        size_t const size_n
    ){
        return size_n*(size_n - 1)/2
            - (size_n - ind_i)*( (size_n - ind_i)-1 )/2 + ind_j ;
    }

    pow_spec(size_t const l_max,size_t const num_fields) throw()
    :m_l_max(l_max)
    ,m_num_fields(num_fields){
        BOOST_ASSERT_MSG(
            l_max <= BLACKPEARL_MAX_LMAX,
            "l_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );
        size_t const num_rows = num_fields*(num_fields+1)/2;
        m_pow_specs = real_matrix_type(num_rows,(l_max+1));
    }

    inline real_scalar_type const & operator()(
        size_t const fld_i,
        size_t const fld_j,
        size_t const mtpl_l
    ) const {
        // http://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
        BOOST_ASSERT(fld_j >= fld_i);
        return m_pow_specs(get_tri_index(fld_i,fld_j,m_num_fields),mtpl_l);
    }

    inline real_scalar_type & operator()(
        size_t const fld_i,
        size_t const fld_j,
        size_t const mtpl_l
    ) {
        BOOST_ASSERT(fld_j >= fld_i);
        return m_pow_specs(get_tri_index(fld_i,fld_j,m_num_fields),mtpl_l);
    }

    inline void set_mtpl (
        size_t const mtpl_l,
        real_matrix_type const & c_ell
    ) {
        for(size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
            for(size_t fld_j = fld_i; fld_j < m_num_fields; ++fld_j){
                m_pow_specs(get_tri_index(fld_i,fld_j,m_num_fields), mtpl_l)
                    = c_ell(fld_i,fld_j);
            }
        }
    }

    inline void get_mtpl (
        size_t const mtpl_l,
        real_matrix_type & c_ell
    ) const {
        for(size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
            for(size_t fld_j = fld_i; fld_j < m_num_fields; ++fld_j){
                 c_ell(fld_i,fld_j) = m_pow_specs(
                     get_tri_index(fld_i,fld_j,m_num_fields),
                     mtpl_l
                 );
                 c_ell(fld_j,fld_i) = m_pow_specs(
                     get_tri_index(fld_i,fld_j,m_num_fields),
                     mtpl_l
                 );
            }
        }
    }

    inline size_t l_max() const {
        return m_l_max;
    }

    inline size_t num_fields() const {
        return m_num_fields;
    }

    inline real_matrix_type const & data() const {
        return m_pow_specs;
    }

    inline real_matrix_type & data() {
        return m_pow_specs;
    }

private:
    size_t m_l_max;
    size_t m_num_fields;
    real_matrix_type m_pow_specs;
};

}}

#endif //BLACKPEARL_CORE_POW_SPEC_HPP
