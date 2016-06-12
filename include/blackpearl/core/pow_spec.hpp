#ifndef BLACKPEARL_CORE_POW_SPEC_HPP
#define BLACKPEARL_CORE_POW_SPEC_HPP

#include <cstddef>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>

namespace blackpearl { namespace core {

template<typename real_scalar_type>
class pow_spec {
public:
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );
    typedef boost::numeric::ublas::matrix<real_scalar_type> real_matrix_type;
    pow_spec(size_t const l_max,size_t const num_fields) throw()
    :m_l_max(l_max)
    ,m_num_fields(num_fields){
        real_matrix_type cov_mat_l(num_fields,num_fields);
        m_pow_specs
            = std::vector<real_matrix_type>(l_max+1,cov_mat_l);
    }

    inline real_scalar_type const & operator()(
        size_t const mtpl_l,
        size_t const fld_i,
        size_t const fld_j
    ) const {
        return m_pow_specs[mtpl_l](fld_i,fld_j);
    }

    inline real_scalar_type & operator()(
        size_t const mtpl_l,
        size_t const fld_i,
        size_t const fld_j
    ) {
        return m_pow_specs[mtpl_l](fld_i,fld_j);
    }

    inline real_matrix_type const & operator() (
        size_t const mtpl_l
    ) const {
        return m_pow_specs[mtpl_l];
    }

    inline real_matrix_type & operator() (
        size_t const mtpl_l
    ) {
        return m_pow_specs[mtpl_l];
    }

    inline size_t l_max() const {
        return m_l_max;
    }

    inline size_t num_fields() const {
        return m_num_fields;
    }
private:
    size_t m_l_max;
    size_t m_num_fields;
    std::vector<real_matrix_type> m_pow_specs;
};

}}

#endif //BLACKPEARL_CORE_POW_SPEC_HPP
