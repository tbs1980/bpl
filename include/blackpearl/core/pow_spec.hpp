#ifndef BLACKPEARL_CORE_POW_SPEC_HPP
#define BLACKPEARL_CORE_POW_SPEC_HPP

#include <cstddef>
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
    pow_spec(size_t const l_max) throw()
    :m_l_max(l_max){

    }
private:
    size_t m_l_max;
};

}}

#endif //BLACKPEARL_CORE_POW_SPEC_HPP
