#ifndef BLACKPEARL_SPH_HRM_COEFFS_HPP
#define BLACKPEARL_SPH_HRM_COEFFS_HPP

#include <complex>
#include <type_traits>
#include <cstddef>
#include <exception>
#include <sstream>
#include <string>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "../config.hpp"

namespace blackpearl{ namespace core {

template<real_scalar_type>
class sph_hrm_coeffs
{
public:
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );
    typedef std::complex<real_scalar_type> complex_scalar_type;
    typedef boost::numeric::ublas::matrix<
        complex_scalar_type> complex_matrix_type;

    static size_t num_sph_hrm_coeffs(size_t const l_max,size_t const m_max) {
        return (m_max+size_type(1))*(l_max+size_type(2))/size_type(2)
            +(l_max+size_type(1))*(l_max-m_max);
    }

    sph_hrm_coeffs(
        std::size_t const l_max,
        std::size_t const m_max,
        std::size_t const num_fields
    ) throw()
    :m_hrm_coeffs(){

    }
private:

    complex_matrix_type m_hrm_coeffs;
};

}}

#endif //BLACKPEARL_SPH_HRM_COEFFS_HPP