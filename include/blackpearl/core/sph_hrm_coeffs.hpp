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

template<class real_scalar_type>
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
        return (m_max+size_t(1))*(l_max+size_t(2))/size_t(2)
            +(l_max+size_t(1))*(l_max-m_max);
    }

    sph_hrm_coeffs(
        std::size_t const l_max,
        std::size_t const m_max,
        std::size_t const num_fields
    ) throw()
    :m_l_max(l_max)
    ,m_m_max(m_max)
    ,m_num_fields(num_fields){
        BOOST_ASSERT_MSG(
            l_max <= BLACKPEARL_MAX_LMAX,
            "l_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_max <= BLACKPEARL_MAX_LMAX,
            "m_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );

        if( m_max > l_max ){
            std::stringstream msg;
            msg << "m_max = "
                << m_max
                << " should be less than or equal to"
                << " l_max = "
                << l_max
                <<".";
            throw std::length_error(msg.str());
        }

        size_t const num_coeffs = num_sph_hrm_coeffs(m_l_max,m_max);
        std::cout<<"num_coeffs = "<<num_coeffs<<std::endl;
        m_hrm_coeffs = complex_matrix_type(
            num_coeffs,
            num_fields
        );

        m_num_shp_hrm_coeffs = num_sph_hrm_coeffs(l_max,m_max);

    }

    inline size_t num_sph_hrm_coeffs() const {
        return m_num_shp_hrm_coeffs;
    }
private:
    size_t m_l_max;
    size_t m_m_max;
    size_t m_num_fields;
    size_t m_num_shp_hrm_coeffs;
    complex_matrix_type m_hrm_coeffs;
};

}}

#endif //BLACKPEARL_SPH_HRM_COEFFS_HPP
