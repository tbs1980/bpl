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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
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
    typedef boost::numeric::ublas::matrix_row<complex_matrix_type>
        complex_matrix_row_type;
    typedef boost::numeric::ublas::vector<complex_scalar_type>
        complex_vector_type;

    static size_t num_sph_hrm_coeffs(size_t const l_max,size_t const m_max) {
        return (m_max+size_t(1))*(l_max+size_t(2))/size_t(2)
            +(l_max+size_t(1))*(l_max-m_max);
    }

    static size_t num_real_indep_coeffs(
        size_t const num_fields,
        size_t const l_max,
        size_t const m_max
    ) {
        return 2*num_sph_hrm_coeffs(l_max,m_max)*num_fields;
    }

    sph_hrm_coeffs(
        std::size_t const num_fields,
        std::size_t const l_max,
        std::size_t const m_max
    ) throw()
    :m_num_fields(num_fields)
    ,m_l_max(l_max)
    ,m_m_max(m_max){
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            l_max <= BLACKPEARL_MAX_LMAX,
            "l_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_max <= BLACKPEARL_MAX_LMAX,
            "m_max too big. Please modify the config.hpp and recompile."
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

        m_num_shp_hrm_coeffs = num_sph_hrm_coeffs(m_l_max,m_max);
        m_hrm_coeffs = complex_matrix_type(
            num_fields,
            m_num_shp_hrm_coeffs
        );
    }

    inline complex_scalar_type & operator() (
        size_t const field,
        size_t const mtpl_l,
        size_t const mtpl_m
    ){
        return m_hrm_coeffs(
            field,
            ((mtpl_m*(2*m_l_max+1-mtpl_m))>>1 ) + mtpl_l
        );
    }

    inline complex_scalar_type const & operator() (
        size_t const field,
        size_t const mtpl_l,
        size_t const mtpl_m
    ) const {
        return m_hrm_coeffs(
            field,
            ((mtpl_m*(2*m_l_max+1-mtpl_m))>>1 ) + mtpl_l
        );
    }

    inline size_t l_max() const {
        return m_l_max;
    }

    inline size_t m_max() const {
        return m_l_max;
    }

    inline size_t num_fields() const {
        return m_num_fields;
    }

    inline size_t num_sph_hrm_coeffs() const {
        return m_num_shp_hrm_coeffs;
    }

    inline void set_row(
        size_t const mtpl_l,
        size_t const mtpl_m,
        complex_vector_type const & shc_row
    ){
        for(size_t fld_i=0;fld_i<m_num_fields;++fld_i){
            m_hrm_coeffs(
                fld_i,
                ((mtpl_m*(2*m_l_max+1-mtpl_m))>>1 ) + mtpl_l
            ) = shc_row(fld_i);
        }
    }

    inline void get_row(
        size_t const mtpl_l,
        size_t const mtpl_m,
        complex_vector_type & shc_row
    ) const {
        for(size_t fld_i=0;fld_i<m_num_fields;++fld_i){
            shc_row(fld_i) = m_hrm_coeffs(
                fld_i,
                ((mtpl_m*(2*m_l_max+1-mtpl_m))>>1 ) + mtpl_l
            );
        }
    }

    inline complex_matrix_type const & data() const {
        return m_hrm_coeffs;
    }

    inline complex_matrix_type & data() {
        return m_hrm_coeffs;
    }
private:
    size_t m_num_fields;
    size_t m_l_max;
    size_t m_m_max;
    size_t m_num_shp_hrm_coeffs;
    complex_matrix_type m_hrm_coeffs;
};

}}

#endif //BLACKPEARL_SPH_HRM_COEFFS_HPP
