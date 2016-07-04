#ifndef BLACKPEARL_CORE_SPH_PREC_MAT_HPP
#define BLACKPEARL_CORE_SPH_PREC_MAT_HPP

#include <cstddef>
#include <exception>
#include <sstream>
#include <string>
#include <algorithm>
#include <type_traits>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "../config.hpp"

namespace blackpearl{ namespace core {

template<typename real_scalar_type>
class sph_diag_prec_mat{
public:
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );
    typedef boost::numeric::ublas::matrix<real_scalar_type> real_matrix_type;

    sph_diag_prec_mat()
    : m_num_fields(0)
    , m_num_pixels(0)
    , m_prec_mat(0,0) {

    }

    sph_diag_prec_mat(
        size_t const num_fields,
        size_t const num_pixels
    )
    : m_num_fields(num_fields)
    , m_num_pixels(num_pixels){
        BOOST_ASSERT_MSG(
            num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            num_pixels <= BLACKPEARL_MAX_NUM_PIXELS,
            "num_pixels too big. Please modify the config.hpp and recompile."
        );
        if( num_pixels == std::size_t(0) ){
            std::stringstream msg;
            msg << "The number of pixels"
                << " in the data should be greater than zero.";
            throw std::length_error(msg.str());
        }
        if( num_fields == std::size_t(0) ){
            std::stringstream msg;
            msg << "The number of fields"
                << " in the data should be greater than zero.";
            throw std::length_error(msg.str());
        }

        m_prec_mat = real_matrix_type(num_fields,num_pixels);
    }

    inline real_scalar_type & operator()(
        size_t const field,
        size_t const pix
    ){
        BOOST_ASSERT(field < m_num_fields);
        BOOST_ASSERT(pix < m_num_pixels);
        return m_prec_mat(field,pix);
    }

    inline real_scalar_type const & operator()(
        size_t const field,
        size_t const pix
    ) const{
        BOOST_ASSERT(field < m_num_fields);
        BOOST_ASSERT(pix < m_num_pixels);
        return m_prec_mat(field,pix);
    }

    inline size_t num_fields() const {
        return m_num_fields;
    }

    inline size_t num_pixels() const{
        return m_num_pixels;
    }

    inline real_matrix_type const & data() const {
        return m_prec_mat;
    }

    inline real_matrix_type & data() {
        return m_prec_mat;
    }

private:
    size_t m_num_fields;
    size_t m_num_pixels;
    real_matrix_type m_prec_mat;
};

}}

#endif //BLACKPEARL_CORE_SPH_PREC_MAT_HPP
