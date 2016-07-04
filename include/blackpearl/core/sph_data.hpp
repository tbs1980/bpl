#ifndef BLACKPEARL_CORE_SPH_DATA_HPP
#define BLACKPEARL_CORE_SPH_DATA_HPP

#include <cstddef>
#include <exception>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "../config.hpp"

namespace blackpearl { namespace core {

template<class real_scalar_type>
class sph_data{
public:
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );

    typedef boost::numeric::ublas::matrix<real_scalar_type> real_matrix_type;

    sph_data()
    : m_spins(0)
    , m_num_fields(0)
    , m_num_pixels(0)
    , m_data(0,0) {

    }

    sph_data(
        std::vector<std::size_t> spins,
        std::size_t const num_pixels
    )
    : m_spins(spins)
    , m_num_fields( spins.size() )
    , m_num_pixels(num_pixels) {
        BOOST_ASSERT_MSG(
            num_pixels <= BLACKPEARL_MAX_NUM_PIXELS,
            "num_pixels too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "m_num_fields too big. Please modify the config.hpp and recompile."
        );

        if( num_pixels == std::size_t(0) ){
            std::stringstream msg;
            msg << "The number of pixels"
                << " in the data should be greater than zero.";
            throw std::length_error(msg.str());
        }
        if( m_num_fields == std::size_t(0) ){
            std::stringstream msg;
            msg << "The number of fields"
                << " in the data should be greater than zero.";
            throw std::length_error(msg.str());
        }

        for(std::size_t i=0;i<m_num_fields;++i){
            if( m_spins[i] != 0 and m_spins[i] !=2 ){
                std::stringstream msg;
                msg << "input spins should be either zero or two.";
                throw std::length_error(msg.str());
            }
        }
        std::sort(spins.begin(),spins.end());
        for(std::size_t i=0;i<m_num_fields;++i){
            if( spins[i] != m_spins[i] ){
                std::stringstream msg;
                msg << "input spins should be sorted in ascending order.";
                throw std::length_error(msg.str());
            }
        }
        std::size_t nstf(0);
        for(std::size_t i=0;i<m_num_fields;++i){
            if( m_spins[i] == 2){
                ++nstf;
            }
        }
        if( nstf % size_t(2) != size_t(0) ){
            std::stringstream msg;
            msg << "the number of spin-2 fields should be an even number.";
            throw std::length_error(msg.str());
        }

        m_num_spin_two_fields = nstf;
        m_num_spin_zero_fields = m_num_fields - nstf;

        m_data = real_matrix_type(m_num_fields,m_num_pixels);
    }

    inline real_scalar_type & operator()(
        size_t const field,
        size_t const pix
    ){
        BOOST_ASSERT(field < m_num_fields);
        BOOST_ASSERT(pix < m_num_pixels);
        return m_data(field,pix);
    }

    inline real_scalar_type const & operator()(
        size_t const field,
        size_t const pix
    ) const{
        BOOST_ASSERT(field < m_num_fields);
        BOOST_ASSERT(pix < m_num_pixels);
        return m_data(field,pix);
    }

    inline size_t num_fields() const {
        return m_num_fields;
    }

    inline size_t num_spin_zero_fields() const {
        return m_num_spin_zero_fields;
    }

    inline size_t num_spin_two_fields() const {
        return m_num_spin_two_fields;
    }

    inline size_t num_pixels() const{
        return m_num_pixels;
    }

    inline real_matrix_type const & data() const {
        return m_data;
    }

    inline real_matrix_type & data() {
        return m_data;
    }

    inline std::vector<std::size_t> const & spins() const {
        return m_spins;
    }

private:
    std::vector<std::size_t> m_spins;
    std::size_t m_num_fields;
    std::size_t m_num_pixels;
    real_matrix_type m_data;
    std::size_t m_num_spin_zero_fields;
    std::size_t m_num_spin_two_fields;
};

}}

#endif //BLACKPEARL_CORE_SPH_DATA_HPP
