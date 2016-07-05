#ifndef BLACKPEARL_CORE_POW_SPEC_HPP
#define BLACKPEARL_CORE_POW_SPEC_HPP

#include <cstddef>
#include <vector>
#include <iostream>
#include <iomanip>
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

    inline static std::size_t num_power_coeffs(std::size_t const l_max){
        return std::size_t(l_max+1);
    }

    inline static std::size_t num_real_indep_coeffs(
        std::size_t const num_fields,
        std::size_t const l_max
    ) {
        return num_power_coeffs(l_max)*num_fields*(num_fields+1)/2;
    }

    inline static std::size_t get_tri_index(
        std::size_t const ind_i,
        std::size_t const ind_j,
        std::size_t const size_n
    ){
        return size_n*(size_n - 1)/2
            - (size_n - ind_i)*( (size_n - ind_i)-1 )/2 + ind_j ;
    }

    pow_spec(std::size_t const num_fields,std::size_t const l_max) throw()
    :m_num_fields(num_fields)
    ,m_l_max(l_max){
        BOOST_ASSERT_MSG(
            l_max <= BLACKPEARL_MAX_LMAX,
            "l_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );
        std::size_t const num_rows = num_fields*(num_fields+1)/2;
        m_pow_specs = real_matrix_type(num_rows,(l_max+1));
    }

    inline real_scalar_type const & operator()(
        std::size_t const fld_i,
        std::size_t const fld_j,
        std::size_t const mtpl_l
    ) const {
        // http://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
        BOOST_ASSERT(fld_i < m_num_fields);
        BOOST_ASSERT(fld_j < m_num_fields);
        BOOST_ASSERT(mtpl_l <= m_l_max);
        if( fld_j >= fld_i){
            return m_pow_specs(get_tri_index(fld_i,fld_j,m_num_fields),mtpl_l);
        }else {
            return m_pow_specs(get_tri_index(fld_j,fld_i,m_num_fields),mtpl_l);
        }
    }

    inline real_scalar_type & operator()(
        std::size_t const fld_i,
        std::size_t const fld_j,
        std::size_t const mtpl_l
    ) {
        BOOST_ASSERT(fld_i < m_num_fields);
        BOOST_ASSERT(fld_j < m_num_fields);
        BOOST_ASSERT(mtpl_l <= m_l_max);
        if( fld_j >= fld_i){
            return m_pow_specs(get_tri_index(fld_i,fld_j,m_num_fields),mtpl_l);
        }else {
            return m_pow_specs(get_tri_index(fld_j,fld_i,m_num_fields),mtpl_l);
        }
    }

    inline void set_mtpl (
        std::size_t const mtpl_l,
        real_matrix_type const & c_ell
    ) {
        BOOST_ASSERT(mtpl_l <= m_l_max);
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
            for(std::size_t fld_j = fld_i; fld_j < m_num_fields; ++fld_j){
                m_pow_specs(get_tri_index(fld_i,fld_j,m_num_fields), mtpl_l)
                    = c_ell(fld_i,fld_j);
            }
        }
    }

    inline void get_mtpl (
        std::size_t const mtpl_l,
        real_matrix_type & c_ell
    ) const {
        BOOST_ASSERT(mtpl_l <= m_l_max);
        for(std::size_t fld_i = 0; fld_i < m_num_fields; ++fld_i){
            for(std::size_t fld_j = fld_i; fld_j < m_num_fields; ++fld_j){
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

    inline std::size_t l_max() const {
        return m_l_max;
    }

    inline std::size_t num_fields() const {
        return m_num_fields;
    }

    inline real_matrix_type const & data() const {
        return m_pow_specs;
    }

    inline real_matrix_type & data() {
        return m_pow_specs;
    }

    void print_pow_spec() const{
        for(std::size_t mtpl_l =0; mtpl_l<=m_l_max; ++mtpl_l){
            std::cout<< "multipole = " << mtpl_l << std::endl;
            for(std::size_t fld_i =0; fld_i < m_num_fields; ++fld_i){
                for(std::size_t fld_j = 0; fld_j < m_num_fields; ++fld_j){
                    std::cout << std::scientific << std::setprecision(4)
                        << this->operator()(fld_i,fld_j,mtpl_l)<< "    ";
                }
                std::cout<< std::endl;
            }
            std::cout<< std::endl;
        }
    }

private:
    std::size_t m_num_fields;
    std::size_t m_l_max;
    real_matrix_type m_pow_specs;
};

}}

#endif //BLACKPEARL_CORE_POW_SPEC_HPP
