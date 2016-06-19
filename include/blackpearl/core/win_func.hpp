#ifndef BLACKPEARL_CORE_WIN_FUNC_HPP
#define BLACKPEARL_CORE_WIN_FUNC_HPP

#include <type_traits>
#include <boost/numeric/ublas/matrix.hpp>
#include "../config.hpp"

namespace blackpearl{ namespace core {

template<typename real_scalar_type>
class win_func{
public:
    static_assert(
        std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type"
    );
    typedef boost::numeric::ublas::matrix<real_scalar_type> real_matrix_type;
    win_func(size_t const num_fields,size_t const l_max) throw()
    :m_num_fields(num_fields)
    ,m_l_max(l_max){
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            l_max <= BLACKPEARL_MAX_LMAX,
            "l_max too big. Please modify the config.hpp and recompile."
        );
        m_win_func = real_matrix_type(num_fields,l_max+1);
    }
    inline real_scalar_type const & operator()(
        size_t const fld_i,
        size_t const mtpl_l
    ) const {
        return m_win_func(fld_i, mtpl_l);
    }

    inline real_scalar_type & operator()(
        size_t const fld_i,
        size_t const mtpl_l
    ) {
        return m_win_func(fld_i, mtpl_l);
    }
    inline size_t l_max() const {
        return m_l_max;
    }

    inline size_t num_fields() const {
        return m_num_fields;
    }
private:
    size_t m_num_fields;
    size_t m_l_max;
    real_matrix_type m_win_func;
};

}}

#endif //BLACKPEARL_CORE_WIN_FUNC_HPP
