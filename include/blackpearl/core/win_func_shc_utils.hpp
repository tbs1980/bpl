#ifndef BLACKPEARL_CORE_WIN_FUC_SHC_UTILS_HPP
#define BLACKPEARL_CORE_WIN_FUC_SHC_UTILS_HPP

#include <cstddef>
#include <boost/assert.hpp>
#include "sph_hrm_coeffs.hpp"
#include "win_func.hpp"

namespace blackpearl{ namespace core {

template<typename real_scalar_type>
void apply_win_func(
    blackpearl::core::win_func<real_scalar_type> const & w_func,
    blackpearl::core::sph_hrm_coeffs<real_scalar_type> & shc
){
    BOOST_ASSERT_MSG(
        w_func.num_fields() == shc.num_fields(),
        "Window function and the alms should have the same dimensionality"
    );
    BOOST_ASSERT_MSG(
        w_func.l_max() == shc.l_max(),
        "Window function and the alms should have the same l_max"
    );
    std::size_t const l_max = shc.l_max();
    std::size_t const m_max = shc.l_max();
    for(std::size_t fld_i = 0; fld_i < w_func.num_fields(); ++fld_i){
        for(std::size_t mtpl_m = 0; mtpl_m <= m_max; ++mtpl_m) {
            for(std::size_t mtpl_l = mtpl_m; mtpl_l <= l_max; ++mtpl_l) {
                shc(fld_i,mtpl_l,mtpl_m) *= w_func(fld_i,mtpl_l);
            }
        }
    }

}

}}

#endif //BLACKPEARL_CORE_WIN_FUC_SHC_UTILS_HPP
