#define BOOST_TEST_MODULE power spectrum
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/pow_spec.hpp>

template<typename real_scalar_type>
void test_pow_spec(){
    using namespace blackpearl::core;
    std::size_t l_max = 2048;
    std::size_t num_fields = 6;
    pow_spec<real_scalar_type> ps(num_fields,l_max);

    for(size_t fld_i =0; fld_i < num_fields; ++fld_i){
        for(size_t fld_j = fld_i; fld_j < num_fields; ++fld_j){
            for(size_t mtpl_l = 0; mtpl_l <= l_max; ++mtpl_l){
                ps(fld_i,fld_j,mtpl_l) = 0;
                if( fld_i == fld_j ){
                    ps(fld_i,fld_j,mtpl_l) = 1;
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(shp_data){
    test_pow_spec<float>();
    test_pow_spec<double>();
}
