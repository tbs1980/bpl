#define BOOST_TEST_MODULE window function
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/win_func.hpp>

template<typename real_scalar_type>
void test_win_func(){
    using namespace blackpearl::core;
    std::size_t l_max = 2048;
    std::size_t num_fields = 6;
    win_func<real_scalar_type> ps(l_max,num_fields);
}

BOOST_AUTO_TEST_CASE(win_func){
    test_win_func<float>();
    test_win_func<double>();
}
