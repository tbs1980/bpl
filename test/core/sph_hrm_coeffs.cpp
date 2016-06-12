#define BOOST_TEST_MODULE shperical harmonic coefficients
#define BOOST_TEST_DYN_LINK
#include <cmath>
#include <random>
#include <vector>
#include <chrono>
#include <complex>
#include <boost/test/unit_test.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

template<typename real_scalar_type>
void test_shp_hrm_coeffs(){
    using namespace blackpearl::core;
    size_t l_max(6143);
    size_t m_max(l_max);
    size_t num_fields(6);
    sph_hrm_coeffs<real_scalar_type> sh(l_max,m_max,num_fields);
}

template<typename real_scalar_type>
void test_shp_hrm_coeffs_row(){
    using namespace boost::numeric::ublas;
    using namespace blackpearl::core;
    typedef std::complex<real_scalar_type> complex_scalar_type;
    size_t l_max(32);
    size_t m_max(l_max);
    size_t num_fields(6);
    sph_hrm_coeffs<real_scalar_type> sh(l_max,m_max,num_fields);

    vector<complex_scalar_type> v(num_fields);
    sh.get_row(0,0,v);
    sh.set_row(1,0,v);
}




BOOST_AUTO_TEST_CASE(shp_hrm_coeffs){
    test_shp_hrm_coeffs<float>();
    test_shp_hrm_coeffs<double>();
}

BOOST_AUTO_TEST_CASE(shp_hrm_coeffs_row){
    test_shp_hrm_coeffs_row<float>();
    test_shp_hrm_coeffs_row<double>();
}
