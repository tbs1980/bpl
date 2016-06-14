#include <random>
#include <cmath>
#include <complex>
#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <sharp_lowlevel.h>
#include <sharp_geomhelpers.h>
#include <sharp_almhelpers.h>
template<typename T> struct cxxjobhelper__ {};

template<> struct cxxjobhelper__<double>
  { enum {val=SHARP_DP}; };

template<> struct cxxjobhelper__<float>
  { enum {val=0}; };


void test_sht_data_layout(){
    using namespace boost::numeric::ublas;
    size_t const num_sides = 2048;
    size_t const num_pixels = 12*num_sides*num_sides;
    size_t const num_fields = 2;
    size_t const l_max = 2*num_sides;
    size_t const m_max = l_max;
    size_t const num_alms =
        (m_max+size_t(1))*(l_max+size_t(2))/size_t(2)
        +(l_max+size_t(1))*(l_max-m_max);

    matrix<double> maps(num_fields,num_pixels);
    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> uni_real;
    for(size_t pix_i = 0; pix_i < num_pixels; ++pix_i){
        for(size_t fld_i = 0; fld_i < num_fields; ++fld_i){
            maps(fld_i,pix_i) = uni_real(rng);
        }
        //maps(1,pix_i) = 0;
        //std::cout<<pix_i<<"\t"<<maps(0,pix_i)<<"\t"<<maps(1,pix_i)<<std::endl;
    }
    matrix<std::complex<double>> alms(num_fields,num_alms);
    for(size_t cf_i = 0; cf_i < num_alms; ++cf_i){
        for(size_t fld_i = 0; fld_i < num_fields-1; ++fld_i){
            alms(fld_i,cf_i) = std::complex<double>(0,0);
        }
    }
    sharp_alm_info * m_p_alm_info;
    sharp_geom_info * m_p_geom_info;
    int const stride = 1;
    sharp_make_healpix_geom_info (
        (int) num_sides,
        stride,
        &m_p_geom_info
    );
    sharp_make_triangular_alm_info (
        (int) l_max,
        (int) m_max,
        stride,
        &m_p_alm_info
    );
    typedef const double* const_real_scalar_pntr_type;
    const_real_scalar_pntr_type * p_ptr_m
        = new const_real_scalar_pntr_type[num_fields];
    p_ptr_m[0] = & maps(0,0);
    typedef const std::complex<double>* const_complex_scalar_pntr_type;
    const_complex_scalar_pntr_type * p_ptr_a
        = new const_complex_scalar_pntr_type[num_fields];
    p_ptr_a[0] = & alms(0,0);
    for(size_t i=1;i<num_fields;++i){
        p_ptr_m[i] = p_ptr_m[i-1] + num_pixels;
        p_ptr_a[i] = p_ptr_a[i-1] + num_alms;
    }
    int const add_output = 0;
    int const sharp_nv = 0;
    int num_parallel_transforms = (int) num_fields;
    int spin = 0;
    double wall_time = 0;
    unsigned long long op_count = 0;
    sharp_execute(SHARP_MAP2ALM,
        spin,
        add_output,
        (void*) & p_ptr_a[0],
        (void*) & p_ptr_m[0],
        (const sharp_geom_info*) m_p_geom_info,
        (const sharp_alm_info*) m_p_alm_info,
        num_parallel_transforms,
        cxxjobhelper__<double>::val,
        sharp_nv,
        &wall_time,
        &op_count
    );
    // for(size_t cf_i = 0; cf_i < num_alms; ++cf_i){
    //     std::cout<<cf_i<<"\t"<<alms(0,cf_i)<<"\t"<<alms(1,cf_i)<<std::endl;
    // }
    std::cout<<"wall time = "<<wall_time<<std::endl;
    std::cout<<"op count = "<<double(op_count)/1e9<<std::endl;

    delete[] p_ptr_a;
    delete[] p_ptr_m;

    if( m_p_geom_info ) {
        sharp_destroy_geom_info(m_p_geom_info);
    }
    if( m_p_alm_info ) {
        sharp_destroy_alm_info(m_p_alm_info);
    }

}

int main(){
    test_sht_data_layout();
    return 0;
}
