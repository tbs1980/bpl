#ifndef BLACKPEARL_CORE_SHT_HPP
#define BLACKPEARL_CORE_SHT_HPP

#include <memory>
#include <cmath>
#include <complex>
#include <exception>
#include <sstream>
#include <string>
#include <vector>
#include <sharp_cxx.h>
#include "sph_data.hpp"
#include "sph_hrm_coeffs.hpp"

namespace blackpearl{ namespace core {


template<class real_scalar_type>
class sht
{
public:
    static_assert(std::is_floating_point<real_scalar_type>::value,
        "The real_scalar_type should be a floating point type");

    typedef std::complex<real_scalar_type> complex_scalar_type;

    inline static void * convert_to_void_ptr (real_scalar_type * ptr) {
        return reinterpret_cast<void *>(ptr);
    }

    inline static void * convert_to_void_ptr (real_scalar_type const * ptr) {
        return const_cast<void *>( reinterpret_cast<const void *>(ptr) );
    }

    static bool is_power_of_2(size_t const x){
        return (x != 0) && ((x & (x - 1)) == 0);
    }

    sht(
        size_t const l_max,
        size_t const m_max,
        size_t const num_pixels,
        size_t const num_fields
    ) throw()
    :m_l_max(l_max)
    ,m_m_max(m_max)
    ,m_num_pixels(num_pixels)
    ,m_num_fields(num_fields){
        if(num_pixels % size_t(12) !=  size_t(0)){
            std::stringstream msg;
            msg << "number of pixels = "
                << num_pixels
                << " should be a multiple of  12.";
            throw std::invalid_argument(msg.str());
        }
        size_t const num_sides = (size_t) std::sqrt(double(num_pixels)/double(12));
        if(is_power_of_2(num_sides) == false){
            std::stringstream msg;
            msg << "number of sides = "
                << num_sides
                << " should be a power of 2.";
            throw std::invalid_argument(msg.str());
        }

        sharp_make_healpix_geom_info (
            (int) num_sides,
            int(1),
            &m_p_geom_info
        );
        sharp_make_triangular_alm_info (
            (int) l_max,
            (int) m_max,
            int(1),
            &m_p_alm_info
        );

    }

    ~sht(){
        if( m_p_geom_info ) sharp_destroy_geom_info(m_p_geom_info);
        if( m_p_alm_info ) sharp_destroy_alm_info(m_p_alm_info);
    }

    void synthesise(
        sph_hrm_coeffs<real_scalar_type> const & alms,
        shp_data<real_scalar_type> & maps
    ) const {
        size_t const num_alms = alms.num_sph_hrm_coeffs();
        size_t const num_pixels = maps.num_pixels();
        size_t const num_fields = maps.num_fields();
        size_t const num_fields_spin_0 = maps.num_spin_zero_fields();
        size_t const num_fields_spin_2 = maps.num_spin_two_fields();

        typedef const complex_scalar_type* const_complex_scalar_pntr_type;
        const_complex_scalar_pntr_type * p_ptr_a
            = new const_complex_scalar_pntr_type[num_fields];
        p_ptr_a[0] = & alms(0,0,0);

        typedef real_scalar_type* real_scalar_pntr_type;
        real_scalar_pntr_type * p_ptr_m = new real_scalar_pntr_type[num_fields];
        p_ptr_m[0] = & maps(0,0);

        for(size_t i=1;i<num_fields;++i){
            p_ptr_a[i] = p_ptr_a[i-1] + num_alms;
            p_ptr_m[i] = p_ptr_m[i-1] + num_pixels;
        }

        int num_parallel_transforms = (int) num_fields_spin_0;
        int spin = 0;
        if(num_parallel_transforms >0)
        {
            sharp_execute(SHARP_ALM2MAP,
                spin,
                (void*) & p_ptr_a[0],
                (void*) & p_ptr_m[0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                cxxjobhelper__<real_scalar_type>::val,
                0,
                0
            );
        }
        num_parallel_transforms = (int) num_fields_spin_2/2;
        spin = 2;
        if(num_parallel_transforms > 0)
        {
            sharp_execute(SHARP_ALM2MAP,
                spin,
                (void*) & p_ptr_a[num_fields_spin_0],
                (void*) & p_ptr_m[num_fields_spin_0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                cxxjobhelper__<real_scalar_type>::val,
                0,
                0
            );
        }

        delete[] p_ptr_a;
        delete[] p_ptr_m;
    }

    void analyse(
        shp_data<real_scalar_type> const & maps,
        sph_hrm_coeffs<real_scalar_type> & alms
    ) const {
        size_t const num_alms = alms.num_sph_hrm_coeffs();
        size_t const num_pixels = maps.num_pixels();
        size_t const num_fields = maps.num_fields();
        size_t const num_fields_spin_0 = maps.num_spin_zero_fields();
        size_t const num_fields_spin_2 = maps.num_spin_two_fields();

        typedef const real_scalar_type* const_real_scalar_pntr_type;
        const_real_scalar_pntr_type * p_ptr_m
            = new const_real_scalar_pntr_type[num_fields];
        p_ptr_m[0] = & maps(0,0);

        typedef complex_scalar_type* complex_scalar_pntr_type;
        complex_scalar_pntr_type * p_ptr_a
            = new complex_scalar_pntr_type[num_fields];
        p_ptr_a[0] = & alms(0,0,0);

        for(size_t i=1;i<num_fields;++i){
            p_ptr_m[i] = p_ptr_m[i-1] + num_pixels;
            p_ptr_a[i] = p_ptr_a[i-1] + num_alms;
        }

        int num_parallel_transforms = (int) num_fields_spin_0;
        int spin = 0;
        if( num_fields_spin_0 > 0)
        {
            sharp_execute(SHARP_MAP2ALM,
                spin,
                (void*) & p_ptr_a[0],
                (void*) & p_ptr_m[0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                cxxjobhelper__<real_scalar_type>::val,
                0,
                0
            );
        }
        num_parallel_transforms = (int) num_fields_spin_2/2;
        spin = 2;
        if( num_fields_spin_2 > 0)
        {
            sharp_execute(SHARP_MAP2ALM,
                spin,
                (void*) & p_ptr_a[num_fields_spin_0],
                (void*) & p_ptr_m[num_fields_spin_0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                cxxjobhelper__<real_scalar_type>::val,
                0,
                0
            );
        }

        delete[] p_ptr_a;
        delete[] p_ptr_m;
    }

private:
    size_t m_l_max;
    size_t m_m_max;
    size_t m_num_pixels;
    size_t m_num_fields;
    sharp_alm_info * m_p_alm_info;
    sharp_geom_info * m_p_geom_info;
};

}}

#endif //BLACKPEARL_CORE_SHT_HPP
