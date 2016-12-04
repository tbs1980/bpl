#ifndef BLACKPEARL_CORE_SHT_HPP
#define BLACKPEARL_CORE_SHT_HPP

#include <memory>
#include <cmath>
#include <complex>
#include <exception>
#include <sstream>
#include <string>
#include <vector>
#include <sharp_lowlevel.h>
#include <sharp_geomhelpers.h>
#include <sharp_almhelpers.h>

#include "sph_data.hpp"
#include "sph_hrm_coeffs.hpp"

namespace blackpearl{ namespace core {

template<typename T> struct cxxjobhelper__ {};

template<> struct cxxjobhelper__<double>
  { enum {val=SHARP_DP}; };

template<> struct cxxjobhelper__<float>
  { enum {val=100}; };

template<class real_scalar_type>
class sht {
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

    static bool is_power_of_2(std::size_t const x){
        return (x != 0) && ((x & (x - 1)) == 0);
    }

    sht(
        std::size_t const num_fields,
        std::size_t const l_max,
        std::size_t const m_max,
        std::size_t const num_pixels
    ) throw()
    :m_num_fields(num_fields)
    ,m_l_max(l_max)
    ,m_m_max(m_max)
    ,m_num_pixels(num_pixels)
    {
        BOOST_ASSERT_MSG(
            m_num_fields <= BLACKPEARL_MAX_NUM_FIELDS,
            "num_fields too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            l_max <= BLACKPEARL_MAX_LMAX,
            "l_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            m_max <= BLACKPEARL_MAX_LMAX,
            "m_max too big. Please modify the config.hpp and recompile."
        );
        BOOST_ASSERT_MSG(
            num_pixels <= BLACKPEARL_MAX_NUM_PIXELS,
            "num_pixels too big. Please modify the config.hpp and recompile."
        );

        if(num_pixels % std::size_t(12) !=  std::size_t(0)){
            std::stringstream msg;
            msg << "number of pixels = "
                << num_pixels
                << " should be a multiple of  12.";
            throw std::invalid_argument(msg.str());
        }
        std::size_t const num_sides = (std::size_t) std::sqrt(
            real_scalar_type(num_pixels)/real_scalar_type(12)
        );
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
        sph_data<real_scalar_type> & maps
    ) const {
        std::size_t const num_alms = alms.num_sph_hrm_coeffs();
        std::size_t const num_pixels = maps.num_pixels();
        std::size_t const num_fields = maps.num_fields();
        std::size_t const num_fields_s_0 = maps.num_spin_zero_fields();
        std::size_t const num_fields_s_2 = maps.num_spin_two_fields();

        typedef const complex_scalar_type* const_complex_scalar_pntr_type;
        const_complex_scalar_pntr_type * p_ptr_a
            = new const_complex_scalar_pntr_type[num_fields];
        p_ptr_a[0] = & alms(0,0,0);

        typedef real_scalar_type* real_scalar_pntr_type;
        real_scalar_pntr_type * p_ptr_m = new real_scalar_pntr_type[num_fields];
        p_ptr_m[0] = & maps(0,0);

        for(std::size_t i=1;i<num_fields;++i){
            p_ptr_a[i] = p_ptr_a[i-1] + num_alms;
            p_ptr_m[i] = p_ptr_m[i-1] + num_pixels;
        }
        // int const add_output = 0;
        // int const sharp_nv = 0;
        int num_parallel_transforms = (int) num_fields_s_0;
        int spin = 0;
        int flags = (int) cxxjobhelper__<real_scalar_type>::val;
        if(num_parallel_transforms >0)
        {
            sharp_execute(SHARP_ALM2MAP,
                spin,
                /*add_output,*/
                (void*) & p_ptr_a[0],
                (void*) & p_ptr_m[0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                flags,
                /*sharp_nv,*/
                NULL,
                NULL
            );
        }
        num_parallel_transforms = (int) num_fields_s_2/2;
        spin = 2;
        if(num_parallel_transforms > 0)
        {
            sharp_execute(SHARP_ALM2MAP,
                spin,
                /*add_output,*/
                (void*) & p_ptr_a[num_fields_s_0],
                (void*) & p_ptr_m[num_fields_s_0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                flags,
                /*sharp_nv,*/
                NULL,
                NULL
            );
        }

        delete[] p_ptr_a;
        delete[] p_ptr_m;
    }

    void analyse(
        sph_data<real_scalar_type> const & maps,
        sph_hrm_coeffs<real_scalar_type> & alms
    ) const {
        std::size_t const num_alms = alms.num_sph_hrm_coeffs();
        std::size_t const num_pixels = maps.num_pixels();
        std::size_t const num_fields = maps.num_fields();
        std::size_t const num_fields_s_0 = maps.num_spin_zero_fields();
        std::size_t const num_fields_s_2 = maps.num_spin_two_fields();
        BOOST_ASSERT(num_pixels == m_num_pixels);
        BOOST_ASSERT(alms.l_max() == m_l_max);
        BOOST_ASSERT(alms.m_max() == m_m_max);

        typedef const real_scalar_type* const_real_scalar_pntr_type;
        const_real_scalar_pntr_type * p_ptr_m
            = new const_real_scalar_pntr_type[num_fields];
        p_ptr_m[0] = & maps(0,0);

        typedef complex_scalar_type* complex_scalar_pntr_type;
        complex_scalar_pntr_type * p_ptr_a
            = new complex_scalar_pntr_type[num_fields];
        p_ptr_a[0] = & alms(0,0,0);

        for(std::size_t i=1;i<num_fields;++i){
            p_ptr_m[i] = p_ptr_m[i-1] + num_pixels;
            p_ptr_a[i] = p_ptr_a[i-1] + num_alms;
        }
        // int const add_output = 0;
        // int const sharp_nv = 0;
        int num_parallel_transforms = (int) num_fields_s_0;
        int spin = 0;
        int flags = (int) cxxjobhelper__<real_scalar_type>::val;
        if( num_fields_s_0 > 0)
        {
            sharp_execute(SHARP_MAP2ALM,
                spin,
                /*add_output,*/
                (void*) & p_ptr_a[0],
                (void*) & p_ptr_m[0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                flags,
                /*sharp_nv,*/
                NULL,
                NULL
            );
        }
        num_parallel_transforms = (int) num_fields_s_2/2;
        spin = 2;
        if( num_fields_s_2 > 0)
        {
            sharp_execute(SHARP_MAP2ALM,
                spin,
                /*add_output,*/
                (void*) & p_ptr_a[num_fields_s_0],
                (void*) & p_ptr_m[num_fields_s_0],
                (const sharp_geom_info*) m_p_geom_info,
                (const sharp_alm_info*) m_p_alm_info,
                num_parallel_transforms,
                flags,
                /*sharp_nv,*/
                NULL,
                NULL
            );
        }

        delete[] p_ptr_a;
        delete[] p_ptr_m;
    }

private:
    std::size_t m_num_fields;
    std::size_t m_l_max;
    std::size_t m_m_max;
    std::size_t m_num_pixels;
    sharp_alm_info * m_p_alm_info;
    sharp_geom_info * m_p_geom_info;
};

}}

#endif //BLACKPEARL_CORE_SHT_HPP
