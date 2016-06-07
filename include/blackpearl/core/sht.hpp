#ifndef BLACKPEARL_CORE_SHT_HPP
#define BLACKPEARL_CORE_SHT_HPP

#include <memory>
#include <cmath>
#include <exception>
#include <sstream>
#include <string>
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
        shp_data<real_scalar_type> & map
    ){

    }

    void analyse(
        shp_data<real_scalar_type> const & map,
        sph_hrm_coeffs<real_scalar_type> & alms
    ){

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