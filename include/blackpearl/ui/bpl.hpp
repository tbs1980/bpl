#ifndef BLACKPEARL_UI_BPL_HPP
#define BLACKPEARL_UI_BPL_HPP

#include <string>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <exception>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/tokenizer.hpp>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <arr.h>
#include <datatypes.h>
#include <blackpearl/core/sph_data.hpp>
#include <blackpearl/core/sph_prec_mat.hpp>
#include <blackpearl/core/win_func.hpp>
#include <blackpearl/utils/csv_io.hpp>
#include <blackpearl/core/sph_hrm_coeffs.hpp>
#include <blackpearl/core/pow_spec.hpp>
#include <blackpearl/core/shc_ps_utils.hpp>
#include <blackpearl/core/sht.hpp>

namespace blackpearl { namespace ui {

class bpl {
public:
    typedef double real_scalar_t;
    bpl(){

    }

    void parse_input( std::string const & file_name ) {
        using namespace blackpearl::core;
        using namespace blackpearl::utils;
        using namespace boost::property_tree;
        using namespace boost::numeric::ublas;

        ini_parser::read_ini(file_name , m_ptree);

        Healpix_Map<real_scalar_t> hpix_data;
        std::string data_file_name = m_ptree.get<std::string>("input.data");
        std::cout<< "--> Input data = " << data_file_name << std::endl;
        read_Healpix_map_from_fits(
            data_file_name,
            hpix_data
        );
        std::size_t const num_pixels_data = (std::size_t) hpix_data.Npix();
        std::cout<<"--> Number of pixels in data = "
            << num_pixels_data << std::endl;
        std::vector<std::size_t> spins(1,0);
        spins[0] = m_ptree.get<std::size_t>("input.spins");
        std::cout<<"--> Spins of data = " << spins[0] << std::endl;
        std::size_t const num_fields = spins.size();
        m_data = sph_data<real_scalar_t>(spins,num_pixels_data);
        real_scalar_t const * p_hpix_data_begin = hpix_data.Map().begin();
        real_scalar_t const * p_hpix_data_end = hpix_data.Map().end();
        real_scalar_t * p_data = &m_data(0,0);
        std::copy(
            p_hpix_data_begin,
            p_hpix_data_end,
            p_data
        );

        Healpix_Map<real_scalar_t> hpix_n_inv;
        std::string n_inv_file_name
            = m_ptree.get<std::string>("input.inv_noise");
        std::cout<< "--> Input inverse noise = "
            << n_inv_file_name << std::endl;
        read_Healpix_map_from_fits (
            n_inv_file_name,
            hpix_n_inv
        );
        std::size_t const num_pixels_n_inv = (std::size_t) hpix_n_inv.Npix();
        std::cout<<"--> Number of pixels in inverse noise = "
            << num_pixels_n_inv << std::endl;
        if( num_pixels_n_inv != num_pixels_data ){
            std::stringstream msg;
            msg << "The number of pixels"
                << " in the data does not match that from inverse noise.";
            throw std::length_error(msg.str());
        }
        m_prec_mat
            = sph_diag_prec_mat<real_scalar_t>( num_fields , num_pixels_n_inv);
        real_scalar_t const * p_hpix_n_inv_begin = hpix_n_inv.Map().begin();
        real_scalar_t const * p_hpix_n_inv_end = hpix_n_inv.Map().end();
        real_scalar_t * p_n_inv = &m_prec_mat(0,0);
        std::copy(
            p_hpix_n_inv_begin,
            p_hpix_n_inv_end,
            p_n_inv
        );

        std::string pix_win_file_name
            = m_ptree.get<std::string>("input.pix_win");
        std::cout<<"--> Input pixel window = " << pix_win_file_name<< std::endl;
        std::vector< std::vector<real_scalar_t> > p_win_mat
            = read_from_csv_file<real_scalar_t>(pix_win_file_name);
        if(p_win_mat[0].size()-1 != num_fields ){
            std::stringstream msg;
            msg << "The number of fileds in pixel window file "
                << p_win_mat[0].size()
                << " does not match data specification "
                << num_fields << ".";
            throw std::length_error(msg.str());
        }
        std::size_t const l_max = p_win_mat.size() - 1;
        std::size_t const m_max = l_max;
        m_win_func = win_func<real_scalar_t>(num_fields,l_max);
        for(std::size_t mtpl_l = 0; mtpl_l <= l_max; ++ mtpl_l){
            for(std::size_t fld_i = 0; fld_i < num_fields; ++fld_i){
                if( p_win_mat[mtpl_l][fld_i+1] < 0. or
                        p_win_mat[mtpl_l][fld_i+1] > 1.
                 ){
                    std::stringstream msg;
                    msg << "The pixel window for ell = "
                        << mtpl_l
                        << " should be between [0,1]. We got ("
                        << p_win_mat[mtpl_l][fld_i+1] << ").";
                    throw std::domain_error(msg.str());
                }
                else {
                    m_win_func(fld_i,mtpl_l) = p_win_mat[mtpl_l][fld_i+1];
                }
            }
        }

        sph_hrm_coeffs<real_scalar_t> alms(num_fields,l_max,m_max);
        sht<real_scalar_t> sh_trans(num_fields,l_max,m_max,num_pixels_data);
        sh_trans.analyse(m_data,alms);
        pow_spec<real_scalar_t> sigma =  extract_pow_spec(alms);
        std::size_t num_cls = sigma.num_real_indep_coeffs(num_fields,l_max);
        std::size_t num_alms
            = alms.num_real_indep_coeffs(num_fields,l_max,m_max);
        m_start_q0 = vector<real_scalar_t>(num_cls+num_alms);
        convert_to_real_vector<real_scalar_t>(alms, sigma, m_start_q0);
    }

    void run_sampler() {

    }

    void write_chains() {

    }

    void write_stats() {

    }
private:
    boost::property_tree::ptree m_ptree;
    blackpearl::core::sph_data<real_scalar_t> m_data;
    blackpearl::core::sph_diag_prec_mat<real_scalar_t> m_prec_mat;
    blackpearl::core::win_func<real_scalar_t> m_win_func;
    boost::numeric::ublas::vector<real_scalar_t> m_start_q0;
};

}}

#endif //BLACKPEARL_UI_BPL_HPP
