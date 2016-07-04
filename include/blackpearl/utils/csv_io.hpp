#ifndef BLACKPEARL_UTILS_CSV_IO_HPP
#define BLACKPEARL_UTILS_CSV_IO_HPP

#include <cstddef>
#include <string>
#include <fstream>
#include <exception>
#include <sstream>
#include <vector>
#include <boost/assert.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace blackpearl{ namespace utils {

template<class real_scalar_t>
std::vector< std::vector<real_scalar_t> > read_from_csv_file(
    std::string const & file_name
)  {
    static_assert(
        std::is_floating_point<real_scalar_t>::value,
        "The real_scalar_t should be a floating point type"
    );
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
    typedef boost::char_separator<char> char_separator_t;

    char_separator_t sep(",");
    std::vector< std::vector<real_scalar_t> > mat;
    std::size_t num_rows(0);
    std::size_t num_cols(0);

    std::ifstream in_file;
    in_file.open( file_name.c_str(),std::ios::in );
    if( in_file.is_open() ) {
        while( !in_file.eof() ) {
            std::string line;
            std::getline(in_file,line);
            tokenizer_t tokens(line, sep);
            std::vector<real_scalar_t> row_entry;
            for(auto iter = tokens.begin(); iter != tokens.end(); ++iter ) {
                row_entry.push_back(
                    boost::lexical_cast<real_scalar_t>(*iter)
                );
            }
            if( num_rows == 0 ) {
                num_cols = row_entry.size();
            }
            else {
                if(row_entry.size() == num_cols) {
                    mat.push_back(row_entry);
                }
                else if(row_entry.size() == 0) {
                    break;
                }
                else {
                    std::stringstream msg;
                    msg << "Row  " << (num_rows + 1)
                        << " has " << row_entry.size() << " columns "
                        << "which is not equal to previous rows ("
                        << num_cols << ").";
                    throw std::runtime_error(msg.str());
                }
            }
            ++ num_rows;
        }
        in_file.close();
    }
    else {
        std::stringstream msg;
        msg << "The file "
            << file_name
            << " cannot be opened.";
        throw std::runtime_error(msg.str());
    }
    return mat;
}

}}
#endif //BLACKPEARL_UTILS_CSV_IO_HPP
