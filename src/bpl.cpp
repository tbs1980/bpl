#include <iostream>
#include <boost/program_options.hpp>
#include <blackpearl/ui/bpl.hpp>

int main(int argc, char** argv){
    using namespace boost::program_options;
    try {
        options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("input-file,f", value< std::string >(), "input file");
        positional_options_description pos_opt;
        pos_opt.add("input-file", -1);
        variables_map var_map;
        store(
            command_line_parser(
                argc, argv
            ).options(desc).positional(pos_opt).run(),
            var_map
        );

        if (var_map.count("help")) {
            std::cout << desc << "\n";
            return 0;
        }
        else if (var_map.count("input-file")) {
            using namespace blackpearl::ui;
            bpl bpl_ui;
            bpl_ui.parse_input( var_map["input-file"].as< std::string >() );
            bpl_ui.run_sampler();
        }
        else {
            std::cout << ":( No arguments passed. "
                << "Please type Cat2Map --help" << std::endl;
        }
    }
    catch(std::exception& except) {
        std::cerr << "==>> Stopping execution: "<< except.what() << std::endl;
    }
    return 0;
}
