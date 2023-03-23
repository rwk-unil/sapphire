#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sys/mman.h>

#include "bcf_traversal.hpp"
#include "CLI11.hpp"
#include "hts.h"
#include "fifo.hpp"
#include "het_info.hpp"
#include "het_info_loader.hpp"
#include "var_info.hpp"
#include "sample_info.hpp"
#include "synced_bcf_reader.h"
#include "vcf.h"

int main(int argc, char**argv) {
    CLI::App app{"Binary file merger utility app"};
    std::string bin_fname = "-";
    std::string bin_ofname = "-";
    app.add_option("-b,--input", bin_fname, "Binary file prefix");
    app.add_option("-o,--output", bin_ofname, "Merged binary file name (output)");
    bool verbose = false;
    app.add_flag("--verbose", verbose, "Be more verbose");
    bool more = false;
    app.add_flag("--more", more, "Be even more verbose");
    verbose = verbose || more;

    CLI11_PARSE(app, argc, argv);

    if (bin_fname.compare("-") == 0) {
        std::cerr << "Requires binary filename prefix\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bin_ofname.compare("-") == 0) {
        std::cerr << "Requires output binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bin_fname.compare(bin_ofname) == 0) {
        std::cerr << "Input and Output files must have different names !" << std::endl;
        exit(-1);
    }

    uint32_t suffix_number = 0;
    std::vector<std::string> filenames;
    auto file = bin_fname + "_" + std::to_string(suffix_number);
    while (fs::exists(file)) {
        filenames.push_back(file);
        // Next file
        file = bin_fname + "_" + std::to_string(++suffix_number);
    }

    auto num_samples = HetInfoMemoryMapMerger::get_num_samples_from_filenames(filenames);

    if (verbose) {
        std::cout << "Total number of samples : " << num_samples << std::endl;
    }

    // Merge file in smaller scope (destructor closes file)
    {
        HetInfoMemoryMapMerger himmm(bin_ofname, num_samples);
        for (const auto& file : filenames) {
            himmm.merge(file);
        }
    }
    HetInfoMemoryMap himm(bin_ofname);
    if (!himm.integrity_check_pass()) {
        std::cerr << "Generated file " << bin_ofname << " has problems" << std::endl;
    }

    if (more) {
        himm.show_info();
    }

    return 0;
}
