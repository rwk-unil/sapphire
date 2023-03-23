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
    CLI::App app{"Binary file splitter utility app"};
    std::string bin_fname = "-";
    std::string bin_ofname = "-";
    app.add_option("-b,--input", bin_fname, "Original binary file name");
    app.add_option("-o,--output", bin_ofname, "Subsampled binary file name (output)");
    std::string samples_fname = "-";
    app.add_option("-S,--samples", samples_fname, "Sample names (text file) as appear in full BCF");
    std::string sub_fname = "-";
    app.add_option("-l,--sample-list", sub_fname, "Unordered sub sample list (text file)");
    bool verbose = false;
    app.add_flag("--verbose", verbose, "Be more verbose");
    bool more = false;
    app.add_flag("--more", more, "Be even more verbose");
    verbose = verbose || more;

    CLI11_PARSE(app, argc, argv);

    if (bin_fname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
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

    if (samples_fname.compare("-") == 0) {
        std::cerr << "Requires samples filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (sub_fname.compare("-") == 0) {
        std::cerr << "Requires sub samples filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    SampleInfoLoader all_sil(samples_fname);
    SampleInfoLoader sub_sil(sub_fname);

    std::vector<uint32_t> idx_to_extract;

    if (verbose) {
        std::cout << "Generating the idx extract list" << std::endl;
        std::cout << sub_sil.sample_names.size() << " IDs requested" << std::endl;
    }

    // Generate the idx to extract list
    for (size_t i = 0; i < sub_sil.sample_names.size(); ++i) {
        auto result = std::find(all_sil.sample_names.begin(), all_sil.sample_names.end(),
                                sub_sil.sample_names[i]);
        if (result != all_sil.sample_names.end()) {
            if (more) {
                std::cerr << i << " : Found sample " << sub_sil.sample_names[i] << " in the original file" << std::endl;
            }
            uint32_t original_idx = result - all_sil.sample_names.begin();
            if (original_idx > all_sil.sample_names.size()) {
                std::cerr << "Something went terribly wrong" << std::endl;
                exit(-1);
            }
            idx_to_extract.push_back(original_idx);
        } else {
            if (verbose) {
                std::cerr << "Sample " << sub_sil.sample_names[i] << " is not in the original file" << std::endl;
            }
        }
    }

    if (verbose) {
        std::cout << "Sorting the idx extract list" << std::endl;
    }

    // Sort the list (it could be that the sub sample file name is shuffled)
    std::sort(idx_to_extract.begin(), idx_to_extract.end());

    // Extract the binary file

    if (verbose) {
        std::cout << "Memory mapping the binary file" << std::endl;
    }

    HetInfoMemoryMap himm(bin_fname);

    if (verbose) {
        std::cout << "Writing the sub binary file" << std::endl;
    }

    himm.write_sub_file(idx_to_extract, bin_ofname);

    if (verbose) {
        std::cout << "Extracted binary data for " << idx_to_extract.size() << " samples to file : " << bin_ofname << std::endl;
    }

    return 0;
}
