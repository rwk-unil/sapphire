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

void write_to_file(uint32_t i, uint32_t size, const std::string& bin_ofname, const HetInfoMemoryMap& himm) {
    std::vector<uint32_t> ids_to_extract(size);
    std::iota(ids_to_extract.begin(), ids_to_extract.end(), i * size);
    auto nth_bin_ofname = bin_ofname + "_" + std::to_string(i);
    himm.write_sub_file(ids_to_extract, nth_bin_ofname);
}

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
    uint32_t split_size = 0;
    app.add_option("-n,--split-size", split_size, "Split in subfiles of this size, if 0 (default) use sub sample list");
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

    if (!split_size && (samples_fname.compare("-") == 0)) {
        std::cerr << "Requires samples filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (!split_size && (sub_fname.compare("-") == 0)) {
        std::cerr << "Requires sub samples filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (split_size) {
        HetInfoMemoryMap himm(bin_fname);
        auto full_chunks = himm.num_samples / split_size;
        auto last_chunk_size = himm.num_samples % split_size;

        if (verbose) {
            std::cout << "Splitting into " << (last_chunk_size ? full_chunks : full_chunks+1) << " chunks of size " << split_size << std::endl;
        }

        for (uint32_t i = 0; i < full_chunks; ++i) {
            /*
            std::vector<uint32_t> ids_to_extract(split_size);
            std::iota(ids_to_extract.begin(), ids_to_extract.end(), i * split_size);
            auto nth_bin_ofname = bin_ofname + "_" + std::to_string(i);
            himm.write_sub_file(ids_to_extract, nth_bin_ofname);
            */
            write_to_file(i, split_size, bin_ofname, himm);
            if (more) {
                std::cout << "Splitted chunk " << i << std::endl;
            }
        }
        if (last_chunk_size) {
            /*
            std::vector<uint32_t> ids_to_extract(last_chunk_size);
            std::iota(ids_to_extract.begin(), ids_to_extract.end(), full_chunks * split_size);
            auto nth_bin_ofname = bin_ofname + "_" + std::to_string(full_chunks);
            himm.write_sub_file(ids_to_extract, nth_bin_ofname);
            */
           write_to_file(full_chunks, last_chunk_size, bin_ofname, himm);
        }
    } else {
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
    }

    return 0;
}
