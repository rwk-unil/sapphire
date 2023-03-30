#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sys/mman.h>

#include "CLI11.hpp"
#include "het_info.hpp"
#include "het_info_loader.hpp"
#include "var_info.hpp"
#include "sample_info.hpp"

int main(int argc, char**argv) {
    CLI::App app{"Binary file diff utility app"};
    //std::string vcf_fname = "-";
    //app.add_option("-f,--vcf-file", vcf_fname, "Variant file name (needed for extra info)");
    std::string bin1_fname = "-";
    std::string bin2_fname = "-";
    app.add_option("-a,--binary1", bin1_fname, "First binary file name");
    app.add_option("-b,--binary2", bin2_fname, "Second binary file name");
    //std::string samples_fname = "-";
    //app.add_option("-S,--samples", samples_fname, "Sample names (text file) as appear in full BCF");
    bool extra = false;
    app.add_flag("--extra-info", extra, "Output extra information");
    //bool more = false;
    //app.add_flag("--more", more, "Output more information");

    CLI11_PARSE(app, argc, argv);

    if (bin1_fname.compare("-") == 0) {
        std::cerr << "Requires first binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bin2_fname.compare("-") == 0) {
        std::cerr << "Requires second binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    HetInfoMemoryMap himm_1(bin1_fname);
    HetInfoMemoryMap himm_2(bin2_fname);

    std::vector<uint32_t> idx1;
    std::vector<uint32_t> idx2;

    std::cout << "First file : " << bin1_fname << std::endl;
    std::cout << "Second file : " << bin2_fname << std::endl;
    std::cout << "The first file has " << himm_1.num_samples << " samples" << std::endl;
    std::cout << "The second file has " << himm_2.num_samples << " samples" << std::endl;

    for (uint32_t i = 0; i < himm_1.num_samples; ++i) {
        idx1.push_back(himm_1.get_orig_idx_of_nth(i));
    }

    for (uint32_t i = 0; i < himm_2.num_samples; ++i) {
        idx2.push_back(himm_2.get_orig_idx_of_nth(i));
    }

    std::sort(idx1.begin(), idx1.end());
    std::sort(idx2.begin(), idx2.end());
    std::vector<uint32_t> common_idx;
    std::set_intersection(idx1.begin(), idx1.end(), idx2.begin(), idx2.end(),
                          std::back_inserter(common_idx));

    std::cout << "There are " << common_idx.size() << " samples in common" << std::endl;

    auto map1 = himm_1.get_orig_idx_to_nth_map();
    auto map2 = himm_2.get_orig_idx_to_nth_map();

    size_t vector_diffs = 0;
    size_t diffs = 0;
    size_t commons = 0;
    for (auto idx : common_idx) {
        auto his1 = himm_1.get_het_info_for_nth(map1.at(idx));
        auto his2 = himm_2.get_het_info_for_nth(map2.at(idx));
        if (his1.size() != his2.size()) {
            vector_diffs++;
            continue;
        } else {
            for (size_t i = 0; i < his1.size(); ++i) {
                if (his1[i] != his2[i]) {
                    diffs++;
                    if (extra) {
                        std::cout << his1[i].to_string() << " " << his1[i].vcf_line << std::endl;
                        std::cout << his2[i].to_string() << " " << his2[i].vcf_line << std::endl;
                    }
                } else {
                    commons++;
                }
            }
        }
    }

    std::cout << "There are " << vector_diffs << " samples that don't have the same amount of variants" << std::endl;
    std::cout << "There are " << diffs << " variants that differ in common vectors" << std::endl;
    std::cout << "There are " << commons << " variants that are the same in common vectors" << std::endl;

    return 0;
}
