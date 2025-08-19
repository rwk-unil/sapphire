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
#include "het_info.hpp"
#include "het_info_loader.hpp"
#include "var_info.hpp"
#include "synced_bcf_reader.h"
#include "vcf.h"

int main(int argc, char**argv) {
    CLI::App app{"Binary file analysis utility app"};
    std::string filename = "-";
    app.add_option("-f,--file", filename, "Input VCF/BCF file name");
    std::string bfname = "-";
    app.add_option("-b,--binary", bfname, "Binary file name");
    bool compute_stats = false;
    app.add_flag("--stats", compute_stats, "Compute statistics");

    CLI11_PARSE(app, argc, argv);

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bfname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    std::cout << "Loading variants ..." << std::endl;
    VarInfoLoader vars(filename);
    std::cout << "Num VCF lines : " << vars.vars.size() << std::endl;

    HetInfoMemoryMap himm(bfname);

    if (compute_stats) {
        const float PP_THRESHOLD = 0.99;
        const size_t DIST_THRESHOLD = 750;
        auto low_pp_count = himm.count_all_vars_below_threshold(PP_THRESHOLD);
        auto low_pp_snps = himm.count_all_snps_below_threshold(PP_THRESHOLD, vars.vars);
        auto solvable_low_pp_snps = himm.count_all_snps_below_threshold_that_can_be_linked_within(PP_THRESHOLD, vars.vars, DIST_THRESHOLD);
        auto maybe_solvable_low_pp_snps = himm.count_all_snps_below_threshold_that_can_be_linked_within(PP_THRESHOLD, vars.vars, DIST_THRESHOLD, false);
        auto rephased_snps = himm.count_rephased_variants();

        std::cout << "For " << himm.num_samples << " samples :" << std::endl;
        std::cout << "There are " << low_pp_count << " low PP variants" << std::endl;
        std::cout << "of which " << low_pp_snps << " are SNPs : "
                  << low_pp_snps * 100.0 / low_pp_count << "%" << std::endl;
        std::cout << "of which " << solvable_low_pp_snps << " have a SNP neighbor within " << DIST_THRESHOLD << " base pairs : "
                  << solvable_low_pp_snps * 100.0 / low_pp_snps << "%" << std::endl;
        std::cout << "of which " << solvable_low_pp_snps << " have a neighbor within " << DIST_THRESHOLD << " base pairs : "
                  << maybe_solvable_low_pp_snps * 100.0 / low_pp_snps << "%" << std::endl;
        std::cout << "Number of already rephased SNPs : " << rephased_snps << std::endl;

        return 0;
    }
}
