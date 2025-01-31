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

std::vector<std::string> read_first_column(const std::string& filename) {
    std::vector<std::string> firstColumn;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        throw "Error opening file";
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string firstElement;
        iss >> firstElement;

        if (!firstElement.empty()) {
            firstColumn.push_back(firstElement);
        }
    }

    file.close();
    return firstColumn;
}

std::string cram_filename(const std::string& sample_name, const std::string& cram_path, const std::string& project_id) {
    std::string cram_file(cram_path);
    cram_file += "/" + sample_name.substr(0, 2);
    cram_file += "/" + sample_name + "_" + project_id + "_0_0.cram";
    return cram_file;
}

int main(int argc, char**argv) {
    CLI::App app{"Sample list"};
    std::string filename = "-";
    app.add_option("-f,--file", filename, "Sample file name (extract with bcftools query -l <file.vcf|bcf>)");
    std::string bfname = "-";
    app.add_option("-b,--binary", bfname, "Binary file name");
    std::string cram_path;
    app.add_option("-p,--cram-path", cram_path, "CRAM files path");
    std::string project_id;
    app.add_option("-I,--project-id", project_id, "UKB Project ID");
    bool verbose = false;
    app.add_flag("-v,--verbose", verbose, "Verbose (human output)");

    CLI11_PARSE(app, argc, argv);

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bfname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    auto samples = read_first_column(filename);

    HetInfoMemoryMap himm(bfname);

    if (verbose) {
        std::cout << "There are " << samples.size() << " samples" << std::endl;
        std::cout << "There are " << himm.num_samples << " samples in the binary file" << std::endl;
    }

    for (uint32_t himm_idx = 0; himm_idx < himm.num_samples; ++himm_idx) {
        // Because the himm is subsampled we need the original index wrt sample list
        uint32_t orig_idx = himm.get_orig_idx_of_nth(himm_idx);
        if (verbose) {
            std::cout << "Sample number : " << himm_idx << " : " << samples[orig_idx] << std::endl;
        }
        std::cout << cram_filename(samples[orig_idx], cram_path, project_id) << std::endl;
        std::cout << cram_filename(samples[orig_idx], cram_path, project_id) + ".crai" << std::endl;
    }

    return 0;

}
