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
#include "sample_info.hpp"
#include "synced_bcf_reader.h"
#include "vcf.h"

int main(int argc, char**argv) {
    CLI::App app{"Binary file merger utility (vertical) app"};
    std::string bin_fname = "-";
    std::string bin_ofname = "-";
    app.add_option("-b,--input", bin_fname, "Binary file prefix");
    app.add_option("-o,--output", bin_ofname, "Merged binary file name (output)");
    bool verbose = false;
    //app.add_flag("-v,--verbose", verbose, "Be more verbose");
    bool more = false;
    //app.add_flag("--more", more, "Be even more verbose");
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

    uint32_t suffix_number = 0;
    std::vector<std::string> filenames;
    auto file = bin_fname + "_" + std::to_string(suffix_number);
    while (fs::exists(file)) {
        filenames.push_back(file);
        // Next file
        file = bin_fname + "_" + std::to_string(++suffix_number);
    }

    // Create a memory map for each file
    std::vector<std::unique_ptr<HetInfoMemoryMap> > himms;
    uint32_t num_samples = 0;

    /* Memory map all the files */
    for (auto& filename : filenames) {
        himms.emplace_back(std::make_unique<HetInfoMemoryMap>(filename));
        auto& himm_p = himms.back();
        if (!himm_p->integrity_check_pass()) {
            std::cerr << "File " << filename << " doesn't pass integrity checks" << std::endl;
        }
        if (!num_samples) {
            num_samples = himm_p->num_samples;
        }

        if (num_samples != himm_p->num_samples) {
            std::cerr << "File " << filename << " has " << himm_p->num_samples << " samples vs expected " << num_samples << " samples" << std::endl;
            exit(-1);
        }
    }

    std::fstream ofs(bin_ofname, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open file " << bin_ofname << std::endl;
        exit(-1);
    }

    const uint32_t endianness = 0xaabbccdd;
    // Write endianness
    ofs.write(reinterpret_cast<const char*>(&endianness), sizeof(uint32_t));
    // Write number of samples
    ofs.write(reinterpret_cast<const char*>(&num_samples), sizeof(uint32_t));
    // Write offset table
    uint64_t dummy_offset = 0xdeadc0dedeadc0de;
    auto table_seek = ofs.tellp();
    for (size_t i = 0; i < num_samples; ++i) {
        ofs.write(reinterpret_cast<const char*>(&dummy_offset), sizeof(uint64_t));
    }

    std::vector<uint64_t> offset_table(num_samples);

    /* For all samples concatenate vertically */
    for (uint32_t i = 0; i < num_samples; ++i) {
        std::vector<HetInfo> his;
        offset_table[i] = ofs.tellp();
        /* For all files memory maps */
        for (auto& h : himms) {
            auto h_his = h->get_het_info_for_nth(i);
            if (his.size()) {
                auto last_line = his.back().vcf_line;
                /* Conditionally append Het Information wrt overlap */
                for (auto hi : h_his) {
                    if (hi.vcf_line > last_line) {
                        his.push_back(hi);
                    }
                }
            } else {
                /* Append first Het Information */
                his.insert(his.end(), h_his.begin(), h_his.end());
            }
        }
        SampleBlock::write_to_stream(ofs, his, i);
    }

    // Rewrite the offset table
    ofs.seekp(table_seek);
    for (size_t i = 0; i < num_samples; ++i) {
        ofs.write(reinterpret_cast<const char*>(&offset_table[i]), sizeof(decltype(offset_table.front())));
    }

    std::cout << "Done writing file " << bin_ofname << std::endl;

    ofs.close();

    return 0;
}
