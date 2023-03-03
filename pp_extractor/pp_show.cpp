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
#include "synced_bcf_reader.h"
#include "vcf.h"

class HetInfoLoader {
public:
    HetInfoLoader(std::string bfname, size_t sample) {
        std::ifstream ifs(bfname, std::ios_base::binary | std::ios_base::in);
        if (!ifs.is_open()) {
            std::cerr << "Failed to open binary file " << bfname << std::endl;
            throw "Cannot open file";
        }
        uint32_t endianness = 0;
        ifs.read(reinterpret_cast<char *>(&endianness), sizeof(endianness));
        if (endianness != 0xaabbccdd) {
            std::cerr << "File " << bfname << " has wrong endianness" << std::endl;
            throw "Wrong endianness in file";
        }
        uint32_t num_samples = 0;
        ifs.read(reinterpret_cast<char *>(&num_samples), sizeof(num_samples));
        if (num_samples == 0) {
            std::cerr << "File " << bfname << " has no samples" << std::endl;
            throw "File has no samples";
        }

        std::vector<uint64_t> offset_table(num_samples);
        for(size_t i = 0; i < offset_table.size(); ++i) {
            ifs.read(reinterpret_cast<char *>(&offset_table[i]), sizeof(uint64_t));
        }

        std::cout << "Sample " << sample << " offset : 0x" << std::hex << offset_table[sample] << std::dec << std::endl;

        ifs.seekg(offset_table[sample]);

        uint32_t mark = 0;
        ifs.read(reinterpret_cast<char *>(&mark), sizeof(mark));
        if (mark != 0xd00dc0de) {
            std::cerr << "File " << bfname << " has wrong mark" << std::endl;
            throw "Wrong mark in file";
        }

        uint32_t id = 0;
        ifs.read(reinterpret_cast<char *>(&id), sizeof(id));
        std::cout << "Id : " << id << std::endl;

        uint32_t size = 0;
        ifs.read(reinterpret_cast<char *>(&size), sizeof(size));

        his.resize(size);
        for(auto& hi : his) {
            hi.from_stream(ifs);
        }
    }

    std::vector<HetInfo> his;
};

class PrintTraversal : public BcfTraversal {
public:
    PrintTraversal(HetInfoLoader& hil) : line_counter(0), entry(0), hil_ref(hil) {}

    virtual void handle_bcf_file_reader() override {
    }

    virtual void handle_bcf_line() override {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;

        if (entry >= hil_ref.his.size()) {
            stop = true;
            return;
        }

        // This is all very slow code, but it is just for testing

        int* ac = NULL;
        if (hil_ref.his[entry].vcf_line == line_counter) {
            std::string result(bcf_hdr_id2name(header, line->rid));
            result += "\t" + std::to_string(line->pos+1);
            result += "\t" + std::string(line->d.id);
            result += "\t" + std::string(line->d.allele[0]);
            result += "\t" + std::string(line->d.allele[1]);
            result += "\t" + std::string(".");
            result += "\t" + std::string(".");
            int nac = 0;
            int ret = bcf_get_info_int32(header, line, "AC", &ac, &nac);
            result += "\tAC=" + std::to_string(*ac);
            result += "\t" + hil_ref.his[entry].to_string();

            std::cout << result << std::endl;
            entry++;
        }
        if (ac) { free(ac); }

        line_counter++;
    }
protected:
    size_t line_counter;
    size_t entry;
    HetInfoLoader& hil_ref;
};

int main(int argc, char**argv) {
    CLI::App app{"PP Show utility app"};
    std::string filename = "-";
    app.add_option("-f,--file", filename, "Input file name");
    std::string bfname = "-";
    app.add_option("-b,--binary", bfname, "Binary file name");
    std::string ofname = "-";
    app.add_option("-o,--output", ofname, "Output file name");
    size_t sample = -1;
    app.add_option("-s,--sample", sample, "Sample index");

    CLI11_PARSE(app, argc, argv);

    if (!(sample+1)) {
        std::cerr << "Please provide sample indice" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bfname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }


#if 1

    VarInfoLoader vars(filename);
    std::cout << "Num VCF lines : " << vars.vars.size() << std::endl;
    //for (auto & v : vars.vars) {
    //    std::cout << v.to_string() << std::endl;
    //}

    HetInfoMemoryMap himm(bfname);

if (0 /* Show global info */) {
    const float PP_THRESHOLD = 0.99;
    const size_t DIST_THRESHOLD = 750;
    auto low_pp_count = himm.count_all_vars_below_threshold(PP_THRESHOLD);
    auto low_pp_snps = himm.count_all_snps_below_threshold(PP_THRESHOLD, vars.vars);
    auto solvable_low_pp_snps = himm.count_all_snps_below_threshold_that_can_be_linked_within(PP_THRESHOLD, vars.vars, DIST_THRESHOLD);
    auto maybe_solvable_low_pp_snps = himm.count_all_snps_below_threshold_that_can_be_linked_within(PP_THRESHOLD, vars.vars, DIST_THRESHOLD, false);

    std::cout << "For " << himm.num_samples << " samples :" << std::endl;
    std::cout << "There are " << low_pp_count << " low PP variants" << std::endl;
    std::cout << "of which " << low_pp_snps << " are SNPs : "
              << low_pp_snps * 100.0 / low_pp_count << "%" << std::endl;
    std::cout << "of which " << solvable_low_pp_snps << " have a SNP neighbor within " << DIST_THRESHOLD << " base pairs : "
              << solvable_low_pp_snps * 100.0 / low_pp_snps << "%" << std::endl;
    std::cout << "of which " << solvable_low_pp_snps << " have a neighbor within " << DIST_THRESHOLD << " base pairs : "
              << maybe_solvable_low_pp_snps * 100.0 / low_pp_snps << "%" << std::endl;

    return 0;
}
    // Per sample info

    //himm.print_positions(sample);
    std::vector<HetInfo> v;
    himm.fill_het_info(v, sample);
    std::cout << "---" << std::endl;
    for (auto hi : v) {
        std::cout << vars.vars[hi.vcf_line].to_string() << "\t";
        std::cout << hi.to_string() << std::endl;
    }

#else
    std::fstream ifs(bfname, std::ios_base::binary | std::ios_base::in);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open binary file " << bfname << std::endl;
        exit(-1);
    }

    uint32_t endianness = 0;
    ifs.read(reinterpret_cast<char *>(&endianness), sizeof(endianness));
    if (endianness != 0xaabbccdd) {
        std::cerr << "File " << bfname << " has wrong endianness" << std::endl;
        exit(-1);
    }
    uint32_t num_samples = 0;
    ifs.read(reinterpret_cast<char *>(&num_samples), sizeof(num_samples));
    if (num_samples == 0) {
        std::cerr << "File " << bfname << " has no samples" << std::endl;
        exit(-1);
    }

    std::vector<uint64_t> offset_table(num_samples);
    for(size_t i = 0; i < offset_table.size(); ++i) {
        ifs.read(reinterpret_cast<char *>(&offset_table[i]), sizeof(uint64_t));
    }

    std::cout << "Sample " << sample << " offset : " << offset_table[sample] << std::endl;

    ifs.close();

    HetInfoLoader hil(bfname, sample);
    //for (auto& hi : hil.his) {
    //    std::cout << hi.to_string() << std::endl;
    //}

    PrintTraversal pt(hil);
    pt.traverse_no_unpack_no_destroy(filename);
    pt.destroy();

    return 0;

    HetInfoMemoryMap himm(bfname);
    //himm.print_positions(sample);
    std::vector<HetInfo> v;
    himm.fill_het_info(v, sample);
    std::cout << "---" << std::endl;
    for (auto hi : v) {
        std::cout << hi.to_string() << std::endl;
    }
#endif

    return 0;
}
