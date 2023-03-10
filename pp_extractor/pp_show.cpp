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

    if (sample == -1) {
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

    std::cout << "Loading variants ..." << std::endl;
    VarInfoLoader vars(filename);
    std::cout << "Num VCF lines : " << vars.vars.size() << std::endl;
    //for (auto& v : vars.vars) {
    //    std::cout << v.to_string() << std::endl;
    //}

    HetInfoMemoryMap himm(bfname);

    // Per sample info

    //himm.print_positions(sample);
    std::vector<HetInfo> v;
    himm.fill_het_info(v, sample);
    std::cout << "---" << std::endl;
    for (auto hi : v) {
        std::cout << vars.vars[hi.vcf_line].to_string() << "\t";
        std::cout << hi.to_string() << std::endl;
    }

    return 0;
}
