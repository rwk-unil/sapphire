#include <iostream>
#include "CLI11.hpp"
#include "bcf_traversal.hpp"
#include "time.hpp"

class BcfLineAdder : protected BcfTransformer {
public:
    BcfLineAdder() {}
    virtual ~BcfLineAdder() {}

    void add_lines(const std::string& ifname, const std::string& ofname) {
        transform(ifname, ofname);
    }

protected:
    void transform_header() override {
        line_num = 0;
        bcf_hdr_append(hdr, "##INFO=<ID=LINE,Number=1,Type=Integer,Description=\"Record number (line) in original BCF (0 based)\">");
    }

    void transform_record() override {
        bcf_update_info_int32(hdr, rec, "LINE", &line_num, 1);
        line_num++;
    }

    int32_t line_num;
};

class GlobalAppOptions {
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name");
        app.add_option("-o,--output", ofname, "Output file name");
        app.add_option("-p,--progress", progress, "Number of VCF lines to show progress");
        app.add_flag("-v,--verbose", verbose, "Will show progress and other messages");
        app.add_flag("--show-number", show_number, "Shows the number extracted");
    }

    CLI::App app{"PP Extractor app"};
    std::string filename = "-";
    std::string ofname = "-";
    size_t progress = 0;
    bool verbose = false;
    bool show_number = false;
};

GlobalAppOptions global_app_options;

int main(int argc, char**argv) {
    auto begin_time = std::chrono::steady_clock::now();

    CLI::App& app = global_app_options.app;
    auto& filename = global_app_options.filename;
    auto& ofname = global_app_options.ofname;
    CLI11_PARSE(app, argc, argv);

    if (filename.compare("-") == 0) {
        std::cerr << "No filename given, will read from stdin\n";
    }

    if (ofname.compare("-") == 0) {
        std::cerr << "Requires output filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    std::cout << "Processing...\n" << std::endl;

    BcfLineAdder bla;
    bla.add_lines(filename, ofname);

    std::cout << "Done !" << std::endl;

    printElapsedTime(begin_time, std::chrono::steady_clock::now());

    return 0;
}