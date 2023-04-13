#include <iostream>
#include "CLI11.hpp"
#include "fifo.hpp"
#include "extractors.hpp"
#include "time.hpp"

class GlobalAppOptions {
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name");
        app.add_option("-o,--output", ofname, "Output file name");
        app.add_option("-s,--start", start, "Starting sample position");
        app.add_option("-e,--end", end, "End sample position (excluded)");
        app.add_option("-p,--progress", progress, "Number of VCF lines to show progress");
        app.add_flag("--pp-from-maf", pp_from_maf, "Genrate PP score from minor allele frequency (MAF)");
        app.add_option("--fifo-size", fifo_size, "FIFO size (number of hets extracted centered on het of interest)");
        app.add_flag("-v,--verbose", verbose, "Will show progress and other messages");
    }

    CLI::App app{"PP Extractor app"};
    std::string filename = "-";
    std::string ofname = "-";
    size_t start = 0;
    size_t end = -1;
    size_t progress = 0;
    bool pp_from_maf = false;
    size_t fifo_size = 5;
    bool verbose = false;
};

GlobalAppOptions global_app_options;

int main(int argc, char**argv) {
    auto begin_time = std::chrono::steady_clock::now();

    CLI::App& app = global_app_options.app;
    auto& start = global_app_options.start;
    auto& end = global_app_options.end;
    auto& filename = global_app_options.filename;
    auto& ofname = global_app_options.ofname;
    CLI11_PARSE(app, argc, argv);

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (ofname.compare("-") == 0) {
        std::cerr << "Requires output filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (!(global_app_options.fifo_size & 1)) {
        std::cerr << "Warning the FIFO size must be odd\n";
        global_app_options.fifo_size++;
        std::cerr << "FIFO size updated to " << global_app_options.fifo_size << std::endl;
    }

    std::cout << "Extracting...\n" << std::endl;

    PPExtractTraversal ppet(start, end, global_app_options.fifo_size,
                            global_app_options.pp_from_maf);
    ppet.set_progress(global_app_options.progress);
    ppet.traverse_no_destroy(filename);
    ppet.finalize();
    ppet.write_to_file(ofname);

    std::cout << "Done !" << std::endl;

    printElapsedTime(begin_time, std::chrono::steady_clock::now());

    return 0;
}