#include <iostream>
#include <cmath>
#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"
#include "bcf_traversal.hpp"
#include "CLI11.hpp"
#include "het_info.hpp"
#include "het_info_loader.hpp"
#include "time.hpp"

constexpr size_t PLOIDY_2 = 2;

class GlobalAppOptions {
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name");
        app.add_option("-o,--output", ofname, "Output file name");
        app.add_option("-b,--binary-file", bfname, "Binary file name");
        app.add_flag("-v,--verbose", verbose, "Will show progress and other messages");
    }

    CLI::App app{"PP Update app"};
    std::string filename = "-";
    std::string ofname = "-";
    std::string bfname = "-";
    bool verbose = false;
};

GlobalAppOptions global_app_options;

class VCFLineWork {
public:
    VCFLineWork() : vcf_line_num(0) {}
    VCFLineWork(size_t vcf_line_num) : vcf_line_num(vcf_line_num) {}
    VCFLineWork(HetInfo& hi, size_t id) : vcf_line_num(hi.vcf_line) {
        updated_data[id] = hi;
    }

    size_t vcf_line_num;
    std::map<size_t, HetInfo> updated_data;
};

void insert_in_work(std::map<size_t, VCFLineWork>& work, HetInfo& hi, size_t id) {
    if (auto work_line = work.find(hi.vcf_line); work_line != work.end()) {
        work_line->second.updated_data[id] = hi;
    } else {
        work[hi.vcf_line] = VCFLineWork(hi, id);
    }
}

void fill_work_from_himm(std::map<size_t, VCFLineWork>& work, HetInfoMemoryMap& himm) {
    // For all samples
    for (size_t i = 0; i < himm.num_samples; ++i) {
        // Get Het variants
        std::vector<HetInfo> his;
        himm.fill_het_info(his, i);

        // For all het variants
        for (auto& hi : his) {
            // If it has been rephased (>1.0)
            if (!std::isnan(hi.pp) && hi.pp > 1.0) {
                // Insert in work
                insert_in_work(work, hi, i);
            }
        }
    }
}

class PPUpdateTransformer : public BcfTransformer {
public:
    PPUpdateTransformer(std::map<size_t, VCFLineWork>& work) : work(work), line_counter(0),
        pp_arr(NULL), pp_arr_size(0) {}
    virtual ~PPUpdateTransformer() {}

    virtual void handle_bcf_line() override {
        auto line = bcf_fri.line;
        auto header = hdr;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY_2) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        // If there is work to do
        if (auto work_line = work.find(line_counter); work_line != work.end()) {
            bool has_pp = false;
            int res = bcf_get_format_float(header, line, "PP", &pp_arr, &pp_arr_size);
            // There are PP values
            if (res > 0) {
                has_pp = true;
            }
            if (!has_pp) {
                std::cerr << "Cannot update PP, it is missing from line " << line_counter << std::endl;
            } else {
                if (work_line->second.vcf_line_num != line_counter) {
                    std::cerr << "Work line is different from line counter !" << std::endl;
                    std::cerr << "Something went wrong in the machinery" << std::endl;
                } else {
                    for (auto todo : work_line->second.updated_data) {
                        auto idx = todo.first;
                        auto new_pp = todo.second.pp;
                        pp_arr[idx] = new_pp;
                    }
                }
                // Update the record
                bcf_update_format_float(hdr, line, "PP", pp_arr, pp_arr_size);
            }
        }

        int ret = bcf_write1(fp, hdr, line);
        if (ret) {
            std::cerr << "Failed to write record" << std::endl;
            exit(-1);
        }

        line_counter++;
    }

protected:
    std::map<size_t, VCFLineWork>& work;
    size_t line_counter;
    float* pp_arr;
    int pp_arr_size;
};

int main(int argc, char**argv) {
    auto start_time = std::chrono::steady_clock::now();

    CLI::App& app = global_app_options.app;
    auto& filename = global_app_options.filename;
    auto& ofname = global_app_options.ofname;
    auto& bfname = global_app_options.bfname;
    CLI11_PARSE(app, argc, argv);

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (ofname.compare("-") == 0) {
        std::cerr << "Requires output filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bfname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    std::cout << "Generating workload ..." << std::endl;

    HetInfoMemoryMap himm(bfname);
    std::map<size_t, VCFLineWork> work;
    fill_work_from_himm(work, himm);

    std::cout << "Updating...\n" << std::endl;

    PPUpdateTransformer pput(work);
    pput.transform(filename, ofname);

    std::cout << "Done !" << std::endl;
    printElapsedTime(start_time, std::chrono::steady_clock::now());

    return 0;
}