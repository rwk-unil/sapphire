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
#include "git_rev.h"

constexpr size_t PLOIDY_2 = 2;

class GlobalAppOptions {
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name");
        app.add_option("-o,--output", ofname, "Output file name");
        app.add_option("-b,--binary-file", bfname, "Binary file name");
        app.add_option("--main-var-vcf", main_var_vcf, "Main var VCF if input file is split VCF");
        app.add_flag("-v,--verbose", verbose, "Will show progress and other messages");
        app.add_flag("--no-pp", nopp, "Don't write/update the PP field");
    }

    CLI::App app{"PP Update app"};
    std::string filename = "-";
    std::string ofname = "-";
    std::string bfname = "-";
    std::string main_var_vcf = "";
    bool verbose = false;
    bool nopp = false;
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

        // This is for split VCF/BCFs
        if (search_line_value) /* [[unlikely]] */ {
            // Load global variant file (in scope to release it after, takes some time)
            VarInfoLoader vil(search_in_file);
            std::string contig = bcf_hdr_id2name(header, line->rid);
            uint32_t pos1 = line->pos;
            std::string ref = std::string(line->d.allele[0]);
            std::string alt = std::string(line->d.allele[1]);
            line_counter = vil.find_vcf_line(contig, pos1, ref, alt);
            if (line_counter == -1) {
                std::cerr << "Could not find VCF record in original VCF variant file " << search_in_file << std::endl;
                throw "vcf line counter error";
            }
            // Do not search next record (for performance), records should be contiguous
            search_line_value = false;
            /* To search for all records use get_vcf_line_map() or get_vcf_line_umap()
             * to generate a map once and do look-ups in the map */
        }

        // If there is work to do
        if (auto work_line = work.find(line_counter); work_line != work.end()) {
            bool has_pp = false;
            int res = bcf_get_format_float(header, line, "PP", &pp_arr, &pp_arr_size);
            // There are PP values
            if (res > 0) {
                has_pp = true;
            }
            if (!has_pp && !global_app_options.nopp) {
                std::cerr << "Cannot update PP, it is missing from line " << line_counter << std::endl;
                errors++;
            } else {
                if (work_line->second.vcf_line_num != line_counter) {
                    std::cerr << "Work line is different from line counter !" << std::endl;
                    std::cerr << "Something went wrong in the machinery" << std::endl;
                    errors++;
                } else {
                    for (auto todo : work_line->second.updated_data) {
                        auto idx = todo.first;
                        // If the PP is defined and smaller than 1.0 or bigger than 1.9 (rephased), it was a rephase target
                        if (!std::isnan(todo.second.pp) and ((todo.second.pp < 1.0) or (todo.second.pp > 1.9))) {
                            rephase_targets++;
                            if (todo.second.pp > 1.9) {
                                number_with_pir++;
                            }
                        }
                        if (!global_app_options.nopp) {
                            auto new_pp = todo.second.pp;
                            if (pp_arr[idx] != new_pp) {
                                updated_pp++;
                            }
                            pp_arr[idx] = new_pp;
                        }
                        // Also update GT
                        if (bcf_fri.gt_arr[idx * PLOIDY_2] != todo.second.a0) {
                            updated_gts++;
                        }
                        bcf_fri.gt_arr[idx * PLOIDY_2] = todo.second.a0;
                        bcf_fri.gt_arr[idx * PLOIDY_2 + 1] = todo.second.a1;
                    }
                }
                // Update the record
                if (!global_app_options.nopp) {
                    bcf_update_format_float(hdr, line, "PP", pp_arr, pp_arr_size);
                }
                bcf_update_genotypes(hdr, line, bcf_fri.gt_arr, bcf_fri.size_gt_arr);
            }
        }

        int ret = bcf_write1(fp, hdr, line);
        if (ret) {
            std::cerr << "Failed to write record" << std::endl;
            exit(-1);
        }

        line_counter++;
    }

    void set_search_line_counter(const std::string& filename) {
        search_line_value = true;
        search_in_file = filename;
        if (!fs::exists(search_in_file)) {
            std::cerr << "File : " << search_in_file << " does not exist !" << std::endl;
            throw "File does not exist";
        }
    }

    void print_stats() const {
        std::cout << "Rephase targets              : " << rephase_targets << std::endl;
        std::cout << "With phase informative reads : " << number_with_pir << std::endl;
        std::cout << "Updated PP entries           : " << updated_pp << std::endl;
        std::cout << "Rephased GTs                 : " << updated_gts << std::endl;
        std::cout << "Errors                       : " << errors << std::endl;
    }

protected:
    std::map<size_t, VCFLineWork>& work;
    size_t line_counter;
    float* pp_arr;
    int pp_arr_size;

    bool search_line_value = false;
    std::string search_in_file;

    size_t rephase_targets = 0;
    size_t number_with_pir = 0;
    size_t updated_pp = 0;
    size_t updated_gts = 0;
    size_t errors = 0;
};

int main(int argc, char**argv) {
    auto start_time = std::chrono::steady_clock::now();

    CLI::App& app = global_app_options.app;
    auto& filename = global_app_options.filename;
    auto& ofname = global_app_options.ofname;
    auto& bfname = global_app_options.bfname;
    CLI11_PARSE(app, argc, argv);

    std::cout << "PP Update git rev : " << std::hex << GIT_REVISION << std::dec << std::endl;

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

    std::cout << "Generating workload..." << std::endl;

    HetInfoMemoryMap himm(bfname);
    std::map<size_t, VCFLineWork> work;
    fill_work_from_himm(work, himm);

    std::cout << "Updating...\n" << std::endl;

    PPUpdateTransformer pput(work);

    if (global_app_options.main_var_vcf != "") {
        pput.set_search_line_counter(global_app_options.main_var_vcf);
    }

    pput.transform(filename, ofname);

    std::cout << "Done !" << std::endl;
    printElapsedTime(start_time, std::chrono::steady_clock::now());
    pput.print_stats();

    return 0;
}