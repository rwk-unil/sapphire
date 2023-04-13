#include <iostream>
#include <cmath>
#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"
#include "bcf_traversal.hpp"
#include "CLI11.hpp"
#include "fifo.hpp"
#include "het_info.hpp"
#include "time.hpp"

constexpr size_t PLOIDY_2 = 2;

class GlobalAppOptions {
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name");
        app.add_option("-o,--output", ofname, "Output file name");
        app.add_option("-s,--start", start, "Starting sample position");
        app.add_option("-e,--end", end, "End sample position (excluded)");
        app.add_option("-p,--progress", progress, "Number of VCF lines to show progress");
        app.add_option("--fifo-size", fifo_size, "FIFO size (number of hets extracted centered on het of interest)");
        app.add_flag("-v,--verbose", verbose, "Will show progress and other messages");
    }

    CLI::App app{"PP Extractor app"};
    std::string filename = "-";
    std::string ofname = "-";
    size_t start = 0;
    size_t end = -1;
    size_t progress = 0;
    size_t fifo_size = 5;
    bool verbose = false;
};

GlobalAppOptions global_app_options;

class PPPred {
public:
    PPPred(float pp_threshold) : pp_threshold(pp_threshold) {}

    bool operator()(HetInfo hi) const { return !(std::isnan(hi.pp)) && hi.pp < pp_threshold; }

protected:
    float pp_threshold;
};

class SampleBlock {
public:
    static void write_to_stream(std::fstream& ofs, const std::vector<HetInfo>& his, uint32_t id) {
        const uint32_t mark = 0xd00dc0de;
        const uint32_t size = his.size();
        // Write a mark
        ofs.write(reinterpret_cast<const char*>(&mark), sizeof(uint32_t));
        // Write the ID
        ofs.write(reinterpret_cast<const char*>(&id), sizeof(uint32_t));
        // Write the size
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
        // Write the Het Infos
        for (auto& hi : his) {
            hi.to_stream(ofs);
        }
    }
};

class PPExtractTraversal : public BcfTraversal {
public:
    PPExtractTraversal(size_t start_id, size_t stop_id, size_t fifo_size) :
        FIFO_SIZE(fifo_size),
        PP_THRESHOLD(0.99),
        pp_arr(NULL),
        pp_arr_size(0),
        start_id(start_id),
        stop_id(stop_id),
        line_counter(0),
        print_counter(0),
        pred(PP_THRESHOLD) {
    }


    virtual void handle_bcf_file_reader() override {
        number_of_het_sites.clear();
        number_of_low_pp_sites.clear();
        number_of_non_snp.clear();
        number_of_het_sites.resize(bcf_fri.n_samples, 0);
        number_of_low_pp_sites.resize(bcf_fri.n_samples, 0);
        number_of_snp_low_pp_sites.resize(bcf_fri.n_samples, 0);
        number_of_non_snp.resize(bcf_fri.n_samples, 0);

        /* Handle garbage input and limit to number of samples */
        if (start_id > bcf_fri.n_samples) {
            start_id = bcf_fri.n_samples;
        }
        if (stop_id > bcf_fri.n_samples) {
            stop_id = bcf_fri.n_samples;
        }
        if (stop_id < start_id) {
            stop_id = start_id;
        }
        std::cout << "Start ID : " << start_id << " Stop ID : " << stop_id << std::endl;
        fifos.resize(stop_id-start_id, GenericKeepFifo<HetInfo, PPPred>(FIFO_SIZE, PPPred(PP_THRESHOLD)));
    }

    virtual void handle_bcf_line() override {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;

        bool has_pp = false;
        int res = bcf_get_format_float(header, line, "PP", &pp_arr, &pp_arr_size);
        // There are PP values
        if (res > 0) {
            has_pp = true;
        }

        bool non_snp = false;
        if (strlen(line->d.allele[0]) > 1 || strlen(line->d.allele[1]) > 1) {
            non_snp = true;
        }

        // Extract heterozygous sites and PP
        for (size_t i = start_id; i < stop_id; ++i) {
            int encoded_a0 = bcf_fri.gt_arr[i*PLOIDY_2];
            int encoded_a1 = bcf_fri.gt_arr[i*PLOIDY_2+1];
            int a0 = bcf_gt_allele(encoded_a0);
            int a1 = bcf_gt_allele(encoded_a1);

            if (a0 != a1) {
                number_of_het_sites[i]++;
                float pp = NAN; // Assume perfect phasing if there is no PP field
                if (has_pp) {
                    pp = pp_arr[i];
                }

                HetInfo hi(line_counter, encoded_a0, encoded_a1, pp);

                if (pred(hi)) {
                    number_of_low_pp_sites[i]++;
                    if (!non_snp) {
                        number_of_snp_low_pp_sites[i]++;
                    }
                }
                if (non_snp) {
                    number_of_non_snp[i]++;
                }

                fifos[i-start_id].insert(hi);
            }
        }

        line_counter++;
        if (global_app_options.progress) {
            if (++print_counter == global_app_options.progress) {
                print_counter = 0;
                printf("\033[A\033[2K");
                std::cout << "Handled " << bcf_fri.line_num << " VCF entries (lines)" << std::endl;
            }
        }
    }

    void finalize() {
        // Finalize FIFOs
        for (auto& f : fifos) {
            f.finalize();
        }
    }

    void write_to_file(std::string filename) {
        std::fstream ofs(filename, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        if (!ofs.is_open()) {
            std::cerr << "Cannot open file " << filename << std::endl;
        }

        const uint32_t endianness = 0xaabbccdd;
        // Write endianness
        ofs.write(reinterpret_cast<const char*>(&endianness), sizeof(uint32_t));
        // Write number of samples
        uint32_t num_samples = stop_id-start_id;
        ofs.write(reinterpret_cast<const char*>(&num_samples), sizeof(uint32_t));
        // Write offset table
        uint64_t dummy_offset = 0xdeadc0dedeadc0de;
        auto table_seek = ofs.tellp();
        for (size_t i = start_id; i < stop_id; ++i) {
            ofs.write(reinterpret_cast<const char*>(&dummy_offset), sizeof(uint64_t));
        }

        // Write the data for all the samples
        std::vector<uint64_t> offset_table(stop_id-start_id);
        for (size_t i = start_id; i < stop_id; ++i) {
            const size_t idx = i-start_id;
            offset_table[idx] = ofs.tellp();
            SampleBlock::write_to_stream(ofs, fifos[idx].get_kept_items_ref(), i);
        }

        // Rewrite the offset table
        ofs.seekp(table_seek);
        for (size_t i = start_id; i < stop_id; ++i) {
            const size_t idx = i-start_id;
            ofs.write(reinterpret_cast<const char*>(&offset_table[idx]), sizeof(decltype(offset_table.front())));
        }

        std::cout << "Done writing file " << filename << std::endl;

        ofs.close();
    }

    const size_t FIFO_SIZE;
    const float PP_THRESHOLD;
    float *pp_arr;
    int pp_arr_size;
    size_t start_id;
    size_t stop_id;
    size_t line_counter;
    size_t print_counter;
    const PPPred pred;
    std::vector<uint32_t> number_of_het_sites;
    std::vector<uint32_t> number_of_low_pp_sites;
    std::vector<uint32_t> number_of_snp_low_pp_sites;
    std::vector<uint32_t> number_of_non_snp;
    std::vector<GenericKeepFifo<HetInfo, PPPred> > fifos;
};

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

    PPExtractTraversal ppet(start, end, global_app_options.fifo_size);
    ppet.traverse_no_destroy(filename);
    ppet.finalize();
    ppet.write_to_file(ofname);

    std::cout << "Done !" << std::endl;

    printElapsedTime(begin_time, std::chrono::steady_clock::now());

    return 0;
}