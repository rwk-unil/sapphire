#include <iostream>
#include <cmath>
#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"
#include "../include/bcf_traversal.hpp"
#include "../include/CLI11.hpp"
#include "include/fifo.hpp"
#include "include/het_info.hpp"

constexpr size_t PLOIDY_2 = 2;

class GlobalAppOptions
{
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", filename, "Input file name");
        app.add_option("-o,--output", ofname, "Output file name");
        app.add_option("-s,--start", start, "Starting sample position");
        app.add_option("-e,--end", end, "End sample position (excluded)");
        app.add_flag("-v,--verbose", verbose, "Will show progress and other messages");
    }

    CLI::App app{"PP Extractor app"};
    std::string filename = "-";
    std::string ofname = "-";
    size_t start = 0;
    size_t end = -1;
    bool verbose = false;
};

GlobalAppOptions global_app_options;

class Het_Fifo {
public:
    Het_Fifo(const size_t size, const float pp_threshold) :
        size(size),
        mid(size/2),
        pp_threshold(pp_threshold) {

    }

    void insert(size_t entry, float pp_e) {
        entries.push_back(entry);
        pp.push_back(pp_e);
        kept.push_back(false);

        if (entries.size() == size) {
            // When we reach size we need to check for small pp in front
            for (size_t i = 0; i <= mid; ++i) {
                if (pp[i] < pp_threshold) {
                    keep();
                    break;
                }
            }
        } else if (entries.size() > size) {
            // When size is passed, pop front
            entries.pop_front();
            pp.pop_front();
            kept.pop_front();

            // If small PP
            if (pp[mid] < pp_threshold) {
                // Keep the information
                keep();
            }
        }
    }

    void finalize() {
        // Search for small PP at the end
        for (size_t i = ((entries.size() < size) ? 0 : (mid + 1)); i < entries.size(); ++i) {
            if (pp[i] < pp_threshold) {
                keep();
                break;
            }
        }
    }

    const std::vector<uint32_t>& get_pos_ref() const {
        return positions;
    }

private:
    void keep() {
        for (size_t i = 0; i < entries.size(); ++i) {
            if (!kept[i]) {
                positions.push_back(entries[i]);
                kept[i] = true;
            }
        }
    }

protected:
    const size_t size;
    const size_t mid;
    const float pp_threshold;
    /// @note uint32_t to save space (a whole human chromosome will not have more positions than uint32_t can represent)
    std::deque<uint32_t> entries;
    std::deque<float> pp;
    std::deque<bool> kept;
    std::vector<uint32_t> positions;
};

class PPPred {
public:
    PPPred(float pp_threshold) : pp_threshold(pp_threshold) {}

    bool operator()(HetInfo hi) const { return !(std::isnan(hi.pp)) && hi.pp < pp_threshold; }

protected:
    const float pp_threshold;
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
    PPExtractTraversal(size_t start_id, size_t stop_id) :
        FIFO_SIZE(5),
        PP_THRESHOLD(0.99),
        pp_arr(NULL),
        pp_arr_size(0),
        start_id(start_id),
        stop_id(stop_id),
        line_counter(0),
        print_counter(0),
        pred(PP_THRESHOLD),
        fifos(stop_id-start_id, GenericKeepFifo<HetInfo, PPPred>(FIFO_SIZE, pred)) {
    }


    virtual void handle_bcf_file_reader() override {
        bool need_resize = false;

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
            need_resize = true;
        }
        if (stop_id > bcf_fri.n_samples) {
            stop_id = bcf_fri.n_samples;
            need_resize = true;
        }
        if (stop_id < start_id) {
            stop_id = start_id;
            need_resize = true;
        }

        if (need_resize) {
            fifos.resize(stop_id-start_id, GenericKeepFifo<HetInfo, PPPred>(FIFO_SIZE, PPPred(PP_THRESHOLD)));
        }
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
        if (global_app_options.verbose) {
            if (++print_counter == 1000) {
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

class ScanTraversal : public BcfTraversal {
public:
    ScanTraversal(size_t start_id, size_t stop_id) :
        FIFO_SIZE(5),
        PP_THRESHOLD(0.99),
        pp_arr(NULL),
        pp_arr_size(0),
        start_id(start_id),
        stop_id(stop_id),
        line_counter(0),
        print_counter(0) {
    }

    virtual ~ScanTraversal() {
        if (pp_arr) {
            free(pp_arr);
        }
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

        fifos.resize(stop_id-start_id, Het_Fifo(FIFO_SIZE, PP_THRESHOLD));
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
            int a0 = bcf_gt_allele(bcf_fri.gt_arr[i*PLOIDY_2]);
            int a1 = bcf_gt_allele(bcf_fri.gt_arr[i*PLOIDY_2+1]);

            if (a0 != a1) {
                number_of_het_sites[i]++;
                float pp = 1.0; // Assume perfect phasing if there is no PP field
                if (has_pp) {
                    pp = pp_arr[i];
                }
                if (pp < PP_THRESHOLD) {
                    number_of_low_pp_sites[i]++;
                    if (!non_snp) {
                        number_of_snp_low_pp_sites[i]++;
                    }
                }
                if (non_snp) {
                    number_of_non_snp[i]++;
                }

                fifos[i-start_id].insert(line_counter, pp);
            }
        }

        line_counter++;
        if (++print_counter == 1000) {
            print_counter = 0;
            printf("\033[A\033[2K");
            std::cout << "Handled " << bcf_fri.line_num << " VCF entries (lines)" << std::endl;
        }
    }

    void finalize() {
        // Finalize FIFOs
        for (auto & f : fifos) {
            f.finalize();
        }
        // Generate extraction map
        for (size_t i = start_id; i < stop_id; ++i) {
            auto& pos = fifos[i-start_id].get_pos_ref();
            for (const auto& p : pos) {
                extractions[p].push_back(i);
            }
        }
    }

    const std::map<uint32_t, std::vector<uint32_t> >& get_extractions_ref() const {
        return extractions;
    }

protected:
    const size_t FIFO_SIZE;
    const float PP_THRESHOLD;
    float *pp_arr;
    int pp_arr_size;
    size_t start_id;
    size_t stop_id;
    size_t line_counter;
    size_t print_counter;
    std::vector<uint32_t> number_of_het_sites;
    std::vector<uint32_t> number_of_low_pp_sites;
    std::vector<uint32_t> number_of_snp_low_pp_sites;
    std::vector<uint32_t> number_of_non_snp;
    std::vector<Het_Fifo> fifos;
    std::map<uint32_t, std::vector<uint32_t> > extractions;
};

class ExtractTraversal : public BcfTraversal {
public:
    ExtractTraversal(const std::map<uint32_t, std::vector<uint32_t> >& extractions) :
        line_counter(0),
        print_counter(0),
        pp_arr(NULL),
        pp_arr_size(0),
        extractions(extractions) {
    }

    virtual ~ExtractTraversal() {
        if (pp_arr) {
            free(pp_arr);
        }
    }

    virtual void handle_bcf_file_reader() override {
        fp.clear();
        fp.resize(bcf_fri.n_samples, NULL);
        hdr.clear();
        hdr.resize(bcf_fri.n_samples, NULL);
    }

    virtual void handle_bcf_line() override {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;

        // If there is something to extract
        if (extractions.find(line_counter) != extractions.end()) {

            // Unpack the line and get genotypes
            bcf_unpack(bcf_fri.line, BCF_UN_STR);
            bcf_fri.ngt = bcf_get_genotypes(header, line, &(bcf_fri.gt_arr), &(bcf_fri.size_gt_arr));
            line_max_ploidy = bcf_fri.ngt / bcf_fri.n_samples;

            bool has_pp = false;
            int res = bcf_get_format_float(header, line, "PP", &pp_arr, &pp_arr_size);
            // There are PP values
            if (res > 0) {
                has_pp = true;
            }

            for (const auto& id : extractions.at(line_counter)) {
                // If the file is not yet opened, open it
                if (fp[id] == NULL) {
                    bcf_hdr_t *original_hdr = header;

                    /// @todo This will crash with segfault if file already exists... Not sure ?

                    std::string ofname = std::string("flanked_low_pp") + std::string(original_hdr->samples[id]) + ".bcf";
                    fp[id] = hts_open(ofname.c_str(), "wz");

                    if (!fp[id]) {
                        std::cerr << "Cannot create BCF file " << ofname << std::endl;
                        throw "Cannot create BCF file";
                    }

                    if (hdr[id] == NULL) {
                        char* const samples[1] = {original_hdr->samples[id]};
                        int imap[1];
                        hdr[id] = bcf_hdr_subset(original_hdr, 1, samples, imap);
                        if (!hdr[id]) {
                            std::cerr << "Cannot duplicate VCF header" << std::endl;
                            throw "Cannot duplicate VCF header";
                        }
                    }

                    // Write the header
                    int ret = bcf_hdr_write(fp[id], hdr[id]);
                    if (ret < 0) {
                        std::cerr << "Failed to write header to file " << ofname << std::endl;
                        throw "Failed to remove samples";
                    }
                }

                // Write record
                bcf1_t *rec = NULL;
                rec = bcf_dup(bcf_fri.line);
                bcf_unpack(rec, BCF_UN_ALL);
                rec->n_sample = 1;

                if (has_pp) {
                    int ret = bcf_update_format_float(hdr[id], rec, "PP", &pp_arr[id] , 1);
                    if (ret) {
                        std::cerr << "Could not update PP in record for id " << id << std::endl;
                        throw "Could not update PP";
                    }
                }
                int ret = bcf_update_genotypes(hdr[id], rec, &bcf_fri.gt_arr[id], 1 * PLOIDY_2);
                if (ret) {
                    std::cerr << "Could not update genotypes in record for id " << id << std::endl;
                    throw "Could not update genotypes";
                }
                ret = bcf_write(fp[id], hdr[id], rec);
                if (ret) {
                    std::cerr << "Could not write record for id " << id << std::endl;
                    throw "Failed to write record";
                }
                bcf_destroy(rec);
            }
        }

        line_counter++;
        if (++print_counter == 1000) {
            print_counter = 0;
            printf("\033[A\033[2K");
            std::cout << "Handled " << bcf_fri.line_num << " VCF entries (lines)" << std::endl;
        }
    }

    void finalize() {
        for (auto& f : fp) {
            if (f) {
                hts_close(f);
            }
        }
        fp.clear();
        for (auto& h : hdr) {
            if (h) {
                bcf_hdr_destroy(h);
            }
        }
        hdr.clear();
    }

protected:
    size_t line_counter;
    size_t print_counter;
    float *pp_arr;
    int pp_arr_size;
    std::vector<bcf_hdr_t *> hdr;
    std::vector<htsFile *> fp;
    const std::map<uint32_t, std::vector<uint32_t> >& extractions;
};

int main(int argc, char**argv) {
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

#if 0

    std::cout << "Scanning...\n" << std::endl;

    ScanTraversal ut(start, end);
    ut.traverse_no_destroy(filename);
    ut.finalize();

    const auto& extractions = ut.get_extractions_ref();
    std::cout << "Number of sites to extract : " << extractions.size() << std::endl;

    std::cout << "Extracting...\n" << std::endl;

    ExtractTraversal et(ut.get_extractions_ref());
    et.traverse_no_unpack_no_destroy(filename);
    et.finalize();

    ut.destroy();
    et.destroy();
#else

    std::cout << "Extracting...\n" << std::endl;

    PPExtractTraversal ppet(start, end);
    ppet.traverse_no_destroy(filename);
    ppet.finalize();
    ppet.write_to_file(ofname);

#endif

    std::cout << "Done !" << std::endl;

    return 0;
}