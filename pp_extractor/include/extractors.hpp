#ifndef __EXTRACTORS_HPP__
#define __EXTRACTORS_HPP__

#include <cmath>

#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"
#include "bcf_traversal.hpp"
#include "het_info.hpp"
#include "var_info.hpp"
#include "fs.hpp"

constexpr size_t PLOIDY_2 = 2;

class PPPred {
public:
    PPPred(float pp_threshold) : pp_threshold(pp_threshold) {}

    bool operator()(HetInfo hi) const { return !(std::isnan(hi.pp)) && hi.pp < pp_threshold; }

protected:
    float pp_threshold;
};

class PPExtractTraversal : public BcfTraversal {
public:
    PPExtractTraversal(size_t start_id, size_t stop_id, size_t fifo_size, bool pp_from_maf, bool pp_from_af) :
        FIFO_SIZE(fifo_size),
        PP_THRESHOLD(0.99),
        MAF_THRESHOLD(0.001),
        pp_arr(NULL),
        pp_arr_size(0),
        start_id(start_id),
        stop_id(stop_id),
        line_counter(0),
        print_counter(0),
        pred(PP_THRESHOLD),
        progress(0),
        pp_from_maf(pp_from_maf),
        pp_from_af(pp_from_af),
        extract_acan(pp_from_maf), // MAF requires AC/AN
        search_line_value(false) {
        if (pp_from_maf) {
            std::cout << "The PP score will be generated from MAF" << std::endl;
        }
        if (pp_from_af) {
            std::cout << "The PP score will be generated from AF" << std::endl;
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

        // We need AC (Allele Count) and AN (Allele Number) for PP from MAF
        int AC = 0;
        int* pAC = NULL;
        int nAC = 0;
        int AN = 0;
        int* pAN = NULL;
        int nAN = 0;
        float AF = 0.0;
        float* pAF = NULL;
        int nAF = 0;
        float synthetic_pp = 0.0;

        if (extract_acan) {
            res = bcf_get_info_int32(header, line, "AC", &pAC, &nAC);
            // If not in the VCF then compute it
            if (res < 0) {
                for (size_t i = start_id; i < stop_id; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[i*PLOIDY_2]) == 1) {
                        AC++;
                    }
                    if (bcf_gt_allele(bcf_fri.gt_arr[i*PLOIDY_2+1]) == 1) {
                        AC++;
                    }
                }
            } else {
                AC = *pAC;
            }
            // If not in the VCF then compute it
            res = bcf_get_info_int32(header, line, "AN", &pAN, &nAN);
            if (res < 0) {
                AN = (stop_id - start_id) * PLOIDY_2;
            } else {
                AN = *pAN;
            }
            float maf = float(AC)/AN;
            synthetic_pp = (maf > MAF_THRESHOLD) ? NAN : 0.5 + maf / 2.0;
            //std::cout << "AC : " << AC << "\tAN : " << AN << "\tsynth PP " << synthetic_pp << std::endl;
        }

        if (pp_from_af) {
            res = bcf_get_info_float(header, line, "AF", &pAF, &nAF);
            if (res < 0) {
                std::cerr << "Could not get AF info field ! PP set to NaN" << std::endl;
                synthetic_pp = NAN;
            } else {
                AF = *pAF;
                if (AF > 0.5) {
                    // Complement to get the minor allele frequency
                    AF = 1.0 - AF;
                }
                synthetic_pp = (AF > MAF_THRESHOLD) ? NAN : 0.5 + AF / 2.0;
                //std::cout << "AF: " << AF << " Synth PP : " << synthetic_pp << std::endl;
            }
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
                } else if (pp_from_maf || pp_from_af) {
                    pp = synthetic_pp;
                }

                if (AC == 1 && pp >= PP_THRESHOLD) {
                    /* Edge case for old version of SHAPEIT5 that would score
                       singletons phased with only one of their parents with
                       a PP of 1.0, In about 95% of cases the singleton comes
                       from the other parent if not observed in the known parent
                       but it could also be a de novo mutation on the known
                       parent haplotype, these cases should not be scored 1.0 */
                       pp = 0.97; /* arbitrary value */
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
        if (progress) {
            if (++print_counter == progress) {
                print_counter = 0;
                printf("\033[A\033[2K");
                std::cout << "Handled " << bcf_fri.line_num << " VCF entries (lines)" << std::endl;
            }
        }

        // Free memory
        if (pAN) {
            free(pAN);
            pAN = NULL;
        }
        if (pAC) {
            free(pAC);
            pAC = NULL;
        }
        if (pAF) {
            free(pAF);
            pAF = NULL;
        }
    }

    void set_progress(const size_t progress) {
        this->progress = progress;
    }

    void set_search_line_counter(const std::string& filename) {
        search_line_value = true;
        search_in_file = filename;
        if (!fs::exists(search_in_file)) {
            std::cerr << "File : " << search_in_file << " does not exist !" << std::endl;
            throw "File does not exist";
        }
    }

    void set_extract_acan() {
        extract_acan = true;
    }

    void set_maf_threshold(const float maf_threhsold) {
        std::cout << "Setting MAF threshold to " << maf_threhsold << std::endl;
        MAF_THRESHOLD = maf_threhsold;
    }

    void finalize() {
        // Finalize FIFOs
        for (auto& f : fifos) {
            f.finalize();
        }
    }

    void show_info() {
        size_t total_kept = 0;
        size_t total_kept_pred = 0;
        for (auto& f : fifos) {
            auto kept_items = f.get_kept_items_ref();
            total_kept += kept_items.size();
            total_kept_pred += f.get_number_kept_with_pred();
        }

        std::cout << "Extracted a total of " << total_kept << " genotypes" << std::endl;
        std::cout << "From which a total of " << total_kept_pred << " were selected given the predicate" << std::endl;
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
    float MAF_THRESHOLD;
    float *pp_arr;
    int pp_arr_size;
    size_t start_id;
    size_t stop_id;
    size_t line_counter;
    size_t print_counter;
    const PPPred pred;
    size_t progress;
    bool pp_from_maf;
    bool pp_from_af;
    bool extract_acan;
    std::vector<uint32_t> number_of_het_sites;
    std::vector<uint32_t> number_of_low_pp_sites;
    std::vector<uint32_t> number_of_snp_low_pp_sites;
    std::vector<uint32_t> number_of_non_snp;
    std::vector<GenericKeepFifo<HetInfo, PPPred> > fifos;
    bool search_line_value;
    std::string search_in_file;
};

#endif /* __EXTRACTORS_HPP__ */