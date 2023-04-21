#include <iostream>
#include <string>
#include <mutex>
#include <condition_variable>
#include <thread>

#include "CLI11.hpp"
#include "bcf_traversal.hpp"
#include "hts.h"
#include "sam.h"
#include "vcf.h"
#include "het_info_loader.hpp"
#include "sample_info.hpp"
#include "time.hpp"

#define DIST(x,y) (std::max((x),(y))-std::min((x),(y)))

#define DEBUG_SHOW_PILEUP 0

inline char getBase (int code) {
    switch (code) {
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        case 15: return 'N';
        default: return -1;
	}
}

class GlobalAppOptions {
public:
    GlobalAppOptions() {
        app.add_option("-f,--file", var_filename, "Input variant file name");
        app.add_option("-b,--binary-file", bin_filename, "Input-Output het binary file name");
        app.add_option("-S,--sample-file", sample_filename, "Sample list file name");
        app.add_option("-l,--sample-list", sample_list_filename, "Sample to use list file name");
        app.add_option("-I,--project-id", project_id, "UKB Project ID");
        app.add_option("-p,--cram-path", cram_path, "CRAM files path");
        app.add_option("-s,--start", start, "Starting sample position");
        app.add_option("-e,--end", end, "End sample position (excluded)");
        app.add_option("-t,--num-threads", n_threads, "Number of threads, default is 1, set to 0 for auto");
        app.add_flag("-v,--verbose", verbose, "Verbose mode, display more messages");
    }

    CLI::App app{"Ultralight Phase Caller"};
    std::string project_id = "XXXXX";
    std::string cram_path = "/mnt/project/Bulk/Whole genome sequences/Whole genome CRAM files";
    std::string var_filename = "-";
    std::string bin_filename = "-";
    std::string sample_filename = "-";
    std::string sample_list_filename = "-";
    size_t start = 0;
    size_t end = -1;
    size_t n_threads = 1;
    bool verbose = false;
};

GlobalAppOptions global_app_options;

class Hetp {
public:
    Hetp() {}
    Hetp(const VarInfo *var_info) :
        var_info(var_info),
        a0_reads_p(new std::set<std::string>),
        a1_reads_p(new std::set<std::string>) {}
    Hetp(float *pp, int *gt, const VarInfo *var_info) :
        Hetp(var_info) {
        pp_arr = pp;
        gt_arr = gt;
    }
    Hetp(uint32_t* ptr, const std::vector<VarInfo>& vi) : Hetp((float*)ptr+3, (int*)ptr+1, &vi[*ptr]) {}

    bool is_snp() {
        return var_info->snp;
    }

    char get_allele0() const {
        return bcf_gt_allele(gt_arr[0]) ? var_info->alt[0] : var_info->ref[0];
    }

    char get_allele1() const {
        return bcf_gt_allele(gt_arr[1]) ? var_info->alt[0] : var_info->ref[0];
    }

    // Very inefficient, but used only for debug
    std::string to_string() const {
        //std::string result(bcf_hdr_id2name(hdr, rec->rid)); // segfault ...
        std::string result(var_info->to_string());
        result += "\t" + std::to_string(bcf_gt_allele(gt_arr[0])) + "|" + std::to_string(bcf_gt_allele(gt_arr[1])) + ":" + std::to_string(get_pp());

        return result;
    }

    void reverse_phase() {
        int a0 = bcf_gt_allele(gt_arr[0]);
        int a1 = bcf_gt_allele(gt_arr[1]);

        gt_arr[0] = bcf_gt_unphased(a1); // First allele is always unphased per BCF standard
        gt_arr[1] = bcf_gt_phased(a0);
        reversed = true;
        // The reads that associate to that allele are now swapped
        a0_reads_p.swap(a1_reads_p);
    }

    float get_pp() const {
        /// @note NaN is when PP is not given (e.g., common variants)
        return std::isnan(pp_arr[0]) ? 1.0 : pp_arr[0];
    }

    void set_validated_pp(size_t number_of_reads) {
        pp_arr[0] += number_of_reads+1;
    }

    void set_bad_pp(size_t number_of_reads) {
        // Set a negative value
        pp_arr[0] = -number_of_reads - pp_arr[0];
    }

    const VarInfo *var_info;
    std::unique_ptr<std::set<std::string> > a0_reads_p;
    std::unique_ptr<std::set<std::string> > a1_reads_p;

protected:
    bool reversed = false;
    float *pp_arr = NULL;
    int *gt_arr = NULL;
};

class Het : public Hetp {
public:
    Het() {}
    Het(bcf1_t *rec, bcf_hdr_t *hdr) :
        Hetp(new VarInfo(rec, hdr)),
        var_info_up(var_info)
    {
        int res = bcf_get_format_float(hdr, rec, "PP", &pp_arr, &pp_arr_size);
        if (res < 0) {
            std::cerr << "Could not extract PP for pos : " << rec->pos << std::endl;
        }
        res = bcf_get_genotypes(hdr, rec, &gt_arr, &gt_arr_size);
        if (res < 0) {
            std::cerr << "Could not extract GT for pos : " << rec->pos << std::endl;
        }
    }

    ~Het() {
        if (gt_arr) {
            free(gt_arr);
            gt_arr = NULL;
        }
        if (pp_arr) {
            free(pp_arr);
            pp_arr = NULL;
        }
    }

    std::unique_ptr<const VarInfo> var_info_up;

private:
    int pp_arr_size = 0;
    int gt_arr_size = 0;
};

class HetTrio {
public:
    HetTrio(HetTrio* prev, Hetp* self, HetTrio* next) :
        prev(prev),
        self(self),
        next(next) {
    }
    HetTrio *prev = NULL;
    Hetp *self = NULL;
    HetTrio *next = NULL;

    bool has_prev() const {
        return !!prev;
    }

    bool has_next() const {
        return !!next;
    }

    size_t distance_to_prev() const {
        return DIST(prev->self->var_info->pos1, self->var_info->pos1);
    }

    size_t distance_to_next() const {
        return DIST(self->var_info->pos1, next->self->var_info->pos1);
    }
};

void het_trio_list_from_hets(std::vector<std::unique_ptr<HetTrio> >& het_trios, std::vector<std::unique_ptr<Hetp> >& hets) {
    // Create the HetTrio linked list
    HetTrio* prev = NULL;
    for (size_t i = 0; i < hets.size(); ++i) {
        // Filter out non SNPs
        if (!hets[i]->is_snp()) {
            continue;
        }
        het_trios.push_back(std::make_unique<HetTrio>(prev, hets[i].get(), (HetTrio*)NULL));
        if (prev) {
            prev->next = het_trios.back().get();
        }
        prev = het_trios.back().get();
    }
}

void het_trio_list_from_hets_filter_distance(std::vector<std::unique_ptr<HetTrio> >& het_trios, std::vector<std::unique_ptr<Hetp> >& hets, const size_t MAX_DISTANCE) {
    // Create the HetTrio linked list
    HetTrio* prev = NULL;
    size_t dist_back;
    size_t dist_forward;
    for (size_t i = 0; i < hets.size(); ++i) {
        // Filter out non SNPs
        if (!hets[i]->is_snp()) {
            continue;
        }
        if (i) {
            dist_back = DIST(hets[i-1]->var_info->pos1, hets[i]->var_info->pos1);
        } else {
            dist_back = std::numeric_limits<size_t>::max();
        }
        if (i+1 < hets.size()) {
            dist_forward = DIST(hets[i+1]->var_info->pos1, hets[i]->var_info->pos1);
        } else {
            dist_forward = std::numeric_limits<size_t>::max();
        }
        auto closest = std::min(dist_back, dist_forward);
        // Filter out SNPs that have no close neighbors
        if (closest > MAX_DISTANCE) {
            continue;
        }

        het_trios.push_back(std::make_unique<HetTrio>(prev, hets[i].get(), (HetTrio*)NULL));
        if (prev) {
            prev->next = het_trios.back().get();
        }
        prev = het_trios.back().get();
    }
}


// This function modifies the linked-list but not the vector size
HetTrio* filter_het_trios(std::vector<std::unique_ptr<HetTrio> >& het_trios, size_t distance_threshold = 1000) {
    HetTrio* het = het_trios.front().get();
    HetTrio* first = NULL;

    size_t dist_back;
    size_t dist_forward;

    while (het) {
        if (het->prev) {
            dist_back = DIST(het->prev->self->var_info->pos1, het->self->var_info->pos1);
        } else {
            dist_back = std::numeric_limits<size_t>::max();
        }
        if (het->next) {
            dist_forward = DIST(het->next->self->var_info->pos1, het->self->var_info->pos1);
        } else {
            dist_forward = std::numeric_limits<size_t>::max();
        }

        auto closest = std::min(dist_back, dist_forward);

        if (closest > distance_threshold) {
            // Remove het from linked list
            if (het->prev) {
                het->prev->next = het->next;
            }
            if (het->next) {
                het->next->prev = het->prev;
            }
            het->self = NULL;
        } else {
            if (!first) {
                first = het;
            }
        }

        het = het->next;
    }

    return first;
}

class DataCaller {
public:
    class DataCallerError {
    public:
        DataCallerError(std::string what) : what(what) {
            std::cerr << what << std::endl;
        }
        std::string what;
    };

    htsFile * fp;				// File handler
    sam_hdr_t * hdr;			// File header
    hts_idx_t * idx;			// Index handler
    hts_itr_t * iter;			// NULL if a region not specified
    int min_baseQ;
    int min_mapQ;				// mapQ filter
    bool opened;

    DataCaller (int _min_baseQ = 30 /** @todo */,int _min_mapQ = 50) :
        fp(NULL),
        hdr(NULL),
        idx(NULL),
        iter(NULL),
        min_baseQ(_min_baseQ),
        min_mapQ(_min_mapQ),
        opened(false)
    {
    }

    ~DataCaller() {
        close();
    }

    void open (std::string cram_file) {
        fp = hts_open(cram_file.c_str(), "r");
        if (!fp) {
            std::string error("Cannot open ");
            error += cram_file;
            throw DataCallerError(error);
        }
        idx = sam_index_load(fp, std::string(cram_file + ".crai").c_str());
        if (!idx) {
            throw DataCallerError(std::string("Failed to load index file"));
        }
        hdr = sam_hdr_read(fp);
        if (!hdr) {
            std::string error("Failed to read header from file ");
            error += cram_file;
            throw DataCallerError(error);
        }
        opened = true;
    }

    bool isOpened() { return opened; }

    void close() {
        if (hdr) {
            bam_hdr_destroy(hdr);
            hdr = NULL;
        }
        if (idx) {
            hts_idx_destroy(idx);
            idx = NULL;
        }
        if (fp) {
            hts_close(fp);
            fp = NULL;
        }
        if (iter) {
            hts_itr_destroy(iter);
            iter = NULL;
        }
        opened = false;
    }

    void jump(std::string chr, int start, int end) {
        static thread_local std::string extra("");
        std::string region = extra + chr + ":" + std::to_string(start) + "-" + std::to_string(end);
        if (iter) {
            /* Destroy old iterator (otherwise memory leak) */
            sam_itr_destroy(iter);
            iter = NULL;
        }
        /// @todo here querys (as string) is used, queryi should be more efficient (skips the to string and back)
        iter = sam_itr_querys(idx, hdr, region.c_str());
        if (!iter) {
            // 1KGP3 aligns on GRCh38, the CRAMs name the contigs with "chr" e.g., "chr17" but the VCF calls have "17"
            // Try if adding "chr" solves the problem
            std::string test("chr");
            region = test + region;
            iter = sam_itr_querys(idx, hdr, region.c_str());
            if (!iter) {
                std::cerr << "Could not jump to region [" << region << "]" << std::endl;
            } else {
                // If the extra "chr" solved the problem, always add it
                extra = test;
            }
        }
    }

    int begin() {
        return iter->beg;
    }

    int end() {
        return iter->end;
    }

    void pileup_reads(const bam_pileup1_t * v_plp, int n_plp, Hetp* het) {
        n_bases_total++;
        for (int i = 0 ; i < n_plp ; ++i) {
            const bam_pileup1_t *p = v_plp + i;
            if (p->is_del || p->is_refskip || p->indel == 1) {
                n_bases_indel++;
                continue;
            } else {
                char base = getBase(bam_seqi(bam_get_seq(p->b), p->qpos));
                char qual = (char)bam_get_qual(p->b)[p->qpos];

                if (qual < min_baseQ) {
                    n_bases_lowqual ++;
                    continue;
                }

                char a0 = het->get_allele0();
                char a1 = het->get_allele1();

                if constexpr (DEBUG_SHOW_PILEUP) {
                    std::cout << "Read name : " << bam_get_qname(p->b) << " position : " << p->qpos << " base : " << base << std::endl;
                }

                if (base != a0 && base != a1) {
                    n_bases_mismatch ++;
                    continue;
                }

                if (base == a0) {
                    // Read that has a0
                    het->a0_reads_p->insert(std::string(bam_get_qname(p->b)));
                } else if (base == a1) {
                    // Read that has a1
                    het->a1_reads_p->insert(std::string(bam_get_qname(p->b)));
                }
            }
        }
        if constexpr (DEBUG_SHOW_PILEUP) {
            std::cout << " --- " << std::endl;
        }
    }

    size_t n_bases_indel = 0;
    size_t n_bases_total = 0;
    size_t n_bases_lowqual = 0;
    size_t n_bases_mismatch = 0;
};

static int pileup_filter(void *data, bam1_t *b) {
    DataCaller *aux = (DataCaller*) data;
    int ret;
    while (1) {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if (ret < 0) break;
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
        if (b->core.flag & BAM_FPAIRED) {
            if (!(b->core.flag & BAM_FPROPER_PAIR)) continue;
            if (b->core.flag & BAM_FMUNMAP) continue;
            if ((b->core.flag & BAM_FREVERSE) == (b->core.flag & BAM_FMREVERSE)) continue;
        }
        if ((int)b->core.qual < aux->min_mapQ) continue;
        break;
    }
    return ret;
}

class HetTraversal : public BcfTraversal {
    public:
    HetTraversal(std::vector<std::unique_ptr<Hetp> >& hets) :
        hets(hets)
    {}

    virtual void handle_bcf_line() override {
        hets.push_back(std::make_unique<Het>(bcf_fri.line, bcf_fri.sr->readers[0].header));
    }

    std::vector<std::unique_ptr<Hetp> >& hets;
};

class HetInfoPtrContainerExt : HetInfoMemoryMap::HetInfoPtrContainer {
public:
    HetInfoPtrContainerExt (HetInfoMemoryMap& parent, size_t sample_idx, const std::vector<VarInfo>& vi) :
        HetInfoMemoryMap::HetInfoPtrContainer(parent, sample_idx), vi(vi) {}

    void fill_het_info_ext(std::vector<std::unique_ptr<Hetp> >& v) {
        v.clear();
        for (size_t i = 0; i < size; ++i) {
            v.emplace_back(std::make_unique<Hetp>(start_pos+i*Iterator_type::skip(), vi));
        }
    }

    const std::vector<VarInfo>& vi;
};

class Evidence {
public:
    Evidence() {
        std::cerr << "Warning empty evidence created !" << std::endl;
    }

    Evidence(Hetp* self, Hetp* other) :
        self(self), other(other)
    {
        // Check strand 0 concordance
        for (auto& s0r : *self->a0_reads_p) {
            if (other->a0_reads_p->find(s0r) != other->a0_reads_p->end()) {
                correct_phase_pir++;
            }
            if (other->a1_reads_p->find(s0r) != other->a1_reads_p->end()) {
                reverse_phase_pir++;
            }
        }
        // Check strand 1 concordance
        for (auto& s1r : *self->a1_reads_p) {
            if (other->a0_reads_p->find(s1r) != other->a0_reads_p->end()) {
                reverse_phase_pir++;
            }
            if (other->a1_reads_p->find(s1r) != other->a1_reads_p->end()) {
                correct_phase_pir++;
            }
        }
    }

    Hetp* self = NULL;
    Hetp* other = NULL;
    size_t correct_phase_pir = 0;
    size_t reverse_phase_pir = 0;
};

class Rephaser {
public:
    /* Rephase if PP below (strict) < PP_THRESHOLD */
    Rephaser() : PP_THRESHOLD(1.0) {}

protected:

#if 0
    // Validates the phase, this requires high PP neighbors
    inline void validated_phase(HetTrio* het, HetTrio* other_het, size_t& correct_phase_pir, size_t& reverse_phase_pir) {
        if (other_het && other_het->self->get_pp() > OTHER_PP_THRESHOLD) {
            check_phase(het, other_het, correct_phase_pir, reverse_phase_pir);
        }
    }
#endif

    inline void gather_evidence(HetTrio* het, std::vector<Evidence>& evidence, size_t max_dist, size_t max_steps) {
        look_back(het, evidence, max_dist, max_steps);
        look_ahead(het, evidence, max_dist, max_steps);
    }

    inline void look_back(HetTrio* het, std::vector<Evidence>& evidence, size_t max_dist, size_t max_steps) {
        HetTrio *prev = het->prev;
        while (prev && max_steps--) {
            auto distance = DIST(prev->self->var_info->pos1, het->self->var_info->pos1);
            if (distance > max_dist) {
                return;
            } else {
                evidence.push_back(Evidence(het->self, prev->self));
            }
            prev = prev->prev;
        }
    }

    inline void look_ahead(HetTrio* het, std::vector<Evidence>& evidence, size_t max_dist, size_t max_steps) {
        HetTrio *next = het->next;
        while (next && max_steps--) {
            auto distance = DIST(next->self->var_info->pos1, het->self->var_info->pos1);
            if (distance > max_dist) {
                return;
            } else {
                evidence.push_back(Evidence(het->self, next->self));
            }
            next = next->next;
        }
    }

public:
    class RephaserStatistics {
    public:
        size_t rephase_tries = 0;
        size_t rephase_success = 0;
        size_t rephase_switch = 0;
        size_t rephase_mixed = 0;
        size_t no_reads = 0;
        size_t num_hets = 0;
    };

    void rephase(std::vector<std::unique_ptr<HetTrio> >& het_trios, const std::string& cram_file) {
        DataCaller dc;
        dc.open(cram_file);
        if (!dc.isOpened()) {
            std::string error("Cannot open data for file ");
            error += cram_file;
            throw DataCaller::DataCallerError(error);
        }

        stats.num_hets = het_trios.size();

        HetTrio* current_het = het_trios.front().get();

        if (global_app_options.verbose) {
            std::cout << "There are " << stats.num_hets << " het sites after initial filter (e.g., non SNP)" << std::endl;
            std::cout << "Starting piling up reads for those hets..." << std::endl;
        }

        size_t het_counter = 0;

        /* PILE UP */
        while(current_het) {
            std::string stid = current_het->self->var_info->contig;
            int target_tid = sam_hdr_name2tid(dc.hdr, stid.c_str());

            /* The iterator is really important for performance */
            /** @todo Not sure about the boundaries around the iterator though ... this could be reduced*/
            /* This will create an iterator that is used for the pileup instead of going through all the reads */
            dc.jump(stid, current_het->self->var_info->pos1 - 300, current_het->self->var_info->pos1 + 300);

            /* Init the pileup with the pileup function, it will use an iterator instead of reading all recs from file */
            const bam_pileup1_t *v_plp;
            int n_plp(0), curr_tid(0), curr_pos(0);
            bam_plp_t s_plp = bam_plp_init(pileup_filter, (void*)&dc);

            while ((v_plp = bam_plp_auto(s_plp, &curr_tid, &curr_pos, &n_plp)) != 0) {
                // The position in the VCF/BCF is 1 based not 0 based
                if (curr_tid == target_tid && curr_pos == (int)current_het->self->var_info->pos1) {
                    dc.pileup_reads(v_plp, n_plp, current_het->self);
                    //current_het = current_het->next;
                    break;
                }
                if (!current_het) {
                    break;
                }
            }
            bam_plp_reset(s_plp);
            bam_plp_destroy(s_plp);

            current_het = current_het->next;

            het_counter++;
            if (global_app_options.verbose) {
                printf("\033[A\033[2K");
                std::cout << "Handled " << het_counter << " het variants" << std::endl;
            }
        }

        /* REPHASE/VALIDATE */
        for (auto& h : het_trios) {
            if (!h->self) continue; // Filtered out
            // Sanity check
            if (!h->self->a0_reads_p->size() && !h->self->a1_reads_p->size()) {
                if (global_app_options.verbose) {
                    std::cerr << "No reads mapped to " << h->self->var_info->contig << ":" << h->self->var_info->pos1 << std::endl;
                }
                stats.no_reads++;
            }

            if (0 /* validate */) {
                std::vector<Evidence> evidence(0);

                gather_evidence(h.get(), evidence, MAX_DISTANCE, 1 /* single step */);

                if (global_app_options.verbose) {
                    std::cout << "Evidence vector is of size : " << evidence.size() << std::endl;

                    for (auto& e : evidence) {
                        std::cout << "Evidence relative to : " << e.other->to_string() << std::endl;
                        std::cout << "Correct phase PIRs : " << e.correct_phase_pir << std::endl;
                        std::cout << "Reverse phase PIRs : " << e.reverse_phase_pir << std::endl;

                        if (e.correct_phase_pir && e.reverse_phase_pir) {
                            std::cerr << "Warning ! " << e.correct_phase_pir << " reads confirm the phase and " << e.reverse_phase_pir << " reads say the phase is wrong" << std::endl;
                            stats.rephase_mixed++;
                        }
                    }
                }

                // Cumulate over evidence, this could be another algorithm
                size_t correct_phase_pir = 0;
                size_t reverse_phase_pir = 0;
                for (auto& e : evidence) {
                    correct_phase_pir += e.correct_phase_pir;
                    reverse_phase_pir += e.reverse_phase_pir;
                }
                // We need at least to have seen some reads
                if (correct_phase_pir || reverse_phase_pir) {
                    stats.rephase_success++;
                    if (correct_phase_pir > reverse_phase_pir) {
                        // Phase is correct
                        h->self->set_validated_pp(correct_phase_pir);
                    } else {
                        // Phase is incorrect
                        stats.rephase_switch++;
                        h->self->set_bad_pp(reverse_phase_pir);
                    }
                }
            } else /* rephase */ {
                // Requires to be rephased
                if (h->self->get_pp() < PP_THRESHOLD) {
                    stats.rephase_tries++;
                    if (global_app_options.verbose) {
                        std::cout << h->self->to_string() << " requires work" << std::endl;
                    }

                    std::vector<Evidence> evidence(0);

                    gather_evidence(h.get(), evidence, MAX_DISTANCE, 100 /* max steps */);

                    if (global_app_options.verbose) {
                        std::cout << "Evidence vector is of size : " << evidence.size() << std::endl;

                        for (auto& e : evidence) {
                            std::cout << "Evidence relative to : " << e.other->to_string() << std::endl;
                            std::cout << "Correct phase PIRs : " << e.correct_phase_pir << std::endl;
                            std::cout << "Reverse phase PIRs : " << e.reverse_phase_pir << std::endl;

                            if (e.correct_phase_pir && e.reverse_phase_pir) {
                                std::cerr << "Warning ! " << e.correct_phase_pir << " reads confirm the phase and " << e.reverse_phase_pir << " reads say the phase is wrong" << std::endl;
                                stats.rephase_mixed++;
                            }
                        }
                    }

                    // Cumulate over evidence, this could be another algorithm
                    /// @todo check the PP score of the other variants
                    size_t correct_phase_pir = 0;
                    size_t reverse_phase_pir = 0;
                    for (auto& e : evidence) {
                        correct_phase_pir += e.correct_phase_pir;
                        reverse_phase_pir += e.reverse_phase_pir;
                    }
                    // We need at least to have seen some reads
                    if (correct_phase_pir || reverse_phase_pir) {
                        stats.rephase_success++;
                        if (correct_phase_pir > reverse_phase_pir) {
                            // Phase is correct
                            h->self->set_validated_pp(correct_phase_pir);
                        } else {
                            // Phase is incorrect
                            h->self->reverse_phase();
                            h->self->set_validated_pp(reverse_phase_pir);
                        }
                        if (global_app_options.verbose) {
                            std::cout << "This is the read validated entry :" << std::endl;
                            std::cout << h->self->to_string() << std::endl;
                            std::cout << "---" << std::endl;
                        }
                    }
                }
            }
        }

        dc.close();

        if (global_app_options.verbose) {
            std::cout << "Tried to rephase " << stats.rephase_tries << " het sites, succeeded with " << stats.rephase_success << std::endl;
            std::cout << "On these " << stats.rephase_success << ", " << stats.rephase_switch << " needed flipping" << std::endl;
        }
    }

    const float PP_THRESHOLD;
    const float OTHER_PP_THRESHOLD = 0.9;
    const size_t MAX_DISTANCE = 1000;
    RephaserStatistics stats;
};

void rephase_sample(const std::vector<VarInfo>& vi, HetInfoMemoryMap& himm, const std::string& cram_file, size_t himm_sample_idx) {
    std::vector<std::unique_ptr<Hetp> > hets;
    std::vector<std::unique_ptr<HetTrio> > het_trios;

    // Get hets from memory map
    HetInfoPtrContainerExt hipce(himm, himm_sample_idx, vi);
    // Hets created from the memory map will directy edit the file on rephase
    hipce.fill_het_info_ext(hets);
    if (global_app_options.verbose) {
        std::cout << "Loaded " << hets.size() << " genotypes from file" << std::endl;
    }
    //het_trio_list_from_hets(het_trios, hets);
    het_trio_list_from_hets_filter_distance(het_trios, hets, 1000);

    try {
        Rephaser r;
        r.rephase(het_trios, cram_file);
    } catch (DataCaller::DataCallerError e) {
        return;
    }
}

class PhaseCaller {
public:
    PhaseCaller(std::string& vcf_filename, std::string& bin_filename, std::string& sample_filename, std::string& samples_to_do_filename,
                size_t n_threads) :
        threads(n_threads, NULL),
        active_threads(n_threads, false),
        samples_to_do(samples_to_do_filename),
        sil(sample_filename),
        vil(vcf_filename),
        himm(bin_filename, PROT_READ | PROT_WRITE)
    {
    }

    PhaseCaller(std::string& vcf_filename, std::string& bin_filename, std::string& sample_filename, size_t n_threads) :
        threads(n_threads, NULL),
        active_threads(n_threads, false),
        samples_to_do("-"),
        sil(sample_filename),
        vil(vcf_filename),
        himm(bin_filename, PROT_READ | PROT_WRITE)
    {
    }

    ~PhaseCaller() {
        std::unique_lock<std::mutex> lk(mutex);
        // Final cleanup
        for (size_t i = 0; i < threads.size(); ++i) {
            if (threads[i]) {
                threads[i]->join();
                delete threads[i];
                threads[i] = NULL;
            }
        }
    }

    void rephase_orchestrator(size_t start_id, size_t stop_id) {
        for (size_t i = start_id; i < stop_id; ++i) {
            thread_fun(0, i, i);
        }
    }

private:
    inline size_t find_free(const std::vector<bool> v) {
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i] == false) return i;
        }
        return v.size();
    }

    inline std::string cram_filename(const std::string& sample_name) const {
        std::string cram_file(global_app_options.cram_path);
        cram_file += "/" + sample_name.substr(0, 2);
        cram_file += "/" + sample_name + "_" + global_app_options.project_id + "_0_0.cram";
        return cram_file;
    }

    /*****************/
    /* MAIN FUNCTION */
    /*****************/
    std::function<void(size_t, size_t, size_t)> thread_fun = [this](size_t thread_idx, size_t sample_idx, size_t himm_idx){
        // Get sample
        auto sample = sil.samples[sample_idx];
        const std::string& sample_name = sample.name;
        // Don't try withdrawn samples
        if (sample_name[0] == 'W') {
            /// @todo remove this is UKB specific
            std::lock_guard lk(mutex);
            std::cerr << "Withdrawn sample " << sample_name << " will not rephase because sequencing data is not available" << std::endl;
        } else {
            std::string cram_file = sample.cram_file_path;
            if (cram_file.length() == 0) {
                // Generate the corresponding cram file path
                std::string cram_file = cram_filename(sample_name);
                std::cout << "Sample idx: " << sample_idx << " name: " << sample_name << " cram path: " << cram_file << std::endl;
            } else {
                std::cout << sample.to_string() << std::endl;
            }

            if (!fs::exists(cram_file)) {
                std::lock_guard lk(mutex);
                std::cerr << "Cannot find file " << cram_file << " skipping ..." << std::endl;
            } else {
                rephase_sample(vil.vars, himm, cram_file, himm_idx);
            }
        }
        {
            std::lock_guard lk(mutex);
            std::cout << "Thread " << thread_idx << " finished" << std::endl;
            active_threads[thread_idx] = false;
        }
        cv.notify_all();
    };

public:
    void rephase_orchestrator_multi_thread(size_t start_id, size_t stop_id) {
        for (size_t i = start_id; i < stop_id; ++i) {
            std::unique_lock<std::mutex> lk(mutex);
            size_t ti = find_free(active_threads);
            cv.wait(lk, [&]{ti = find_free(active_threads); return ti < active_threads.size(); });

            if (threads[ti]) {
                // If a thread was launched but finished
                threads[ti]->join();
                std::cout << "Joined thread " << ti << std::endl;
                delete threads[ti];
                threads[ti] = NULL;
            }

            std::cout << "Launching thread " << ti << std::endl;
            active_threads[ti] = true;
            threads[ti] = new std::thread(thread_fun, ti, i, i);
        }

        // Final cleanup
        for (size_t i = 0; i < threads.size(); ++i) {
            if (threads[i]) {
                threads[i]->join();
                delete threads[i];
                threads[i] = NULL;
            }
        }
    }

    void rephase_orchestrator_multi_thread() {
        if (samples_to_do.sample_names.size()) {
            rephase_orchestrator_multi_thread_with_list();
        } else {
            /** @todo remove this message */
            std::cout << "Running version compatible with splitted bins" << std::endl;
            rephase_orchestrator_multi_thread_without_list();
        }
    }

    /* This one is a special case for subsampled binary files */
    void rephase_orchestrator_multi_thread_without_list() {
        for (uint32_t himm_idx = 0; himm_idx < himm.num_samples; ++himm_idx) {
            // Because the himm is subsampled we need the original index wrt sample list
            uint32_t orig_idx = himm.get_orig_idx_of_nth(himm_idx);
            std::unique_lock<std::mutex> lk(mutex);
            size_t ti = find_free(active_threads);
            cv.wait(lk, [&]{ti = find_free(active_threads); return ti < active_threads.size(); });

            if (threads[ti]) {
                // If a thread was launched but finished
                threads[ti]->join();
                std::cout << "Joined thread " << ti << std::endl;
                delete threads[ti];
                threads[ti] = NULL;
            }

            std::cout << "Launching thread " << ti << std::endl;
            active_threads[ti] = true;
            threads[ti] = new std::thread(thread_fun, ti, orig_idx, himm_idx);
        }

        // Final cleanup
        for (size_t i = 0; i < threads.size(); ++i) {
            if (threads[i]) {
                threads[i]->join();
                delete threads[i];
                threads[i] = NULL;
            }
        }
    }

    void rephase_orchestrator_multi_thread_with_list() {
        for (size_t i = 0; i < sil.sample_names.size(); ++i) {
            if (std::find(samples_to_do.sample_names.begin(), samples_to_do.sample_names.end(),
                sil.sample_names[i]) != samples_to_do.sample_names.end()) {
                std::unique_lock<std::mutex> lk(mutex);
                size_t ti = find_free(active_threads);
                cv.wait(lk, [&]{ti = find_free(active_threads); return ti < active_threads.size(); });

                if (threads[ti]) {
                    // If a thread was launched but finished
                    threads[ti]->join();
                    std::cout << "Joined thread " << ti << std::endl;
                    delete threads[ti];
                    threads[ti] = NULL;
                }

                std::cout << "Launching thread " << ti << std::endl;
                active_threads[ti] = true;
                threads[ti] = new std::thread(thread_fun, ti, i, i);
            }
        }

        // Final cleanup
        for (size_t i = 0; i < threads.size(); ++i) {
            if (threads[i]) {
                threads[i]->join();
                delete threads[i];
                threads[i] = NULL;
            }
        }
    }

    std::vector<std::thread*> threads;
    std::vector<bool> active_threads;
    std::mutex mutex;
    std::condition_variable cv;

    SampleInfoLoader samples_to_do;
    SampleInfoLoader sil;
    VarInfoLoader vil;
    // We memory map as to save RAM space and update the file in place
    HetInfoMemoryMap himm;
};

int main(int argc, char**argv) {
    auto start_time = std::chrono::steady_clock::now();

    auto& opt = global_app_options;
    auto& app = global_app_options.app;
    CLI11_PARSE(app, argc, argv);

    if (opt.var_filename.compare("-") == 0) {
        std::cerr << "Requires variant VCF/BCF file" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }
    if (opt.bin_filename.compare("-") == 0) {
        std::cerr << "Requires het binary file" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }
    if (opt.sample_filename.compare("-") == 0) {
        std::cerr << "Requires sample file name" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }
    if (opt.n_threads == 0) {
        opt.n_threads = std::thread::hardware_concurrency();
        std::cerr << "Setting number of threads to " << opt.n_threads << std::endl;
    }

    PhaseCaller pc(opt.var_filename, opt.bin_filename, opt.sample_filename, opt.sample_list_filename, opt.n_threads);
    pc.rephase_orchestrator_multi_thread();
    printElapsedTime(start_time, std::chrono::steady_clock::now());

    return 0;
}