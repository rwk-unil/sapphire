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
    std::string cram_path = "/mnt/Bulk/Whole Genome Sequences/Whole genome CRAM files";
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
        return pp_arr[0];
    }

    void set_validated_pp(size_t number_of_reads) {
        pp_arr[0] += number_of_reads+1;
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

    bool has_prex() const {
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

    void jump(std::string& chr, int start, int end) {
        std::string region = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
        if (iter) {
            /* Destroy old iterator (otherwise memory leak) */
            sam_itr_destroy(iter);
            iter = NULL;
        }
        /// @todo here querys (as string) is used, queryi should be more efficient (skips the to string and back)
        iter = sam_itr_querys(idx, hdr, region.c_str());
        if (!iter) {
            std::cerr << "Could not jump to region [" << region << "]" << std::endl;
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

void rephase(std::vector<std::unique_ptr<HetTrio> >& het_trios, const std::string& cram_file) {
    DataCaller dc;
    dc.open(cram_file);
    if (!dc.isOpened()) {
        std::cerr << "Cannot open data for file " << cram_file << std::endl;
        exit(-1);
    }

    /* Do some pileup */

    HetTrio* current_het = het_trios.front().get();

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
            /// @todo This should not even happen... plus if the iterator it NULL => segfault
            //if (curr_pos < dc.begin() || curr_pos >= dc.end()) continue;

            // The position in the VCF/BCF is 1 based not 0 based
            if (curr_tid == target_tid && curr_pos == current_het->self->var_info->pos1) {
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
    }

    for (auto& h : het_trios) {
        if (h->self->get_pp() < 1.0) {
            std::cout << h->self->to_string() << " requires work" << std::endl;

            size_t correct_phase_pir = 0;
            size_t reverse_phase_pir = 0;

            auto& strand_0_reads = *h->self->a0_reads_p.get();
            auto& strand_1_reads = *h->self->a1_reads_p.get();

            /**
             * @todo
             *  - Check distances not to waste time searching in the sets of reads
             *  - Check prev and next chain if within distance
             */

            // For all the reads in both strands check the reads in prev and next
            for (auto& s0r : strand_0_reads) {
                if (h->prev) {
                    if (h->prev->self->a0_reads_p->find(s0r) != h->prev->self->a0_reads_p->end()) {
                        correct_phase_pir++;
                    }
                    if (h->prev->self->a1_reads_p->find(s0r) != h->prev->self->a1_reads_p->end()) {
                        reverse_phase_pir++;
                    }
                }
                if (h->next) {
                    if (h->next->self->a0_reads_p->find(s0r) != h->next->self->a0_reads_p->end()) {
                        correct_phase_pir++;
                    }
                    if (h->next->self->a1_reads_p->find(s0r) != h->next->self->a1_reads_p->end()) {
                        reverse_phase_pir++;
                    }
                }
            }
            for (auto& s1r : strand_1_reads) {
                if (h->prev) {
                    if (h->prev->self->a0_reads_p->find(s1r) != h->prev->self->a0_reads_p->end()) {
                        reverse_phase_pir++;
                    }
                    if (h->prev->self->a1_reads_p->find(s1r) != h->prev->self->a1_reads_p->end()) {
                        correct_phase_pir++;
                    }
                }
                if (h->next) {
                    if (h->next->self->a0_reads_p->find(s1r) != h->next->self->a0_reads_p->end()) {
                        reverse_phase_pir++;
                    }
                    if (h->next->self->a1_reads_p->find(s1r) != h->next->self->a1_reads_p->end()) {
                        correct_phase_pir++;
                    }
                }
            }

            std::cout << "Correct phase PIRs : " << correct_phase_pir << std::endl;
            std::cout << "Reverse phase PIRs : " << reverse_phase_pir << std::endl;

            if (correct_phase_pir && reverse_phase_pir) {
                std::cerr << "Warning ! " << correct_phase_pir << " reads confirm the phase and " << reverse_phase_pir << " reads say the phase is wrong" << std::endl;
            }
            // We need at least to have seen some reads
            if (correct_phase_pir || reverse_phase_pir) {
                if (correct_phase_pir > reverse_phase_pir) {
                    // Phase is correct
                    h->self->set_validated_pp(correct_phase_pir);
                } else {
                    // Phase is incorrect
                    h->self->reverse_phase();
                    h->self->set_validated_pp(reverse_phase_pir);
                }
                std::cout << "This is the read validated entry :" << std::endl;
                std::cout << h->self->to_string() << std::endl;
                std::cout << "---" << std::endl;
            }
        }
    }

    dc.close();
}

void rephase_sample(const std::vector<VarInfo>& vi, HetInfoMemoryMap& himm, const std::string& cram_file, size_t sample_idx) {
    std::vector<std::unique_ptr<Hetp> > hets;
    std::vector<std::unique_ptr<HetTrio> > het_trios;

    // Get hets from memory map
    HetInfoPtrContainerExt hipce(himm, sample_idx, vi);
    // Hets created from the memory map will directy edit the file on rephase
    hipce.fill_het_info_ext(hets);
    het_trio_list_from_hets(het_trios, hets);

    try {
        rephase(het_trios, cram_file);
    } catch (DataCaller::DataCallerError e) {
        return;
    }
}

void rephase_example(std::string& vcf_file, std::string& cram_file) {
    std::vector<std::unique_ptr<Hetp> > hets;
    std::vector<std::unique_ptr<HetTrio> > het_trios;
    HetTraversal ht(hets);
    // Will fill the hets vector
    ht.traverse(vcf_file);

    // This filters out the non SNPs
    het_trio_list_from_hets(het_trios, hets); /// @todo handle non SNPs ?

    try {
        rephase(het_trios, cram_file);
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
        himm(bin_filename)
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
            thread_fun(i, 0);
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

    std::function<void(size_t, size_t)> thread_fun = [this](size_t thread_idx, size_t sample_idx){
        // Get sample name
        std::string sample_name = sil.sample_names[sample_idx];
        // Don't try withdrawn samples
        if (sample_name[0] == 'W') {
            std::lock_guard lk(mutex);
            std::cerr << "Withdrawn sample " << sample_name << " will not rephase because sequencing data is not available" << std::endl;
        } else {
            // Generate the corresponding cram file path
            std::string cram_file = cram_filename(sample_name);

            std::cout << "Sample idx: " << sample_idx << " name: " << sample_name << " cram path: " << cram_file << std::endl;
            if (!fs::exists(cram_file)) {
                std::lock_guard lk(mutex);
                std::cerr << "Cannot find file " << cram_file << " skipping ..." << std::endl;
            }
            //rephase_sample(vil.vars, himm, cram_file, sample_idx);
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
            threads[ti] = new std::thread(thread_fun, ti, i);
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
                threads[ti] = new std::thread(thread_fun, ti, i);
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

#if 0
    CLI::App app{"Ultralight phase caller"};
    std::string cram_file = "-";
    std::string vcf_file = "-";
    app.add_option("-f,--file", cram_file, "Input file name");
    app.add_option("--vcf", vcf_file, "VCF input file name");

    CLI11_PARSE(app, argc, argv);

    if (cram_file.compare("-") == 0) {
        std::cerr << "Requires cram_file" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }
    if (vcf_file.compare("-") == 0) {
        std::cerr << "Requires vcf_file" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }

    rephase_example(vcf_file, cram_file);
#else
    auto& opt = global_app_options;
    auto& app = global_app_options.app;
    std::string cram_file = "-";
    app.add_option("-c, --cram", cram_file, "TODO REMOVE");
    CLI11_PARSE(global_app_options.app, argc, argv);

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

#if 0
    std::cout << "Loading sample names" << std::endl;
    SampleInfoLoader sil(opt.sample_filename);
    std::cout << "Loaded " << sil.sample_names.size() << " sample names" << std::endl;

    std::cout << "Loading variants" << std::endl;

    // Load all variants from variant BCF/VCF file once
    VarInfoLoader vil(opt.var_filename);
    std::cout << "Loaded " << vil.vars.size() << " variant entries" << std::endl;

    std::cout << "Memory mapping binary file - This file will be edited !" << std::endl;

    // Memory map the binary het variant file
    HetInfoMemoryMap himm(opt.bin_filename, PROT_READ | PROT_WRITE);
#endif

    PhaseCaller pc(opt.var_filename, opt.bin_filename, opt.sample_filename, opt.sample_list_filename, opt.n_threads);
    pc.rephase_orchestrator_multi_thread(0, 32);
    std::cout << "TEST 2" << std::endl;
    pc.rephase_orchestrator_multi_thread();

    //rephase_orchestrator(vil, himm, sil, 0, 15);
    //rephase_orchestrator_multi_thread(vil, himm, sil, 0, 32, opt.n_threads);
    return 0;

    std::cout << "Rephasing a sample" << std::endl;
    //rephase_sample(vil.vars, himm, cram_file, 465);

#endif


#if 0
    htsFile *file = hts_open(cram_file.c_str(), "r");
    hts_idx_t *idx = sam_index_load(file, std::string(cram_file + ".crai").c_str());

    if (!idx) {
        std::cerr << "Failed to load index file" << std::endl;
        exit(-1);
    }

    if (file->is_cram) {
        std::cout << cram_file << " is a CRAM file" << std::endl;
    }

    sam_hdr_t *hdr = sam_hdr_read(file);
    if (!hdr) {
        std::cerr << "Failed to read header from file " << cram_file << std::endl;
        exit(-1);
    }

    int tid = sam_hdr_name2tid(hdr, "chr1");
    hts_itr_t *it = sam_itr_queryi(idx, tid, 108188080, 108190284);

    if (!it) {
        std::cerr << "Failed to get iterator for query" << std::endl;
        exit(-1);
    }

    int ret;
    bam1_t *r = bam_init1();
    while ((ret = sam_itr_next(file, it, r)) >= 0) {
        std::cout << bam_get_qname(r) << std::endl;
    }

    bam_destroy1(r);
    sam_hdr_destroy(hdr);
    hts_close(file);
#endif

    return 0;
}