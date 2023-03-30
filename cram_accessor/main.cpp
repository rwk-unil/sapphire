#include <iostream>
#include <string>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <filesystem>

namespace fs = std::filesystem;

#include "CLI11.hpp"
#include "hts.h"
#include "sam.h"
#include "vcf.h"

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
        app.add_option("-p,--cram-path", cram_path, "CRAM files path");
        app.add_option("-t,--num-threads", n_threads, "Number of threads, default is 1, set to 0 for auto");
        app.add_flag("-v,--verbose", verbose, "Verbose mode, display more messages");
    }

    CLI::App app{"CRAM Accessor"};
    std::string cram_path = "/mnt/project/Bulk/Whole genome sequences/Whole genome CRAM files/20";
    size_t n_threads = 1;
    bool verbose = false;
};

GlobalAppOptions global_app_options;

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

    void pileup_reads(const bam_pileup1_t * v_plp, int n_plp) {
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

                if constexpr (DEBUG_SHOW_PILEUP) {
                    std::cout << "Read name : " << bam_get_qname(p->b) << " position : " << p->qpos << " base : " << base << std::endl;
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

class Accessor {
public:
    Accessor(size_t num_accesses = 1000) : NUM_ACCESSES(num_accesses) {}

public:

    void access(const std::string& cram_file) {
        DataCaller dc;
        dc.open(cram_file);
        if (!dc.isOpened()) {
            std::string error("Cannot open data for file ");
            error += cram_file;
            throw DataCaller::DataCallerError(error);
        }

        for (size_t i = 0; i < NUM_ACCESSES; ++i) {


            std::string stid("chr20");
            int target_tid = sam_hdr_name2tid(dc.hdr, stid.c_str());

            size_t position = rand() % MAX_REGION_SIZE;

            /* The iterator is really important for performance */
            /** @todo Not sure about the boundaries around the iterator though ... this could be reduced*/
            /* This will create an iterator that is used for the pileup instead of going through all the reads */
            dc.jump(stid, position - 300, position + 300);

            /* Init the pileup with the pileup function, it will use an iterator instead of reading all recs from file */
            const bam_pileup1_t *v_plp;
            int n_plp(0), curr_tid(0), curr_pos(0);
            bam_plp_t s_plp = bam_plp_init(pileup_filter, (void*)&dc);

            while ((v_plp = bam_plp_auto(s_plp, &curr_tid, &curr_pos, &n_plp)) != 0) {
                // The position in the VCF/BCF is 1 based not 0 based
                if (curr_tid == target_tid && curr_pos == position) {
                    dc.pileup_reads(v_plp, n_plp);
                    //current_het = current_het->next;
                    break;
                }
            }
            bam_plp_reset(s_plp);
            bam_plp_destroy(s_plp);
        }

        dc.close();
    }

    const size_t NUM_ACCESSES;
    const size_t MAX_DISTANCE = 1000;
    const size_t MAX_REGION_SIZE = 10000000;
};

class AccessorMain {
public:
    AccessorMain(size_t n_threads) :
        threads(n_threads, NULL),
        active_threads(n_threads, false)
    {
    }

    ~AccessorMain() {
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

private:
    inline size_t find_free(const std::vector<bool> v) {
        for (size_t i = 0; i < v.size(); ++i) {
            if (v[i] == false) return i;
        }
        return v.size();
    }

    /*****************/
    /* MAIN FUNCTION */
    /*****************/
    std::function<void(size_t, std::string)> thread_fun = [this](size_t thread_idx, std::string cram_file){
        std::cout << "CRAM file: " << cram_file << std::endl;
        if (!fs::exists(cram_file)) {
            std::lock_guard lk(mutex);
            std::cerr << "Cannot find file " << cram_file << " skipping ..." << std::endl;
        } else {
                try {
                    Accessor r;
                    r.access(cram_file);
                } catch (DataCaller::DataCallerError e) {
                    return;
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
    void orchestrator_multi_thread() {
        size_t tid = 0;
        for (const auto & entry : fs::directory_iterator(global_app_options.cram_path)) {
            std::cout << entry.path() << std::endl;
            if (entry.path().extension().compare("cram") != 0) {
                // Skip non CRAM files
                continue;
            }

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
            threads[ti] = new std::thread(thread_fun, tid, entry.path());
            tid++;
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
};

int main(int argc, char**argv) {
    auto& opt = global_app_options;
    auto& app = global_app_options.app;
    CLI11_PARSE(app, argc, argv);

    if (opt.n_threads == 0) {
        opt.n_threads = std::thread::hardware_concurrency();
        std::cerr << "Setting number of threads to " << opt.n_threads << std::endl;
    }

    AccessorMain am(opt.n_threads);
    am.orchestrator_multi_thread();

    return 0;
}