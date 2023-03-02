#include <iostream>
#include <fstream>
#include <iterator>
#include <cstddef>
#include <cmath>
#include <fcntl.h>
#include <sys/mman.h>
#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"
#include "../include/bcf_traversal.hpp"
#include "../include/CLI11.hpp"
#include "include/fifo.hpp"
#include "include/het_info.hpp"
#include "include/var_info.hpp"
#include "../include/fs.hpp"
#include "../include/compression.hpp"

class HetInfoLoader {
public:
    HetInfoLoader(std::string bfname, size_t sample) {
        std::ifstream ifs(bfname, std::ios_base::binary | std::ios_base::in);
        if (!ifs.is_open()) {
            std::cerr << "Failed to open binary file " << bfname << std::endl;
            throw "Cannot open file";
        }
        uint32_t endianness = 0;
        ifs.read(reinterpret_cast<char *>(&endianness), sizeof(endianness));
        if (endianness != 0xaabbccdd) {
            std::cerr << "File " << bfname << " has wrong endianness" << std::endl;
            throw "Wrong endianness in file";
        }
        uint32_t num_samples = 0;
        ifs.read(reinterpret_cast<char *>(&num_samples), sizeof(num_samples));
        if (num_samples == 0) {
            std::cerr << "File " << bfname << " has no samples" << std::endl;
            throw "File has no samples";
        }

        std::vector<uint64_t> offset_table(num_samples);
        for(size_t i = 0; i < offset_table.size(); ++i) {
            ifs.read(reinterpret_cast<char *>(&offset_table[i]), sizeof(uint64_t));
        }

        std::cout << "Sample " << sample << " offset : 0x" << std::hex << offset_table[sample] << std::dec << std::endl;

        ifs.seekg(offset_table[sample]);

        uint32_t mark = 0;
        ifs.read(reinterpret_cast<char *>(&mark), sizeof(mark));
        if (mark != 0xd00dc0de) {
            std::cerr << "File " << bfname << " has wrong mark" << std::endl;
            throw "Wrong mark in file";
        }

        uint32_t id = 0;
        ifs.read(reinterpret_cast<char *>(&id), sizeof(id));
        std::cout << "Id : " << id << std::endl;

        uint32_t size = 0;
        ifs.read(reinterpret_cast<char *>(&size), sizeof(size));

        his.resize(size);
        for(auto& hi : his) {
            hi.from_stream(ifs);
        }
    }

    std::vector<HetInfo> his;
};

class HetInfoMemoryMap {
public:
    HetInfoMemoryMap(std::string bfname) : file_size(fs::file_size(bfname)) {
        fd = open(bfname.c_str(), O_RDONLY, 0);
        if (fd < 0) {
            std::cerr << "Failed to open file : " << bfname << std::endl;
            throw "Failed to open file";
        }

        file_mmap_p = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (file_mmap_p == NULL) {
            std::cerr << "Failed to memory map file : " << bfname << std::endl;
            close(fd);
            throw "Failed to mmap file";
        }

        // Test the memory map (first thing is the endianness in the file)
        uint32_t endianness = *(uint32_t*)(file_mmap_p);
        if (endianness != ENDIANNESS) {
            std::cerr << "Bad endianness in memory map" << std::endl;
            throw "Bad endianness";
        }

        num_samples = *((uint32_t*)file_mmap_p+1);
        offset_table = ((uint64_t*)file_mmap_p+1);
    }

    ~HetInfoMemoryMap() {
        if (file_mmap_p) {
            munmap(file_mmap_p, file_size);
            file_mmap_p = NULL;
        }
        if (fd > 0) {
            close(fd);
            fd = 0;
        }
    }

    /// @note inspired by https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
    template <typename T, size_t SKIP = 4>
    class Iterator
    {
    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;  // or also value_type*
        using reference         = T&;  // or also value_type&

        Iterator(pointer ptr) : m_ptr(ptr) {}

        reference operator*() const { return *m_ptr; }
        pointer operator->() { return m_ptr; }

        // Prefix increment - jump 4 because there is pos, int a0, int a1, float pp
        Iterator& operator++() { m_ptr += SKIP; return *this; }

        // Postfix increment
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }

        friend bool operator== (const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; };
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; };

        static constexpr size_t skip() { return SKIP; }

    protected:
        pointer m_ptr;
    };

    class PositionContainer {
        using Iterator_type = Iterator<uint32_t, 4>;
    public:
        PositionContainer (HetInfoMemoryMap& parent, size_t sample) : parent(parent), sample(sample) {
            uint32_t mark = *(uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample]);
            if (mark != 0xd00dc0de) {
                std::cerr << "Wrong mark on sample index" << sample << std::endl;
                throw "Wrong mark";
            }
            size = *(((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 2);
            start_pos = ((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 3;
        }
        Iterator_type begin() { return Iterator(start_pos); }
        Iterator_type end()   { return Iterator(start_pos+size*Iterator_type::skip()); }

    protected:
        HetInfoMemoryMap& parent;
        uint32_t *start_pos;
        size_t sample;
        size_t size;
    };

    class HetInfoPtrContainer {
        using Iterator_type = Iterator<uint32_t, 4>;
    public:
        HetInfoPtrContainer (HetInfoMemoryMap& parent, size_t sample) : parent(parent), sample(sample) {
            uint32_t mark = *(uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample]);
            if (mark != 0xd00dc0de) {
                std::cerr << "Wrong mark on sample index" << sample << std::endl;
                throw "Wrong mark";
            }
            sample_id = *(((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 1);
            size = *(((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 2);
            start_pos = ((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 3;
        }

        void fill_het_info(std::vector<HetInfo>& v) {
            v.clear();
            for (size_t i = 0; i < size; ++i) {
                v.push_back(HetInfo(start_pos+i*Iterator_type::skip()));
            }
        }

    protected:
        HetInfoMemoryMap& parent;
        uint32_t *start_pos;
        size_t sample;
        uint32_t sample_id;
        size_t size;
    };

    void fill_het_info(std::vector<HetInfo>& v, size_t sample) {
        HetInfoPtrContainer hipc(*this, sample);
        hipc.fill_het_info(v);
    }

    void print_positions(size_t sample) {
        PositionContainer p(*this, sample);
        for (auto it = p.begin(); it != p.end(); it++) {
            std::cout << "Pos : " << *it << std::endl;
        }
        for (auto pos : p) {
            std::cout << "Pos : " << pos << std::endl;
        }
    }

    size_t count_all_vars_below_threshold(float threshold) {
        size_t counter = 0;
        // For each sample
        for (size_t i = 0; i < num_samples; ++i) {
            std::vector<HetInfo> v;
            // For each het site in dataset
            fill_het_info(v, i);
            for (auto hi : v) {
                if (!std::isnan(hi.pp) && hi.pp < threshold) {
                    counter++;
                }
            }
        }
        return counter;
    }

    size_t count_all_snps_below_threshold(float threshold, const std::vector<VarInfo>& vi) {
        size_t counter = 0;
        // For each sample
        for (size_t i = 0; i < num_samples; ++i) {
            std::vector<HetInfo> v;
            // For each het site in dataset
            fill_het_info(v, i);
            for (auto hi : v) {
                if (vi[hi.position].snp && !std::isnan(hi.pp) && hi.pp < threshold) {
                    counter++;
                }
            }
        }
        return counter;
    }

    size_t count_all_snps_below_threshold_that_can_be_linked_within(float threshold, const std::vector<VarInfo>& vi, size_t dist, bool snp_nei_req = true) {
        size_t counter = 0;
        // For each sample
        for (size_t i = 0; i < num_samples; ++i) {
            std::vector<HetInfo> v;
            // For each het site in dataset
            fill_het_info(v, i);
            for (size_t j = 0; j < v.size(); ++j) {
                auto& hi = v[j];
                // If low PP and SNP
                if (vi[hi.position].snp && !std::isnan(hi.pp) && hi.pp < threshold) {
                    // Check the neighbors
                    size_t lower_bound = j < 2 ? 0 : j-2;
                    size_t upper_bound = j < v.size()-2 ? j+2 : v.size();
                    // For each neighbor
                    for (size_t k = lower_bound; k <= upper_bound; ++k) {
                        // Don't check variant with itself
                        if (k == j) continue;
                        auto& nei = v[k];
                        // If the neighbor is also a SNP and within distance
                        if ((!snp_nei_req || vi[nei.position].snp) && vi[hi.position].distance(vi[nei.position]) < dist) {
                            counter++;
                            break; // Don't count twice
                        }
                    }
                }
            }
        }
        return counter;
    }

    int fd;
    size_t file_size;
    void *file_mmap_p;
    uint32_t num_samples;
    uint64_t *offset_table;
};

class PrintTraversal : public BcfTraversal {
public:
    PrintTraversal(HetInfoLoader& hil) : line_counter(0), entry(0), hil_ref(hil) {}

    virtual void handle_bcf_file_reader() override {
    }

    virtual void handle_bcf_line() override {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;

        if (entry >= hil_ref.his.size()) {
            stop = true;
            return;
        }

        // This is all very slow code, but it is just for testing

        int* ac = NULL;
        if (hil_ref.his[entry].position == line_counter) {
            std::string result(bcf_hdr_id2name(header, line->rid));
            result += "\t" + std::to_string(line->pos+1);
            result += "\t" + std::string(line->d.id);
            result += "\t" + std::string(line->d.allele[0]);
            result += "\t" + std::string(line->d.allele[1]);
            result += "\t" + std::string(".");
            result += "\t" + std::string(".");
            int nac = 0;
            int ret = bcf_get_info_int32(header, line, "AC", &ac, &nac);
            result += "\tAC=" + std::to_string(*ac);
            result += "\t" + hil_ref.his[entry].to_string();

            std::cout << result << std::endl;
            entry++;
        }
        if (ac) { free(ac); }

        line_counter++;
    }
protected:
    size_t line_counter;
    size_t entry;
    HetInfoLoader& hil_ref;
};

int main(int argc, char**argv) {
    CLI::App app{"PP Show utility app"};
    std::string filename = "-";
    app.add_option("-f,--file", filename, "Input file name");
    std::string bfname = "-";
    app.add_option("-b,--binary", bfname, "Binary file name");
    std::string ofname = "-";
    app.add_option("-o,--output", ofname, "Output file name");
    size_t sample = -1;
    app.add_option("-s,--sample", sample, "Sample index");

    CLI11_PARSE(app, argc, argv);

    if (!(sample+1)) {
        std::cerr << "Please provide sample indice" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bfname.compare("-") == 0) {
        std::cerr << "Requires binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }


#if 1

    Vars vars(filename);
    std::cout << "Num VCF lines : " << vars.vars.size() << std::endl;
    //for (auto & v : vars.vars) {
    //    std::cout << v.to_string() << std::endl;
    //}

    HetInfoMemoryMap himm(bfname);

    const float PP_THRESHOLD = 0.99;
    const size_t DIST_THRESHOLD = 750;
    auto low_pp_count = himm.count_all_vars_below_threshold(PP_THRESHOLD);
    auto low_pp_snps = himm.count_all_snps_below_threshold(PP_THRESHOLD, vars.vars);
    auto solvable_low_pp_snps = himm.count_all_snps_below_threshold_that_can_be_linked_within(PP_THRESHOLD, vars.vars, DIST_THRESHOLD);
    auto maybe_solvable_low_pp_snps = himm.count_all_snps_below_threshold_that_can_be_linked_within(PP_THRESHOLD, vars.vars, DIST_THRESHOLD, false);

    std::cout << "For " << himm.num_samples << " samples :" << std::endl;
    std::cout << "There are " << low_pp_count << " low PP variants" << std::endl;
    std::cout << "of which " << low_pp_snps << " are SNPs : "
              << low_pp_snps * 100.0 / low_pp_count << "%" << std::endl;
    std::cout << "of which " << solvable_low_pp_snps << " have a SNP neighbor within " << DIST_THRESHOLD << " base pairs : "
              << solvable_low_pp_snps * 100.0 / low_pp_snps << "%" << std::endl;
    std::cout << "of which " << solvable_low_pp_snps << " have a neighbor within " << DIST_THRESHOLD << " base pairs : "
              << maybe_solvable_low_pp_snps * 100.0 / low_pp_snps << "%" << std::endl;

    return 0;

    // Per sample info

    //himm.print_positions(sample);
    std::vector<HetInfo> v;
    himm.fill_het_info(v, sample);
    std::cout << "---" << std::endl;
    for (auto hi : v) {
        std::cout << vars.vars[hi.position].to_string() << "\t";
        std::cout << hi.to_string() << std::endl;
    }

#else
    std::fstream ifs(bfname, std::ios_base::binary | std::ios_base::in);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open binary file " << bfname << std::endl;
        exit(-1);
    }

    uint32_t endianness = 0;
    ifs.read(reinterpret_cast<char *>(&endianness), sizeof(endianness));
    if (endianness != 0xaabbccdd) {
        std::cerr << "File " << bfname << " has wrong endianness" << std::endl;
        exit(-1);
    }
    uint32_t num_samples = 0;
    ifs.read(reinterpret_cast<char *>(&num_samples), sizeof(num_samples));
    if (num_samples == 0) {
        std::cerr << "File " << bfname << " has no samples" << std::endl;
        exit(-1);
    }

    std::vector<uint64_t> offset_table(num_samples);
    for(size_t i = 0; i < offset_table.size(); ++i) {
        ifs.read(reinterpret_cast<char *>(&offset_table[i]), sizeof(uint64_t));
    }

    std::cout << "Sample " << sample << " offset : " << offset_table[sample] << std::endl;

    ifs.close();

    HetInfoLoader hil(bfname, sample);
    //for (auto& hi : hil.his) {
    //    std::cout << hi.to_string() << std::endl;
    //}

    PrintTraversal pt(hil);
    pt.traverse_no_unpack_no_destroy(filename);
    pt.destroy();

    return 0;

    HetInfoMemoryMap himm(bfname);
    //himm.print_positions(sample);
    std::vector<HetInfo> v;
    himm.fill_het_info(v, sample);
    std::cout << "---" << std::endl;
    for (auto hi : v) {
        std::cout << hi.to_string() << std::endl;
    }
#endif

    return 0;
}
