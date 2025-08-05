#ifndef __HET_INFO_LOADER_HPP__
#define __HET_INFO_LOADER_HPP__

#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sys/mman.h>

#include "fs.hpp"
#include "het_info.hpp"
#include "var_info.hpp"

const uint32_t ENDIANNESS = 0xaabbccdd;

class HetInfoMemoryMap {
public:
    HetInfoMemoryMap(std::string bfname) : HetInfoMemoryMap(bfname, PROT_READ) {}
    HetInfoMemoryMap(std::string bfname, int mmflags) : file_size(fs::file_size(bfname)) {
        fd = open(bfname.c_str(), (mmflags & PROT_WRITE) ? O_RDWR : O_RDONLY, 0);
        if (fd < 0) {
            std::cerr << "Failed to open file : " << bfname << std::endl;
            throw "Failed to open file";
        }

        file_mmap_p = mmap(NULL, file_size, mmflags, MAP_SHARED, fd, 0);
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
            msync(file_mmap_p, file_size, MS_SYNC);
            munmap(file_mmap_p, file_size);
            file_mmap_p = NULL;
        }
        if (fd > 0) {
            close(fd);
            fd = 0;
        }
    }

    uint32_t *get_ptr_on_nth(uint32_t n) const {
        uint32_t *start = ((uint32_t*)(((char*)file_mmap_p) + offset_table[n]));
        if (*start != 0xd00dc0de) {
            std::cerr << "Something is wrong, mark not found for idx " << n << std::endl;
            return nullptr;
        } else {
            return start;
        }
    }

    std::vector<HetInfo> get_het_info_for_nth(uint32_t n) const {
        std::vector<HetInfo> his;
        auto start = get_ptr_on_nth(n);
        auto size = *(start + 2); /* Skip Mark, id */
        auto p = start + 3; /* Skip Mark, id, size */
        for (auto todo = size; todo > 0; todo--) {
            his.push_back(HetInfo(p));
            p += 4; /* Size of HetInfo */
        }
        return his;
    }

    uint32_t get_size_of_nth(uint32_t n) const {
        const auto start = get_ptr_on_nth(n);
        return *(start+2) * sizeof(uint32_t) * 4 /* Size of HetInfo */ + 3 * sizeof(uint32_t); /* Mark, id, size */
    }

    uint32_t get_orig_idx_of_nth(uint32_t n) const {
        const auto start = get_ptr_on_nth(n);
        return *(start+1);
    }

    std::map<uint32_t, uint32_t> get_orig_idx_to_nth_map() const {
        std::map<uint32_t, uint32_t> map;
        for (uint32_t i = 0; i < num_samples; ++i) {
            map.insert({get_orig_idx_of_nth(i), i});
        }
        return map;
    }

    bool integrity_check_pass() const {
        size_t computed_size = (num_samples + 1) * sizeof(uint64_t);

        bool pass = true;
        for (size_t i = 0; i < num_samples-1; ++i) {
            computed_size += get_size_of_nth(i);
            auto p1 = get_ptr_on_nth(i);
            auto p2 = get_ptr_on_nth(i+1);
            if (*p1 != 0xd00dc0de) pass = false;
            if (*p2 != 0xd00dc0de) pass = false;
            if (((uint32_t*)((uint8_t*)p1 + get_size_of_nth(i))) != p2) {
                std::cerr << "Offset " << i << "and next one don't match up !" << std::endl;
                pass = false;
            }
        }
        computed_size += get_size_of_nth(num_samples-1);
        if (computed_size != file_size) {
            std::cerr << "File has different size than it should be" << std::endl;
            pass = false;
        }
        return pass;
    }

    void show_info() const {
        std::cout << "File size : " << file_size << std::endl;
        std::cout << "Number of samples : " << num_samples << std::endl;
        for (size_t i = 0; i < num_samples; ++i) {
            std::cout << "Size of sample " << i << " : " << get_size_of_nth(i) << std::endl;
        }
        bool pass = integrity_check_pass();
        std::cout << "File passes integrity check : " << (pass ? "YES" : "NO") << std::endl;
    }

    void write_sub_file(const std::vector<uint32_t> ids_to_extract, const std::string& filename) const {
        std::fstream ofs(filename, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        if (!ofs.is_open()) {
            std::cerr << "Cannot open file " << filename << std::endl;
            throw "Cannot open file";
        }

        const uint32_t endianness = 0xaabbccdd;
        // Write endianness
        ofs.write(reinterpret_cast<const char*>(&endianness), sizeof(uint32_t));
        // Write number of samples
        uint32_t n = ids_to_extract.size();
        ofs.write(reinterpret_cast<const char*>(&n), sizeof(uint32_t));
        // Write offset table
        std::vector<uint64_t> new_offset_table(n, 0xdeadc0dedeadc0de);
        uint64_t dummy_offset = 0xdeadc0dedeadc0de;
        auto table_seek = ofs.tellp();
        for (uint32_t i = 0; i < n; ++i) {
            ofs.write(reinterpret_cast<const char*>(&dummy_offset), sizeof(uint64_t));
        }

        for (uint32_t i = 0; i < n; ++i) {
            uint32_t *start = ((uint32_t*)(((char*)file_mmap_p) + offset_table[ids_to_extract[i]]));
            if (*start != 0xd00dc0de) {
                std::cerr << "Something is wrong, mark not found for idx " << i << std::endl;
            } else {
                uint32_t id = *(start+1);
                uint32_t size = *(start+2);
                if (id != ids_to_extract[i]) {
                    std::cerr << "Trying to extract id " << ids_to_extract[i] << " but found id " << id << std::endl;
                }
                new_offset_table[i] = ofs.tellp();
                ofs.write(reinterpret_cast<const char*>(start), size * sizeof(uint32_t) * 4 /* Size of HetInfo */ + 3 * sizeof(uint32_t) /* Mark, id, size */);
            }
        }

        // Rewrite the offset table
        ofs.seekp(table_seek);
        for (const auto& offset : new_offset_table) {
            ofs.write(reinterpret_cast<const char*>(&offset), sizeof(decltype(offset)));
        }

        ofs.close();
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
    public:
        using Iterator_type = Iterator<uint32_t, 4>;
        HetInfoPtrContainer (HetInfoMemoryMap& parent, size_t sample) : parent(parent), sample(sample),
            sample_id(*(((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 1)),
            size(*(((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 2)),
            start_pos(((uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample])) + 3)
        {
            uint32_t mark = *(uint32_t*)(((char*)parent.file_mmap_p) + parent.offset_table[sample]);
            if (mark != 0xd00dc0de) {
                std::cerr << "Wrong mark on sample index" << sample << std::endl;
                throw "Wrong mark";
            }
        }

        void fill_het_info(std::vector<HetInfo>& v) {
            v.clear();
            for (size_t i = 0; i < size; ++i) {
                v.push_back(HetInfo(start_pos+i*Iterator_type::skip()));
            }
        }

        /* Switch phase of all genotypes, this is for switching vertically split
         * files that require to be switched during ligation, this should be
         * extremely fast */
        void switch_phase() {
            /* The reason we use "int" here is to be coherent with HTSLIB */
            static_assert(sizeof(int) == sizeof(uint32_t), "int should be of same size than uint32_t");
            for (size_t i = 0; i < size; ++i) {
                int *gt_arr = (int*)start_pos+i*Iterator_type::skip()+1;
                int a0 = bcf_gt_allele(gt_arr[0]);
                int a1 = bcf_gt_allele(gt_arr[1]);

                /* First allele is always unphased per VCF/BCF standard.
                 * Therefore, we can't just switch the values directly */
                gt_arr[0] = bcf_gt_unphased(a1);
                gt_arr[1] = bcf_gt_phased(a0);
            }
        }

    protected:
        HetInfoMemoryMap& parent;
    public:
        const size_t sample;
        const uint32_t sample_id;
        const size_t size;
        uint32_t* const start_pos;
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
                if (vi[hi.vcf_line].snp && !std::isnan(hi.pp) && hi.pp < threshold) {
                    counter++;
                }
            }
        }
        return counter;
    }

    size_t count_all_snps_with_pir(float threshold, const std::vector<VarInfo>& vi) {
        size_t counter = 0;
        // For each sample
        for (size_t i = 0; i < num_samples; ++i) {
            std::vector<HetInfo> v;
            // For each het site in dataset
            fill_het_info(v, i);
            for (auto hi : v) {
                if (vi[hi.vcf_line].snp && !std::isnan(hi.pp) && hi.pp < threshold) {
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
                if (vi[hi.vcf_line].snp && !std::isnan(hi.pp) && hi.pp < threshold) {
                    // Check the neighbors
                    size_t lower_bound = j < 2 ? 0 : j-2;
                    size_t upper_bound = j < v.size()-2 ? j+2 : v.size();
                    // For each neighbor
                    for (size_t k = lower_bound; k < upper_bound; ++k) {
                        // Don't check variant with itself
                        if (k == j) continue;
                        auto& nei = v[k];
                        // If the neighbor is also a SNP and within distance
                        if ((!snp_nei_req || vi[nei.vcf_line].snp) && vi[hi.vcf_line].distance(vi[nei.vcf_line]) < dist) {
                            counter++;
                            break; // Don't count twice
                        }
                    }
                }
            }
        }
        return counter;
    }

    size_t count_rephased_variants() {
        size_t counter = 0;
        // For each sample
        for (size_t i = 0; i < num_samples; ++i) {
            std::vector<HetInfo> v;
            // For each het site in dataset
            fill_het_info(v, i);
            for (size_t j = 0; j < v.size(); ++j) {
                auto& hi = v[j];
                // If SNP and PP above 1 (means rephase with sequencing reads)
                if (!std::isnan(hi.pp) && hi.pp > 1.0) {
                    counter++;
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

class HetInfoMemoryMapMerger {
public:

    static uint32_t get_num_samples_from_filenames(const std::vector<std::string>& filenames) {
        uint32_t total = 0;
        for (const auto& f : filenames) {
            HetInfoMemoryMap himm(f);
            total += himm.num_samples;
        }
        return total;
    }

    HetInfoMemoryMapMerger(std::string ofname, uint32_t num_samples) :
        current_sample(0),
        offset_table(num_samples, 0xdeadc0dedeadc0de),
        ofs(ofname, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc)
    {
        if (!ofs.is_open()) {
            std::cerr << "Cannot open file " << ofname << std::endl;
        }

        const uint32_t endianness = 0xaabbccdd;
        // Write endianness
        ofs.write(reinterpret_cast<const char*>(&endianness), sizeof(uint32_t));
        // Write number of samples
        ofs.write(reinterpret_cast<const char*>(&num_samples), sizeof(uint32_t));
        // Write offset table
        table_seek = ofs.tellp();
        for (const auto& offset : offset_table) {
            ofs.write(reinterpret_cast<const char*>(&offset), sizeof(decltype(offset)));
        }
    }

    void merge(const std::string& filename) {
        HetInfoMemoryMap himm(filename);
        if (!himm.integrity_check_pass()) {
            std::cerr << "File " << filename << " doesn't pass integrity checks" << std::endl;
        }
        for (size_t i = 0; i < himm.num_samples; ++i) {
            if (current_sample > offset_table.size()) {
                std::cerr << "Number of samples is not enough for merge" << std::endl;
                return;
            }
            offset_table[current_sample] = ofs.tellp();
            ofs.write(reinterpret_cast<const char*>(himm.get_ptr_on_nth(i)), himm.get_size_of_nth(i));
            current_sample++;
        }
    }

    virtual ~HetInfoMemoryMapMerger() {
        /* Update the offset table before closing the file */
        ofs.seekp(table_seek);
        for (const auto& offset : offset_table) {
            ofs.write(reinterpret_cast<const char*>(&offset), sizeof(decltype(offset)));
        }
        ofs.close();
    }

protected:
    std::streampos table_seek;
    uint32_t current_sample;
    std::vector<uint64_t> offset_table;
    std::fstream ofs;
};

#endif /* __HET_INFO_LOADER_HPP__ */
