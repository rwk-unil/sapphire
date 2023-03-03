#ifndef __HET_INFO_LOADER_HPP__
#define __HET_INFO_LOADER_HPP__

#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sys/mman.h>

#include "compression.hpp"
#include "fs.hpp"
#include "het_info.hpp"
#include "var_info.hpp"

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
                    for (size_t k = lower_bound; k <= upper_bound; ++k) {
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

    int fd;
    size_t file_size;
    void *file_mmap_p;
    uint32_t num_samples;
    uint64_t *offset_table;
};

#endif /* __HET_INFO_LOADER_HPP__ */
