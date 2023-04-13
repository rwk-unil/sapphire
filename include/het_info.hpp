#ifndef __HET_INFO_HPP__
#define __HET_INFO_HPP__

#include <string>
#include "vcf.h"

class HetInfo {
public:
    HetInfo(std::ifstream& ifs) {
        from_stream(ifs);
    }
    HetInfo(uint32_t* ptr) : HetInfo(*ptr, *(ptr+1), *(ptr+2), *(float*)(ptr+3)) {}
    HetInfo() : vcf_line(0), a0(0), a1(0), pp(NAN) {}
    HetInfo(int vcf_line, int a0, int a1, float pp) : vcf_line(vcf_line), a0(a0), a1(a1), pp(pp) {}
    int vcf_line;
    int a0;
    int a1;
    float pp;

    std::string to_string() const {
        std::string output("Position : ");
        output += std::to_string(vcf_line);
        output += std::string(" ") + std::to_string(bcf_gt_allele(a0)) + (bcf_gt_is_phased(a1) ? "|" : "/") + std::to_string(bcf_gt_allele(a1));
        output += std::string(" PP : ") + std::to_string(pp);

        return output;
    }

    void to_stream(std::fstream& ofs) const {
        ofs.write(reinterpret_cast<const char*>(&vcf_line), sizeof(vcf_line));
        ofs.write(reinterpret_cast<const char*>(&a0), sizeof(a0));
        ofs.write(reinterpret_cast<const char*>(&a1), sizeof(a1));
        ofs.write(reinterpret_cast<const char*>(&pp), sizeof(pp));
    }

    void from_stream(std::ifstream& ifs) {
        ifs.read(reinterpret_cast<char*>(&vcf_line), sizeof(vcf_line));
        ifs.read(reinterpret_cast<char*>(&a0), sizeof(a0));
        ifs.read(reinterpret_cast<char*>(&a1), sizeof(a1));
        ifs.read(reinterpret_cast<char*>(&pp), sizeof(pp));
    }
};

constexpr bool operator==(const HetInfo& lhs, const HetInfo& rhs) {
    return lhs.vcf_line == rhs.vcf_line &&
           lhs.a0 == rhs.a0 && lhs.a1 == rhs.a1 &&
           ((std::isnan(lhs.pp) && std::isnan(rhs.pp)) || (lhs.pp == rhs.pp));
}

constexpr bool operator!=(const HetInfo& lhs, const HetInfo& rhs) {
    return !(lhs == rhs);
}

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

#endif /* __HET_INFO_HPP__ */
