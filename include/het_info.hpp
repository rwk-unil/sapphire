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
    HetInfo() : position(0), a0(0), a1(0), pp(NAN) {}
    HetInfo(int position, int a0, int a1, float pp) : position(position), a0(a0), a1(a1), pp(pp) {}
    int position;
    int a0;
    int a1;
    float pp;

    std::string to_string() const {
        std::string output("Position : ");
        output += std::to_string(position);
        output += std::string(" ") + std::to_string(bcf_gt_allele(a0)) + (bcf_gt_is_phased(a0) ? "|" : "/") + std::to_string(bcf_gt_allele(a1));
        output += std::string(" PP : ") + std::to_string(pp);

        return output;
    }

    void to_stream(std::fstream& ofs) const {
        ofs.write(reinterpret_cast<const char*>(&position), sizeof(position));
        ofs.write(reinterpret_cast<const char*>(&a0), sizeof(a0));
        ofs.write(reinterpret_cast<const char*>(&a1), sizeof(a1));
        ofs.write(reinterpret_cast<const char*>(&pp), sizeof(pp));
    }

    void from_stream(std::ifstream& ifs) {
        ifs.read(reinterpret_cast<char*>(&position), sizeof(position));
        ifs.read(reinterpret_cast<char*>(&a0), sizeof(a0));
        ifs.read(reinterpret_cast<char*>(&a1), sizeof(a1));
        ifs.read(reinterpret_cast<char*>(&pp), sizeof(pp));
    }
};

#endif /* __HET_INFO_HPP__ */
