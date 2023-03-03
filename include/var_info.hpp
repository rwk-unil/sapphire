#ifndef __VAR_INFO_HPP__
#define __VAR_INFO_HPP__

#include <bcf_traversal.hpp>
#include <xcf.hpp>

#ifndef DIST
#define DIST(x,y) (std::max((x),(y))-std::min((x),(y)))
#endif

class VarInfo {
public:
    VarInfo(bcf1_t *line, bcf_hdr_t *hdr) {
        /// @todo This constructor is a bit slow because of the repeated alloc/free
        contig = bcf_hdr_id2name(hdr, line->rid);
        pos1 = line->pos; // Is 1 based and not 0 based
        id = std::string(line->d.id);
        ref = std::string(line->d.allele[0]);
        alt = std::string(line->d.allele[1]);
        snp = ref.length() == 1 && alt.length() == 1;
        int *acp = NULL;
        int nac = 0;
        bcf_unpack(line, BCF_UN_INFO);
        int ret = bcf_get_info_int32(hdr, line, "AC", &acp, &nac);
        if (ret < 0) {
            std::cerr << "Warning AC cannot be retrieved at chr : " << contig << " pos1 : " << pos1 << std::endl;
            ac = 0;
        } else {
            ac = *acp;
        }
        #if 0 /* Should be the same for all, so not interesting to extract */
        ret = bcf_get_info_int32(hdr, line, "AN", &acp, &nac);
        if (ret < 0) {
            std::cerr << "Warning AN cannot be retrieved at chr : " << contig << " pos1 : " << pos1 << std::endl;
            an = 0;
        } else {
            an = *acp;
        }
        #endif
        if (acp) { free(acp); }
    }

    VarInfo(const bcf_file_reader_info_t& bcf_fri) : VarInfo(bcf_fri.line, bcf_fri.sr->readers[0].header) {}

    std::string to_string() {
        std::string result(contig);
        result += "\t" + std::to_string(pos1+1); // Is 1 based and not 0 based
        result += "\t" + id;
        result += "\t" + ref;
        result += "\t" + alt;
        result += "\t.";
        result += "\t.";
        result += "\tAC=" + std::to_string(ac);

        return result;
    }

    inline size_t distance(const VarInfo& other) const {
        return DIST(pos1, other.pos1);
    }

    std::string contig;
    uint32_t pos1;
    std::string id;
    std::string ref;
    std::string alt;
    bool snp;
    uint32_t ac;
    uint32_t an;
};

class VarInfoTraversal : public BcfTraversal {
public:
    VarInfoTraversal(std::vector<VarInfo>& vars) : vars(vars) {}

    virtual void handle_bcf_file_reader() override {
    }

    virtual void handle_bcf_line() override {
        vars.push_back(VarInfo(bcf_fri));
    }
protected:
    std::vector<VarInfo>& vars;
};

class Vars {
public:
    Vars(std::string vcf_file) {
        VarInfoTraversal vit(vars);
        vit.traverse_no_unpack_no_destroy(vcf_file);
        vit.destroy();
    }

    std::vector<VarInfo> vars;
};

#endif /* __VAR_INFO_HPP__ */
