#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sys/mman.h>

#include "bcf_traversal.hpp"
#include "CLI11.hpp"
#include "hts.h"
#include "het_info.hpp"
#include "het_info_loader.hpp"
#include "var_info.hpp"
#include "sample_info.hpp"
#include "synced_bcf_reader.h"
#include "vcf.h"

class Data {
public:
    Data(size_t vcf_line, float pp, bool switched, size_t num_pir, int a0, int a1) :
        vcf_line(vcf_line), pp(pp), switched(switched), num_pir(num_pir),
        a0(bcf_gt_allele(a0)), a1(bcf_gt_allele(a1))
    {}

    std::string to_string() const {
        std::stringstream ss;
        ss << pp << "," << switched << "," << num_pir;
        return ss.str();
    }

    static std::string csv_header() {
        return std::string("PP, Switched, Num PIR");
    }

    size_t vcf_line;
    float pp;
    bool switched;
    size_t num_pir;
    int a0;
    int a1;
};

class BinDiffExtractor {
public:
    BinDiffExtractor(std::string& bin1_fname, std::string& bin2_fname,
        std::string& vcf_fname, std::string& samples_fname,
        std::string& sub_fname) :
        bin1_fname(bin1_fname), bin2_fname(bin2_fname),
        vcf_fname(vcf_fname), samples_fname(samples_fname),
        sub_fname(sub_fname),
        himm_original(bin1_fname),
        himm_rephased(bin2_fname) {
        if (sub_fname.compare("-") == 0) {
            ids.resize(himm_original.num_samples);
            std::iota(ids.begin(), ids.end(), 0);
        } else {
            ids = SampleInfoLoader::ids_from_files(samples_fname, sub_fname);
        }
    }

    void fill_internal_data_from_files() {
        /* For all samples */
        for (const auto& sample_idx : ids) {
            data.push_back({});

            HetInfoMemoryMap::HetInfoPtrContainer orig_hipc(himm_original, sample_idx);
            HetInfoMemoryMap::HetInfoPtrContainer reph_hipc(himm_rephased, sample_idx);

            std::vector<HetInfo> orig_hi;
            std::vector<HetInfo> reph_hi;
            orig_hipc.fill_het_info(orig_hi);
            reph_hipc.fill_het_info(reph_hi);

            if (orig_hi.size() != reph_hi.size()) {
                std::cerr << "Binary files do not match !" << std::endl;
                exit(-1);
            }

            for (size_t i = 0; i < orig_hi.size(); ++i) {
                if (!std::isnan(orig_hi.at(i).pp) && orig_hi.at(i).pp < 0.99) {
                    data.back().emplace_back(Data(
                        orig_hi.at(i).vcf_line,
                        orig_hi.at(i).pp,
                        orig_hi.at(i).a0 != reph_hi.at(i).a0,
                        (reph_hi.at(i).pp > 1.0) ? (size_t)(reph_hi.at(i).pp - 1.0) : 0,
                        reph_hi.at(i).a0,
                        reph_hi.at(i).a1
                    ));
                }
            }
        }
    }

    void print_csv(bool extra, bool more, bool ac) {
        if (extra) {
            if (vcf_fname.compare("-") == 0) {
                std::cerr << "Requires variant VCF/BCF filename\n";
                exit(-1);
            }
            if (samples_fname.compare("-") == 0) {
                std::cerr << "Requires samples filename\n";
                exit(-1);
            }

            VarInfoLoader vil(vcf_fname);
            SampleInfoLoader sil_full(samples_fname);

            std::cout << "Sample name, " <<  Data::csv_header() << ", is SNP" << (ac ? ", AC" : "") << (more ? ", VCF line" : "") << std::endl;

            for (size_t i = 0; i < ids.size(); ++i) {
                for (auto& e : data[i]) {
                    std::cout << sil_full.sample_names[ids[i]] << "," << e.to_string() << "," << vil.vars[e.vcf_line].snp <<
                    (ac ? std::string(",") + std::to_string(vil.vars[e.vcf_line].ac) : "") <<
                    (more ? std::string(",") + vil.vars[e.vcf_line].to_string() : "") << std::endl;
                }
            }
        } else {
            std::cout << Data::csv_header() << std::endl;

            for (auto& v : data) {
                for (auto& e : v) {
                    std::cout << e.to_string() << std::endl;
                }
            }
        }
    }

    void print_csv_filter_ac(size_t ac_threshold) {
        constexpr bool ac = true;
        constexpr bool more = true;
        if (vcf_fname.compare("-") == 0) {
            std::cerr << "Requires variant VCF/BCF filename\n";
            exit(-1);
        }
        if (samples_fname.compare("-") == 0) {
            std::cerr << "Requires samples filename\n";
            exit(-1);
        }

        VarInfoLoader vil(vcf_fname);
        SampleInfoLoader sil_full(samples_fname);

        std::cout << "Sample name, " <<  Data::csv_header() << ", is SNP" << (ac ? ", AC" : "") << ", GT" << (more ? ", VCF line" : "") << std::endl;

        for (size_t i = 0; i < ids.size(); ++i) {
            for (auto& e : data[i]) {
                if (vil.vars[e.vcf_line].ac <= ac_threshold) {
                    std::cout << sil_full.sample_names[ids[i]] << "," << e.to_string() << "," << vil.vars[e.vcf_line].snp <<
                    (ac ? std::string(",") + std::to_string(vil.vars[e.vcf_line].ac) : "") <<
                    std::string(",") + std::to_string(e.a0) + "|" + std::to_string(e.a1) <<
                    (more ? std::string(",") + vil.vars[e.vcf_line].to_string() : "") << std::endl;
                }
            }
        }
    }

    std::string bin1_fname;
    std::string bin2_fname;
    std::string vcf_fname;
    std::string samples_fname;
    std::string sub_fname;
    HetInfoMemoryMap himm_original;
    HetInfoMemoryMap himm_rephased;
    std::vector<std::vector<Data> > data;
    std::vector<size_t> ids;
};

int main(int argc, char**argv) {
    CLI::App app{"Binary file diff utility app"};
    std::string vcf_fname = "-";
    app.add_option("-f,--vcf-file", vcf_fname, "Variant file name (needed for extra info)");
    std::string bin1_fname = "-";
    std::string bin2_fname = "-";
    app.add_option("-a,--binary1", bin1_fname, "Original binary file name");
    app.add_option("-b,--binary2", bin2_fname, "Rephased binary file name");
    std::string samples_fname = "-";
    app.add_option("-S,--samples", samples_fname, "Sample names (text file) as appear in full BCF");
    std::string sub_fname = "-";
    app.add_option("-l,--sample-list", sub_fname, "Unordered sub sample list (text file)");
    bool extra = false;
    app.add_flag("--extra-info", extra, "Output extra information");
    bool more = false;
    app.add_flag("--more", more, "Output more information");
    bool ac = false;
    app.add_flag("--ac", ac, "Added allele count field");
    size_t ac_threshold = 0;
    app.add_option("--ac-threshold", ac_threshold, "AC threshold filtering value (0 means no filtering)");

    CLI11_PARSE(app, argc, argv);

    if (bin1_fname.compare("-") == 0) {
        std::cerr << "Requires original binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (bin2_fname.compare("-") == 0) {
        std::cerr << "Requires rephased binary filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    BinDiffExtractor bde(bin1_fname, bin2_fname, vcf_fname, samples_fname, sub_fname);

    bde.fill_internal_data_from_files();

    if (ac_threshold) {
        bde.print_csv_filter_ac(ac_threshold);
    } else {
        bde.print_csv(extra, more, ac);
    }

    return 0;
}
