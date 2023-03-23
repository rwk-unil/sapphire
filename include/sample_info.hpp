#ifndef __SAMPLE_INFO_HPP__
#define __SAMPLE_INFO_HPP__

#include <fstream>
#include <string>

class SampleInfoLoader {
public:
    SampleInfoLoader(const std::string& sample_filename) {
        if (sample_filename.compare("-") == 0) return;
        std::ifstream file(sample_filename);
        std::string line;
        while (std::getline(file, line)) {
            if (line.length()) {
                sample_names.push_back(line);
            }
        }
    }

    static std::vector<size_t> ids_from_files(std::string& samples_fname, std::string& sub_fname) {
        SampleInfoLoader sil_full(samples_fname);
        SampleInfoLoader sil_sub(sub_fname);

        std::vector<size_t> ids;

        for (auto& s : sil_sub.sample_names) {
            auto it = std::find(sil_full.sample_names.begin(), sil_full.sample_names.end(), s);
            if (it == sil_full.sample_names.end()) {
                std::cerr << "Sample with name " << s << " not found in " << samples_fname << std::endl;
            } else {
                ids.push_back(it-sil_full.sample_names.begin());
            }
        }

        return ids;
    }

    std::vector<std::string> sample_names;
};

#endif /* __SAMPLE_INFO_HPP__ */