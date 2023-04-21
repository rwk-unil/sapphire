#ifndef __SAMPLE_INFO_HPP__
#define __SAMPLE_INFO_HPP__

#include <fstream>
#include <string>
#include <sstream>

class SampleInfo {
public:
    SampleInfo(size_t index, const std::string& name) :
        original_index(index),
        name(name),
        cram_file_path("")
    {}

    SampleInfo(size_t index, const std::string& name, const std::string& cram_file_path) :
        original_index(index),
        name(name),
        cram_file_path(cram_file_path)
    {}

    std::string to_string() const {
        std::stringstream ss;
        ss << "Sample idx: " << original_index << " name: " << name << " cram path: " << cram_file_path;
        return ss.str();
    }

    size_t original_index = -1;
    std::string name;
    std::string cram_file_path;
};

std::vector<std::string> tokenize(const std::string& s, const std::string& delim) {
    size_t pos_start = 0, pos_end, delim_len = delim.length();
    std::string token;
    std::vector<std::string> result;

    while ((pos_end = s.find(delim, pos_start)) != std::string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        result.push_back(token);
    }
    token = s.substr(pos_start);
    if (token.length()) {
        result.push_back(token);
    }

    return result;
}

class SampleInfoLoader {
public:
    SampleInfoLoader(const std::string& sample_filename) {
        if (sample_filename.compare("-") == 0) return;
        std::ifstream file(sample_filename);
        std::string line;
        size_t i = 0;
        //std::cout << "File : " << sample_filename << std::endl;
        while (std::getline(file, line)) {
            //std::cout << "line : " << line << std::endl;
            if (line.length()) {
                auto tokens = tokenize(line, ",");
                //std::cout << "Tokens size : " << tokens.size() << std::endl;
                if (tokens.size() == 1) {
                    //std::cout << "Pushing line : " << line << std::endl;
                    sample_names.push_back(line); // Same as tokens [0]
                    samples.push_back(SampleInfo(i, line));
                } else if (tokens.size() > 1) {
                    sample_names.push_back(tokens[1]);
                    const std::string& cram_path = (tokens.size() > 2) ? tokens[2] : std::string("");
                    samples.push_back(SampleInfo(std::stoi(tokens[0]), tokens[1], cram_path));
                } else {
                    std::cerr << "Could not load line : '" << line << "' from file" << sample_filename << std::endl;
                    /// @todo throw
                }
                i++;
            }
        }
        //std::cout << "Samples size : " << samples.size() << std::endl;
        //std::cout << "Found samples : " << std::endl;
        for (auto& s : samples) std::cout << s.to_string() << std::endl;
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
    std::vector<SampleInfo> samples;
};

#endif /* __SAMPLE_INFO_HPP__ */