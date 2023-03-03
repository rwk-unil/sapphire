#ifndef __SAMPLE_INFO_HPP__
#define __SAMPLE_INFO_HPP__

#include <fstream>
#include <string>

class SampleInfoLoader {
public:
    SampleInfoLoader(const std::string& sample_filename) {
        std::ifstream file(sample_filename);
        std::string line;
        while (std::getline(file, line)) {
            if (line.length()) {
                sample_names.push_back(line);
            }
        }
    }

    std::vector<std::string> sample_names;
};

#endif /* __SAMPLE_INFO_HPP__ */