#pragma once
#include <string>
#include <fstream>
#include <stdexcept>

// reads a FASTA file and returns just the sequence as a plain string
// skips header lines (they start with '>') and any whitespace
inline std::string read_fasta(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    std::string sequence, line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;  // skip headers
        for (char c : line)
            if (c != '\r' && c != ' ') sequence += c;
    }
    return sequence;
}
