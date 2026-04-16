#include "Utils.h"

bool Utils::isValidDNA(char c) {
    char upperC = std::toupper(c);
    return upperC == 'A' || upperC == 'C' || upperC == 'G' || upperC == 'T';
}

std::string Utils::toUpper(const std::string& str) {
    std::string result;
    for (char c : str) {
        result += std::toupper(c);
    }
    return result;
}