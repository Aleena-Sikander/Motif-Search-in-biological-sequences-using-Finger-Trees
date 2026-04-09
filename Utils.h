#ifndef UTILS_H
#define UTILS_H

#include <string>

class Utils {
public:
    static bool isValidDNA(char c);
    static std::string toUpper(const std::string& str);
};

#endif
