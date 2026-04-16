#ifndef UTILS_H
#define UTILS_H

#include <string>

class Utils {
public:
    static bool isValidDNA(char c); //to check if this letter is a valid DNA char or not
    static std::string toUpper(const std::string& str);
};

#endif
