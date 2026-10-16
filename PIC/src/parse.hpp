/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
Value parsing helpers for the key=value input files.
*/
#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include "vector3.hpp"

// Trims leading and trailing whitespace (spaces, tabs, newlines).
inline std::string trim(const std::string &s) {
    const std::string::size_type first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const std::string::size_type last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

// True if only whitespace is left in the stream (no trailing garbage).
inline bool fullyConsumed(std::istringstream &in) {
    return (in >> std::ws).eof();
}

// Parses one real number. Throws invalid_argument on bad input.
inline double parseDouble(const std::string &s, const std::string &key) {
    std::istringstream in(s);
    double value;
    if (!(in >> value) || !fullyConsumed(in)) {
        throw std::invalid_argument("bad value for '" + key +
                                    "', expected a number");
    }
    return value;
}

// Parses one integer. Throws invalid_argument on bad input.
inline int parseInt(const std::string &s, const std::string &key) {
    std::istringstream in(s);
    int value;
    if (!(in >> value) || !fullyConsumed(in)) {
        throw std::invalid_argument("bad value for '" + key +
                                    "', expected an integer");
    }
    return value;
}

// Parses "x y z" into a Vector3. Throws invalid_argument on bad input.
inline Vector3 parseVector3(const std::string &s, const std::string &key) {
    std::istringstream in(s);
    Vector3 v;
    if (!(in >> v.x >> v.y >> v.z) || !fullyConsumed(in)) {
        throw std::invalid_argument("bad value for '" + key +
                                    "', expected three numbers");
    }
    return v;
}
