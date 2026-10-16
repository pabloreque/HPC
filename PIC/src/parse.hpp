/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
Parsing helpers shared by the input files (key=value input, profile file).
*/
#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include "vector3.hpp"

// Throws invalid_argument with 'message' if 'condition' is false.
inline void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

// Formats "path:lineNo: message" for errors tied to one input line.
inline std::string lineError(const std::string &path, int lineNo,
                             const std::string &message) {
    return path + ":" + std::to_string(lineNo) + ": " + message;
}

// Trims leading and trailing whitespace (spaces, tabs, newlines).
inline std::string trim(const std::string &text) {
    const std::string::size_type first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const std::string::size_type last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

// True if only whitespace is left in the stream (no trailing garbage).
inline bool fullyConsumed(std::istringstream &in) {
    return (in >> std::ws).eof();
}

// Parses one real number. Throws invalid_argument on bad input.
inline double parseDouble(const std::string &text, const std::string &key) {
    std::istringstream in(text);
    double value;
    if (!(in >> value) || !fullyConsumed(in)) {
        throw std::invalid_argument("bad value for '" + key +
                                    "', expected a number");
    }
    return value;
}

// Parses one integer. Throws invalid_argument on bad input.
inline int parseInt(const std::string &text, const std::string &key) {
    std::istringstream in(text);
    int value;
    if (!(in >> value) || !fullyConsumed(in)) {
        throw std::invalid_argument("bad value for '" + key +
                                    "', expected an integer");
    }
    return value;
}

// Parses "x y z" into a Vector3. Throws invalid_argument on bad input.
inline Vector3 parseVector3(const std::string &text, const std::string &key) {
    std::istringstream in(text);
    Vector3 value;
    if (!(in >> value.x >> value.y >> value.z) || !fullyConsumed(in)) {
        throw std::invalid_argument("bad value for '" + key +
                                    "', expected three numbers");
    }
    return value;
}

// Drops a trailing '#' comment and the surrounding whitespace.
inline std::string stripComment(const std::string &line) {
    return trim(line.substr(0, line.find('#')));
}
