/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
.txt input parsing: initial conditions for the generic driver.

Reads key=value files (see examples/): dt, steps, E, B, then one
particle after another, each with its own q, m, x, v lines (a q=/m=
line applies to the next particle defined by its x=/v= pair, so each
particle can be a different species). Unknown keys and missing or
invalid values fail with a descriptive error.
*/
#pragma once

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "boris_pusher.hpp"
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

// Parses "x y z" into a Vector3, failing clearly on bad input.
inline Vector3 parseVector3(const std::string &s, const std::string &key) {
    std::istringstream in(s);
    Vector3 v;
    if (!(in >> v.x >> v.y >> v.z)) {
        throw std::runtime_error("bad value for '" + key + "', expected three numbers");
    }
    return v;
}

struct Input {
    double dt = 0.0;
    int numSteps = 0;
    Vector3 E;
    Vector3 B;
    std::vector<Particle> particles;
};

struct KeyLine {
    std::string key;
    std::string value;
};

// Reads the next meaningful "key=value" line: strips '#' comments,
// skips blank lines, and fails on lines without '=' or with an empty
// value. Returns false at EOF. lineNo tracks physical lines for errors.
inline bool nextKeyLine(std::ifstream &file, const std::string &path, int &lineNo, KeyLine &out) {
    std::string line;
    while (std::getline(file, line)) {
        ++lineNo;
        // Strip comments.
        line = trim(line.substr(0, line.find('#')));
        if (line.empty()) {
            continue;
        }
        const std::string::size_type eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": expected key=value");
        }
        out.key = trim(line.substr(0, eq));
        out.value = trim(line.substr(eq + 1));
        if (out.value.empty()) {
            throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": empty value for '" + out.key + "'");
        }
        return true;
    }
    return false;
}

// Accumulates particles and global settings line by line: a q=/m=
// line applies to the next particle, and each particle needs its own
// pair -- carrying values over silently would mask a forgotten line
// with wrong physics instead of failing loudly.
struct InputBuilder {
    Input input;
    std::string path;
    int lineNo = 0;
    bool hasDt = false, hasSteps = false;
    bool hasE = false, hasB = false;
    Species pendingSpecies{0.0, 0.0};
    bool hasPendingCharge = false, hasPendingMass = false;
    bool lastParticleOpen = false;

    [[noreturn]] void fail(const std::string &message) const {
        throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": " + message);
    }

    // Opens a new particle from its x= line, carrying the pending q=/m=.
    void openParticle(const std::string &value) {
        if (lastParticleOpen) {
            fail("each particle needs one v= line before the next x=");
        }
        if (!hasPendingCharge || !hasPendingMass) {
            fail("each particle needs its own q= and m= lines before x=");
        }
        input.particles.push_back(Particle{parseVector3(value, "x"), Vector3{}, pendingSpecies});
        hasPendingCharge = false;
        hasPendingMass = false;
        lastParticleOpen = true;
    }

    // Closes the open particle with its v= line.
    void closeParticle(const std::string &value) {
        if (!lastParticleOpen) {
            fail("v= without a preceding x= for the same particle");
        }
        input.particles.back().velocity = parseVector3(value, "v");
        lastParticleOpen = false;
    }

    void applyLine(const std::string &key, const std::string &value) {
        if (key == "dt") {
            input.dt = std::stod(value);
            hasDt = true;
        } else if (key == "steps") {
            input.numSteps = std::stoi(value);
            hasSteps = true;
        } else if (key == "E") {
            input.E = parseVector3(value, key);
            hasE = true;
        } else if (key == "B") {
            input.B = parseVector3(value, key);
            hasB = true;
        } else if (key == "q") {
            pendingSpecies.charge = std::stod(value);
            hasPendingCharge = true;
        } else if (key == "m") {
            pendingSpecies.mass = std::stod(value);
            if (pendingSpecies.mass <= 0.0) {
                fail("m must be > 0");
            }
            hasPendingMass = true;
        } else if (key == "x") {
            openParticle(value);
        } else if (key == "v") {
            closeParticle(value);
        } else {
            fail("unknown key '" + key + "'");
        }
    }
};

// Checks required keys and value ranges once all lines are consumed.
inline void validateInput(const InputBuilder &b) {
    for (const char *k : {"dt", "steps", "E", "B"}) {
        const bool ok = (std::string(k) == "dt" && b.hasDt) || (std::string(k) == "steps" && b.hasSteps) ||
                        (std::string(k) == "E" && b.hasE) || (std::string(k) == "B" && b.hasB);
        if (!ok) {
            throw std::runtime_error(b.path + ": missing required key '" + std::string(k) + "'");
        }
    }
    if (b.input.dt <= 0.0) {
        throw std::runtime_error(b.path + ": dt must be > 0");
    }
    if (b.input.numSteps < 0) {
        throw std::runtime_error(b.path + ": steps must be >= 0");
    }
    if (b.input.particles.empty() || b.lastParticleOpen) {
        throw std::runtime_error(b.path + ": need at least one particle with x= and v= lines");
    }
}

// Reads the whole .txt into a validated Input.
inline Input parseInputFile(const std::string &path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open input file '" + path + "'");
    }
    InputBuilder builder;
    builder.path = path;
    KeyLine raw;
    while (nextKeyLine(file, path, builder.lineNo, raw)) {
        builder.applyLine(raw.key, raw.value);
    }
    validateInput(builder);
    return builder.input;
}
