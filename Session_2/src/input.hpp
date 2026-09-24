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

inline std::string trim(const std::string &s) {
    const std::string::size_type first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const std::string::size_type last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

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

inline Input parseInputFile(const std::string &path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open input file '" + path + "'");
    }
    Input input;
    bool hasDt = false, hasSteps = false;
    bool hasE = false, hasB = false;
    // Species of the particle currently being defined: a q=/m= line
    // applies to the next particle, and each particle needs its own
    // pair -- carrying values over silently would mask a forgotten line
    // with wrong physics instead of failing loudly.
    Species pendingSpecies{0.0, 0.0};
    bool hasPendingCharge = false, hasPendingMass = false;
    bool lastParticleOpen = false;

    std::string line;
    int lineNo = 0;
    auto fail = [&](const std::string &message) {
        throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": " + message);
    };
    while (std::getline(file, line)) {
        ++lineNo;
        // Strip comments.
        line = trim(line.substr(0, line.find('#')));
        if (line.empty()) {
            continue;
        }
        const std::string::size_type eq = line.find('=');
        if (eq == std::string::npos) {
            fail("expected key=value");
        }
        const std::string key = trim(line.substr(0, eq));
        const std::string value = trim(line.substr(eq + 1));
        auto needValue = [&](const std::string &k) {
            if (value.empty()) {
                fail("empty value for '" + k + "'");
            }
        };
        if (key == "dt") {
            needValue(key);
            input.dt = std::stod(value);
            hasDt = true;
        } else if (key == "steps") {
            needValue(key);
            input.numSteps = std::stoi(value);
            hasSteps = true;
        } else if (key == "E") {
            needValue(key);
            input.E = parseVector3(value, key);
            hasE = true;
        } else if (key == "B") {
            needValue(key);
            input.B = parseVector3(value, key);
            hasB = true;
        } else if (key == "q") {
            needValue(key);
            pendingSpecies.charge = std::stod(value);
            hasPendingCharge = true;
        } else if (key == "m") {
            needValue(key);
            pendingSpecies.mass = std::stod(value);
            if (pendingSpecies.mass <= 0.0) {
                fail("m must be > 0");
            }
            hasPendingMass = true;
        } else if (key == "x") {
            needValue(key);
            if (lastParticleOpen) {
                fail("each particle needs one v= line before the next x=");
            }
            if (!hasPendingCharge || !hasPendingMass) {
                fail("each particle needs its own q= and m= lines before x=");
            }
            input.particles.push_back(Particle{parseVector3(value, key), Vector3{}, pendingSpecies});
            hasPendingCharge = false;
            hasPendingMass = false;
            lastParticleOpen = true;
        } else if (key == "v") {
            needValue(key);
            if (!lastParticleOpen) {
                fail("v= without a preceding x= for the same particle");
            }
            input.particles.back().velocity = parseVector3(value, key);
            lastParticleOpen = false;
        } else {
            fail("unknown key '" + key + "'");
        }
    }

    for (const char *k : {"dt", "steps", "E", "B"}) {
        const bool ok = (std::string(k) == "dt" && hasDt) || (std::string(k) == "steps" && hasSteps) ||
                        (std::string(k) == "E" && hasE) || (std::string(k) == "B" && hasB);
        if (!ok) {
            throw std::runtime_error(path + ": missing required key '" + std::string(k) + "'");
        }
    }
    if (input.dt <= 0.0) {
        throw std::runtime_error(path + ": dt must be > 0");
    }
    if (input.numSteps < 0) {
        throw std::runtime_error(path + ": steps must be >= 0");
    }
    if (input.particles.empty() || lastParticleOpen) {
        throw std::runtime_error(path + ": need at least one particle with x= and v= lines");
    }
    return input;
}
