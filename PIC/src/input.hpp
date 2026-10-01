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
#include <utility>
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

// Formats "path:lineNo: message" for errors tied to one input line.
inline std::string lineError(const std::string &path, int lineNo,
                             const std::string &message) {
    return path + ":" + std::to_string(lineNo) + ": " + message;
}

// Splits "key = value" into trimmed key and value. Fails on lines without
// '=' or with an empty value.
inline KeyLine splitKeyValue(const std::string &line, const std::string &path,
                             int lineNo) {
    const std::string::size_type eq = line.find('=');
    if (eq == std::string::npos) {
        throw std::runtime_error(lineError(path, lineNo, "expected key=value"));
    }
    KeyLine keyLine{trim(line.substr(0, eq)), trim(line.substr(eq + 1))};
    if (keyLine.value.empty()) {
        throw std::runtime_error(
            lineError(path, lineNo, "empty value for '" + keyLine.key + "'"));
    }
    return keyLine;
}

// Reads the next meaningful "key=value" line: strips '#' comments and
// skips blank lines. Returns false at EOF. lineNo tracks physical lines
// for errors.
inline bool nextKeyLine(std::ifstream &file, const std::string &path,
                        int &lineNo, KeyLine &out) {
    std::string line;
    while (std::getline(file, line)) {
        ++lineNo;
        // Strip comments.
        line = trim(line.substr(0, line.find('#')));
        if (!line.empty()) {
            out = splitKeyValue(line, path, lineNo);
            return true;
        }
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

    // Throws runtime_error prefixed with the current path:lineNo.
    [[noreturn]] void fail(const std::string &message) const {
        throw std::runtime_error(lineError(path, lineNo, message));
    }

    // Opens a new particle from its x= line, carrying the pending q=/m=.
    void openParticle(const std::string &value) {
        if (lastParticleOpen) {
            fail("each particle needs one v= line before the next x=");
        }
        if (!hasPendingCharge || !hasPendingMass) {
            fail("each particle needs its own q= and m= lines before x=");
        }
        input.particles.push_back(
            Particle{parseVector3(value, "x"), Vector3{}, pendingSpecies});
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

    // Handles dt, steps, E, B. Returns false if the key is not one of them.
    bool applyGlobalKey(const std::string &key, const std::string &value) {
        if (key == "dt") {
            input.dt = parseDouble(value, key);
            hasDt = true;
        } else if (key == "steps") {
            input.numSteps = parseInt(value, key);
            hasSteps = true;
        } else if (key == "E") {
            input.E = parseVector3(value, key);
            hasE = true;
        } else if (key == "B") {
            input.B = parseVector3(value, key);
            hasB = true;
        } else {
            return false;
        }
        return true;
    }

    // Stores the pending mass, rejecting non-positive values.
    void applyMass(const std::string &value) {
        pendingSpecies.mass = parseDouble(value, "m");
        if (pendingSpecies.mass <= 0.0) {
            fail("m must be > 0");
        }
        hasPendingMass = true;
    }

    // Handles q, m, x, v. Returns false if the key is not one of them.
    bool applyParticleKey(const std::string &key, const std::string &value) {
        if (key == "q") {
            pendingSpecies.charge = parseDouble(value, key);
            hasPendingCharge = true;
        } else if (key == "m") {
            applyMass(value);
        } else if (key == "x") {
            openParticle(value);
        } else if (key == "v") {
            closeParticle(value);
        } else {
            return false;
        }
        return true;
    }

    // Dispatches one key=value line; fails on unknown keys and adds
    // path:lineNo to value parsing errors.
    void applyLine(const std::string &key, const std::string &value) {
        try {
            if (!applyGlobalKey(key, value) && !applyParticleKey(key, value)) {
                fail("unknown key '" + key + "'");
            }
        } catch (const std::invalid_argument &e) {
            fail(e.what());
        }
    }
};

// Throws if any of dt, steps, E, B never appeared in the file.
inline void checkRequiredKeys(const InputBuilder &b) {
    const std::pair<const char *, bool> required[] = {
        {"dt", b.hasDt}, {"steps", b.hasSteps}, {"E", b.hasE}, {"B", b.hasB}};
    for (const auto &[key, present] : required) {
        if (!present) {
            throw std::runtime_error(b.path + ": missing required key '" + key +
                                     "'");
        }
    }
}

// Throws if dt or steps are out of range, or no complete particle exists.
inline void checkValues(const InputBuilder &b) {
    if (b.input.dt <= 0.0) {
        throw std::runtime_error(b.path + ": dt must be > 0");
    }
    if (b.input.numSteps < 0) {
        throw std::runtime_error(b.path + ": steps must be >= 0");
    }
    if (b.input.particles.empty() || b.lastParticleOpen) {
        throw std::runtime_error(
            b.path + ": need at least one particle with x= and v= lines");
    }
}

// Checks required keys and value ranges once all lines are consumed.
inline void validateInput(const InputBuilder &b) {
    checkRequiredKeys(b);
    checkValues(b);
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
