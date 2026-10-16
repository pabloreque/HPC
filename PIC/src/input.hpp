/*
Pablo Requeijo, September 24 2026 (Session 2); October 1 2026 (Session 3).
HPC coursework: PIC.
.txt input parsing: initial conditions for the generic driver.

Reads key=value files (see examples/) in one of two modes:

- Test-particle mode (Session 2): dt, steps, E, B, then one particle
  after another, each with its own q, m, x, v lines (a q=/m= line applies
  to the next particle defined by its x=/v= pair, so each particle can be
  a different species). Optional length= makes x periodic in [0, length).
- Hybrid mode (Session 3), selected by an nx= line: dt, steps, B (initial
  uniform field), length, plus the keys of hybrid_input.hpp. E is computed
  from Ohm's law and particles are loaded from ppc=, so E= and explicit
  particles are rejected. A profile_file= line is read here into the
  initial moment profiles (format in profile_file.hpp).

Unknown keys and missing or invalid values fail with a descriptive error.
*/
#pragma once

#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "boris_pusher.hpp"
#include "hybrid_input.hpp"
#include "parse.hpp"
#include "profile_file.hpp"
#include "vector3.hpp"

struct Input {
    double dt = 0.0;
    int numSteps = 0;
    Vector3 E;
    Vector3 B;
    std::vector<Particle3D> particles;
    // Periodic box length along x; unset = unbounded (test-particle mode).
    std::optional<double> length;
    // Hybrid PIC settings, present only in hybrid mode (nx= line).
    std::optional<HybridSettings> hybrid;
};

struct KeyLine {
    std::string key;
    std::string value;
};

// Splits "key = value" into trimmed key and value. Fails on lines without
// '=' or with an empty value.
inline KeyLine splitKeyValue(const std::string &line, const std::string &path,
                             int lineNo) {
    const std::string::size_type equalsPosition = line.find('=');
    if (equalsPosition == std::string::npos) {
        throw std::runtime_error(lineError(path, lineNo, "expected key=value"));
    }
    KeyLine keyLine{trim(line.substr(0, equalsPosition)),
                    trim(line.substr(equalsPosition + 1))};
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
        line = stripComment(line);
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
    HybridSettings hybrid;
    bool hasNx = false, hasPpc = false;
    // Uniform moments set explicitly, which a profile column must not
    // also set.
    bool hasDensity = false, hasVth = false;
    // First hybrid-only key seen, reported if nx= is missing.
    std::string firstHybridKey;
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
            Particle3D{parseVector3(value, "x"), Vector3{}, pendingSpecies});
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

    // Handles dt, steps, E, B, length. Returns false for any other key.
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
        } else if (key == "length") {
            input.length = parseDouble(value, key);
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

    // Handles the hybrid-mode keys. Returns false for any other key.
    bool applyHybridLine(const std::string &key, const std::string &value) {
        if (!applyHybridKey(hybrid, key, value)) {
            return false;
        }
        hasNx = hasNx || key == "nx";
        hasPpc = hasPpc || key == "ppc";
        hasDensity = hasDensity || key == "density";
        hasVth = hasVth || key == "vth";
        if (firstHybridKey.empty()) {
            firstHybridKey = key;
        }
        return true;
    }

    // Dispatches one key=value line; fails on unknown keys and adds
    // path:lineNo to value parsing errors.
    void applyLine(const std::string &key, const std::string &value) {
        try {
            if (!applyGlobalKey(key, value) && !applyParticleKey(key, value) &&
                !applyHybridLine(key, value)) {
                fail("unknown key '" + key + "'");
            }
        } catch (const std::invalid_argument &error) {
            fail(error.what());
        }
    }
};

// Throws runtime_error "path: message" for whole-file errors.
[[noreturn]] inline void failFile(const InputBuilder &builder,
                                  const std::string &message) {
    throw std::runtime_error(builder.path + ": " + message);
}

// Throws if a key required by the file's mode never appeared: dt, steps,
// B always; E in test-particle mode; length and ppc in hybrid mode.
inline void checkRequiredKeys(const InputBuilder &builder) {
    std::vector<std::pair<const char *, bool>> required = {
        {"dt", builder.hasDt},
        {"steps", builder.hasSteps},
        {"B", builder.hasB}};
    if (builder.hasNx) {
        required.push_back({"length", builder.input.length.has_value()});
        required.push_back({"ppc", builder.hasPpc});
    } else {
        required.push_back({"E", builder.hasE});
    }
    for (const auto &[key, present] : required) {
        if (!present) {
            failFile(builder,
                     std::string("missing required key '") + key + "'");
        }
    }
}

// Throws if dt, steps or length are out of range.
inline void checkValues(const InputBuilder &builder) {
    if (builder.input.dt <= 0.0) {
        failFile(builder, "dt must be > 0");
    }
    if (builder.input.numSteps < 0) {
        failFile(builder, "steps must be >= 0");
    }
    if (builder.input.length && *builder.input.length <= 0.0) {
        failFile(builder, "length must be > 0");
    }
}

// Test-particle mode: needs complete particles and no hybrid-only keys.
inline void checkTestParticleMode(const InputBuilder &builder) {
    if (!builder.firstHybridKey.empty()) {
        failFile(builder,
                 "'" + builder.firstHybridKey + "' needs nx= (hybrid mode)");
    }
    if (builder.input.particles.empty() || builder.lastParticleOpen) {
        failFile(builder, "need at least one particle with x= and v= lines");
    }
}

// Hybrid mode: no E= or explicit particles, and valid hybrid settings.
inline void checkHybridMode(const InputBuilder &builder) {
    if (builder.hasE) {
        failFile(builder,
                 "E= is not allowed in hybrid mode (E comes from Ohm's "
                 "law)");
    }
    if (!builder.input.particles.empty() || builder.hasPendingCharge ||
        builder.hasPendingMass) {
        failFile(builder, "explicit particles are not allowed in hybrid mode "
                          "(use ppc=)");
    }
    try {
        validateHybridSettings(builder.hybrid, builder.input.B);
    } catch (const std::invalid_argument &error) {
        failFile(builder, error.what());
    }
}

// Checks required keys and value ranges once all lines are consumed.
inline void validateInput(const InputBuilder &builder) {
    checkRequiredKeys(builder);
    checkValues(builder);
    if (builder.hasNx) {
        checkHybridMode(builder);
    } else {
        checkTestParticleMode(builder);
    }
}

// Resolves 'path' relative to the directory of the input file
// 'inputPath'; an absolute 'path' is returned unchanged.
inline std::string resolveRelativeTo(const std::string &inputPath,
                                     const std::string &path) {
    return (std::filesystem::path(inputPath).parent_path() / path).string();
}

// Throws if a moment is given both by a uniform key and by a column of
// the profile file.
inline void checkProfileConflicts(const InputBuilder &builder,
                                  const NodeProfiles &profiles) {
    if (builder.hasDensity && !profiles.density.empty()) {
        failFile(builder, "density= conflicts with the n column of the "
                          "profile file");
    }
    if (builder.hasVth && !profiles.thermalSpeed.empty()) {
        failFile(builder, "vth= conflicts with the vth column of the "
                          "profile file");
    }
}

// The validated hybrid settings, with the profile file (if any) read
// into population.profiles.
inline HybridSettings loadHybridSettings(const InputBuilder &builder) {
    HybridSettings settings = builder.hybrid;
    if (!settings.profileFile.empty()) {
        settings.population.profiles = readProfileFile(
            resolveRelativeTo(builder.path, settings.profileFile), settings.nx);
        checkProfileConflicts(builder, settings.population.profiles);
    }
    return settings;
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
    if (builder.hasNx) {
        builder.input.hybrid = loadHybridSettings(builder);
    }
    return builder.input;
}
