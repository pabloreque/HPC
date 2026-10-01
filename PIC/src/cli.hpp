/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Command-line parsing for the generic driver.
*/
#pragma once

#include <filesystem>
#include <stdexcept> // std::runtime_error for invalid CLI arguments
#include <string>

struct RunConfig {
    std::string inputPath;
    std::string outputPath = "results.h5";
    std::string groupName = "/trajectory";
};

// True if arg ends in ".h5".
inline bool hasH5Suffix(const std::string &arg) {
    return std::filesystem::path(arg).extension() == ".h5";
}

// Positionals are the plain arguments after the input file (not options like
// --group): at most [output.h5] and [group].
constexpr int kMaxPositionals = 2;

// Throws runtime_error if the number of positionals exceeds the limit.
inline void checkPositionalLimit(const std::string &arg, int count) {
    if (count > kMaxPositionals) {
        throw std::runtime_error("error: unexpected argument '" + arg + "'");
    }
}

// Stores one positional: output path if it ends in .h5, group otherwise.
inline void applyPositional(RunConfig &config, const std::string &arg,
                            int &positionalCount) {
    ++positionalCount;
    checkPositionalLimit(arg, positionalCount);
    if (hasH5Suffix(arg)) {
        config.outputPath = arg;
    } else {
        config.groupName = arg;
    }
}

// Reads the value following "--group" at argv[i]. Returns 1, the extra
// argv entry consumed.
inline int applyGroupOption(RunConfig &config, int i, int argc, char *argv[]) {
    if (i + 1 >= argc) {
        throw std::runtime_error("error: --group needs a value");
    }
    config.groupName = argv[i + 1];
    return 1;
}

// Consumes one option or positional at argv[i]. Returns extra argv entries
// consumed (1 for "--group <value>", 0 otherwise).
inline int applyArg(RunConfig &config, int i, int argc, char *argv[],
                    int &positionalCount) {
    const std::string arg = argv[i];
    if (arg == "--group") {
        return applyGroupOption(config, i, argc, argv);
    }
    if (arg.rfind("--", 0) == 0) {
        throw std::runtime_error("error: unknown option '" + arg + "'");
    }
    applyPositional(config, arg, positionalCount);
    return 0;
}

// Throws runtime_error if the group name is empty (from "" on the CLI).
inline void checkGroupName(const std::string &groupName) {
    if (groupName.empty()) {
        throw std::runtime_error("error: group name cannot be empty");
    }
}

// Builds the usage text shown when the input file is missing.
inline std::string usageMessage(const std::string &program) {
    return "Usage: " + program + " <input.txt> [output.h5] [group]\n" +
           "       " + program + " <input.txt> [output.h5] --group <group>";
}

// Parses argv into a RunConfig. Throws runtime_error carrying the complete
// message to print: bare usage text when args are missing, "error: ..."
// lines for bad options. Defaults: output "results.h5", group "/trajectory".
inline RunConfig parseArgs(int argc, char *argv[]) {
    if (argc < 2) {
        throw std::runtime_error(usageMessage(argv[0]));
    }
    RunConfig config;
    config.inputPath = argv[1];
    int positionalCount = 0;
    for (int i = 2; i < argc; ++i) {
        i += applyArg(config, i, argc, argv, positionalCount);
    }
    checkGroupName(config.groupName);
    return config;
}
