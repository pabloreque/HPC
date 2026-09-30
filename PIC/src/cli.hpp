/*
Pablo Requeijo, September 24 2026.
Session 2, HPC coursework: Boris pusher.
Command-line parsing for the generic driver.
*/
#pragma once

#include <stdexcept>
#include <string>

struct RunConfig {
    std::string inputPath;
    std::string outputPath = "results.h5";
    std::string groupName = "/trajectory";
};

// A positional is the output file, unless it lacks the .h5 suffix,
// in which case it is the group name.
inline void applyPositional(RunConfig &config, const std::string &arg) {
    if (arg.size() >= 3 && arg.substr(arg.size() - 3) == ".h5") {
        config.outputPath = arg;
    } else {
        config.groupName = arg;
    }
}

// Consumes one option or positional at argv[i]. Returns extra argv entries
// consumed (1 for "--group <value>", 0 otherwise).
inline int applyArg(RunConfig &config, int i, int argc, char *argv[]) {
    const std::string arg = argv[i];
    if (arg == "--group") {
        if (i + 1 >= argc) {
            throw std::runtime_error("error: --group needs a value");
        }
        config.groupName = argv[i + 1];
        return 1;
    }
    if (arg.rfind("--", 0) == 0) {
        throw std::runtime_error("error: unknown option '" + arg + "'");
    }
    if (config.outputPath == "results.h5" && (i == 2)) {
        // First positional after input is the output file... unless the
        // caller passed only a group; disambiguate by extension.
        applyPositional(config, arg);
        return 0;
    }
    if (config.groupName == "/trajectory") {
        config.groupName = arg;
        return 0;
    }
    throw std::runtime_error("error: unexpected argument '" + arg + "'");
}

// Parses argv into a RunConfig. Throws runtime_error carrying the complete
// message to print: bare usage text when args are missing, "error: ..."
// lines for bad options. Defaults: output "results.h5", group "/trajectory".
inline RunConfig parseArgs(int argc, char *argv[]) {
    if (argc < 2) {
        throw std::runtime_error(std::string("Usage: ") + argv[0] +
                                 " <input.txt> [output.h5] [group]\n       " + argv[0] +
                                 " <input.txt> [output.h5] --group <group>");
    }
    RunConfig config{argv[1], "results.h5", "/trajectory"};
    for (int i = 2; i < argc; ++i) {
        i += applyArg(config, i, argc, argv);
    }
    // Edge case: "./prog input.txt results.h5" leaves group default -- correct.
    // "./prog input.txt /cyclotron" is detected above by missing .h5 suffix.
    return config;
}
