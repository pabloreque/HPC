/*
Pablo Requeijo, October 1 2026.
Sessions 2-3, HPC coursework: PIC.
Catch2 tests of the .txt input parser (src/input.hpp), in test-particle
and hybrid mode.

Each test writes a small input file to the system temp directory and
checks that parseInputFile either reads it correctly or fails with a
"path:line: message" error pointing at the offending line.
*/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "input.hpp"

using Catch::Matchers::ContainsSubstring;

// A minimal valid input, one key per line (line numbers start at 1).
const std::vector<std::string> kValidLines = {"dt=0.1",  "steps=10", "E=0 0 0",
                                              "B=0 0 1", "q=1",      "m=1",
                                              "x=0 0 0", "v=1 0 0"};

// Returns kValidLines with the line for 'key' replaced by 'line'
// (an empty 'line' removes it), joined into one file content.
std::string inputWith(const std::string &key, const std::string &line) {
    std::string content;
    for (const std::string &valid : kValidLines) {
        const bool matches = valid.rfind(key + "=", 0) == 0;
        const std::string &chosen = matches ? line : valid;
        if (!chosen.empty()) {
            content += chosen + "\n";
        }
    }
    return content;
}

// Writes 'content' to a temp file named 'name' and returns its path.
// Each test case uses its own name so ctest can run them in parallel.
std::string writeTempInput(const std::string &name,
                           const std::string &content) {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / name;
    std::ofstream(path) << content;
    return path.string();
}

struct BadInputCase {
    std::string key;
    std::string line;
    std::string expectedError;
};

// Checks that every case makes parseInputFile throw the expected error.
void checkAllFail(const std::string &name,
                  const std::vector<BadInputCase> &cases) {
    for (const BadInputCase &c : cases) {
        CAPTURE(c.line);
        const std::string path = writeTempInput(name, inputWith(c.key, c.line));
        CHECK_THROWS_WITH(parseInputFile(path),
                          ContainsSubstring(c.expectedError));
    }
}

TEST_CASE("Input parser reads a valid file, skipping comments and blanks",
          "[input]") {
    const std::string path = writeTempInput(
        "pic_input_valid.txt",
        "# header comment\n\n" + inputWith("E", "E = 0.5 0 0  # inline"));
    const Input input = parseInputFile(path);

    REQUIRE(input.dt == 0.1);
    REQUIRE(input.numSteps == 10);
    REQUIRE(input.E.x == 0.5);
    REQUIRE(input.B.z == 1.0);
    REQUIRE(input.particles.size() == 1);
    REQUIRE(input.particles[0].velocity.x == 1.0);
    REQUIRE(input.particles[0].species.mass == 1.0);
}

TEST_CASE("Input parser rejects malformed numbers with path:line", "[input]") {
    checkAllFail(
        "pic_input_numbers.txt",
        {{"dt", "dt=abc", ":1: bad value for 'dt', expected a number"},
         {"dt", "dt=0.1abc", ":1: bad value for 'dt', expected a number"},
         {"steps", "steps=1.5",
          ":2: bad value for 'steps', expected an integer"},
         {"B", "B=0 0 1 7", ":4: bad value for 'B', expected three numbers"},
         {"m", "m=-1", ":6: m must be > 0"}});
}

TEST_CASE("Input parser rejects malformed lines and unknown keys", "[input]") {
    checkAllFail("pic_input_lines.txt",
                 {{"E", "E", ":3: expected key=value"},
                  {"E", "E=", ":3: empty value for 'E'"},
                  {"E", "foo=1", ":3: unknown key 'foo'"}});
}

TEST_CASE("Input parser reports missing required keys", "[input]") {
    checkAllFail("pic_input_missing.txt",
                 {{"dt", "", "missing required key 'dt'"},
                  {"B", "", "missing required key 'B'"}});
}

// A minimal valid hybrid-mode input (one key per line).
const std::vector<std::string> kHybridLines = {
    "dt=0.002", "steps=10", "B=1 0 0", "length=6.28", "nx=64", "ppc=4"};

// Joins kHybridLines plus 'extra' lines into one file content.
std::string hybridInputWith(const std::vector<std::string> &extra) {
    std::string content;
    for (const std::string &line : kHybridLines) {
        content += line + "\n";
    }
    for (const std::string &line : extra) {
        content += line + "\n";
    }
    return content;
}

TEST_CASE("Input parser reads a periodic length in test-particle mode",
          "[input][periodic]") {
    const std::string path = writeTempInput(
        "pic_input_length.txt", inputWith("dt", "dt=0.1\nlength=4"));
    const Input input = parseInputFile(path);
    REQUIRE(input.length.has_value());
    REQUIRE(*input.length == 4.0);
    REQUIRE_FALSE(input.hybrid.has_value());
}

TEST_CASE("Input parser reads a hybrid-mode file with defaults",
          "[input][hybrid]") {
    const std::string path = writeTempInput(
        "pic_input_hybrid.txt",
        hybridInputWith({"Te=0.1", "nu=0.001", "wave_mode=1",
                         "wave_amplitude=0.01", "wave_branch=L"}));
    const Input input = parseInputFile(path);
    REQUIRE(input.hybrid.has_value());
    const HybridSettings &s = *input.hybrid;
    REQUIRE(s.nx == 64);
    REQUIRE(*input.length == 6.28);
    REQUIRE(s.population.ppc == 4);
    REQUIRE(s.population.density == 1.0);
    REQUIRE(s.iterations == 2);
    REQUIRE(s.electrons.temperature == 0.1);
    REQUIRE(s.electrons.resistivity == 0.0);
    REQUIRE(s.electrons.hyperResistivity == 0.001);
    REQUIRE(s.wave.branch == WaveBranch::L);
    REQUIRE(input.particles.empty());
}

TEST_CASE("Input parser rejects inconsistent hybrid-mode files",
          "[input][hybrid]") {
    const std::vector<std::pair<std::vector<std::string>, std::string>> cases =
        {{{"E=0 0 0"}, "E= is not allowed in hybrid mode"},
         {{"q=1", "m=1", "x=0 0 0", "v=0 0 0"},
          "explicit particles are not allowed in hybrid mode"},
         {{"iterations=0"}, "iterations must be >= 1"},
         {{"nu=-0.1"}, "nu must be >= 0"},
         {{"wave_branch=X"}, "expected R or L"},
         {{"wave_amplitude=0.1"}, "wave_amplitude needs wave_mode > 0"},
         {{"wave_mode=1", "density=2"}, "a wave needs density=1"}};
    for (const auto &[extra, expected] : cases) {
        CAPTURE(expected);
        const std::string path =
            writeTempInput("pic_input_hybrid_bad.txt", hybridInputWith(extra));
        CHECK_THROWS_WITH(parseInputFile(path), ContainsSubstring(expected));
    }
}

TEST_CASE("Input parser requires nx= for hybrid-only keys", "[input][hybrid]") {
    checkAllFail("pic_input_hybrid_key.txt",
                 {{"dt", "dt=0.1\nppc=4", "'ppc' needs nx= (hybrid mode)"}});
    const std::string path = writeTempInput("pic_input_hybrid_missing.txt",
                                            "dt=0.1\nsteps=1\nB=1 0 0\nnx=8\n");
    CHECK_THROWS_WITH(parseInputFile(path),
                      ContainsSubstring("missing required key 'length'"));
}

// Content of a density profile file: a comment, a blank line, then
// 'values' one per line.
std::string densityProfileContent(const std::vector<std::string> &values) {
    std::string content = "# density on the primal nodes\n\n";
    for (const std::string &value : values) {
        content += value + "\n";
    }
    return content;
}

// 'count' density values 1, 1.01, 1.02, ... as text.
std::vector<std::string> rampValues(int count) {
    std::vector<std::string> values;
    for (int i = 0; i < count; ++i) {
        values.push_back(std::to_string(1.0 + 0.01 * i));
    }
    return values;
}

TEST_CASE("Input parser reads a density profile file next to the input",
          "[input][hybrid][density]") {
    writeTempInput("pic_density_ramp.txt",
                   densityProfileContent(rampValues(64)));
    const std::string path =
        writeTempInput("pic_input_density.txt",
                       hybridInputWith({"density_file=pic_density_ramp.txt"}));
    const Input input = parseInputFile(path);
    const std::vector<double> &profile =
        input.hybrid->population.densityProfile;
    REQUIRE(profile.size() == 64);
    REQUIRE(profile.front() == 1.0);
    REQUIRE(profile.back() == 1.63);
}

TEST_CASE("Input parser rejects bad density profiles",
          "[input][hybrid][density]") {
    std::vector<std::string> withZero = rampValues(64);
    withZero[5] = "0";
    std::vector<std::string> withText = rampValues(64);
    withText[0] = "abc";
    const std::string fileLine = "density_file=pic_density_bad.txt";
    const std::vector<std::tuple<std::vector<std::string>,
                                 std::vector<std::string>, std::string>>
        cases = {
            {rampValues(63),
             {fileLine},
             "expected nx=64 density values, got 63"},
            {withZero,
             {fileLine},
             "pic_density_bad.txt:8: density must be > 0"},
            {withText, {fileLine}, "pic_density_bad.txt:3: bad value"},
            {rampValues(64), {fileLine, "density=2"}, "mutually exclusive"},
            {rampValues(64),
             {fileLine, "wave_mode=1"},
             "a wave needs a uniform density"},
            {rampValues(64),
             {"density_file=pic_density_missing.txt"},
             "cannot open density file"}};
    for (const auto &[values, extra, expected] : cases) {
        CAPTURE(expected);
        writeTempInput("pic_density_bad.txt", densityProfileContent(values));
        const std::string path =
            writeTempInput("pic_input_density_bad.txt", hybridInputWith(extra));
        CHECK_THROWS_WITH(parseInputFile(path), ContainsSubstring(expected));
    }
}
