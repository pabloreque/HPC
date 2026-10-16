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

// Content of a profile file: a comment, a blank line, the header, then
// the data rows.
std::string profileContent(const std::string &header,
                           const std::vector<std::string> &rows) {
    std::string content = "# initial moments on the primal nodes\n\n";
    content += header + "\n";
    for (const std::string &row : rows) {
        content += row + "\n";
    }
    return content;
}

// 'count' rows "n uy" with n = 1, 1.01, 1.02, ... and uy = 0.5.
std::vector<std::string> rampRows(int count) {
    std::vector<std::string> rows;
    for (int i = 0; i < count; ++i) {
        rows.push_back(std::to_string(1.0 + 0.01 * i) + " 0.5");
    }
    return rows;
}

TEST_CASE("Input parser reads a profile file next to the input",
          "[input][hybrid][profile]") {
    writeTempInput("pic_profile_ramp.txt",
                   profileContent("n uy", rampRows(64)));
    const std::string path = writeTempInput(
        "pic_input_profile.txt",
        hybridInputWith({"profile_file=pic_profile_ramp.txt", "vth=0.2"}));
    const Input input = parseInputFile(path);
    const PopulationSettings &population = input.hybrid->population;
    const NodeProfiles &profiles = population.profiles;
    REQUIRE(profiles.density.size() == 64);
    REQUIRE(profiles.density.front() == 1.0);
    REQUIRE(profiles.density.back() == 1.63);
    REQUIRE(profiles.bulkVelocity.size() == 64);
    REQUIRE(profiles.bulkVelocity[3].x == 0.0);
    REQUIRE(profiles.bulkVelocity[3].y == 0.5);
    REQUIRE(profiles.bulkVelocity[3].z == 0.0);
    // No vth column: the uniform vth= applies.
    REQUIRE(profiles.thermalSpeed.empty());
    REQUIRE(population.vth == 0.2);
}

struct BadProfileCase {
    std::string header;
    std::vector<std::string> rows;
    std::vector<std::string> extraLines;
    std::string expectedError;
};

TEST_CASE("Input parser rejects bad profile files",
          "[input][hybrid][profile]") {
    std::vector<std::string> withZero = rampRows(64);
    withZero[5] = "0 0.5";
    std::vector<std::string> withText = rampRows(64);
    withText[0] = "abc 0.5";
    std::vector<std::string> shortRow = rampRows(64);
    shortRow[2] = "1.0";
    const std::string fileLine = "profile_file=pic_profile_bad.txt";
    const std::vector<BadProfileCase> cases = {
        {"n uy", rampRows(63), {fileLine}, "expected nx=64 rows, got 63"},
        {"n uy", withZero, {fileLine}, "pic_profile_bad.txt:9: n must be > 0"},
        {"n uy", withText, {fileLine}, "pic_profile_bad.txt:4: bad value"},
        {"n uy", shortRow, {fileLine}, "pic_profile_bad.txt:6: expected 2"},
        {"vth", {"-0.1"}, {fileLine}, "vth must be >= 0"},
        {"n T", rampRows(64), {fileLine}, "unknown profile column 'T'"},
        {"n n", rampRows(64), {fileLine}, "column 'n' appears twice"},
        {"n uy", rampRows(64), {fileLine, "density=2"}, "density= conflicts"},
        {"vth",
         std::vector<std::string>(64, "0.1"),
         {fileLine, "vth=0.1"},
         "vth= conflicts"},
        {"n uy",
         rampRows(64),
         {fileLine, "wave_mode=1"},
         "a wave needs uniform moments"},
        {"n uy",
         rampRows(64),
         {"profile_file=pic_profile_missing.txt"},
         "cannot open profile file"}};
    for (const BadProfileCase &c : cases) {
        CAPTURE(c.expectedError);
        writeTempInput("pic_profile_bad.txt", profileContent(c.header, c.rows));
        const std::string path = writeTempInput("pic_input_profile_bad.txt",
                                                hybridInputWith(c.extraLines));
        CHECK_THROWS_WITH(parseInputFile(path),
                          ContainsSubstring(c.expectedError));
    }
}
