/*
Pablo Requeijo, October 1 2026.
Session 2, HPC coursework: Boris pusher.
Catch2 tests of the .txt input parser (src/input.hpp).

Each test writes a small input file to the system temp directory and
checks that parseInputFile either reads it correctly or fails with a
"path:line: message" error pointing at the offending line.
*/

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <string>
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
