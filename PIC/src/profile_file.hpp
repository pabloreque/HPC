/*
Pablo Requeijo, October 16 2026.
Session 3, HPC coursework: hybrid PIC.
Initial moment profile file (profile_file= in hybrid mode).

Plain text; '#' comments and blank lines are skipped. The first line is
a header naming the columns, any subset of n, ux, uy, uz, vth in any
order, followed by one row per primal node x_i = i*dx (exactly nx rows):

  n     uy    vth
  1.10  0.30  0.10
  ...

Moments left out stay uniform (density=, u = 0, vth=). n must be > 0 and
vth >= 0. The result is a NodeProfiles (loader.hpp), the same arrays a
Python front end could pass directly.
*/
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "loader.hpp"
#include "parse.hpp"
#include "vector3.hpp"

constexpr std::array<const char *, 5> kProfileColumns = {"n", "ux", "uy", "uz",
                                                         "vth"};

// The columns of a profile file, named by its header, in file order.
struct ProfileTable {
    std::vector<std::string> names;
    std::vector<std::vector<double>> columns;
};

// Splits a line into whitespace-separated words.
inline std::vector<std::string> splitWords(const std::string &line) {
    std::istringstream in(line);
    std::vector<std::string> words;
    std::string word;
    while (in >> word) {
        words.push_back(word);
    }
    return words;
}

// Throws invalid_argument unless every header name is a known column
// that appears only once.
inline void checkProfileHeader(const std::vector<std::string> &names) {
    for (const std::string &name : names) {
        require(std::find(kProfileColumns.begin(), kProfileColumns.end(),
                          name) != kProfileColumns.end(),
                "unknown profile column '" + name +
                    "' (expected n, ux, uy, uz, vth)");
        require(std::count(names.begin(), names.end(), name) == 1,
                "profile column '" + name + "' appears twice");
    }
}

// Parses one value of column 'name': n must be > 0 and vth >= 0.
inline double parseProfileValue(const std::string &text,
                                const std::string &name) {
    const double value = parseDouble(text, name);
    if (name == "n") {
        require(value > 0.0, "n must be > 0");
    } else if (name == "vth") {
        require(value >= 0.0, "vth must be >= 0");
    }
    return value;
}

// Appends one data row, one value per header column, to the table.
inline void appendProfileRow(ProfileTable &table,
                             const std::vector<std::string> &words) {
    require(words.size() == table.names.size(),
            "expected " + std::to_string(table.names.size()) +
                " values (one per header column), got " +
                std::to_string(words.size()));
    for (std::size_t j = 0; j < words.size(); ++j) {
        table.columns[j].push_back(parseProfileValue(words[j], table.names[j]));
    }
}

// Adds one non-blank line to the table: the header first, then data rows.
inline void readProfileLine(ProfileTable &table, const std::string &line) {
    const std::vector<std::string> words = splitWords(line);
    if (table.names.empty()) {
        checkProfileHeader(words);
        table.names = words;
        table.columns.resize(words.size());
    } else {
        appendProfileRow(table, words);
    }
}

// Number of data rows read so far.
inline std::size_t rowCount(const ProfileTable &table) {
    return table.columns.empty() ? 0 : table.columns.front().size();
}

// Reads a profile file into its table, with path:line errors. Throws
// unless it holds exactly nx data rows.
inline ProfileTable readProfileTable(const std::string &path, int nx) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open profile file '" + path + "'");
    }
    ProfileTable table;
    std::string line;
    for (int lineNo = 1; std::getline(file, line); ++lineNo) {
        line = stripComment(line);
        try {
            if (!line.empty()) {
                readProfileLine(table, line);
            }
        } catch (const std::invalid_argument &error) {
            throw std::runtime_error(lineError(path, lineNo, error.what()));
        }
    }
    if (rowCount(table) != static_cast<std::size_t>(nx)) {
        throw std::runtime_error(path + ": expected nx=" + std::to_string(nx) +
                                 " rows, got " +
                                 std::to_string(rowCount(table)));
    }
    return table;
}

// The column 'name' of the table, or an empty vector if it is absent.
inline std::vector<double> profileColumn(const ProfileTable &table,
                                         const std::string &name) {
    for (std::size_t j = 0; j < table.names.size(); ++j) {
        if (table.names[j] == name) {
            return table.columns[j];
        }
    }
    return {};
}

// The bulk velocity at every node from the ux, uy, uz columns (an absent
// one is 0); empty if none of them is given.
inline std::vector<Vector3> bulkVelocityRows(const ProfileTable &table) {
    const std::vector<double> ux = profileColumn(table, "ux");
    const std::vector<double> uy = profileColumn(table, "uy");
    const std::vector<double> uz = profileColumn(table, "uz");
    if (ux.empty() && uy.empty() && uz.empty()) {
        return {};
    }
    std::vector<Vector3> rows(rowCount(table));
    for (std::size_t i = 0; i < rows.size(); ++i) {
        rows[i] = Vector3{ux.empty() ? 0.0 : ux[i], uy.empty() ? 0.0 : uy[i],
                          uz.empty() ? 0.0 : uz[i]};
    }
    return rows;
}

// Reads a profile file with nx rows into the initial node profiles.
inline NodeProfiles readProfileFile(const std::string &path, int nx) {
    const ProfileTable table = readProfileTable(path, nx);
    return NodeProfiles{profileColumn(table, "n"), bulkVelocityRows(table),
                        profileColumn(table, "vth")};
}
