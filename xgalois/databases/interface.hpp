#ifndef XGALOIS_DATABASES_INTERFACE_HPP_
#define XGALOIS_DATABASES_INTERFACE_HPP_

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/databases/embedded_conway_polys.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace databases {

class DatabaseInterface {
 public:
  DatabaseInterface() = default;
  virtual ~DatabaseInterface() = default;

 protected:
  virtual const std::string& get_db_path() const = 0;

  bool file_exists() const {
    std::ifstream f(get_db_path());
    return f.good();
  }

  std::vector<std::string> parse_csv_line(const std::string& line) const {
    std::vector<std::string> segments;
    std::string current_segment;
    bool in_quotes = false;

    for (char c : line) {
      if (c == '"') {
        in_quotes = !in_quotes;
      } else if (c == ',' && !in_quotes) {
        segments.push_back(current_segment);
        current_segment.clear();
      } else {
        current_segment += c;
      }
    }

    if (!current_segment.empty() || !segments.empty()) {
      segments.push_back(current_segment);
    }
    return segments;
  }
};

struct PrimeFactorsResult {
  std::vector<uint64_t> factors;
  std::vector<int> multiplicities;
  int composite;
};

class PrimeFactorsDatabase : public DatabaseInterface {
 private:
  static inline const std::string kEmptyPath = "";

 protected:
  const std::string& get_db_path() const override { return kEmptyPath; }

 public:
  PrimeFactorsDatabase() = default;

  // Computes exact prime factorization on the fly with zero disk I/O.
  PrimeFactorsResult fetch(uint64_t n) const {
    if (n < 2) {
      throw std::runtime_error(
          "Number not found in prime factors database (must be >= 2): " +
          std::to_string(n));
    }

    auto pairs = xg::utils::PrimeFactorizeWithMultiplicity(n);
    PrimeFactorsResult result;
    result.factors.reserve(pairs.size());
    result.multiplicities.reserve(pairs.size());

    for (const auto& [prime, multiplicity] : pairs) {
      result.factors.push_back(prime);
      result.multiplicities.push_back(static_cast<int>(multiplicity));
    }

    result.composite = (pairs.size() == 1 && pairs[0].second == 1) ? 0 : 1;
    return result;
  }
};

struct IrreduciblePolyResult {
  std::vector<int> nonzero_degrees;
  std::vector<int> nonzero_coeffs;
};

using ConwayPolyResult = IrreduciblePolyResult;

namespace detail {

struct PairHash {
  size_t operator()(const std::pair<uint64_t, uint64_t>& key) const noexcept {
    uint64_t h1 = key.first * 0x9e3779b97f4a7c15ULL;
    uint64_t h2 = key.second * 0xbf58476d1ce4e5b9ULL;
    h1 ^= (h1 >> 30);
    h2 ^= (h2 >> 27);
    return static_cast<size_t>(
        h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2)));
  }
};

using PolyTable =
    std::unordered_map<std::pair<uint64_t, uint64_t>,
                       std::vector<IrreduciblePolyResult>, PairHash>;

inline std::string Trim(std::string_view sv) {
  size_t start = 0;
  while (start < sv.size() &&
         std::isspace(static_cast<unsigned char>(sv[start])) != 0) {
    ++start;
  }
  size_t end = sv.size();
  while (end > start &&
         std::isspace(static_cast<unsigned char>(sv[end - 1])) != 0) {
    --end;
  }
  return std::string(sv.substr(start, end - start));
}

inline std::vector<std::string> ParseCsvLine(std::string_view line) {
  std::vector<std::string> segments;
  std::string current_segment;
  bool in_quotes = false;

  for (char c : line) {
    if (c == '"') {
      in_quotes = !in_quotes;
    } else if (c == ',' && !in_quotes) {
      segments.push_back(Trim(current_segment));
      current_segment.clear();
    } else {
      current_segment += c;
    }
  }

  if (!current_segment.empty() || !segments.empty()) {
    segments.push_back(Trim(current_segment));
  }
  return segments;
}

inline bool ParseUint64Strict(const std::string& token, uint64_t& out) {
  std::string trimmed = Trim(token);
  if (trimmed.empty()) {
    return false;
  }
  for (char c : trimmed) {
    if (std::isdigit(static_cast<unsigned char>(c)) == 0) {
      return false;
    }
  }
  try {
    size_t parsed_len = 0;
    uint64_t val = std::stoull(trimmed, &parsed_len, 10);
    if (parsed_len != trimmed.size()) {
      return false;
    }
    out = val;
    return true;
  } catch (const std::invalid_argument&) {
    return false;
  } catch (const std::out_of_range&) {
    return false;
  }
}

inline bool ParseIntList(const std::string& list_str, std::vector<int>& out) {
  out.clear();
  if (list_str.empty()) {
    return false;
  }
  std::stringstream ss(list_str);
  std::string item;
  while (std::getline(ss, item, ';')) {
    uint64_t val = 0;
    if (!ParseUint64Strict(item, val) ||
        val > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
      out.clear();
      return false;
    }
    out.push_back(static_cast<int>(val));
  }
  return !out.empty();
}

inline std::string ResolveDbPath(const std::string& filename,
                                 const std::string& explicit_path = "") {
  if (!explicit_path.empty()) {
    return explicit_path;
  }

  namespace fs = std::filesystem;
  std::error_code ec;

  // 1. Environment variable override XGALOIS_DB_DIR.
  if (const char* env_dir = std::getenv("XGALOIS_DB_DIR");
      env_dir != nullptr && env_dir[0] != '\0') {
    fs::path candidate = fs::path(env_dir) / filename;
    if (fs::exists(candidate, ec) && !ec) {
      return candidate.string();
    }
  }

  // 2. Bazel runfiles (TEST_SRCDIR and RUNFILES_DIR).
  for (const char* env_name : {"TEST_SRCDIR", "RUNFILES_DIR"}) {
    if (const char* base = std::getenv(env_name);
        base != nullptr && base[0] != '\0') {
      fs::path base_path(base);
      const fs::path candidates[] = {
          base_path / "_main" / "xgalois" / "databases" / filename,
          base_path / "__main__" / "xgalois" / "databases" / filename,
          base_path / "xgalois" / "databases" / filename,
      };
      for (const auto& candidate : candidates) {
        if (fs::exists(candidate, ec) && !ec) {
          return candidate.string();
        }
      }
    }
  }

  // 3. Relative fallback from workspace root.
  return "xgalois/databases/" + filename;
}

inline PolyTable LoadPolyDatabaseFile(const std::string& filepath) {
  std::ifstream db_file(filepath);
  if (!db_file.is_open()) {
    throw std::runtime_error("Database file not found: " + filepath);
  }

  PolyTable table;
  std::string line;
  while (std::getline(db_file, line)) {
    // Strip // inline comments BEFORE CSV splitting.
    size_t comment_pos = line.find("//");
    if (comment_pos != std::string::npos) {
      line.erase(comment_pos);
    }

    std::string trimmed = Trim(line);
    if (trimmed.empty()) {
      continue;
    }

    std::vector<std::string> segments = ParseCsvLine(trimmed);
    if (segments.size() < 4) {
      continue;
    }

    uint64_t characteristic = 0;
    uint64_t degree = 0;
    if (!ParseUint64Strict(segments[0], characteristic) ||
        !ParseUint64Strict(segments[1], degree)) {
      continue;
    }

    IrreduciblePolyResult result;
    if (!ParseIntList(segments[2], result.nonzero_degrees) ||
        !ParseIntList(segments[3], result.nonzero_coeffs)) {
      continue;
    }

    if (result.nonzero_degrees.size() != result.nonzero_coeffs.size()) {
      continue;
    }

    table[{characteristic, degree}].push_back(std::move(result));
  }

  return table;
}

struct CachedPolyTable {
  std::atomic<bool> initialized{false};
  std::mutex init_mutex;
  PolyTable table;

  const PolyTable& Get(const std::string& filename,
                       const std::string& explicit_path = "") {
    if (!initialized.load(std::memory_order_acquire)) {
      std::lock_guard<std::mutex> lock(init_mutex);
      if (!initialized.load(std::memory_order_relaxed)) {
        std::string resolved = ResolveDbPath(filename, explicit_path);
        table = LoadPolyDatabaseFile(resolved);
        initialized.store(true, std::memory_order_release);
      }
    }
    return table;
  }
};

inline CachedPolyTable& DefaultIrreduciblePolyCache() {
  static CachedPolyTable cache;
  return cache;
}

inline CachedPolyTable& DefaultConwayPolyCache() {
  static CachedPolyTable cache;
  return cache;
}

}  // namespace detail

class IrreduciblePolyDatabase : public DatabaseInterface {
 private:
  static inline const std::string kDefaultFilename = "irreducible_polys.db";
  static inline const std::string kDefaultDbPath =
      "xgalois/databases/irreducible_polys.db";

  std::string custom_db_path_;
  std::shared_ptr<detail::CachedPolyTable> custom_cache_;

  const detail::PolyTable& GetTable() const {
    if (!custom_db_path_.empty()) {
      return custom_cache_->Get(kDefaultFilename, custom_db_path_);
    }
    return detail::DefaultIrreduciblePolyCache().Get(kDefaultFilename);
  }

 protected:
  const std::string& get_db_path() const override {
    return custom_db_path_.empty() ? kDefaultDbPath : custom_db_path_;
  }

 public:
  IrreduciblePolyDatabase() = default;

  explicit IrreduciblePolyDatabase(std::string db_path)
      : custom_db_path_(std::move(db_path)),
        custom_cache_(std::make_shared<detail::CachedPolyTable>()) {}

  IrreduciblePolyResult fetch(int characteristic, int degree) const {
    const detail::PolyTable& table = GetTable();
    if (characteristic > 0 && degree > 0) {
      auto it = table.find({static_cast<uint64_t>(characteristic),
                            static_cast<uint64_t>(degree)});
      if (it != table.end() && !it->second.empty()) {
        return it->second.front();
      }
    }
    throw std::runtime_error(
        "Irreducible polynomial not found for characteristic " +
        std::to_string(characteristic) + ", degree " + std::to_string(degree));
  }

  std::vector<IrreduciblePolyResult> fetch_all(int characteristic,
                                               int degree) const {
    const detail::PolyTable& table = GetTable();
    if (characteristic > 0 && degree > 0) {
      auto it = table.find({static_cast<uint64_t>(characteristic),
                            static_cast<uint64_t>(degree)});
      if (it != table.end()) {
        return it->second;
      }
    }
    return {};
  }
};

class ConwayPolyDatabase : public DatabaseInterface {
 private:
  static inline const std::string kDefaultFilename = "conway_polys.db";
  static inline const std::string kDefaultDbPath =
      "xgalois/databases/conway_polys.db";

  std::string custom_db_path_;
  std::shared_ptr<detail::CachedPolyTable> custom_cache_;

  const detail::PolyTable& GetTable() const {
    if (!custom_db_path_.empty()) {
      return custom_cache_->Get(kDefaultFilename, custom_db_path_);
    }
    return detail::DefaultConwayPolyCache().Get(kDefaultFilename);
  }

 protected:
  const std::string& get_db_path() const override {
    return custom_db_path_.empty() ? kDefaultDbPath : custom_db_path_;
  }

 public:
  ConwayPolyDatabase() = default;

  explicit ConwayPolyDatabase(std::string db_path)
      : custom_db_path_(std::move(db_path)),
        custom_cache_(std::make_shared<detail::CachedPolyTable>()) {}

  ConwayPolyResult fetch(int characteristic, int degree) const {
    if (characteristic <= 0 || degree <= 0) {
      if (!custom_db_path_.empty()) {
        // Trigger file validation for explicit custom path.
        (void)GetTable();
      }
      throw std::runtime_error(
          "Conway polynomial not found for characteristic " +
          std::to_string(characteristic) + ", degree " +
          std::to_string(degree));
    }

    // Fast path: O(log N) lookup in embedded constexpr table with zero disk
    // I/O (only when using the default database).
    if (custom_db_path_.empty()) {
      if (auto embedded = detail::LookupEmbeddedConwayPoly(
              static_cast<uint32_t>(characteristic),
              static_cast<uint32_t>(degree));
          embedded.has_value()) {
        ConwayPolyResult result;
        result.nonzero_degrees.reserve(embedded->size());
        result.nonzero_coeffs.reserve(embedded->size());
        for (const auto& term : *embedded) {
          result.nonzero_degrees.push_back(static_cast<int>(term.degree));
          result.nonzero_coeffs.push_back(static_cast<int>(term.coeff));
        }
        return result;
      }
    }

    const detail::PolyTable& table = GetTable();
    auto it = table.find(
        {static_cast<uint64_t>(characteristic), static_cast<uint64_t>(degree)});
    if (it != table.end() && !it->second.empty()) {
      return it->second.front();
    }

    throw std::runtime_error("Conway polynomial not found for characteristic " +
                             std::to_string(characteristic) + ", degree " +
                             std::to_string(degree));
  }
};

inline std::string ConvertPolyResultToString(const std::vector<int>& degrees,
                                             const std::vector<int>& coeffs,
                                             int characteristic) {
  if (degrees.size() != coeffs.size()) {
    throw std::invalid_argument(
        "Mismatch between degrees and coefficients vectors");
  }
  if (characteristic < 2) {
    throw std::invalid_argument("Characteristic must be at least 2");
  }

  if (degrees.empty()) {
    return "0";
  }

  std::string result;
  bool first_term = true;

  std::vector<std::pair<int, int>> degree_coeff_pairs;
  degree_coeff_pairs.reserve(degrees.size());
  for (size_t i = 0; i < degrees.size(); ++i) {
    degree_coeff_pairs.emplace_back(degrees[i], coeffs[i]);
  }

  std::sort(degree_coeff_pairs.begin(), degree_coeff_pairs.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  for (const auto& [degree, raw_coeff] : degree_coeff_pairs) {
    int coeff = raw_coeff % characteristic;
    if (coeff < 0) {
      coeff += characteristic;
    }
    if (coeff == 0) {
      continue;
    }

    if (!first_term) {
      result += "+";
    }

    if (characteristic == 2) {
      if (degree == 0) {
        result += "1";
      } else if (degree == 1) {
        result += "x";
      } else {
        result += "x^" + std::to_string(degree);
      }
    } else {
      if (coeff != 1 || degree == 0) {
        result += std::to_string(coeff);
      }

      if (degree > 0) {
        if (coeff != 1) {
          result += "*";
        }
        if (degree == 1) {
          result += "x";
        } else {
          result += "x^" + std::to_string(degree);
        }
      }
    }

    first_term = false;
  }

  return result.empty() ? "0" : result;
}

inline std::string GetIrreduciblePolynomial(int characteristic, int degree) {
  try {
    databases::IrreduciblePolyDatabase irreducible_db;
    auto result = irreducible_db.fetch(characteristic, degree);

    return ConvertPolyResultToString(result.nonzero_degrees,
                                     result.nonzero_coeffs, characteristic);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error("No irreducible polynomial found for GF(" +
                             std::to_string(characteristic) + "^" +
                             std::to_string(degree) + "): " + e.what());
  }
}

inline std::string GetConwayPolynomial(int characteristic, int degree) {
  try {
    databases::ConwayPolyDatabase conway_db;
    auto result = conway_db.fetch(characteristic, degree);

    return ConvertPolyResultToString(result.nonzero_degrees,
                                     result.nonzero_coeffs, characteristic);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error("No Conway polynomial found for GF(" +
                             std::to_string(characteristic) + "^" +
                             std::to_string(degree) + "): " + e.what());
  }
}

}  // namespace databases
}  // namespace xg

#endif  // XGALOIS_DATABASES_INTERFACE_HPP_
