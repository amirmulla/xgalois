#ifndef XGALOIS_DATABASES_INTERFACE_HPP
#define XGALOIS_DATABASES_INTERFACE_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace xg {
namespace databases {

class DatabaseInterface {
 public:
  DatabaseInterface() = default;
  virtual ~DatabaseInterface() = default;

 protected:
  virtual const std::string& get_db_path() const = 0;

  bool file_exists() const {
    std::ifstream f(get_db_path().c_str());
    return f.good();
  }

  std::vector<std::string> parse_csv_line(const std::string& line) const {
    std::vector<std::string> segments;
    std::string current_segment;
    bool in_quotes = false;

    for (size_t i = 0; i < line.length(); ++i) {
      char c = line[i];
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
  static inline const std::string db_path_ =
      "xgalois/databases/prime_factors.db";

 protected:
  const std::string& get_db_path() const override { return db_path_; }

 public:
  PrimeFactorsDatabase() = default;

  PrimeFactorsResult fetch(uint64_t n) const {
    if (!file_exists()) {
      throw std::runtime_error("Database file not found: " + get_db_path());
    }

    std::ifstream db_file(get_db_path());
    std::string line;
    while (std::getline(db_file, line)) {

      if (line.empty() || line.substr(0, 2) == "//") {
        continue;
      }

      std::stringstream ss(line);
      std::string segment;
      std::vector<std::string> segments;

      while (std::getline(ss, segment,
                          ',')) {
        segments.push_back(segment);
      }

      if (segments.empty()) continue;

      uint64_t current_n;
      try {
        current_n = std::stoll(segments[0]);
      } catch (const std::invalid_argument& ia) {
        std::cerr << "Skipping malformed line (number): " << line << '\n';
        continue;
      }

      if (current_n == n) {
        PrimeFactorsResult result;
        if (segments.size() < 3) {
          std::cerr << "Skipping malformed line (segment count): " << line
                    << '\n';
          continue;
        }

        std::stringstream factors_ss(segments[1]);
        std::string factor_pair_str;
        while (std::getline(factors_ss, factor_pair_str,
                            ';')) {
          std::stringstream factor_pair_ss(factor_pair_str);
          std::string factor_str, mult_str;
          if (std::getline(factor_pair_ss, factor_str, ':') &&
              std::getline(factor_pair_ss, mult_str,
                           ':')) {
            try {
              result.factors.push_back(std::stoll(factor_str));
              result.multiplicities.push_back(std::stoi(mult_str));
            } catch (const std::invalid_argument& ia) {
              std::cerr << "Skipping malformed factor/multiplicity: "
                        << factor_pair_str << '\n';
              result.factors.clear();
              result.multiplicities.clear();
              break;
            }
          }
        }

        try {
          result.composite = std::stoi(segments[2]);
        } catch (const std::invalid_argument& ia) {
          std::cerr << "Skipping malformed line (composite flag): " << line
                    << '\n';
          result.composite = -1;
        }
        return result;
      }
    }
    throw std::runtime_error("Number not found in prime factors database: " +
                             std::to_string(n));
  }
};

struct IrreduciblePolyResult {
  std::vector<int> nonzero_degrees;
  std::vector<int> nonzero_coeffs;
};

class IrreduciblePolyDatabase : public DatabaseInterface {
 private:
  static inline const std::string db_path_ =
      "xgalois/databases/irreducible_polys.db";

 protected:
  const std::string& get_db_path() const override { return db_path_; }

 public:
  IrreduciblePolyDatabase() = default;

  IrreduciblePolyResult fetch(int characteristic, int degree) const {
    if (!file_exists()) {
      throw std::runtime_error("Database file not found: " + get_db_path());
    }

    std::ifstream db_file(get_db_path());
    std::string line;
    while (std::getline(db_file, line)) {

      if (line.empty() || line.substr(0, 2) == "//") {
        continue;
      }

      std::vector<std::string> segments = parse_csv_line(line);

      if (segments.size() < 4) continue;

      int current_char, current_deg;
      try {
        current_char = std::stoi(segments[0]);
        current_deg = std::stoi(segments[1]);
      } catch (const std::invalid_argument& ia) {
        std::cerr << "Skipping malformed line (char/deg): " << line << '\n';
        continue;
      }

      if (current_char == characteristic && current_deg == degree) {
        IrreduciblePolyResult result;
        std::stringstream degrees_ss(segments[2]);
        std::string deg_str;
        while (std::getline(degrees_ss, deg_str,
                            ';')) {
          try {
            result.nonzero_degrees.push_back(std::stoi(deg_str));
          } catch (const std::invalid_argument& ia) {
            std::cerr << "Skipping malformed degree: " << deg_str << '\n';
            result.nonzero_degrees.clear();
            break;
          }
        }
        if (result.nonzero_degrees.empty() && !segments[2].empty()) continue;

        std::stringstream coeffs_ss(segments[3]);
        std::string coeff_str;
        while (std::getline(coeffs_ss, coeff_str,
                            ';')) {
          try {
            result.nonzero_coeffs.push_back(std::stoi(coeff_str));
          } catch (const std::invalid_argument& ia) {
            std::cerr << "Skipping malformed coefficient: " << coeff_str
                      << '\n';
            result.nonzero_coeffs.clear();
            break;
          }
        }
        if (result.nonzero_coeffs.empty() && !segments[3].empty()) continue;

        if (result.nonzero_degrees.size() != result.nonzero_coeffs.size()) {
          std::cerr << "Mismatch between degree and coefficient counts: "
                    << line << '\n';
          continue;
        }
        return result;
      }
    }
    throw std::runtime_error(
        "Irreducible polynomial not found for characteristic " +
        std::to_string(characteristic) + ", degree " + std::to_string(degree));
  }
};

using ConwayPolyResult = IrreduciblePolyResult;

class ConwayPolyDatabase : public DatabaseInterface {
 private:
  static inline const std::string db_path_ =
      "xgalois/databases/conway_polys.db";

 protected:
  const std::string& get_db_path() const override { return db_path_; }

 public:
  ConwayPolyDatabase() = default;

  ConwayPolyResult fetch(int characteristic, int degree) const {
    if (!file_exists()) {
      throw std::runtime_error("Database file not found: " + get_db_path());
    }

    std::ifstream db_file(get_db_path());
    std::string line;
    while (std::getline(db_file, line)) {

      if (line.empty() || line.substr(0, 2) == "//") {
        continue;
      }

      std::vector<std::string> segments = parse_csv_line(line);

      if (segments.size() < 4) continue;

      int current_char, current_deg;
      try {
        current_char = std::stoi(segments[0]);
        current_deg = std::stoi(segments[1]);
      } catch (const std::invalid_argument& ia) {
        std::cerr << "Skipping malformed line (char/deg): " << line << '\n';
        continue;
      }

      if (current_char == characteristic && current_deg == degree) {
        ConwayPolyResult result;
        std::stringstream degrees_ss(segments[2]);
        std::string deg_str;
        while (std::getline(degrees_ss, deg_str,
                            ';')) {
          try {
            result.nonzero_degrees.push_back(std::stoi(deg_str));
          } catch (const std::invalid_argument& ia) {
            std::cerr << "Skipping malformed degree: " << deg_str << '\n';
            result.nonzero_degrees.clear();
            break;
          }
        }
        if (result.nonzero_degrees.empty() && !segments[2].empty()) continue;

        std::stringstream coeffs_ss(segments[3]);
        std::string coeff_str;
        while (std::getline(coeffs_ss, coeff_str,
                            ';')) {
          try {
            result.nonzero_coeffs.push_back(std::stoi(coeff_str));
          } catch (const std::invalid_argument& ia) {
            std::cerr << "Skipping malformed coefficient: " << coeff_str
                      << '\n';
            result.nonzero_coeffs.clear();
            break;
          }
        }
        if (result.nonzero_coeffs.empty() && !segments[3].empty()) continue;

        if (result.nonzero_degrees.size() != result.nonzero_coeffs.size()) {
          std::cerr << "Mismatch between degree and coefficient counts: "
                    << line << '\n';
          continue;
        }
        return result;
      }
    }
    throw std::runtime_error("Conway polynomial not found for characteristic " +
                             std::to_string(characteristic) + ", degree " +
                             std::to_string(degree));
  }
};

static std::string ConvertPolyResultToString(const std::vector<int>& degrees,
                                             const std::vector<int>& coeffs,
                                             int characteristic) {
  if (degrees.size() != coeffs.size()) {
    throw std::invalid_argument(
        "Mismatch between degrees and coefficients vectors");
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

  for (const auto& [degree, coeff] : degree_coeff_pairs) {
    if (coeff == 0) continue;

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
        if (coeff != 1) result += "*";
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

}
}

#endif
