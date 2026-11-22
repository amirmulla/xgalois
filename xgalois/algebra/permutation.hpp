#ifndef XGALOIS_ALGEBRA_PERMUTATION_HPP_
#define XGALOIS_ALGEBRA_PERMUTATION_HPP_

#include <algorithm>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace algebra {

// Finite permutation on {0, 1, ..., degree - 1} using left-action composition:
//   (a * b)(i) = a(b(i)).
class Permutation {
 public:
  Permutation() = default;

  static Permutation Identity(uint32_t degree) {
    std::vector<uint32_t> img(degree);
    std::iota(img.begin(), img.end(), uint32_t{0});
    return Permutation(std::move(img), UncheckedTag{});
  }

  static Permutation FromMapping(std::vector<uint32_t> mapping) {
    const std::size_t n = mapping.size();
    std::vector<uint8_t> seen(n, 0);
    for (uint32_t v : mapping) {
      if (v >= n) {
        throw std::invalid_argument(
            "Permutation::FromMapping: image value out of range");
      }
      if (seen[v] != 0) {
        throw std::invalid_argument(
            "Permutation::FromMapping: mapping is not a bijection");
      }
      seen[v] = 1;
    }
    return Permutation(std::move(mapping), UncheckedTag{});
  }

  // Constructs a permutation from a product of cycles in left-action order:
  //   FromCycles(n, {c_0, c_1, ...}) = c_0 * c_1 * ...
  // Each cycle (a_0 a_1 ... a_{k-1}) maps a_0 -> a_1 -> ... -> a_{k-1} -> a_0.
  static Permutation FromCycles(
      uint32_t degree, const std::vector<std::vector<uint32_t>>& cycles) {
    Permutation result = Identity(degree);
    std::vector<uint8_t> seen(degree, 0);
    for (const auto& cycle : cycles) {
      for (uint32_t v : cycle) {
        if (v >= degree) {
          throw std::invalid_argument(
              "Permutation::FromCycles: cycle entry out of range");
        }
        if (seen[v] != 0) {
          throw std::invalid_argument(
              "Permutation::FromCycles: duplicate entry within a single cycle");
        }
        seen[v] = 1;
      }
      for (uint32_t v : cycle) {
        seen[v] = 0;
      }
      if (cycle.size() <= 1) {
        continue;
      }
      std::vector<uint32_t> img(degree);
      std::iota(img.begin(), img.end(), uint32_t{0});
      for (std::size_t i = 0; i < cycle.size(); ++i) {
        img[cycle[i]] = cycle[(i + 1) % cycle.size()];
      }
      result = result.Compose(Permutation(std::move(img), UncheckedTag{}));
    }
    return result;
  }

  uint32_t Degree() const noexcept {
    return static_cast<uint32_t>(img_.size());
  }

  const std::vector<uint32_t>& Mapping() const noexcept { return img_; }

  uint32_t Apply(uint32_t i) const {
    if (i >= img_.size()) {
      throw std::out_of_range("Permutation::Apply: index out of range");
    }
    return img_[i];
  }

  uint32_t operator()(uint32_t i) const { return Apply(i); }

  // Left-action composition: (this * other)(i) = (*this)(other(i)).
  Permutation Compose(const Permutation& other) const {
    if (img_.size() != other.img_.size()) {
      throw std::invalid_argument(
          "Permutation::Compose: permutations must have equal degree");
    }
    const std::size_t n = img_.size();
    std::vector<uint32_t> out(n);
    for (std::size_t i = 0; i < n; ++i) {
      out[i] = img_[other.img_[i]];
    }
    return Permutation(std::move(out), UncheckedTag{});
  }

  Permutation operator*(const Permutation& other) const {
    return Compose(other);
  }

  Permutation& operator*=(const Permutation& other) {
    *this = Compose(other);
    return *this;
  }

  Permutation Inverse() const {
    const std::size_t n = img_.size();
    std::vector<uint32_t> inv(n);
    for (std::size_t i = 0; i < n; ++i) {
      inv[img_[i]] = static_cast<uint32_t>(i);
    }
    return Permutation(std::move(inv), UncheckedTag{});
  }

  // O(n) cycle-based powering supporting arbitrary signed exponents.
  Permutation Pow(int64_t exponent) const {
    if (exponent == 0) {
      return Identity(Degree());
    }
    if (exponent == 1) {
      return *this;
    }
    if (exponent == -1) {
      return Inverse();
    }
    const bool negative = exponent < 0;
    const uint64_t abs_exp =
        negative ? (uint64_t{0} - static_cast<uint64_t>(exponent))
                 : static_cast<uint64_t>(exponent);
    return PowWithSign(abs_exp, negative);
  }

  // O(n) cycle-based powering for unsigned 64-bit exponents.
  Permutation PowUnsigned(uint64_t exponent) const {
    if (exponent == 0) {
      return Identity(Degree());
    }
    if (exponent == 1) {
      return *this;
    }
    return PowWithSign(exponent, false);
  }

  bool IsIdentity() const noexcept {
    for (std::size_t i = 0; i < img_.size(); ++i) {
      if (img_[i] != static_cast<uint32_t>(i)) {
        return false;
      }
    }
    return true;
  }

  // Disjoint cycles of length >= 2, canonically ordered by smallest element.
  std::vector<std::vector<uint32_t>> Cycles() const {
    const std::size_t n = img_.size();
    std::vector<uint8_t> visited(n, 0);
    std::vector<std::vector<uint32_t>> cycles;
    for (std::size_t i = 0; i < n; ++i) {
      if (visited[i] != 0) {
        continue;
      }
      if (img_[i] == static_cast<uint32_t>(i)) {
        visited[i] = 1;
        continue;
      }
      std::vector<uint32_t> cycle;
      uint32_t cur = static_cast<uint32_t>(i);
      while (visited[cur] == 0) {
        visited[cur] = 1;
        cycle.push_back(cur);
        cur = img_[cur];
      }
      cycles.push_back(std::move(cycle));
    }
    return cycles;
  }

  // Integer partition of Degree() listing all cycle lengths (including
  // 1-cycles) sorted in non-increasing order.
  std::vector<uint32_t> CycleType() const {
    const std::size_t n = img_.size();
    std::vector<uint8_t> visited(n, 0);
    std::vector<uint32_t> lengths;
    for (std::size_t i = 0; i < n; ++i) {
      if (visited[i] != 0) {
        continue;
      }
      uint32_t len = 0;
      uint32_t cur = static_cast<uint32_t>(i);
      while (visited[cur] == 0) {
        visited[cur] = 1;
        cur = img_[cur];
        ++len;
      }
      lengths.push_back(len);
    }
    std::sort(lengths.begin(), lengths.end(), std::greater<uint32_t>{});
    return lengths;
  }

  // Exact group-theoretic order computed as the LCM of disjoint cycle lengths.
  BigInt Order() const {
    BigInt lcm(1);
    const std::size_t n = img_.size();
    std::vector<uint8_t> visited(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
      if (visited[i] != 0) {
        continue;
      }
      uint64_t len = 0;
      uint32_t cur = static_cast<uint32_t>(i);
      while (visited[cur] == 0) {
        visited[cur] = 1;
        cur = img_[cur];
        ++len;
      }
      if (len > 1) {
        const BigInt len_big = BigInt::FromUint64(len);
        const BigInt g = BigInt::Gcd(lcm, len_big);
        lcm = (lcm / g) * len_big;
      }
    }
    return lcm;
  }

  // Parity sign: +1 for even permutations, -1 for odd permutations.
  int Sign() const {
    const std::size_t n = img_.size();
    std::vector<uint8_t> visited(n, 0);
    std::size_t parity = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (visited[i] != 0) {
        continue;
      }
      std::size_t len = 0;
      uint32_t cur = static_cast<uint32_t>(i);
      while (visited[cur] == 0) {
        visited[cur] = 1;
        cur = img_[cur];
        ++len;
      }
      if (len > 1) {
        parity ^= (len - 1) & std::size_t{1};
      }
    }
    return (parity == 0) ? 1 : -1;
  }

  bool IsEven() const { return Sign() == 1; }

  std::string ToString() const {
    const auto cycles = Cycles();
    if (cycles.empty()) {
      return "()";
    }
    std::string out;
    for (const auto& cycle : cycles) {
      out.push_back('(');
      for (std::size_t i = 0; i < cycle.size(); ++i) {
        if (i > 0) {
          out.push_back(' ');
        }
        out += std::to_string(cycle[i]);
      }
      out.push_back(')');
    }
    return out;
  }

  bool operator==(const Permutation& other) const noexcept = default;
  std::strong_ordering operator<=>(const Permutation& other) const noexcept =
      default;

  friend std::ostream& operator<<(std::ostream& os, const Permutation& perm) {
    return os << perm.ToString();
  }

 private:
  struct UncheckedTag {};

  Permutation(std::vector<uint32_t> img, UncheckedTag) noexcept
      : img_(std::move(img)) {}

  Permutation PowWithSign(uint64_t abs_exp, bool negative) const {
    const std::size_t n = img_.size();
    std::vector<uint32_t> out(n);
    std::vector<uint8_t> visited(n, 0);
    std::vector<uint32_t> cycle;
    for (std::size_t i = 0; i < n; ++i) {
      if (visited[i] != 0) {
        continue;
      }
      cycle.clear();
      uint32_t cur = static_cast<uint32_t>(i);
      while (visited[cur] == 0) {
        visited[cur] = 1;
        cycle.push_back(cur);
        cur = img_[cur];
      }
      const uint64_t len = static_cast<uint64_t>(cycle.size());
      const uint64_t rem = abs_exp % len;
      const std::size_t shift =
          static_cast<std::size_t>((negative && rem != 0) ? (len - rem) : rem);
      for (std::size_t j = 0; j < cycle.size(); ++j) {
        out[cycle[j]] = cycle[(j + shift) % cycle.size()];
      }
    }
    return Permutation(std::move(out), UncheckedTag{});
  }

  std::vector<uint32_t> img_;
};

}  // namespace algebra

using algebra::Permutation;

}  // namespace xg

namespace std {

template <>
struct hash<xg::algebra::Permutation> {
  size_t operator()(const xg::algebra::Permutation& perm) const noexcept {
    size_t seed = std::hash<uint32_t>{}(perm.Degree());
    for (uint32_t v : perm.Mapping()) {
      const size_t hv = std::hash<uint32_t>{}(v);
      seed ^= hv + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

}  // namespace std

#endif  // XGALOIS_ALGEBRA_PERMUTATION_HPP_
