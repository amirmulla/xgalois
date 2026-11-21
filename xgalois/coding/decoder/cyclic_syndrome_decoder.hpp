#ifndef XGALOIS_CODING_DECODER_CYCLIC_SYNDROME_DECODER_HPP_
#define XGALOIS_CODING_DECODER_CYCLIC_SYNDROME_DECODER_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoder/syndrome_decoder.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace coding {

// Cyclic syndrome-polynomial decoder using Meggitt cyclic-shift lookup:
// computes s(x) = r(x) mod g(x) and cyclically shifts s(x) <- x*s(x) mod g(x)
// while checking a compact table of coset leaders pinned at position n - 1.
template <typename GaloisField>
class CyclicSyndromeDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using polynomial_type = PolynomialDense<GaloisField>;
  using SparseError = std::vector<std::pair<std::size_t, storage_type>>;

  CyclicSyndromeDecoder(const AbstractLinearCode<GaloisField> &code,
                        polynomial_type generator_poly,
                        std::size_t max_errors = 0,
                        std::size_t max_table_entries = 1'000'000)
      : Base(code),
        n_(code.Length()),
        k_(code.Dimension()),
        r_(n_ - k_),
        max_table_entries_(max_table_entries) {
    auto field = this->code_->Field();
    if (!generator_poly.HasField() || generator_poly.Field() != field) {
      throw std::invalid_argument(
          "CyclicSyndromeDecoder: generator polynomial field mismatch");
    }
    if (generator_poly.IsZero() ||
        static_cast<std::size_t>(generator_poly.Degree()) != r_) {
      throw std::invalid_argument(
          "CyclicSyndromeDecoder: generator polynomial degree must equal "
          "n - k");
    }
    generator_poly_ = generator_poly.Monic();

    if (max_errors == 0) {
      const std::size_t d = code.MinimumDistance();
      max_errors_ = (d > 0) ? (d - 1) / 2 : 0;
    } else {
      if (max_errors > n_) {
        throw std::invalid_argument(
            "CyclicSyndromeDecoder: max_errors cannot exceed code length");
      }
      max_errors_ = max_errors;
    }
    if (max_table_entries_ == 0) {
      throw std::invalid_argument(
          "CyclicSyndromeDecoder: max_table_entries must be positive");
    }

    PrecomputeSyndromeBasis();
    BuildMeggittTable();
  }

  template <typename CyclicCodeType>
    requires requires(const CyclicCodeType &c) { c.GeneratorPolynomial(); }
  explicit CyclicSyndromeDecoder(const CyclicCodeType &code,
                                 std::size_t max_errors = 0,
                                 std::size_t max_table_entries = 1'000'000)
      : CyclicSyndromeDecoder(code, code.GeneratorPolynomial(), max_errors,
                              max_table_entries) {}

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "CyclicSyndromeDecoder::DecodeWithResult: field mismatch");
    }
    if (received_word.Dimension() != 1 || received_word.Size() != n_) {
      throw std::invalid_argument(
          "CyclicSyndromeDecoder::DecodeWithResult: received word must be 1D "
          "with length equal to code length");
    }

    if (r_ == 0) {
      auto msg = this->code_->Unencode(received_word);
      return result_type::Success(received_word, std::move(msg), 0);
    }

    std::vector<storage_type> syn = ComputeSyndromeVector(received_word);
    bool all_zero = true;
    for (const auto &v : syn) {
      if (!field->IsZero(v)) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      auto msg = this->code_->Unencode(received_word);
      return result_type::Success(received_word, std::move(msg), 0);
    }

    if (max_errors_ == 0) {
      return result_type::Failure(
          "Syndrome polynomial is non-zero and decoding radius is 0");
    }

    // Meggitt cyclic-shift search: track minimum-weight match across all n
    // shifts (breaking early on a weight-1 match).
    const SparseError *best_err = nullptr;
    std::size_t best_shift = 0;
    std::vector<storage_type> cur_syn = std::move(syn);
    std::vector<storage_type> next_syn(r_, field->AdditiveIdentity());

    for (std::size_t shift = 0; shift < n_; ++shift) {
      auto it = meggitt_table_.find(cur_syn);
      if (it != meggitt_table_.end()) {
        if (best_err == nullptr || it->second.size() < best_err->size()) {
          best_err = &it->second;
          best_shift = shift;
          if (best_err->size() == 1) {
            break;
          }
        }
      }
      if (shift + 1 < n_) {
        ShiftSyndromeByX(*field, cur_syn, next_syn);
        cur_syn.swap(next_syn);
      }
    }

    if (best_err == nullptr) {
      return result_type::Failure(
          "Syndrome polynomial not found in lookup table (uncorrectable "
          "error)");
    }

    codeword_type corrected = received_word;
    for (const auto &[shifted_pos, val] : *best_err) {
      const std::size_t orig_pos = (shifted_pos + n_ - best_shift) % n_;
      corrected.Raw(orig_pos) = field->Sub(corrected.Raw(orig_pos), val);
    }

    if (!this->code_->Contains(corrected)) {
      return result_type::Failure(
          "CyclicSyndromeDecoder: corrected word failed codeword check");
    }

    auto msg = this->code_->Unencode(corrected);
    return result_type::Success(std::move(corrected), std::move(msg),
                                best_err->size());
  }

  std::string Name() const override { return "CyclicSyndrome"; }

  std::size_t DecodingRadius() const override { return max_errors_; }

  std::size_t MaxErrors() const noexcept { return max_errors_; }

  std::size_t TableSize() const noexcept { return meggitt_table_.size(); }

  std::string ToString() const override {
    return "CyclicSyndrome decoder for [" + std::to_string(n_) + ", " +
           std::to_string(k_) +
           "] cyclic code (max errors: " + std::to_string(max_errors_) + ")";
  }

 private:
  void PrecomputeSyndromeBasis() {
    auto field = this->code_->Field();
    const storage_type zero = field->AdditiveIdentity();
    const storage_type one = field->MultiplicativeIdentity();

    neg_g_low_.assign(r_, zero);
    for (std::size_t i = 0; i < r_; ++i) {
      neg_g_low_[i] = field->Neg(generator_poly_.RawCoeff(i));
    }

    x_pow_mod_g_.assign(n_, std::vector<storage_type>(r_, zero));
    if (r_ == 0) {
      return;
    }

    // x^0 mod g(x) = 1.
    x_pow_mod_g_[0][0] = one;
    for (std::size_t j = 1; j < n_; ++j) {
      ShiftSyndromeByX(*field, x_pow_mod_g_[j - 1], x_pow_mod_g_[j]);
    }
  }

  // Computes next(x) = (x * cur(x)) mod g_monic(x) in O(r).
  void ShiftSyndromeByX(const GaloisField &field,
                        const std::vector<storage_type> &cur,
                        std::vector<storage_type> &next) const {
    const storage_type top = cur[r_ - 1];
    next[0] = field.AdditiveIdentity();
    for (std::size_t i = 1; i < r_; ++i) {
      next[i] = cur[i - 1];
    }
    if (!field.IsZero(top)) {
      field.FmaArray(std::span<storage_type>(next.data(), r_),
                     std::span<const storage_type>(neg_g_low_.data(), r_), top);
    }
  }

  std::vector<storage_type> ComputeSyndromeVector(
      const codeword_type &word) const {
    auto field = this->code_->Field();
    const storage_type zero = field->AdditiveIdentity();
    std::vector<storage_type> syn(r_, zero);
    const storage_type *w_data = word.RawData().data();
    std::span<storage_type> syn_span(syn.data(), r_);
    for (std::size_t j = 0; j < n_; ++j) {
      if (!field->IsZero(w_data[j])) {
        field->FmaArray(
            syn_span, std::span<const storage_type>(x_pow_mod_g_[j].data(), r_),
            w_data[j]);
      }
    }
    return syn;
  }

  void BuildMeggittTable() {
    if (r_ == 0 || max_errors_ == 0) {
      return;
    }
    auto field = this->code_->Field();
    const uint64_t q = field->Order();
    if (q == 0 || (q - 1) > max_table_entries_) {
      throw std::runtime_error(
          "CyclicSyndromeDecoder: table would exceed max_table_entries (" +
          std::to_string(max_table_entries_) + ")");
    }

    const storage_type zero = field->AdditiveIdentity();
    const std::vector<storage_type> nonzero_elems =
        detail::EnumerateNonZeroRawElements(*field);

    std::vector<std::vector<storage_type>> syn_stack(
        max_errors_ + 1, std::vector<storage_type>(r_, zero));
    SparseError current_errors;
    current_errors.reserve(max_errors_);

    // Store all error patterns of weight w = 1 .. max_errors_ that have a
    // non-zero error at position n_ - 1.
    const std::span<const storage_type> last_col(x_pow_mod_g_[n_ - 1].data(),
                                                 r_);
    for (std::size_t weight = 1; weight <= max_errors_; ++weight) {
      for (const auto &last_val : nonzero_elems) {
        std::fill(syn_stack[1].begin(), syn_stack[1].end(), zero);
        field->FmaArray(std::span<storage_type>(syn_stack[1].data(), r_),
                        last_col, last_val);
        current_errors.clear();
        current_errors.emplace_back(n_ - 1, last_val);

        if (weight == 1) {
          InsertTableEntry(syn_stack[1], current_errors);
        } else {
          EnumeratePrefix(*field, nonzero_elems, weight - 1, 0, 1,
                          current_errors, syn_stack);
        }
      }
    }
  }

  void InsertTableEntry(const std::vector<storage_type> &syn,
                        const SparseError &err) {
    if (meggitt_table_.find(syn) != meggitt_table_.end()) {
      return;
    }
    if (meggitt_table_.size() >= max_table_entries_) {
      throw std::runtime_error(
          "CyclicSyndromeDecoder: syndrome table exceeded max_table_entries (" +
          std::to_string(max_table_entries_) + ")");
    }
    meggitt_table_.emplace(syn, err);
  }

  void EnumeratePrefix(const GaloisField &field,
                       const std::vector<storage_type> &nonzero_elems,
                       std::size_t remaining, std::size_t start_pos,
                       std::size_t depth, SparseError &current_errors,
                       std::vector<std::vector<storage_type>> &syn_stack) {
    // Choose `remaining` positions strictly before n_ - 1.
    const std::size_t max_pos = (n_ - 1) - remaining;
    const auto &prev_syn = syn_stack[depth];
    auto &next_syn = syn_stack[depth + 1];

    for (std::size_t pos = start_pos; pos <= max_pos; ++pos) {
      const std::span<const storage_type> col_span(x_pow_mod_g_[pos].data(),
                                                   r_);
      for (const auto &val : nonzero_elems) {
        std::copy(prev_syn.begin(), prev_syn.end(), next_syn.begin());
        field.FmaArray(std::span<storage_type>(next_syn.data(), r_), col_span,
                       val);
        current_errors.emplace_back(pos, val);

        if (remaining == 1) {
          InsertTableEntry(next_syn, current_errors);
        } else {
          EnumeratePrefix(field, nonzero_elems, remaining - 1, pos + 1,
                          depth + 1, current_errors, syn_stack);
        }

        current_errors.pop_back();
      }
    }
  }

  polynomial_type generator_poly_;
  std::size_t n_;
  std::size_t k_;
  std::size_t r_;
  std::size_t max_errors_;
  std::size_t max_table_entries_;
  std::vector<storage_type> neg_g_low_;
  std::vector<std::vector<storage_type>> x_pow_mod_g_;
  std::unordered_map<std::vector<storage_type>, SparseError,
                     detail::RawVectorHash<storage_type>>
      meggitt_table_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_CYCLIC_SYNDROME_DECODER_HPP_
