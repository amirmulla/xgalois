#ifndef XGALOIS_CODING_DECODER_SYNDROME_DECODER_HPP_
#define XGALOIS_CODING_DECODER_SYNDROME_DECODER_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

namespace detail {

template <typename StorageType>
struct RawVectorHash {
  std::size_t operator()(const std::vector<StorageType> &vec) const noexcept {
    std::size_t h = std::hash<std::size_t>{}(vec.size());
    for (const auto &elem : vec) {
      const std::size_t he = std::hash<StorageType>{}(elem);
      h ^= he + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return h;
  }
};

template <typename GaloisField>
inline std::vector<typename GaloisField::element_type>
EnumerateNonZeroRawElements(const GaloisField &field) {
  using StorageType = typename GaloisField::element_type;
  const uint64_t q = field.Order();
  if (q == 0) {
    throw std::runtime_error(
        "Field order too large to enumerate non-zero elements");
  }
  std::vector<StorageType> elems;
  elems.reserve(static_cast<std::size_t>(q - 1));
  if constexpr (std::is_integral_v<StorageType>) {
    for (uint64_t i = 0; i < q; ++i) {
      StorageType v = field.SetElementValue(static_cast<StorageType>(i));
      if (!field.IsZero(v)) {
        elems.push_back(v);
      }
    }
  } else {
    for (uint64_t i = 0; i < q; ++i) {
      StorageType v = field.FromInt(i);
      if (!field.IsZero(v)) {
        elems.push_back(std::move(v));
      }
    }
  }
  return elems;
}

}  // namespace detail

template <typename GaloisField>
class SyndromeDecoder : public Decoder<GaloisField> {
 public:
  using Base = Decoder<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using result_type = typename Base::result_type;
  using SparseError = std::vector<std::pair<std::size_t, storage_type>>;

  explicit SyndromeDecoder(const AbstractLinearCode<GaloisField> &code,
                           std::size_t max_errors = 0,
                           std::size_t max_table_entries = 1'000'000)
      : Base(code), linear_code_(&code), max_table_entries_(max_table_entries) {
    const std::size_t n = this->code_->Length();
    if (max_errors == 0) {
      const std::size_t d = linear_code_->MinimumDistance();
      max_errors_ = (d > 0) ? (d - 1) / 2 : 0;
    } else {
      if (max_errors > n) {
        throw std::invalid_argument(
            "SyndromeDecoder: max_errors cannot exceed code length");
      }
      max_errors_ = max_errors;
    }
    if (max_table_entries_ == 0) {
      throw std::invalid_argument(
          "SyndromeDecoder: max_table_entries must be positive");
    }

    BuildSyndromeTable();
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    auto field = this->code_->Field();
    if (!received_word.HasField() || received_word.Field() != field) {
      throw std::invalid_argument(
          "SyndromeDecoder::DecodeWithResult: field mismatch");
    }
    if (received_word.Dimension() != 1 ||
        received_word.Size() != this->code_->Length()) {
      throw std::invalid_argument(
          "SyndromeDecoder::DecodeWithResult: received word must be 1D with "
          "length equal to code length");
    }

    const auto syn = linear_code_->Syndrome(received_word);
    bool all_zero = true;
    for (const auto &v : syn.RawSpan()) {
      if (!field->IsZero(v)) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      auto msg = this->code_->Unencode(received_word);
      return result_type::Success(received_word, std::move(msg), 0);
    }

    std::vector<storage_type> key(syn.RawSpan().begin(), syn.RawSpan().end());
    auto it = syndrome_table_.find(key);
    if (it == syndrome_table_.end()) {
      return result_type::Failure(
          "Syndrome not found in lookup table (uncorrectable error)");
    }

    codeword_type corrected = received_word;
    const SparseError &err = it->second;
    for (const auto &[pos, val] : err) {
      corrected.Raw(pos) = field->Sub(corrected.Raw(pos), val);
    }

    auto msg = this->code_->Unencode(corrected);
    return result_type::Success(std::move(corrected), std::move(msg),
                                err.size());
  }

  std::string Name() const override { return "Syndrome"; }

  std::size_t DecodingRadius() const override { return max_errors_; }

  std::size_t MaxErrors() const noexcept { return max_errors_; }

  std::size_t TableSize() const noexcept { return syndrome_table_.size(); }

  std::string ToString() const override {
    return "Syndrome decoder for [" + std::to_string(this->code_->Length()) +
           ", " + std::to_string(linear_code_->Dimension()) +
           "] linear code (max errors: " + std::to_string(max_errors_) + ")";
  }

 private:
  void BuildSyndromeTable() {
    auto field = this->code_->Field();
    const auto &H = linear_code_->ParityCheckMatrix();
    const std::size_t r = H.Shape(0);
    const std::size_t n = this->code_->Length();
    const storage_type zero = field->AdditiveIdentity();

    std::vector<storage_type> zero_syn(r, zero);
    syndrome_table_.emplace(zero_syn, SparseError{});

    if (r == 0 || max_errors_ == 0) {
      return;
    }

    const uint64_t q = field->Order();
    if (q == 0 || (q - 1) > max_table_entries_) {
      throw std::runtime_error(
          "SyndromeDecoder: syndrome table would exceed max_table_entries (" +
          std::to_string(max_table_entries_) + ")");
    }

    // Check if the full syndrome space size q^r is small enough to stop early.
    std::size_t full_syndrome_space = 0;
    try {
      const uint64_t qr = utils::SafeIntegerPower(q, r);
      if (qr <= max_table_entries_) {
        full_syndrome_space = static_cast<std::size_t>(qr);
      }
    } catch (const std::overflow_error &) {
      full_syndrome_space = 0;
    }

    // Transpose H into contiguous column vectors so each incremental update
    // uses a single O(r) field->FmaArray call.
    std::vector<std::vector<storage_type>> h_cols(
        n, std::vector<storage_type>(r, zero));
    const storage_type *h_data = H.RawData().data();
    for (std::size_t row = 0; row < r; ++row) {
      const storage_type *h_row = h_data + row * n;
      for (std::size_t col = 0; col < n; ++col) {
        h_cols[col][row] = h_row[col];
      }
    }

    const std::vector<storage_type> nonzero_elems =
        detail::EnumerateNonZeroRawElements(*field);

    // Preallocate per-depth syndrome buffers to avoid heap churn in recursion.
    std::vector<std::vector<storage_type>> syn_stack(
        max_errors_ + 1, std::vector<storage_type>(r, zero));
    SparseError current_errors;
    current_errors.reserve(max_errors_);

    for (std::size_t weight = 1; weight <= max_errors_; ++weight) {
      if (full_syndrome_space != 0 &&
          syndrome_table_.size() >= full_syndrome_space) {
        break;
      }
      EnumerateWeight(*field, h_cols, nonzero_elems, n, r, weight, 0, 0,
                      full_syndrome_space, current_errors, syn_stack);
    }
  }

  void EnumerateWeight(const GaloisField &field,
                       const std::vector<std::vector<storage_type>> &h_cols,
                       const std::vector<storage_type> &nonzero_elems,
                       std::size_t n, std::size_t r, std::size_t remaining,
                       std::size_t start_pos, std::size_t depth,
                       std::size_t full_syndrome_space,
                       SparseError &current_errors,
                       std::vector<std::vector<storage_type>> &syn_stack) {
    if (full_syndrome_space != 0 &&
        syndrome_table_.size() >= full_syndrome_space) {
      return;
    }
    const std::size_t max_pos = n - remaining;
    const auto &prev_syn = syn_stack[depth];
    auto &next_syn = syn_stack[depth + 1];

    for (std::size_t pos = start_pos; pos <= max_pos; ++pos) {
      const std::span<const storage_type> col_span(h_cols[pos].data(), r);
      for (const auto &val : nonzero_elems) {
        std::copy(prev_syn.begin(), prev_syn.end(), next_syn.begin());
        field.FmaArray(std::span<storage_type>(next_syn.data(), r), col_span,
                       val);
        current_errors.emplace_back(pos, val);

        if (remaining == 1) {
          if (syndrome_table_.find(next_syn) == syndrome_table_.end()) {
            if (syndrome_table_.size() >= max_table_entries_) {
              throw std::runtime_error(
                  "SyndromeDecoder: syndrome table exceeded max_table_entries "
                  "(" +
                  std::to_string(max_table_entries_) + ")");
            }
            syndrome_table_.emplace(next_syn, current_errors);
            if (full_syndrome_space != 0 &&
                syndrome_table_.size() >= full_syndrome_space) {
              current_errors.pop_back();
              return;
            }
          }
        } else {
          EnumerateWeight(field, h_cols, nonzero_elems, n, r, remaining - 1,
                          pos + 1, depth + 1, full_syndrome_space,
                          current_errors, syn_stack);
          if (full_syndrome_space != 0 &&
              syndrome_table_.size() >= full_syndrome_space) {
            current_errors.pop_back();
            return;
          }
        }

        current_errors.pop_back();
      }
    }
  }

  const AbstractLinearCode<GaloisField> *linear_code_;
  std::size_t max_errors_;
  std::size_t max_table_entries_;
  std::unordered_map<std::vector<storage_type>, SparseError,
                     detail::RawVectorHash<storage_type>>
      syndrome_table_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_DECODER_SYNDROME_DECODER_HPP_
