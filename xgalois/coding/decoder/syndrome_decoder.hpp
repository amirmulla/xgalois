#ifndef XGALOIS_CODING_SYNDROME_DECODER_HPP
#define XGALOIS_CODING_SYNDROME_DECODER_HPP

#include <unordered_map>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename ElementType>
class SyndromeDecoder : public Decoder<ElementType> {
 public:
  using element_type = ElementType;
  using codeword_type = std::vector<ElementType>;
  using message_type = std::vector<ElementType>;
  using matrix_type = xg::linalg::Matrix<ElementType>;
  using vector_type = xg::linalg::Vector<ElementType>;

  // Constructor
  explicit SyndromeDecoder(const AbstractCode<ElementType>* code,
                           size_t max_errors = 1)
      : Decoder<ElementType>(code), max_errors_(max_errors) {
    // Try to cast to linear code
    linear_code_ = dynamic_cast<const AbstractLinearCode<ElementType>*>(code);
    if (!linear_code_) {
      throw std::invalid_argument("SyndromeDecoder requires a linear code");
    }

    // Build syndrome table
    BuildSyndromeTable();
  }

  // Decode to codeword using syndrome decoding
  codeword_type DecodeToCode(
      const codeword_type& received_word) const override {
    if (received_word.size() != this->code_->Length()) {
      throw std::invalid_argument(
          "Received word length must match code length");
    }

    // Compute syndrome
    auto syndrome = linear_code_->Syndrome(received_word);

    // Look up error pattern in syndrome table
    auto error_pattern = LookupErrorPattern(syndrome);

    // Correct the error
    auto field = this->code_->Field();
    codeword_type corrected = received_word;

    for (size_t i = 0; i < corrected.size(); ++i) {
      corrected[i] = field->Sub(corrected[i], error_pattern[i]);
    }

    return corrected;
  }

  // Decode to message
  message_type DecodeToMessage(
      const codeword_type& received_word) const override {
    auto corrected_codeword = DecodeToCode(received_word);

    // Use the encoder to unencode
    auto encoder = this->code_->GetEncoder();
    return encoder->Unencode(corrected_codeword);
  }

  // Get maximum number of errors this decoder can handle
  size_t MaxErrors() const { return max_errors_; }

  std::string ToString() const override {
    return "Syndrome Decoder for [" + std::to_string(this->code_->Length()) +
           ", " + std::to_string(linear_code_->Dimension()) + "] linear code" +
           " (max errors: " + std::to_string(max_errors_) + ")";
  }

 private:
  const AbstractLinearCode<ElementType>* linear_code_;
  size_t max_errors_;

  // Syndrome table: maps syndrome to error pattern
  std::unordered_map<std::vector<ElementType>, std::vector<ElementType>>
      syndrome_table_;

  void BuildSyndromeTable() {
    auto field = this->code_->Field();
    auto parity_check = linear_code_->ParityCheckMatrix();
    size_t n = this->code_->Length();
    size_t q = field->Order();

    // Generate all error patterns up to max_errors
    GenerateErrorPatterns(n, q, max_errors_);
  }

  void GenerateErrorPatterns(size_t n, size_t q, size_t max_weight) {
    auto field = this->code_->Field();
    auto parity_check = linear_code_->ParityCheckMatrix();

    // Generate all error patterns of weight 0 to max_weight
    for (size_t weight = 0; weight <= max_weight; ++weight) {
      GenerateErrorPatternsOfWeight(n, q, weight, parity_check);
    }
  }

  void GenerateErrorPatternsOfWeight(size_t n, size_t q, size_t weight,
                                     const matrix_type& parity_check) {
    if (weight == 0) {
      // Zero error pattern
      std::vector<ElementType> error_pattern(n, ElementType{});
      auto syndrome = ComputeSyndrome(error_pattern, parity_check);
      syndrome_table_[syndrome] = error_pattern;
      return;
    }

    // Generate all combinations of positions for errors
    std::vector<size_t> positions;
    GenerateErrorCombinations(n, q, weight, 0, positions, parity_check);
  }

  void GenerateErrorCombinations(size_t n, size_t q, size_t remaining_weight,
                                 size_t start_pos,
                                 std::vector<size_t>& positions,
                                 const matrix_type& parity_check) {
    if (remaining_weight == 0) {
      // Generate all non-zero field elements for these positions
      GenerateFieldCombinations(n, q, positions, 0, std::vector<ElementType>(),
                                parity_check);
      return;
    }

    for (size_t pos = start_pos; pos <= n - remaining_weight; ++pos) {
      positions.push_back(pos);
      GenerateErrorCombinations(n, q, remaining_weight - 1, pos + 1, positions,
                                parity_check);
      positions.pop_back();
    }
  }

  void GenerateFieldCombinations(size_t n, size_t q,
                                 const std::vector<size_t>& positions,
                                 size_t pos_index,
                                 std::vector<ElementType> field_values,
                                 const matrix_type& parity_check) {
    if (pos_index == positions.size()) {
      // Create error pattern
      std::vector<ElementType> error_pattern(n, ElementType{});
      for (size_t i = 0; i < positions.size(); ++i) {
        error_pattern[positions[i]] = field_values[i];
      }

      auto syndrome = ComputeSyndrome(error_pattern, parity_check);

      // Only store if syndrome is not already in table (coset leader)
      if (syndrome_table_.find(syndrome) == syndrome_table_.end()) {
        syndrome_table_[syndrome] = error_pattern;
      }
      return;
    }

    // Try all non-zero field elements
    auto field = this->code_->Field();
    for (size_t i = 1; i < q; ++i) {
      field_values.push_back(ElementType(i));
      GenerateFieldCombinations(n, q, positions, pos_index + 1, field_values,
                                parity_check);
      field_values.pop_back();
    }
  }

  std::vector<ElementType> ComputeSyndrome(
      const std::vector<ElementType>& error_pattern,
      const matrix_type& parity_check) const {
    vector_type error_vec(error_pattern);
    vector_type syndrome_vec = parity_check * error_vec;

    std::vector<ElementType> syndrome(syndrome_vec.Size());
    for (size_t i = 0; i < syndrome_vec.Size(); ++i) {
      syndrome[i] = syndrome_vec[i];
    }
    return syndrome;
  }

  std::vector<ElementType> LookupErrorPattern(
      const vector_type& syndrome) const {
    std::vector<ElementType> syndrome_vec(syndrome.Size());
    for (size_t i = 0; i < syndrome.Size(); ++i) {
      syndrome_vec[i] = syndrome[i];
    }

    auto it = syndrome_table_.find(syndrome_vec);
    if (it != syndrome_table_.end()) {
      return it->second;
    }

    // If syndrome not found, return zero error pattern (uncorrectable)
    return std::vector<ElementType>(this->code_->Length(), ElementType{});
  }
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_SYNDROME_DECODER_HPP
