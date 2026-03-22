#ifndef XGALOIS_CODING_SYNDROME_DECODER_HPP
#define XGALOIS_CODING_SYNDROME_DECODER_HPP

#include <map>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class SyndromeDecoder : public Decoder<GaloisField> {
 public:
  using element_type = typename Decoder<GaloisField>::element_type;
  using codeword_type = typename Decoder<GaloisField>::codeword_type;
  using message_type = typename Decoder<GaloisField>::message_type;
  using matrix_type = xg::garray<GaloisField>;
  using vector_type = xg::garray<GaloisField>;

  // Constructor
  explicit SyndromeDecoder(const AbstractCode<GaloisField>* code,
                           size_t max_errors = 1)
      : Decoder<GaloisField>(code), max_errors_(max_errors) {
    // Try to cast to linear code
    linear_code_ = dynamic_cast<const AbstractLinearCode<GaloisField>*>(code);
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
      corrected[i] = field->Sub(corrected[i], error_pattern(i));
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
  const AbstractLinearCode<GaloisField>* linear_code_;
  size_t max_errors_;

  // Syndrome table: maps syndrome to error pattern
  std::map<std::vector<element_type>, codeword_type>
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
    auto field = this->code_->Field();
    if (weight == 0) {
      // Zero error pattern
      codeword_type error_pattern =
          xg::linalg::zeros<GaloisField>({n}, field);
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
      GenerateFieldCombinations(n, q, positions, 0, std::vector<element_type>(),
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
                                 std::vector<element_type> field_values,
                                 const matrix_type& parity_check) {
    auto field = this->code_->Field();
    if (pos_index == positions.size()) {
      // Create error pattern
      codeword_type error_pattern =
          xg::linalg::zeros<GaloisField>({n}, field);
      for (size_t i = 0; i < positions.size(); ++i) {
        error_pattern(positions[i]) = field_values[i];
      }

      auto syndrome = ComputeSyndrome(error_pattern, parity_check);

      // Only store if syndrome is not already in table (coset leader)
      if (syndrome_table_.find(syndrome) == syndrome_table_.end()) {
        syndrome_table_[syndrome] = error_pattern;
      }
      return;
    }

    // Try all non-zero field elements
    for (size_t i = 1; i < q; ++i) {
      field_values.push_back(element_type(i, field));
      GenerateFieldCombinations(n, q, positions, pos_index + 1, field_values,
                                parity_check);
      field_values.pop_back();
    }
  }

  std::vector<element_type> ComputeSyndrome(
      const codeword_type& error_pattern,
      const matrix_type& parity_check) const {
    vector_type syndrome_vec = xg::linalg::dot(parity_check, error_pattern);

    std::vector<element_type> syndrome(syndrome_vec.size());
    for (size_t i = 0; i < syndrome_vec.size(); ++i) {
      syndrome[i] = syndrome_vec(i);
    }
    return syndrome;
  }

  codeword_type LookupErrorPattern(
      const vector_type& syndrome) const {
    std::vector<element_type> syndrome_vec(syndrome.size());
    for (size_t i = 0; i < syndrome.size(); ++i) {
      syndrome_vec[i] = syndrome(i);
    }

    auto it = syndrome_table_.find(syndrome_vec);
    if (it != syndrome_table_.end()) {
      return it->second;
    }

    // If syndrome not found, return zero error pattern (uncorrectable)
    return xg::linalg::zeros<GaloisField>({this->code_->Length()},
                                          this->code_->Field());
  }
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_SYNDROME_DECODER_HPP
