#ifndef XGALOIS_CODING_CONVOLUTIONAL_HPP_
#define XGALOIS_CODING_CONVOLUTIONAL_HPP_

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/coding/decoder/decoder.hpp"
#include "xgalois/coding/decoding_result.hpp"
#include "xgalois/coding/encoder/encoder.hpp"
#include "xgalois/coding/encoder/generator_matrix_encoder.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/linalg/linalg.hpp"

namespace xg {
namespace coding {

enum class TerminationMode : uint8_t {
  kZeroTailed,
  kTruncated,
  kTailBiting,
};

enum class BcjrMode : uint8_t {
  kLogMap,
  kMaxLogMap,
};

namespace conv_detail {

constexpr double kMaxFiniteLlr = 1e6;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();
constexpr std::size_t kMaxConvolutionalMemory = 16;

inline double ClampLlr(double x) noexcept {
  if (x > kMaxFiniteLlr) {
    return kMaxFiniteLlr;
  }
  if (x < -kMaxFiniteLlr) {
    return -kMaxFiniteLlr;
  }
  return x;
}

inline void ValidateLlrVector(const std::vector<double> &llrs,
                              std::size_t expected_size, const char *fn_name) {
  if (llrs.size() != expected_size) {
    throw std::invalid_argument(std::string(fn_name) + ": LLR length (" +
                                std::to_string(llrs.size()) +
                                ") does not match expected length (" +
                                std::to_string(expected_size) + ")");
  }
  for (double v : llrs) {
    if (std::isnan(v)) {
      throw std::invalid_argument(std::string(fn_name) +
                                  ": LLR vector must not contain NaN");
    }
  }
}

inline double MaxStar(double a, double b, BcjrMode mode) noexcept {
  if (a == kNegInf) {
    return b;
  }
  if (b == kNegInf) {
    return a;
  }
  const double m = (a > b) ? a : b;
  if (mode == BcjrMode::kMaxLogMap) {
    return m;
  }
  return m + std::log1p(std::exp(-std::abs(a - b)));
}

inline GFArray<GaloisFieldBinary> BuildBinaryArray(
    const std::shared_ptr<GaloisFieldBinary> &field,
    std::span<const uint8_t> bits) {
  using RawArray = GFArray<GaloisFieldBinary>::raw_array_type;
  using ShapeType = GFArray<GaloisFieldBinary>::shape_type;
  RawArray raw = RawArray::from_shape(ShapeType{bits.size()});
  uint8_t *dst = raw.data();
  for (std::size_t i = 0; i < bits.size(); ++i) {
    dst[i] = static_cast<uint8_t>(bits[i] & 1U);
  }
  return GFArray<GaloisFieldBinary>::FromRaw(field, std::move(raw));
}

// Multiply an m x m binary matrix (stored as m row bitmasks) by an m-bit state
// vector over GF(2).
inline uint32_t MatVecMulF2(const std::vector<uint32_t> &mat,
                            uint32_t vec) noexcept {
  uint32_t out = 0;
  for (std::size_t r = 0; r < mat.size(); ++r) {
    if ((static_cast<uint32_t>(std::popcount(mat[r] & vec)) & 1U) != 0U) {
      out |= (uint32_t{1} << r);
    }
  }
  return out;
}

// Multiply two m x m binary matrices (row bitmasks) over GF(2).
inline std::vector<uint32_t> MatMulF2(const std::vector<uint32_t> &a,
                                      const std::vector<uint32_t> &b) {
  const std::size_t m = a.size();
  std::vector<uint32_t> c(m, 0);
  for (std::size_t i = 0; i < m; ++i) {
    uint32_t row = 0;
    for (std::size_t k = 0; k < m; ++k) {
      if (((a[i] >> k) & 1U) != 0U) {
        row ^= b[k];
      }
    }
    c[i] = row;
  }
  return c;
}

// Binary matrix exponentiation A^exp over GF(2) using square-and-multiply.
inline std::vector<uint32_t> MatPowF2(std::vector<uint32_t> base,
                                      std::size_t exp) {
  const std::size_t m = base.size();
  std::vector<uint32_t> res(m, 0);
  for (std::size_t i = 0; i < m; ++i) {
    res[i] = uint32_t{1} << i;
  }
  while (exp > 0) {
    if ((exp & 1U) != 0U) {
      res = MatMulF2(res, base);
    }
    exp >>= 1;
    if (exp > 0) {
      base = MatMulF2(base, base);
    }
  }
  return res;
}

// Invert an m x m binary matrix over GF(2). Returns std::nullopt if singular.
inline std::optional<std::vector<uint32_t>> InvertMatrixF2(
    std::vector<uint32_t> mat) {
  const std::size_t m = mat.size();
  std::vector<uint32_t> inv(m, 0);
  for (std::size_t i = 0; i < m; ++i) {
    inv[i] = uint32_t{1} << i;
  }
  for (std::size_t col = 0; col < m; ++col) {
    std::size_t pivot = col;
    while (pivot < m && ((mat[pivot] >> col) & 1U) == 0U) {
      ++pivot;
    }
    if (pivot == m) {
      return std::nullopt;
    }
    if (pivot != col) {
      std::swap(mat[pivot], mat[col]);
      std::swap(inv[pivot], inv[col]);
    }
    for (std::size_t row = 0; row < m; ++row) {
      if (row != col && ((mat[row] >> col) & 1U) != 0U) {
        mat[row] ^= mat[col];
        inv[row] ^= inv[col];
      }
    }
  }
  return inv;
}

}  // namespace conv_detail

// Flat, cache-friendly trellis transition table for a rate k0/n0 binary
// convolutional encoder (k0 = 1).
struct Trellis {
  struct PrevTransition {
    uint32_t prev_state = 0;
    uint32_t input = 0;
    uint32_t output_mask = 0;
  };

  std::size_t num_states = 1;
  std::size_t num_input_bits = 1;
  std::size_t num_output_bits = 0;
  std::size_t memory = 0;
  bool is_recursive = false;

  // Flat tables indexed by [state * (1 << num_input_bits) + input].
  std::vector<uint32_t> next_state;
  std::vector<uint32_t> output_mask;
  // Incoming transitions grouped by destination state.
  std::vector<std::vector<PrevTransition>> prev_transitions;
  // For each state s, the input bit that drives the shift-register input to 0
  // (always 0 for feedforward; feedback parity for RSC).
  std::vector<uint32_t> tail_input;

  std::size_t NumInputs() const noexcept {
    return std::size_t{1} << num_input_bits;
  }

  uint32_t NextState(uint32_t state, uint32_t input) const noexcept {
    return next_state[static_cast<std::size_t>(state) * NumInputs() + input];
  }

  uint32_t OutputMask(uint32_t state, uint32_t input) const noexcept {
    return output_mask[static_cast<std::size_t>(state) * NumInputs() + input];
  }
};

struct BcjrOutput {
  std::vector<double> app_llrs;
  std::vector<double> extrinsic_llrs;
  GFArray<GaloisFieldBinary> hard_bits;
};

class ConvolutionalCode;
class ConvolutionalEncoder;
class ViterbiDecoder;
class BcjrDecoder;

// Binary rate-1/n0 convolutional code (feedforward or recursive systematic)
// supporting zero-tailed, truncated, and tail-biting termination over GF(2).
//
// Polynomial bit convention (matching standard octal notation such as NASA
// {0133, 0171} and 3GPP LTE {013, 015}):
//   - Bit `memory` (MSB) is the D^0 (current input/feedback) tap.
//   - Bit `memory - i` is the D^i tap for i = 1..memory.
//   - Register state `s` stores (d_{t-1}, ..., d_{t-m}) with d_{t-1} in bit
//     `memory - 1` and d_{t-m} in bit 0.
class ConvolutionalCode : public AbstractLinearCode<GaloisFieldBinary> {
 public:
  using Base = AbstractLinearCode<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using element_type = GaloisFieldElement<GaloisFieldBinary>;
  using storage_type = typename Base::storage_type;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using matrix_type = GFMatrix<GaloisFieldBinary>;
  using syndrome_type = GFArray<GaloisFieldBinary>;

  // Create a feedforward rate-1/n0 convolutional code from generator
  // polynomials (e.g. NASA K=7, R=1/2: {0133, 0171}).
  static ConvolutionalCode CreateFeedforward(
      std::shared_ptr<GaloisFieldBinary> field, std::size_t input_block_len,
      const std::vector<uint32_t> &generators,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    if (!field) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateFeedforward: field cannot be null");
    }
    if (input_block_len == 0) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateFeedforward: input_block_len must be > 0");
    }
    if (generators.empty()) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateFeedforward: generators cannot be empty");
    }
    uint32_t max_poly = 0;
    for (uint32_t g : generators) {
      if (g == 0) {
        throw std::invalid_argument(
            "ConvolutionalCode::CreateFeedforward: generator polynomial cannot "
            "be zero");
      }
      max_poly = std::max(max_poly, g);
    }
    const std::size_t memory =
        static_cast<std::size_t>(std::bit_width(max_poly) - 1);
    if (memory > conv_detail::kMaxConvolutionalMemory) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateFeedforward: memory exceeds supported "
          "maximum (16)");
    }
    Trellis trellis = BuildFeedforwardTrellis(memory, generators);
    return ConvolutionalCode(std::move(field), input_block_len,
                             std::move(trellis), termination, 0, generators);
  }

  static ConvolutionalCode CreateFeedforward(
      std::size_t input_block_len, const std::vector<uint32_t> &generators,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    return CreateFeedforward(std::make_shared<GaloisFieldBinary>(),
                             input_block_len, generators, termination);
  }

  // Create a recursive systematic convolutional (RSC) code where output bit 0
  // is the systematic input and remaining output bits are parity streams
  // generated by `feedforward_polys` with feedback `feedback_poly`
  // (e.g. 3GPP LTE: feedback = 013, feedforward = {015}).
  static ConvolutionalCode CreateRecursiveSystematic(
      std::shared_ptr<GaloisFieldBinary> field, std::size_t input_block_len,
      uint32_t feedback_poly, const std::vector<uint32_t> &feedforward_polys,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    if (!field) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateRecursiveSystematic: field cannot be null");
    }
    if (input_block_len == 0) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateRecursiveSystematic: input_block_len must "
          "be > 0");
    }
    if (feedback_poly == 0) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateRecursiveSystematic: feedback_poly cannot "
          "be zero");
    }
    if (feedforward_polys.empty()) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateRecursiveSystematic: feedforward_polys "
          "cannot be empty");
    }
    uint32_t max_poly = feedback_poly;
    for (uint32_t g : feedforward_polys) {
      if (g == 0) {
        throw std::invalid_argument(
            "ConvolutionalCode::CreateRecursiveSystematic: feedforward "
            "polynomial cannot be zero");
      }
      max_poly = std::max(max_poly, g);
    }
    const std::size_t memory =
        static_cast<std::size_t>(std::bit_width(max_poly) - 1);
    if (memory > conv_detail::kMaxConvolutionalMemory) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateRecursiveSystematic: memory exceeds "
          "supported maximum (16)");
    }
    if (((feedback_poly >> memory) & 1U) == 0U) {
      throw std::invalid_argument(
          "ConvolutionalCode::CreateRecursiveSystematic: feedback_poly must "
          "have degree equal to encoder memory");
    }
    Trellis trellis = BuildRscTrellis(memory, feedback_poly, feedforward_polys);
    return ConvolutionalCode(std::move(field), input_block_len,
                             std::move(trellis), termination, feedback_poly,
                             feedforward_polys);
  }

  static ConvolutionalCode CreateRecursiveSystematic(
      std::size_t input_block_len, uint32_t feedback_poly,
      const std::vector<uint32_t> &feedforward_polys,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    return CreateRecursiveSystematic(std::make_shared<GaloisFieldBinary>(),
                                     input_block_len, feedback_poly,
                                     feedforward_polys, termination);
  }

  // Convenience preset: NASA K=7, R=1/2 feedforward code {0133, 0171} octal.
  static ConvolutionalCode NasaRateHalf(
      std::shared_ptr<GaloisFieldBinary> field, std::size_t input_block_len,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    return CreateFeedforward(std::move(field), input_block_len, {0133U, 0171U},
                             termination);
  }

  static ConvolutionalCode NasaRateHalf(
      std::size_t input_block_len,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    return NasaRateHalf(std::make_shared<GaloisFieldBinary>(), input_block_len,
                        termination);
  }

  // Convenience preset: 3GPP LTE Turbo constituent RSC (feedback=013,
  // feedforward={015} octal, 8 states).
  static ConvolutionalCode LteTurboConstituent(
      std::shared_ptr<GaloisFieldBinary> field, std::size_t input_block_len,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    return CreateRecursiveSystematic(std::move(field), input_block_len, 013U,
                                     {015U}, termination);
  }

  static ConvolutionalCode LteTurboConstituent(
      std::size_t input_block_len,
      TerminationMode termination = TerminationMode::kZeroTailed) {
    return LteTurboConstituent(std::make_shared<GaloisFieldBinary>(),
                               input_block_len, termination);
  }

  // Rule of Five (mutex is not copyable/movable).
  ConvolutionalCode(const ConvolutionalCode &other)
      : Base(other),
        input_block_len_(other.input_block_len_),
        num_steps_(other.num_steps_),
        termination_(other.termination_),
        feedback_poly_(other.feedback_poly_),
        feedforward_polys_(other.feedforward_polys_),
        trellis_(other.trellis_),
        tail_biting_inv_matrix_(other.tail_biting_inv_matrix_) {
    std::lock_guard<std::mutex> lock(other.matrix_mutex_);
    generator_matrix_ = other.generator_matrix_;
    parity_check_matrix_ = other.parity_check_matrix_;
  }

  ConvolutionalCode &operator=(const ConvolutionalCode &other) {
    if (this != &other) {
      Base::operator=(other);
      input_block_len_ = other.input_block_len_;
      num_steps_ = other.num_steps_;
      termination_ = other.termination_;
      feedback_poly_ = other.feedback_poly_;
      feedforward_polys_ = other.feedforward_polys_;
      trellis_ = other.trellis_;
      tail_biting_inv_matrix_ = other.tail_biting_inv_matrix_;
      std::scoped_lock lock(matrix_mutex_, other.matrix_mutex_);
      generator_matrix_ = other.generator_matrix_;
      parity_check_matrix_ = other.parity_check_matrix_;
    }
    return *this;
  }

  ConvolutionalCode(ConvolutionalCode &&other) noexcept
      : Base(std::move(other)),
        input_block_len_(other.input_block_len_),
        num_steps_(other.num_steps_),
        termination_(other.termination_),
        feedback_poly_(other.feedback_poly_),
        feedforward_polys_(std::move(other.feedforward_polys_)),
        trellis_(std::move(other.trellis_)),
        tail_biting_inv_matrix_(std::move(other.tail_biting_inv_matrix_)),
        generator_matrix_(std::move(other.generator_matrix_)),
        parity_check_matrix_(std::move(other.parity_check_matrix_)) {}

  ConvolutionalCode &operator=(ConvolutionalCode &&other) noexcept {
    if (this != &other) {
      Base::operator=(std::move(other));
      input_block_len_ = other.input_block_len_;
      num_steps_ = other.num_steps_;
      termination_ = other.termination_;
      feedback_poly_ = other.feedback_poly_;
      feedforward_polys_ = std::move(other.feedforward_polys_);
      trellis_ = std::move(other.trellis_);
      tail_biting_inv_matrix_ = std::move(other.tail_biting_inv_matrix_);
      generator_matrix_ = std::move(other.generator_matrix_);
      parity_check_matrix_ = std::move(other.parity_check_matrix_);
    }
    return *this;
  }

  ~ConvolutionalCode() override = default;

  const Trellis &GetTrellis() const noexcept { return trellis_; }
  TerminationMode GetTerminationMode() const noexcept { return termination_; }
  std::size_t InputBlockLength() const noexcept { return input_block_len_; }
  std::size_t NumSteps() const noexcept { return num_steps_; }
  std::size_t Memory() const noexcept { return trellis_.memory; }
  std::size_t NumOutputBits() const noexcept {
    return trellis_.num_output_bits;
  }
  bool IsRecursiveSystematic() const noexcept { return trellis_.is_recursive; }
  uint32_t FeedbackPoly() const noexcept { return feedback_poly_; }
  const std::vector<uint32_t> &FeedforwardPolys() const noexcept {
    return feedforward_polys_;
  }

  // Compute the exact circular initial state S_0 = (I - A^K)^{-1} S_zero_init
  // over GF(2) for kTailBiting mode.
  uint32_t ComputeTailBitingInitialState(
      std::span<const uint8_t> message_bits) const {
    if (message_bits.size() != input_block_len_) {
      throw std::invalid_argument(
          "ConvolutionalCode::ComputeTailBitingInitialState: message length "
          "mismatch");
    }
    if (trellis_.memory == 0) {
      return 0;
    }
    uint32_t state = 0;
    for (std::size_t t = 0; t < input_block_len_; ++t) {
      state = trellis_.NextState(state,
                                 static_cast<uint32_t>(message_bits[t] & 1U));
    }
    return conv_detail::MatVecMulF2(tail_biting_inv_matrix_, state);
  }

  // Fast O(K) shift-register encoding on raw binary bits (0/1).
  std::vector<uint8_t> EncodeBits(std::span<const uint8_t> message_bits) const {
    if (message_bits.size() != input_block_len_) {
      throw std::invalid_argument(
          "ConvolutionalCode::EncodeBits: message length mismatch");
    }
    const std::size_t n0 = trellis_.num_output_bits;
    std::vector<uint8_t> out(this->length_, 0);

    uint32_t state = 0;
    if (termination_ == TerminationMode::kTailBiting) {
      state = ComputeTailBitingInitialState(message_bits);
    }

    std::size_t pos = 0;
    for (std::size_t t = 0; t < input_block_len_; ++t) {
      const uint32_t u = static_cast<uint32_t>(message_bits[t] & 1U);
      const uint32_t mask = trellis_.OutputMask(state, u);
      state = trellis_.NextState(state, u);
      for (std::size_t j = 0; j < n0; ++j) {
        out[pos++] = static_cast<uint8_t>((mask >> j) & 1U);
      }
    }

    if (termination_ == TerminationMode::kZeroTailed) {
      for (std::size_t t = 0; t < trellis_.memory; ++t) {
        const uint32_t u = trellis_.tail_input[state];
        const uint32_t mask = trellis_.OutputMask(state, u);
        state = trellis_.NextState(state, u);
        for (std::size_t j = 0; j < n0; ++j) {
          out[pos++] = static_cast<uint8_t>((mask >> j) & 1U);
        }
      }
    }

    return out;
  }

  // Fast O(K) shift-register encoding returning a GFArray<GaloisFieldBinary>.
  codeword_type EncodeFast(const message_type &message) const {
    if (!message.HasField() || message.Field() != this->field_) {
      throw std::invalid_argument(
          "ConvolutionalCode::EncodeFast: field mismatch");
    }
    if (message.Dimension() != 1 || message.Size() != input_block_len_) {
      throw std::invalid_argument(
          "ConvolutionalCode::EncodeFast: message must be 1D with size equal "
          "to Dimension()");
    }
    std::vector<uint8_t> bits = EncodeBits(message.RawSpan());
    return conv_detail::BuildBinaryArray(this->field_, bits);
  }

  // Fast unencoding with exact codeword membership verification.
  inline message_type UnencodeFast(const codeword_type &codeword) const;

  // Lazy O(K * N) generator matrix construction (only built if requested).
  const matrix_type &GeneratorMatrix() const override {
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!generator_matrix_.has_value()) {
      const std::size_t k = input_block_len_;
      const std::size_t n = this->length_;
      using RawArray = matrix_type::raw_array_type;
      using ShapeType = matrix_type::shape_type;
      RawArray g_raw = RawArray::from_shape(ShapeType{k, n});
      uint8_t *g_data = g_raw.data();
      std::vector<uint8_t> basis(k, 0);
      for (std::size_t r = 0; r < k; ++r) {
        basis[r] = 1U;
        std::vector<uint8_t> cw = EncodeBits(basis);
        std::copy(cw.begin(), cw.end(), g_data + r * n);
        basis[r] = 0U;
      }
      generator_matrix_ = matrix_type::FromRaw(this->field_, std::move(g_raw));
    }
    return *generator_matrix_;
  }

  // Lazy parity-check matrix construction via GF(2) nullspace of G.
  const matrix_type &ParityCheckMatrix() const override {
    const matrix_type &g = GeneratorMatrix();
    std::lock_guard<std::mutex> lock(matrix_mutex_);
    if (!parity_check_matrix_.has_value()) {
      parity_check_matrix_ = Base::ComputeParityCheckFromGenerator(g);
    }
    return *parity_check_matrix_;
  }

  // O(K) random codeword sampling without materializing GeneratorMatrix().
  codeword_type RandomCodeword() const override {
    auto msg = message_type::Random(this->field_, {input_block_len_});
    return EncodeFast(msg);
  }

  codeword_type RandomCodeword(std::mt19937_64 &rng) const override {
    auto msg = message_type::Random(this->field_, {input_block_len_}, rng);
    return EncodeFast(msg);
  }

  std::string ToString() const override {
    std::string mode_str = "zero-tailed";
    if (termination_ == TerminationMode::kTruncated) {
      mode_str = "truncated";
    } else if (termination_ == TerminationMode::kTailBiting) {
      mode_str = "tail-biting";
    }
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(this->dimension_) + "] " +
           (trellis_.is_recursive ? "RSC " : "feedforward ") + mode_str +
           " convolutional code (memory=" + std::to_string(trellis_.memory) +
           ", rate=1/" + std::to_string(trellis_.num_output_bits) + ")";
  }

 private:
  static std::size_t ComputeCodeLength(std::size_t input_block_len,
                                       std::size_t memory,
                                       std::size_t num_output_bits,
                                       TerminationMode termination) {
    const std::size_t steps = (termination == TerminationMode::kZeroTailed)
                                  ? (input_block_len + memory)
                                  : input_block_len;
    return steps * num_output_bits;
  }

  static Trellis BuildFeedforwardTrellis(
      std::size_t memory, const std::vector<uint32_t> &generators) {
    Trellis t;
    t.memory = memory;
    t.num_states = std::size_t{1} << memory;
    t.num_input_bits = 1;
    t.num_output_bits = generators.size();
    t.is_recursive = false;
    t.next_state.resize(t.num_states * 2, 0);
    t.output_mask.resize(t.num_states * 2, 0);
    t.prev_transitions.resize(t.num_states);
    t.tail_input.assign(t.num_states, 0U);

    for (uint32_t s = 0; s < static_cast<uint32_t>(t.num_states); ++s) {
      for (uint32_t u = 0; u < 2U; ++u) {
        const uint32_t reg = (u << memory) | s;
        uint32_t out_mask = 0;
        for (std::size_t j = 0; j < generators.size(); ++j) {
          const uint32_t bit =
              static_cast<uint32_t>(std::popcount(reg & generators[j])) & 1U;
          out_mask |= (bit << j);
        }
        const uint32_t ns =
            (memory == 0) ? 0U : ((u << (memory - 1)) | (s >> 1));
        const std::size_t idx = static_cast<std::size_t>(s) * 2 + u;
        t.next_state[idx] = ns;
        t.output_mask[idx] = out_mask;
        t.prev_transitions[ns].push_back({s, u, out_mask});
      }
    }
    return t;
  }

  static Trellis BuildRscTrellis(
      std::size_t memory, uint32_t feedback_poly,
      const std::vector<uint32_t> &feedforward_polys) {
    Trellis t;
    t.memory = memory;
    t.num_states = std::size_t{1} << memory;
    t.num_input_bits = 1;
    t.num_output_bits = 1 + feedforward_polys.size();
    t.is_recursive = true;
    t.next_state.resize(t.num_states * 2, 0);
    t.output_mask.resize(t.num_states * 2, 0);
    t.prev_transitions.resize(t.num_states);
    t.tail_input.resize(t.num_states, 0U);

    const uint32_t state_mask =
        (memory == 0) ? 0U : ((uint32_t{1} << memory) - 1U);
    const uint32_t fb_state_taps = feedback_poly & state_mask;

    for (uint32_t s = 0; s < static_cast<uint32_t>(t.num_states); ++s) {
      const uint32_t fb_bit =
          static_cast<uint32_t>(std::popcount(s & fb_state_taps)) & 1U;
      // Setting u = fb_bit forces d_t = u ^ fb_bit = 0, shifting 0 into state.
      t.tail_input[s] = fb_bit;

      for (uint32_t u = 0; u < 2U; ++u) {
        const uint32_t d_t = u ^ fb_bit;
        const uint32_t reg = (d_t << memory) | s;
        // Output bit 0 is systematic input u.
        uint32_t out_mask = u;
        for (std::size_t j = 0; j < feedforward_polys.size(); ++j) {
          const uint32_t p =
              static_cast<uint32_t>(std::popcount(reg & feedforward_polys[j])) &
              1U;
          out_mask |= (p << (j + 1));
        }
        const uint32_t ns =
            (memory == 0) ? 0U : ((d_t << (memory - 1)) | (s >> 1));
        const std::size_t idx = static_cast<std::size_t>(s) * 2 + u;
        t.next_state[idx] = ns;
        t.output_mask[idx] = out_mask;
        t.prev_transitions[ns].push_back({s, u, out_mask});
      }
    }
    return t;
  }

  ConvolutionalCode(std::shared_ptr<GaloisFieldBinary> field,
                    std::size_t input_block_len, Trellis trellis,
                    TerminationMode termination, uint32_t feedback_poly,
                    std::vector<uint32_t> feedforward_polys)
      : Base(std::move(field),
             ComputeCodeLength(input_block_len, trellis.memory,
                               trellis.num_output_bits, termination),
             input_block_len, Metric::HAMMING, "Convolutional", "ViterbiHard"),
        input_block_len_(input_block_len),
        num_steps_(termination == TerminationMode::kZeroTailed
                       ? (input_block_len + trellis.memory)
                       : input_block_len),
        termination_(termination),
        feedback_poly_(feedback_poly),
        feedforward_polys_(std::move(feedforward_polys)),
        trellis_(std::move(trellis)) {
    if (termination_ == TerminationMode::kTailBiting && trellis_.memory > 0) {
      InitializeTailBitingMatrix();
    }
    RegisterEncodersAndDecoders();
  }

  void InitializeTailBitingMatrix() {
    const std::size_t m = trellis_.memory;
    // Build state transition matrix A over GF(2): s_{t+1} = A * s_t when u_t =
    // 0. Column c of A is NextState(1 << c, 0).
    std::vector<uint32_t> mat_a(m, 0);
    for (std::size_t c = 0; c < m; ++c) {
      const uint32_t next_s = trellis_.NextState(uint32_t{1} << c, 0U);
      for (std::size_t r = 0; r < m; ++r) {
        if (((next_s >> r) & 1U) != 0U) {
          mat_a[r] |= (uint32_t{1} << c);
        }
      }
    }
    std::vector<uint32_t> a_pow =
        conv_detail::MatPowF2(mat_a, input_block_len_);
    // In GF(2), I - A^K = I + A^K.
    for (std::size_t r = 0; r < m; ++r) {
      a_pow[r] ^= (uint32_t{1} << r);
    }
    auto inv = conv_detail::InvertMatrixF2(a_pow);
    if (!inv.has_value()) {
      throw std::invalid_argument(
          "ConvolutionalCode: (I - A^K) is singular over GF(2) for the "
          "requested tail-biting block length K");
    }
    tail_biting_inv_matrix_ = std::move(*inv);
  }

  inline void RegisterEncodersAndDecoders();

  std::size_t input_block_len_;
  std::size_t num_steps_;
  TerminationMode termination_;
  uint32_t feedback_poly_;
  std::vector<uint32_t> feedforward_polys_;
  Trellis trellis_;
  std::vector<uint32_t> tail_biting_inv_matrix_;

  mutable std::mutex matrix_mutex_;
  mutable std::optional<matrix_type> generator_matrix_;
  mutable std::optional<matrix_type> parity_check_matrix_;
};

// Fast O(K) shift-register encoder for ConvolutionalCode.
class ConvolutionalEncoder : public Encoder<GaloisFieldBinary> {
 public:
  using Base = Encoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;

  explicit ConvolutionalEncoder(const ConvolutionalCode &code)
      : Base(code), conv_code_(&code) {}

  codeword_type Encode(const message_type &message) const override {
    return conv_code_->EncodeFast(message);
  }

  message_type Unencode(const codeword_type &codeword) const override {
    return conv_code_->UnencodeFast(codeword);
  }

  std::string Name() const override { return "Convolutional"; }

 private:
  const ConvolutionalCode *conv_code_;
};

// Maximum-likelihood Viterbi decoder supporting hard-decision Hamming distance
// and soft-decision LLR correlation metrics for zero-tailed, truncated, and
// tail-biting convolutional codes.
class ViterbiDecoder : public Decoder<GaloisFieldBinary> {
 public:
  using Base = Decoder<GaloisFieldBinary>;
  using field_type = GaloisFieldBinary;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;
  using result_type = DecodingResult<GaloisFieldBinary>;

  explicit ViterbiDecoder(const ConvolutionalCode &code)
      : Base(code), conv_code_(&code) {}

  result_type DecodeHard(const codeword_type &received) const {
    auto field = this->code_->Field();
    if (!received.HasField() || received.Field() != field) {
      throw std::invalid_argument("ViterbiDecoder::DecodeHard: field mismatch");
    }
    if (received.Dimension() != 1 || received.Size() != conv_code_->Length()) {
      throw std::invalid_argument(
          "ViterbiDecoder::DecodeHard: received length mismatch");
    }
    // Convert hard bits {0, 1} into unit-magnitude LLRs {+1.0, -1.0}.
    // Maximizing sum_b (1 - 2*y_b) * (1 - 2*r_b) is strictly equivalent to
    // minimizing Hamming distance sum_b (y_b ^ r_b).
    const auto span = received.RawSpan();
    std::vector<double> unit_llrs(span.size());
    for (std::size_t i = 0; i < span.size(); ++i) {
      unit_llrs[i] = ((span[i] & 1U) == 0U) ? 1.0 : -1.0;
    }
    return DecodeSoft(unit_llrs);
  }

  result_type DecodeSoft(const std::vector<double> &channel_llrs) const {
    const std::size_t n = conv_code_->Length();
    conv_detail::ValidateLlrVector(channel_llrs, n,
                                   "ViterbiDecoder::DecodeSoft");

    std::vector<double> llrs(n);
    std::vector<uint8_t> hard_rx(n);
    for (std::size_t i = 0; i < n; ++i) {
      llrs[i] = conv_detail::ClampLlr(channel_llrs[i]);
      hard_rx[i] = (channel_llrs[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
    }

    const Trellis &trellis = conv_code_->GetTrellis();
    const TerminationMode term = conv_code_->GetTerminationMode();

    std::vector<uint8_t> decoded_msg;
    if (term == TerminationMode::kTailBiting && trellis.num_states > 1) {
      decoded_msg = DecodeTailBiting(llrs);
    } else {
      decoded_msg =
          DecodeFromZeroInit(llrs, term == TerminationMode::kZeroTailed);
    }

    std::vector<uint8_t> decoded_cw = conv_code_->EncodeBits(decoded_msg);
    std::size_t num_errors = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if ((decoded_cw[i] & 1U) != (hard_rx[i] & 1U)) {
        ++num_errors;
      }
    }

    auto field = this->code_->Field();
    return result_type::Success(
        conv_detail::BuildBinaryArray(field, decoded_cw),
        conv_detail::BuildBinaryArray(field, decoded_msg), num_errors);
  }

  result_type DecodeWithResult(
      const codeword_type &received_word) const override {
    return DecodeHard(received_word);
  }

  std::string Name() const override { return "ViterbiHard"; }
  std::size_t DecodingRadius() const override { return 0; }

 private:
  static double BranchCorrelation(uint32_t out_mask, const double *step_llrs,
                                  std::size_t n0) noexcept {
    double metric = 0.0;
    for (std::size_t j = 0; j < n0; ++j) {
      const double sign = (((out_mask >> j) & 1U) == 0U) ? 1.0 : -1.0;
      metric += sign * step_llrs[j];
    }
    return metric;
  }

  std::vector<uint8_t> DecodeFromZeroInit(const std::vector<double> &llrs,
                                          bool zero_tailed) const {
    const Trellis &trellis = conv_code_->GetTrellis();
    const std::size_t num_states = trellis.num_states;
    const std::size_t n0 = trellis.num_output_bits;
    const std::size_t k = conv_code_->InputBlockLength();
    const std::size_t steps = conv_code_->NumSteps();

    std::vector<double> pm(num_states, conv_detail::kNegInf);
    std::vector<double> next_pm(num_states, conv_detail::kNegInf);
    pm[0] = 0.0;

    std::vector<uint32_t> prev_state(steps * num_states, 0U);
    std::vector<uint8_t> prev_input(steps * num_states, 0U);

    for (std::size_t t = 0; t < steps; ++t) {
      std::fill(next_pm.begin(), next_pm.end(), conv_detail::kNegInf);
      const double *step_llrs = llrs.data() + t * n0;
      const bool is_tail = zero_tailed && (t >= k);

      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        if (pm[s] == conv_detail::kNegInf) {
          continue;
        }
        const uint32_t u_begin = is_tail ? trellis.tail_input[s] : 0U;
        const uint32_t u_end = is_tail ? (u_begin + 1U) : 2U;
        for (uint32_t u = u_begin; u < u_end; ++u) {
          const uint32_t ns = trellis.NextState(s, u);
          const uint32_t mask = trellis.OutputMask(s, u);
          const double cand = pm[s] + BranchCorrelation(mask, step_llrs, n0);
          if (cand > next_pm[ns]) {
            next_pm[ns] = cand;
            prev_state[t * num_states + ns] = s;
            prev_input[t * num_states + ns] = static_cast<uint8_t>(u);
          }
        }
      }
      pm.swap(next_pm);
    }

    uint32_t best_state = 0;
    if (!zero_tailed) {
      double best_metric = conv_detail::kNegInf;
      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        if (pm[s] > best_metric) {
          best_metric = pm[s];
          best_state = s;
        }
      }
    }

    std::vector<uint8_t> msg(k, 0U);
    uint32_t cur = best_state;
    for (std::size_t t = steps; t-- > 0;) {
      const uint8_t u = prev_input[t * num_states + cur];
      cur = prev_state[t * num_states + cur];
      if (t < k) {
        msg[t] = u;
      }
    }
    return msg;
  }

  // Wrap-Around Viterbi Algorithm (WAVA) with exact ML circular fallback.
  std::vector<uint8_t> DecodeTailBiting(const std::vector<double> &llrs) const {
    const Trellis &trellis = conv_code_->GetTrellis();
    const std::size_t num_states = trellis.num_states;
    const std::size_t n0 = trellis.num_output_bits;
    const std::size_t k = conv_code_->InputBlockLength();

    // Fast path: WAVA (up to 4 wrap-around iterations).
    std::vector<double> pm(num_states, 0.0);
    std::vector<double> next_pm(num_states, conv_detail::kNegInf);
    std::vector<uint32_t> prev_state(k * num_states, 0U);
    std::vector<uint8_t> prev_input(k * num_states, 0U);

    constexpr std::size_t kWavaIterations = 4;
    for (std::size_t iter = 0; iter < kWavaIterations; ++iter) {
      for (std::size_t t = 0; t < k; ++t) {
        std::fill(next_pm.begin(), next_pm.end(), conv_detail::kNegInf);
        const double *step_llrs = llrs.data() + t * n0;
        for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
          if (pm[s] == conv_detail::kNegInf) {
            continue;
          }
          for (uint32_t u = 0; u < 2U; ++u) {
            const uint32_t ns = trellis.NextState(s, u);
            const uint32_t mask = trellis.OutputMask(s, u);
            const double cand = pm[s] + BranchCorrelation(mask, step_llrs, n0);
            if (cand > next_pm[ns]) {
              next_pm[ns] = cand;
              prev_state[t * num_states + ns] = s;
              prev_input[t * num_states + ns] = static_cast<uint8_t>(u);
            }
          }
        }
        pm.swap(next_pm);
      }

      // Check whether the globally best path is tail-biting (start == end).
      uint32_t best_end = 0;
      double best_metric = conv_detail::kNegInf;
      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        if (pm[s] > best_metric) {
          best_metric = pm[s];
          best_end = s;
        }
      }

      std::vector<uint8_t> cand_msg(k, 0U);
      uint32_t cur = best_end;
      for (std::size_t t = k; t-- > 0;) {
        cand_msg[t] = prev_input[t * num_states + cur];
        cur = prev_state[t * num_states + cur];
      }
      if (cur == best_end) {
        return cand_msg;
      }

      // Normalize metrics before next wrap-around pass.
      for (double &m : pm) {
        if (m != conv_detail::kNegInf) {
          m -= best_metric;
        }
      }
    }

    // Exact ML fallback: evaluate each candidate starting state s0 and require
    // end state == s0 (guaranteed optimal for all circular paths).
    double global_best_metric = conv_detail::kNegInf;
    std::vector<uint8_t> global_best_msg(k, 0U);
    for (uint32_t s0 = 0; s0 < static_cast<uint32_t>(num_states); ++s0) {
      std::fill(pm.begin(), pm.end(), conv_detail::kNegInf);
      pm[s0] = 0.0;
      for (std::size_t t = 0; t < k; ++t) {
        std::fill(next_pm.begin(), next_pm.end(), conv_detail::kNegInf);
        const double *step_llrs = llrs.data() + t * n0;
        for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
          if (pm[s] == conv_detail::kNegInf) {
            continue;
          }
          for (uint32_t u = 0; u < 2U; ++u) {
            const uint32_t ns = trellis.NextState(s, u);
            const uint32_t mask = trellis.OutputMask(s, u);
            const double cand = pm[s] + BranchCorrelation(mask, step_llrs, n0);
            if (cand > next_pm[ns]) {
              next_pm[ns] = cand;
              prev_state[t * num_states + ns] = s;
              prev_input[t * num_states + ns] = static_cast<uint8_t>(u);
            }
          }
        }
        pm.swap(next_pm);
      }

      if (pm[s0] > global_best_metric) {
        global_best_metric = pm[s0];
        uint32_t cur = s0;
        for (std::size_t t = k; t-- > 0;) {
          global_best_msg[t] = prev_input[t * num_states + cur];
          cur = prev_state[t * num_states + cur];
        }
      }
    }
    return global_best_msg;
  }

  const ConvolutionalCode *conv_code_;
};

// Soft-Input Soft-Output (SISO) BCJR (Bahl-Cocke-Jelinek-Raviv) MAP decoder in
// log-domain supporting exact Log-MAP and Max-Log-MAP.
class BcjrDecoder {
 public:
  explicit BcjrDecoder(const ConvolutionalCode &code) : conv_code_(&code) {}

  BcjrOutput DecodeSoft(const std::vector<double> &channel_llrs,
                        const std::vector<double> &a_priori_llrs = {},
                        BcjrMode mode = BcjrMode::kLogMap) const {
    const std::size_t n = conv_code_->Length();
    const std::size_t k = conv_code_->InputBlockLength();
    const std::size_t steps = conv_code_->NumSteps();
    const Trellis &trellis = conv_code_->GetTrellis();
    const std::size_t num_states = trellis.num_states;
    const std::size_t n0 = trellis.num_output_bits;
    const TerminationMode term = conv_code_->GetTerminationMode();

    conv_detail::ValidateLlrVector(channel_llrs, n, "BcjrDecoder::DecodeSoft");
    if (!a_priori_llrs.empty()) {
      conv_detail::ValidateLlrVector(a_priori_llrs, k,
                                     "BcjrDecoder::DecodeSoft (a_priori)");
    }

    std::vector<double> llrs(n);
    for (std::size_t i = 0; i < n; ++i) {
      llrs[i] = conv_detail::ClampLlr(channel_llrs[i]);
    }
    std::vector<double> prior(k, 0.0);
    if (!a_priori_llrs.empty()) {
      for (std::size_t i = 0; i < k; ++i) {
        prior[i] = conv_detail::ClampLlr(a_priori_llrs[i]);
      }
    }

    // Precompute branch metrics gamma[t][s * 2 + u] = 0.5*(1-2u)*L_a(u_t) +
    // 0.5*sum_j (1-2*c_{t,j})*L_c(r_{t,j}).
    const std::size_t stride = num_states * 2;
    std::vector<double> gamma(steps * stride, conv_detail::kNegInf);
    for (std::size_t t = 0; t < steps; ++t) {
      const double *step_llrs = llrs.data() + t * n0;
      const bool is_tail = (term == TerminationMode::kZeroTailed) && (t >= k);
      const double half_la = is_tail ? 0.0 : (0.5 * prior[t]);

      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        const uint32_t u_begin = is_tail ? trellis.tail_input[s] : 0U;
        const uint32_t u_end = is_tail ? (u_begin + 1U) : 2U;
        for (uint32_t u = u_begin; u < u_end; ++u) {
          const uint32_t mask = trellis.OutputMask(s, u);
          double g = (u == 0U) ? half_la : -half_la;
          for (std::size_t j = 0; j < n0; ++j) {
            const double sign = (((mask >> j) & 1U) == 0U) ? 0.5 : -0.5;
            g += sign * step_llrs[j];
          }
          gamma[t * stride + static_cast<std::size_t>(s) * 2 + u] = g;
        }
      }
    }

    // Forward recursion alpha[t * num_states + s] with per-step normalization.
    std::vector<double> alpha((steps + 1) * num_states, conv_detail::kNegInf);
    if (term == TerminationMode::kTailBiting) {
      // Uniform initialization + 1 wrap-around warm-up pass for circular state
      // distribution.
      std::vector<double> cur_a(num_states, 0.0);
      std::vector<double> nxt_a(num_states, conv_detail::kNegInf);
      for (std::size_t t = 0; t < steps; ++t) {
        std::fill(nxt_a.begin(), nxt_a.end(), conv_detail::kNegInf);
        double max_val = conv_detail::kNegInf;
        for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
          const double a_s = cur_a[s];
          if (a_s == conv_detail::kNegInf) {
            continue;
          }
          for (uint32_t u = 0; u < 2U; ++u) {
            const double g =
                gamma[t * stride + static_cast<std::size_t>(s) * 2 + u];
            if (g == conv_detail::kNegInf) {
              continue;
            }
            const uint32_t ns = trellis.NextState(s, u);
            nxt_a[ns] = conv_detail::MaxStar(nxt_a[ns], a_s + g, mode);
            if (nxt_a[ns] > max_val) {
              max_val = nxt_a[ns];
            }
          }
        }
        if (max_val != conv_detail::kNegInf) {
          for (double &v : nxt_a) {
            if (v != conv_detail::kNegInf) {
              v -= max_val;
            }
          }
        }
        cur_a.swap(nxt_a);
      }
      std::copy(cur_a.begin(), cur_a.end(), alpha.begin());
    } else {
      alpha[0] = 0.0;
    }

    for (std::size_t t = 0; t < steps; ++t) {
      const double *cur = alpha.data() + t * num_states;
      double *nxt = alpha.data() + (t + 1) * num_states;
      double max_val = conv_detail::kNegInf;
      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        const double a_s = cur[s];
        if (a_s == conv_detail::kNegInf) {
          continue;
        }
        for (uint32_t u = 0; u < 2U; ++u) {
          const double g =
              gamma[t * stride + static_cast<std::size_t>(s) * 2 + u];
          if (g == conv_detail::kNegInf) {
            continue;
          }
          const uint32_t ns = trellis.NextState(s, u);
          nxt[ns] = conv_detail::MaxStar(nxt[ns], a_s + g, mode);
          if (nxt[ns] > max_val) {
            max_val = nxt[ns];
          }
        }
      }
      if (max_val != conv_detail::kNegInf) {
        for (std::size_t s = 0; s < num_states; ++s) {
          if (nxt[s] != conv_detail::kNegInf) {
            nxt[s] -= max_val;
          }
        }
      }
    }

    // Backward recursion beta[t * num_states + s] with per-step normalization.
    std::vector<double> beta((steps + 1) * num_states, conv_detail::kNegInf);
    if (term == TerminationMode::kZeroTailed) {
      beta[steps * num_states + 0] = 0.0;
    } else if (term == TerminationMode::kTruncated) {
      std::fill(beta.begin() + static_cast<std::ptrdiff_t>(steps * num_states),
                beta.end(), 0.0);
    } else {
      // Tail-biting: 1 wrap-around warm-up pass for circular backward metrics.
      std::vector<double> cur_b(num_states, 0.0);
      std::vector<double> prv_b(num_states, conv_detail::kNegInf);
      for (std::size_t t = steps; t-- > 0;) {
        std::fill(prv_b.begin(), prv_b.end(), conv_detail::kNegInf);
        double max_val = conv_detail::kNegInf;
        for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
          double acc = conv_detail::kNegInf;
          for (uint32_t u = 0; u < 2U; ++u) {
            const double g =
                gamma[t * stride + static_cast<std::size_t>(s) * 2 + u];
            if (g == conv_detail::kNegInf) {
              continue;
            }
            const uint32_t ns = trellis.NextState(s, u);
            if (cur_b[ns] != conv_detail::kNegInf) {
              acc = conv_detail::MaxStar(acc, g + cur_b[ns], mode);
            }
          }
          prv_b[s] = acc;
          if (acc > max_val) {
            max_val = acc;
          }
        }
        if (max_val != conv_detail::kNegInf) {
          for (double &v : prv_b) {
            if (v != conv_detail::kNegInf) {
              v -= max_val;
            }
          }
        }
        cur_b.swap(prv_b);
      }
      std::copy(cur_b.begin(), cur_b.end(),
                beta.begin() + static_cast<std::ptrdiff_t>(steps * num_states));
    }

    for (std::size_t t = steps; t-- > 0;) {
      const double *nxt = beta.data() + (t + 1) * num_states;
      double *cur = beta.data() + t * num_states;
      double max_val = conv_detail::kNegInf;
      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        double acc = conv_detail::kNegInf;
        for (uint32_t u = 0; u < 2U; ++u) {
          const double g =
              gamma[t * stride + static_cast<std::size_t>(s) * 2 + u];
          if (g == conv_detail::kNegInf) {
            continue;
          }
          const uint32_t ns = trellis.NextState(s, u);
          if (nxt[ns] != conv_detail::kNegInf) {
            acc = conv_detail::MaxStar(acc, g + nxt[ns], mode);
          }
        }
        cur[s] = acc;
        if (acc > max_val) {
          max_val = acc;
        }
      }
      if (max_val != conv_detail::kNegInf) {
        for (std::size_t s = 0; s < num_states; ++s) {
          if (cur[s] != conv_detail::kNegInf) {
            cur[s] -= max_val;
          }
        }
      }
    }

    // Compute APP LLRs, extrinsic LLRs, and hard bit decisions for t = 0..K-1.
    BcjrOutput output;
    output.app_llrs.resize(k, 0.0);
    output.extrinsic_llrs.resize(k, 0.0);
    std::vector<uint8_t> hard(k, 0U);

    for (std::size_t t = 0; t < k; ++t) {
      const double *a_t = alpha.data() + t * num_states;
      const double *b_next = beta.data() + (t + 1) * num_states;
      double acc0 = conv_detail::kNegInf;
      double acc1 = conv_detail::kNegInf;

      for (uint32_t s = 0; s < static_cast<uint32_t>(num_states); ++s) {
        const double a_s = a_t[s];
        if (a_s == conv_detail::kNegInf) {
          continue;
        }
        const uint32_t ns0 = trellis.NextState(s, 0U);
        const double g0 =
            gamma[t * stride + static_cast<std::size_t>(s) * 2 + 0U];
        if (g0 != conv_detail::kNegInf && b_next[ns0] != conv_detail::kNegInf) {
          acc0 = conv_detail::MaxStar(acc0, a_s + g0 + b_next[ns0], mode);
        }

        const uint32_t ns1 = trellis.NextState(s, 1U);
        const double g1 =
            gamma[t * stride + static_cast<std::size_t>(s) * 2 + 1U];
        if (g1 != conv_detail::kNegInf && b_next[ns1] != conv_detail::kNegInf) {
          acc1 = conv_detail::MaxStar(acc1, a_s + g1 + b_next[ns1], mode);
        }
      }

      const double app = acc0 - acc1;
      output.app_llrs[t] = app;
      const double sys_llr = trellis.is_recursive ? llrs[t * n0] : 0.0;
      output.extrinsic_llrs[t] = app - prior[t] - sys_llr;
      hard[t] = (app >= 0.0) ? uint8_t{0} : uint8_t{1};
    }

    output.hard_bits = conv_detail::BuildBinaryArray(conv_code_->Field(), hard);
    return output;
  }

 private:
  const ConvolutionalCode *conv_code_;
};

inline ConvolutionalCode::message_type ConvolutionalCode::UnencodeFast(
    const codeword_type &codeword) const {
  if (!codeword.HasField() || codeword.Field() != this->field_) {
    throw std::invalid_argument(
        "ConvolutionalCode::UnencodeFast: field mismatch");
  }
  if (codeword.Dimension() != 1 || codeword.Size() != this->length_) {
    throw std::invalid_argument(
        "ConvolutionalCode::UnencodeFast: codeword length mismatch");
  }
  const auto span = codeword.RawSpan();
  std::vector<uint8_t> msg_bits(input_block_len_, 0U);
  if (trellis_.is_recursive) {
    const std::size_t n0 = trellis_.num_output_bits;
    for (std::size_t t = 0; t < input_block_len_; ++t) {
      msg_bits[t] = static_cast<uint8_t>(span[t * n0] & 1U);
    }
    const std::vector<uint8_t> reencoded = EncodeBits(msg_bits);
    for (std::size_t i = 0; i < this->length_; ++i) {
      if ((reencoded[i] & 1U) != (span[i] & 1U)) {
        throw std::invalid_argument(
            "ConvolutionalCode::UnencodeFast: input is not a valid codeword");
      }
    }
    return conv_detail::BuildBinaryArray(this->field_, msg_bits);
  }

  ViterbiDecoder viterbi(*this);
  auto res = viterbi.DecodeHard(codeword);
  if (!res.Ok() || res.num_errors != 0) {
    throw std::invalid_argument(
        "ConvolutionalCode::UnencodeFast: input is not a valid codeword");
  }
  return std::move(res.message);
}

inline void ConvolutionalCode::RegisterEncodersAndDecoders() {
  this->RegisterEncoder(
      "Convolutional", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &conv = dynamic_cast<const ConvolutionalCode &>(code);
        return std::make_shared<ConvolutionalEncoder>(conv);
      });
  this->RegisterEncoder(
      "GeneratorMatrix", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &lin =
            dynamic_cast<const AbstractLinearCode<GaloisFieldBinary> &>(code);
        return std::make_shared<GeneratorMatrixEncoder<GaloisFieldBinary>>(lin);
      });
  this->RegisterDecoder(
      "ViterbiHard", [](const AbstractCode<GaloisFieldBinary> &code) {
        const auto &conv = dynamic_cast<const ConvolutionalCode &>(code);
        return std::make_shared<ViterbiDecoder>(conv);
      });
}

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_CONVOLUTIONAL_HPP_
