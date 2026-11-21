#ifndef XGALOIS_CODING_ABSTRACT_LINEAR_CODE_HPP_
#define XGALOIS_CODING_ABSTRACT_LINEAR_CODE_HPP_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/coding/abstract_code.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/utils/math.hpp"

namespace xg {
namespace coding {

template <typename GaloisField>
class AbstractLinearCode : public AbstractCode<GaloisField> {
 public:
  using Base = AbstractCode<GaloisField>;
  using field_type = GaloisField;
  using element_type = typename Base::element_type;
  using storage_type = typename Base::storage_type;
  using codeword_type = typename Base::codeword_type;
  using message_type = typename Base::message_type;
  using matrix_type = GFMatrix<GaloisField>;
  using syndrome_type = GFArray<GaloisField>;

  explicit AbstractLinearCode(std::shared_ptr<GaloisField> field,
                              std::size_t length, std::size_t dimension,
                              Metric metric = Metric::HAMMING,
                              std::string default_encoder = "GeneratorMatrix",
                              std::string default_decoder = "Syndrome")
      : Base(std::move(field), length, metric, std::move(default_encoder),
             std::move(default_decoder)),
        dimension_(dimension) {
    if (dimension_ > length) {
      throw std::invalid_argument(
          "AbstractLinearCode: dimension cannot exceed code length");
    }
  }

  ~AbstractLinearCode() override = default;

  std::size_t Dimension() const override { return dimension_; }

  virtual const matrix_type &GeneratorMatrix() const = 0;
  virtual const matrix_type &ParityCheckMatrix() const = 0;

  // Compute exact parity-check matrix H (shape {(n-k), n}) satisfying G*H^T = 0
  // from any full-row-rank (or rank-deficient) generator matrix G using
  // linalg::Nullspace.
  static matrix_type ComputeParityCheckFromGenerator(
      const matrix_type &generator) {
    if (!generator.HasField() || generator.Dimension() != 2) {
      throw std::invalid_argument(
          "ComputeParityCheckFromGenerator: generator must be a 2D matrix");
    }
    return linalg::Nullspace(generator);
  }

  // Compute exact generator matrix G (shape {k, n}) satisfying G*H^T = 0
  // from any parity-check matrix H using linalg::Nullspace.
  static matrix_type ComputeGeneratorFromParityCheck(
      const matrix_type &parity_check) {
    if (!parity_check.HasField() || parity_check.Dimension() != 2) {
      throw std::invalid_argument(
          "ComputeGeneratorFromParityCheck: parity_check must be a 2D matrix");
    }
    return linalg::Nullspace(parity_check);
  }

  virtual matrix_type SystematicGeneratorMatrix() const {
    return linalg::SystematicForm(GeneratorMatrix()).g_sys;
  }

  bool IsSystematic() const {
    const auto &G = GeneratorMatrix();
    if (G.Shape(0) != dimension_ || G.Shape(1) != this->length_) {
      return false;
    }
    for (std::size_t i = 0; i < dimension_; ++i) {
      for (std::size_t j = 0; j < dimension_; ++j) {
        const auto &val = G.Raw(i, j);
        if (i == j) {
          if (!this->field_->IsOne(val)) return false;
        } else {
          if (!this->field_->IsZero(val)) return false;
        }
      }
    }
    return true;
  }

  syndrome_type Syndrome(const codeword_type &word) const {
    if (!word.HasField() || word.Field() != this->field_) {
      throw std::invalid_argument("Syndrome: field mismatch");
    }
    if (word.Dimension() != 1 || word.Size() != this->length_) {
      throw std::invalid_argument(
          "Syndrome: word must be 1D with size equal to code length");
    }
    const auto &H = ParityCheckMatrix();
    if (H.Shape(0) == 0) {
      return GFArray<GaloisField>::Zeros(this->field_, {0});
    }
    return linalg::Dot(H, word);
  }

  bool Contains(const codeword_type &word) const override {
    if (!word.HasField() || word.Field() != this->field_ ||
        word.Dimension() != 1 || word.Size() != this->length_) {
      return false;
    }
    const auto syn = Syndrome(word);
    for (const auto &v : syn.RawSpan()) {
      if (!this->field_->IsZero(v)) {
        return false;
      }
    }
    return true;
  }

  codeword_type RandomCodeword() const override {
    if (dimension_ == 0) {
      return codeword_type::Zeros(this->field_, {this->length_});
    }
    auto msg = message_type::Random(this->field_, {dimension_});
    return linalg::Dot(msg, GeneratorMatrix());
  }

  codeword_type RandomCodeword(std::mt19937_64 &rng) const override {
    if (dimension_ == 0) {
      return codeword_type::Zeros(this->field_, {this->length_});
    }
    auto msg = message_type::Random(this->field_, {dimension_}, rng);
    return linalg::Dot(msg, GeneratorMatrix());
  }

  codeword_type LinearCombination(const codeword_type &c1,
                                  const element_type &a1,
                                  const codeword_type &c2,
                                  const element_type &a2) const {
    return (c1 * a1) + (c2 * a2);
  }

  // Enumerate all codewords (for small codes with q^k <= max_codewords).
  std::vector<codeword_type> GetCodewords(
      std::size_t max_codewords = 65536) const {
    const uint64_t q = this->field_->Order();
    if (q == 0) {
      throw std::runtime_error(
          "GetCodewords: field order too large for enumeration");
    }
    uint64_t total = 1;
    try {
      total = utils::SafeIntegerPower(q, dimension_);
    } catch (const std::overflow_error &) {
      throw std::runtime_error("GetCodewords: code too large for enumeration");
    }
    if (total > max_codewords) {
      throw std::runtime_error("GetCodewords: code size (" +
                               std::to_string(total) + ") exceeds limit (" +
                               std::to_string(max_codewords) + ")");
    }

    const auto &G = GeneratorMatrix();
    std::vector<storage_type> alphabet = EnumerateFieldElements();
    std::vector<codeword_type> codewords;
    codewords.reserve(static_cast<std::size_t>(total));

    std::vector<std::size_t> digits(dimension_, 0);
    using RawArray = typename message_type::raw_array_type;
    using ShapeType = typename message_type::shape_type;
    auto msg_raw = RawArray::from_shape(ShapeType{dimension_});

    for (uint64_t idx = 0; idx < total; ++idx) {
      for (std::size_t j = 0; j < dimension_; ++j) {
        msg_raw(j) = alphabet[digits[j]];
      }
      auto msg = message_type::FromRaw(this->field_, msg_raw);
      if (dimension_ == 0) {
        codewords.push_back(
            codeword_type::Zeros(this->field_, {this->length_}));
      } else {
        codewords.push_back(linalg::Dot(msg, G));
      }
      // Increment base-q counter.
      for (std::size_t j = 0; j < dimension_; ++j) {
        if (++digits[j] < alphabet.size()) {
          break;
        }
        digits[j] = 0;
      }
    }
    return codewords;
  }

  // Exact minimum distance by exhaustive search over all non-zero codewords
  // (cached; subclasses like GRS/Hamming/Golay override with O(1) formulas).
  std::size_t MinimumDistance() const override {
    if (cached_min_distance_.has_value()) {
      return *cached_min_distance_;
    }
    if (dimension_ == 0) {
      cached_min_distance_ = this->length_;
      return this->length_;
    }
    const auto codewords = GetCodewords();
    std::size_t min_w = std::numeric_limits<std::size_t>::max();
    for (const auto &cw : codewords) {
      const std::size_t w = this->Weight(cw);
      if (w > 0 && w < min_w) {
        min_w = w;
        if (min_w == 1) {
          break;
        }
      }
    }
    cached_min_distance_ = (min_w == std::numeric_limits<std::size_t>::max())
                               ? this->length_
                               : min_w;
    return *cached_min_distance_;
  }

  std::string ToString() const override {
    return "[" + std::to_string(this->length_) + ", " +
           std::to_string(dimension_) + "] linear code over GF(" +
           std::to_string(this->field_->Order()) + ")";
  }

 protected:
  // Enumerate all field elements in internal representation (works for GF(2),
  // GF(p), GF(2^m) including Zech, and GF(p^m) ExtensionElement!).
  std::vector<storage_type> EnumerateFieldElements() const {
    const uint64_t q = this->field_->Order();
    std::vector<storage_type> elems;
    elems.reserve(static_cast<std::size_t>(q));
    if constexpr (std::is_integral_v<storage_type>) {
      for (uint64_t i = 0; i < q; ++i) {
        elems.push_back(
            this->field_->SetElementValue(static_cast<storage_type>(i)));
      }
    } else {
      for (uint64_t i = 0; i < q; ++i) {
        elems.push_back(this->field_->FromInt(i));
      }
    }
    return elems;
  }

  std::size_t dimension_;
  mutable std::optional<std::size_t> cached_min_distance_;
};

}  // namespace coding
}  // namespace xg

#endif  // XGALOIS_CODING_ABSTRACT_LINEAR_CODE_HPP_
