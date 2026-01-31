#ifndef XGALOIS_CODING_GRS_HPP
#define XGALOIS_CODING_GRS_HPP

#include "xgalois/coding/abstract_linear_code.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include <xtensor/containers/xarray.hpp>
#include <xtensor/views/xview.hpp>
#include <xtensor/containers/xadapt.hpp>
#include <memory>
#include <vector>

namespace xg {
namespace coding {

template <typename ElementType>
class GeneralizedReedSolomonCode : public AbstractLinearCode<ElementType> {
public:
    using element_type = ElementType;
    using codeword_type = xt::xarray<ElementType>;
    using message_type = xt::xarray<ElementType>;
    using matrix_type = xt::xarray<ElementType>;
    using vector_type = xt::xarray<ElementType>;
    using polynomial_type = xg::PolynomialDense<ElementType>;

    // Constructor for GRS code
    GeneralizedReedSolomonCode(std::shared_ptr<GaloisFieldBase<ElementType>> field,
                              const xt::xarray<ElementType>& evaluation_points,
                              const xt::xarray<ElementType>& column_multipliers,
                              size_t dimension)
        : AbstractLinearCode<ElementType>(evaluation_points.size(),
                                         dimension,
                                         "GRSEncoder",
                                         "GRSDecoder",
                                         Metric::HAMMING),
          field_(field),
          evaluation_points_(evaluation_points),
          column_multipliers_(column_multipliers) {

        if (evaluation_points.size() != column_multipliers.size()) {
            throw std::invalid_argument("Evaluation points and column multipliers must have same size");
        }

        if (dimension > evaluation_points.size()) {
            throw std::invalid_argument("Dimension cannot exceed number of evaluation points");
        }

        // Check that evaluation points are distinct
        for (size_t i = 0; i < evaluation_points.size(); ++i) {
            for (size_t j = i + 1; j < evaluation_points.size(); ++j) {
                if (evaluation_points(i) == evaluation_points(j)) {
                    throw std::invalid_argument("Evaluation points must be distinct");
                }
            }
        }

        // Check that column multipliers are non-zero
        for (size_t i = 0; i < column_multipliers.size(); ++i) {
            if (column_multipliers(i) == ElementType{}) {
                throw std::invalid_argument("Column multipliers must be non-zero");
            }
        }

        // Build generator and parity check matrices
        BuildMatrices();

        // Register encoders and decoders
        RegisterEncodersAndDecoders();
    }

    // Standard Reed-Solomon constructor (column multipliers all 1)
    static std::unique_ptr<GeneralizedReedSolomonCode<ElementType>>
    StandardReedSolomon(std::shared_ptr<GaloisFieldBase<ElementType>> field,
                       const xt::xarray<ElementType>& evaluation_points,
                       size_t dimension) {
        xt::xarray<ElementType> multipliers = xt::ones<ElementType>({evaluation_points.size()});
        return std::make_unique<GeneralizedReedSolomonCode>(field, evaluation_points,
                                                           multipliers, dimension);
    }

    // Reed-Solomon constructor using consecutive powers
    static std::unique_ptr<GeneralizedReedSolomonCode<ElementType>>
    ReedSolomon(std::shared_ptr<GaloisFieldBase<ElementType>> field,
                size_t length,
                size_t dimension,
                const ElementType& primitive_element) {
        xt::xarray<ElementType> evaluation_points = xt::zeros<ElementType>({length});
        xt::xarray<ElementType> multipliers = xt::ones<ElementType>({length});

        // Use consecutive powers of primitive element
        ElementType power = ElementType(1);
        for (size_t i = 0; i < length; ++i) {
            evaluation_points(i) = power;
            power = field->Mul(power, primitive_element);
        }

        return std::make_unique<GeneralizedReedSolomonCode>(field, evaluation_points,
                                                           multipliers, dimension);
    }

    // Override methods from AbstractLinearCode
    std::shared_ptr<GaloisFieldBase<ElementType>> Field() const override {
        return field_;
    }

    matrix_type GeneratorMatrix() const override {
        return generator_matrix_;
    }

    matrix_type ParityCheckMatrix() const override {
        return parity_check_matrix_;
    }

    size_t MinimumDistance() const override {
        // For GRS codes, minimum distance is n - k + 1 (MDS property)
        return this->length_ - this->dimension_ + 1;
    }

    std::unique_ptr<AbstractLinearCode<ElementType>> DualCode() const override {
        // Dual of GRS code is also GRS
        return ComputeDualCode();
    }

    std::string ToString() const override {
        return "[" + std::to_string(this->length_) + ", " +
               std::to_string(this->dimension_) + ", " +
               std::to_string(MinimumDistance()) + "] GRS code over GF(" +
               std::to_string(field_->Order()) + ")";
    }

    // GRS-specific methods
    const xt::xarray<ElementType>& EvaluationPoints() const {
        return evaluation_points_;
    }

    const xt::xarray<ElementType>& ColumnMultipliers() const {
        return column_multipliers_;
    }

    // Encode polynomial directly
    codeword_type EncodePolynomial(const polynomial_type& message_poly) const {
        if (message_poly.Degree() >= this->dimension_) {
            throw std::invalid_argument("Message polynomial degree must be less than dimension");
        }

        codeword_type codeword = xt::zeros<ElementType>({this->length_});
        for (size_t i = 0; i < this->length_; ++i) {
            ElementType eval = message_poly.Evaluate(evaluation_points_(i));
            codeword(i) = field_->Mul(eval, column_multipliers_(i));
        }

        return codeword;
    }

    // Decode to polynomial (if possible)
    polynomial_type DecodeToPolynomial(const codeword_type& received_word) const {
        // This would use polynomial interpolation
        // For now, use the standard decoder and then interpolate
        auto decoder = this->GetDecoder();
        auto corrected = decoder->DecodeToCode(received_word);

        return InterpolateFromCodeword(corrected);
    }

    // Check if code is Maximum Distance Separable (MDS)
    bool IsMDS() const {
        return MinimumDistance() == this->length_ - this->dimension_ + 1;
    }

    // Systematic encoding positions
    xt::xarray<size_t> SystematicPositions() const {
        xt::xarray<size_t> positions = xt::zeros<size_t>({this->dimension_});
        for (size_t i = 0; i < this->dimension_; ++i) {
            positions(i) = i;
        }
        return positions;
    }

private:
    std::shared_ptr<GaloisFieldBase<ElementType>> field_;
    xt::xarray<ElementType> evaluation_points_;
    xt::xarray<ElementType> column_multipliers_;
    matrix_type generator_matrix_;
    matrix_type parity_check_matrix_;

    void BuildMatrices() {
        BuildGeneratorMatrix();
        BuildParityCheckMatrix();
    }

    void BuildGeneratorMatrix() {
        size_t k = this->dimension_;
        size_t n = this->length_;

        generator_matrix_ = matrix_type(k, n);

        // Generator matrix: G[i][j] = v_j * alpha_j^i
        // where v_j is column multiplier and alpha_j is evaluation point
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ElementType power = field_->Pow(evaluation_points_(j), i);
                ElementType entry = field_->Mul(column_multipliers_(j), power);
                generator_matrix_.Set(i, j, entry);
            }
        }
    }

    void BuildParityCheckMatrix() {
        size_t k = this->dimension_;
        size_t n = this->length_;
        size_t r = n - k;

        parity_check_matrix_ = matrix_type(r, n);

        // Parity check matrix for GRS dual code
        // H[i][j] = w_j * beta_j^i / v_j
        // where w_j are dual multipliers and beta_j are dual evaluation points

        // For standard construction, use same evaluation points
        // and compute dual multipliers
        auto dual_multipliers = ComputeDualMultipliers();

        for (size_t i = 0; i < r; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ElementType power = field_->Pow(evaluation_points_(j), i);
                ElementType entry = field_->Mul(dual_multipliers(j), power);
                parity_check_matrix_.Set(i, j, entry);
            }
        }
    }

    xt::xarray<ElementType> ComputeDualMultipliers() const {
        // Compute dual multipliers for GRS dual code
        xt::xarray<ElementType> dual_multipliers = xt::zeros<ElementType>({this->length_});

        for (size_t j = 0; j < this->length_; ++j) {
            // Compute product of (alpha_j - alpha_i) for all i != j
            ElementType product = ElementType(1);
            for (size_t i = 0; i < this->length_; ++i) {
                if (i != j) {
                    ElementType diff = field_->Sub(evaluation_points_(j), evaluation_points_(i));
                    product = field_->Mul(product, diff);
                }
            }

            // Dual multiplier is 1 / (v_j * product)
            ElementType denominator = field_->Mul(column_multipliers_(j), product);
            dual_multipliers(j) = field_->Inv(denominator);
        }

        return dual_multipliers;
    }

    std::unique_ptr<AbstractLinearCode<ElementType>> ComputeDualCode() const {
        auto dual_multipliers = ComputeDualMultipliers();
        size_t dual_dimension = this->length_ - this->dimension_;

        return std::make_unique<GeneralizedReedSolomonCode>(field_, evaluation_points_,
                                                           dual_multipliers, dual_dimension);
    }

    polynomial_type InterpolateFromCodeword(const codeword_type& codeword) const {
        // Lagrange interpolation to recover polynomial
        xt::xarray<ElementType> values = xt::zeros<ElementType>({this->length_});

        // Remove column multipliers
        for (size_t i = 0; i < this->length_; ++i) {
            values(i) = field_->Div(codeword(i), column_multipliers_(i));
        }

        // Use first k points for interpolation
        xt::xarray<ElementType> points = xt::zeros<ElementType>({this->dimension_});
        xt::xarray<ElementType> vals = xt::zeros<ElementType>({this->dimension_});

        for (size_t i = 0; i < this->dimension_; ++i) {
            points(i) = evaluation_points_(i);
            vals(i) = values(i);
        }

        return LagrangeInterpolation(points, vals);
    }

    polynomial_type LagrangeInterpolation(const xt::xarray<ElementType>& points,
                                         const xt::xarray<ElementType>& values) const {
        if (points.size() != values.size()) {
            throw std::invalid_argument("Points and values must have same size");
        }

        size_t n = points.size();
        xt::xarray<ElementType> result_coeffs = xt::zeros<ElementType>({n});

        for (size_t i = 0; i < n; ++i) {
            // Compute Lagrange basis polynomial L_i(x)
            xt::xarray<ElementType> basis_coeffs = xt::ones<ElementType>({1});

            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    // Multiply by (x - points[j]) / (points[i] - points[j])
                    ElementType denominator = field_->Sub(points(i), points(j));
                    ElementType inv_denom = field_->Inv(denominator);

                    // Multiply basis_coeffs by (x - points[j])
                    xt::xarray<ElementType> new_coeffs = xt::zeros<ElementType>({basis_coeffs.size() + 1});
                    for (size_t k = 0; k < basis_coeffs.size(); ++k) {
                        new_coeffs(k) = field_->Add(new_coeffs(k),
                                                   field_->Mul(basis_coeffs(k),
                                                              field_->Neg(points(j))));
                        new_coeffs(k + 1) = field_->Add(new_coeffs(k + 1), basis_coeffs(k));
                    }

                    // Multiply by 1 / (points[i] - points[j])
                    for (size_t k = 0; k < new_coeffs.size(); ++k) {
                        new_coeffs(k) = field_->Mul(new_coeffs(k), inv_denom);
                    }

                    basis_coeffs = new_coeffs;
                }
            }

            // Add values[i] * L_i(x) to result
            for (size_t k = 0; k < basis_coeffs.size() && k < result_coeffs.size(); ++k) {
                ElementType term = field_->Mul(values(i), basis_coeffs(k));
                result_coeffs(k) = field_->Add(result_coeffs(k), term);
            }
        }

        // Convert xtensor array to std::vector for polynomial constructor
        std::vector<ElementType> coeffs_vec(result_coeffs.size());
        for (size_t i = 0; i < result_coeffs.size(); ++i) {
            coeffs_vec[i] = result_coeffs(i);
        }

        return polynomial_type(coeffs_vec, field_);
    }

    void RegisterEncodersAndDecoders() {
        // Register GRS encoder
        this->RegisterEncoder("GRSEncoder",
            [](const AbstractCode<ElementType>* code) {
                return std::make_unique<GRSEncoder<ElementType>>(code);
            });

        // Register GRS decoder (Peterson-Gorenstein-Zierler)
        this->RegisterDecoder("GRSDecoder",
            [](const AbstractCode<ElementType>* code) {
                return std::make_unique<GRSDecoder<ElementType>>(code);
            });

        // Also register Berlekamp-Massey decoder
        this->RegisterDecoder("BerlekampMassey",
            [](const AbstractCode<ElementType>* code) {
                return std::make_unique<BerlekampMasseyDecoder<ElementType>>(code);
            });
    }
};

// Forward declarations for GRS-specific encoder/decoder
template <typename ElementType> class GRSEncoder;
template <typename ElementType> class GRSDecoder;
template <typename ElementType> class BerlekampMasseyDecoder;

} // namespace coding
} // namespace xg

#endif // XGALOIS_CODING_GRS_HPP
