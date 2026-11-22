#ifndef XGALOIS_POLY_SUBPRODUCT_TREE_HPP_
#define XGALOIS_POLY_SUBPRODUCT_TREE_HPP_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/poly/half_gcd.hpp"
#include "xgalois/poly/poly_dense.hpp"

namespace xg {
namespace poly {

// Balanced subproduct tree for O(M(n) log n) / O(n log^2 n) multipoint
// evaluation and Lagrange interpolation over any FiniteField F.
template <FiniteField F>
class SubproductTree {
 public:
  using field_type = F;
  using storage_type = typename F::element_type;
  using poly_type = PolynomialDense<F>;

  static constexpr std::size_t kEvalDirectCutoff = 16;

  SubproductTree(std::shared_ptr<F> field, std::vector<storage_type> points,
                 std::string var = "x")
      : field_(std::move(field)),
        points_(std::move(points)),
        var_(std::move(var)) {
    if (!field_) {
      throw std::invalid_argument(
          "SubproductTree: field pointer cannot be null");
    }
    if (var_.empty()) {
      throw std::invalid_argument(
          "SubproductTree: variable string cannot be empty");
    }
    if (!points_.empty()) {
      tree_.resize(4 * points_.size());
      BuildTree(1, 0, points_.size());
      PrecomputeWeights();
    }
  }

  std::size_t Size() const noexcept { return points_.size(); }

  bool Empty() const noexcept { return points_.empty(); }

  const std::shared_ptr<F>& Field() const noexcept { return field_; }

  const std::vector<storage_type>& Points() const noexcept { return points_; }

  const poly_type& RootProduct() const {
    if (points_.empty()) {
      throw std::logic_error(
          "SubproductTree::RootProduct: tree has no evaluation points");
    }
    return tree_[1];
  }

  // Evaluates f at all points [alpha_0, ..., alpha_{n-1}] in O(M(n) log n)
  // via a top-down remainder tree.
  std::vector<storage_type> MultipointEvaluate(const poly_type& f) const {
    if (!f.HasField() || f.Field() != field_) {
      throw std::invalid_argument(
          "SubproductTree::MultipointEvaluate: polynomial field mismatch");
    }
    const std::size_t n = points_.size();
    if (n == 0) {
      return {};
    }
    std::vector<storage_type> out(n, field_->AdditiveIdentity());
    if (f.IsZero()) {
      return out;
    }

    const poly_type root_rem =
        (f.Degree() >= tree_[1].Degree()) ? FastRem(f, tree_[1]) : f;
    EvalRec(1, 0, n, root_rem, out);
    return out;
  }

  // Interpolates the unique polynomial f of degree < n such that
  // f(points_[j]) == values[j] in O(M(n) log n) using cached barycentric
  // weights w_j = M'(alpha_j)^(-1) and a bottom-up linear combination tree.
  poly_type Interpolate(std::span<const storage_type> values) const {
    const std::size_t n = points_.size();
    if (values.size() != n) {
      throw std::invalid_argument(
          "SubproductTree::Interpolate: values length must match point count");
    }
    if (n == 0) {
      return poly_type::Zero(field_, var_);
    }
    if (has_duplicate_points_) {
      throw std::domain_error(
          "SubproductTree::Interpolate: evaluation points must be distinct");
    }

    std::vector<storage_type> scaled(n);
    field_->MulArray(std::span<storage_type>(scaled), values,
                     std::span<const storage_type>(weights_));
    return InterpRec(1, 0, n, std::span<const storage_type>(scaled));
  }

 private:
  static poly_type FastRem(const poly_type& a, const poly_type& b) {
    if (a.Degree() < b.Degree()) {
      return a;
    }
    auto [q, r] = detail::FastDivRem(a, b);
    (void)q;
    return std::move(r);
  }

  void BuildTree(std::size_t node, std::size_t l, std::size_t r) {
    if (r - l == 1) {
      std::vector<storage_type> raw = {field_->Neg(points_[l]),
                                       field_->MultiplicativeIdentity()};
      tree_[node] = poly_type::FromRaw(std::move(raw), field_, var_);
      return;
    }
    const std::size_t mid = l + (r - l) / 2;
    BuildTree(2 * node, l, mid);
    BuildTree(2 * node + 1, mid, r);
    tree_[node] = tree_[2 * node] * tree_[2 * node + 1];
  }

  void EvalRec(std::size_t node, std::size_t l, std::size_t r,
               const poly_type& rem, std::vector<storage_type>& out) const {
    const std::size_t len = r - l;
    if (len == 1) {
      out[l] = rem.RawCoeff(0);
      return;
    }
    if (len <= kEvalDirectCutoff) {
      const auto batch =
          rem.EvalBatch(std::span<const storage_type>(points_.data() + l, len));
      for (std::size_t i = 0; i < len; ++i) {
        out[l + i] = batch[i];
      }
      return;
    }
    const std::size_t mid = l + len / 2;
    const poly_type rem_left = FastRem(rem, tree_[2 * node]);
    const poly_type rem_right = FastRem(rem, tree_[2 * node + 1]);
    EvalRec(2 * node, l, mid, rem_left, out);
    EvalRec(2 * node + 1, mid, r, rem_right, out);
  }

  void PrecomputeWeights() {
    const poly_type deriv = tree_[1].Derivative();
    const std::vector<storage_type> deriv_evals = MultipointEvaluate(deriv);
    for (const auto& v : deriv_evals) {
      if (field_->IsZero(v)) {
        has_duplicate_points_ = true;
        return;
      }
    }
    weights_.resize(deriv_evals.size());
    if constexpr (requires(const F& f, std::span<storage_type> dst,
                           std::span<const storage_type> src) {
                    f.BatchInv(dst, src);
                  }) {
      field_->BatchInv(std::span<storage_type>(weights_),
                       std::span<const storage_type>(deriv_evals));
    } else {
      // Portable Montgomery batch inversion fallback.
      const std::size_t n = deriv_evals.size();
      std::vector<storage_type> prefix(n);
      prefix[0] = deriv_evals[0];
      for (std::size_t i = 1; i < n; ++i) {
        prefix[i] = field_->Mul(prefix[i - 1], deriv_evals[i]);
      }
      storage_type inv_acc = field_->Inv(prefix[n - 1]);
      for (std::size_t i = n - 1; i > 0; --i) {
        weights_[i] = field_->Mul(inv_acc, prefix[i - 1]);
        inv_acc = field_->Mul(inv_acc, deriv_evals[i]);
      }
      weights_[0] = inv_acc;
    }
  }

  poly_type InterpRec(std::size_t node, std::size_t l, std::size_t r,
                      std::span<const storage_type> scaled) const {
    if (r - l == 1) {
      return poly_type::FromRaw({scaled[l]}, field_, var_);
    }
    const std::size_t mid = l + (r - l) / 2;
    const poly_type left_poly = InterpRec(2 * node, l, mid, scaled);
    const poly_type right_poly = InterpRec(2 * node + 1, mid, r, scaled);
    return left_poly * tree_[2 * node + 1] + right_poly * tree_[2 * node];
  }

  std::shared_ptr<F> field_;
  std::vector<storage_type> points_;
  std::string var_ = "x";
  std::vector<poly_type> tree_;
  std::vector<storage_type> weights_;
  bool has_duplicate_points_ = false;
};

}  // namespace poly

using poly::SubproductTree;

}  // namespace xg

#endif  // XGALOIS_POLY_SUBPRODUCT_TREE_HPP_
