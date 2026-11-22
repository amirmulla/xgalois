#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "xgalois/algebra/group_algorithms.hpp"
#include "xgalois/algebra/groups.hpp"
#include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
#include "xgalois/coding/gabidulin.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/interleaved_rs.hpp"
#include "xgalois/coding/lrc.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/ring/ideal.hpp"
#include "xgalois/ring/integer_ring.hpp"
#include "xgalois/ring/modular_ring.hpp"
#include "xgalois/ring/quotient_ring.hpp"
#include "xgalois/utils/bigint.hpp"

namespace nb = nanobind;

namespace xg {
namespace python {
namespace {

using PrimeField = GaloisFieldPrime<uint32_t>;
using BinaryExtField = GaloisFieldBinaryExtension<uint32_t>;
using U32Array1D =
    nb::ndarray<uint32_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
using ConstU32Array1D =
    nb::ndarray<const uint32_t, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

std::pair<uint64_t, std::string> BinaryCoeffsToMaskAndPoly(
    const std::vector<uint32_t>& coeffs, uint8_t degree) {
  if (coeffs.size() != static_cast<std::size_t>(degree) + 1) {
    throw std::invalid_argument(
        "Binary extension field modulus must have exactly m + 1 ascending "
        "coefficients");
  }
  uint64_t mask = 0;
  for (std::size_t deg = 0; deg <= static_cast<std::size_t>(degree); ++deg) {
    if (coeffs[deg] > 1U) {
      throw std::invalid_argument(
          "Binary extension field modulus coefficients must be 0 or 1");
    }
    if (coeffs[deg] == 1U) {
      mask |= (uint64_t{1} << deg);
    }
  }
  if (((mask >> degree) & 1ULL) == 0ULL) {
    throw std::invalid_argument(
        "Binary extension field modulus must be monic of degree m");
  }
  std::string poly;
  for (std::size_t deg = static_cast<std::size_t>(degree) + 1; deg-- > 0;) {
    if (((mask >> deg) & 1ULL) == 0ULL) {
      continue;
    }
    if (!poly.empty()) {
      poly += "+";
    }
    if (deg == 0) {
      poly += "1";
    } else if (deg == 1) {
      poly += "x";
    } else {
      poly += "x^" + std::to_string(deg);
    }
  }
  return {mask, poly};
}

template <typename GaloisField>
GFArray<GaloisField> MakeGFArray(const std::shared_ptr<GaloisField>& field,
                                 const std::vector<uint32_t>& values) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;
  RawArray raw = RawArray::from_shape(ShapeType{values.size()});
  for (std::size_t i = 0; i < values.size(); ++i) {
    raw(i) = field->SetElementValue(static_cast<StorageType>(values[i]));
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
}

template <typename GaloisField>
std::vector<uint32_t> GFArrayToVec(const std::shared_ptr<GaloisField>& field,
                                   const GFArray<GaloisField>& vec) {
  const auto span = vec.RawSpan();
  std::vector<uint32_t> out;
  out.reserve(span.size());
  for (const auto& raw_v : span) {
    out.push_back(static_cast<uint32_t>(field->GetElementValue(raw_v)));
  }
  return out;
}

template <typename GaloisField>
GFMatrix<GaloisField> MakeGFMatrix(
    const std::shared_ptr<GaloisField>& field,
    const std::vector<std::vector<uint32_t>>& rows) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;
  if (rows.empty()) {
    throw std::invalid_argument("Matrix cannot be empty");
  }
  const std::size_t num_rows = rows.size();
  const std::size_t num_cols = rows[0].size();
  RawArray raw = RawArray::from_shape(ShapeType{num_rows, num_cols});
  for (std::size_t r = 0; r < num_rows; ++r) {
    if (rows[r].size() != num_cols) {
      throw std::invalid_argument("Ragged 2D matrix input");
    }
    for (std::size_t c = 0; c < num_cols; ++c) {
      raw(r, c) = field->SetElementValue(static_cast<StorageType>(rows[r][c]));
    }
  }
  return GFMatrix<GaloisField>::FromRaw(field, std::move(raw));
}

template <typename GaloisField>
std::vector<std::vector<uint32_t>> GFMatrixToVec(
    const std::shared_ptr<GaloisField>& field,
    const GFMatrix<GaloisField>& mat) {
  const std::size_t rows = mat.Shape(0);
  const std::size_t cols = mat.Shape(1);
  std::vector<std::vector<uint32_t>> out(rows, std::vector<uint32_t>(cols));
  for (std::size_t r = 0; r < rows; ++r) {
    for (std::size_t c = 0; c < cols; ++c) {
      out[r][c] = static_cast<uint32_t>(field->GetElementValue(mat.Raw(r, c)));
    }
  }
  return out;
}

template <typename GaloisField>
std::vector<uint32_t> PolyToVec(const std::shared_ptr<GaloisField>& field,
                                const PolynomialDense<GaloisField>& poly) {
  std::vector<uint32_t> out;
  if (poly.IsZero()) {
    return out;
  }
  const int deg = poly.Degree();
  out.reserve(static_cast<std::size_t>(deg) + 1);
  for (int i = 0; i <= deg; ++i) {
    out.push_back(static_cast<uint32_t>(
        field->GetElementValue(poly.RawCoeff(static_cast<std::size_t>(i)))));
  }
  return out;
}

std::optional<coding::LRCPartitionType> ParsePartitionType(
    const std::optional<std::string>& pt) {
  if (!pt.has_value()) {
    return std::nullopt;
  }
  if (*pt == "multiplicative") {
    return coding::LRCPartitionType::kMultiplicative;
  }
  if (*pt == "additive") {
    return coding::LRCPartitionType::kAdditive;
  }
  throw std::invalid_argument("Unsupported partition_type: '" + *pt + "'");
}

coding::IRSDecoderMode ParseIrsMode(const std::string& mode_str) {
  if (mode_str == "ssb") {
    return coding::IRSDecoderMode::kSsb;
  }
  if (mode_str == "ssb_exact") {
    return coding::IRSDecoderMode::kSsbExact;
  }
  if (mode_str == "bky_hankel") {
    return coding::IRSDecoderMode::kBkyHankel;
  }
  if (mode_str == "independent") {
    return coding::IRSDecoderMode::kIndependent;
  }
  throw std::invalid_argument("Unknown IRS decoder mode: '" + mode_str + "'");
}

// ============================================================================
// Field Wrappers
// ============================================================================

struct PyGFPrime {
  std::shared_ptr<PrimeField> field;

  explicit PyGFPrime(uint32_t prime, bool prime_testing = true)
      : field(std::make_shared<PrimeField>(prime, "int", prime_testing)) {}

  uint64_t Order() const { return field->Order(); }
  uint64_t Characteristic() const { return field->Characteristic(); }
  uint32_t Degree() const { return field->Degree(); }
  uint64_t Modulus() const { return field->Modulus(); }

  uint32_t Add(uint32_t a, uint32_t b) const { return field->Add(a, b); }
  uint32_t Sub(uint32_t a, uint32_t b) const { return field->Sub(a, b); }
  uint32_t Mul(uint32_t a, uint32_t b) const { return field->Mul(a, b); }
  uint32_t Div(uint32_t a, uint32_t b) const { return field->Div(a, b); }
  uint32_t Inv(uint32_t a) const { return field->Inv(a); }
  uint32_t Pow(uint32_t a, uint64_t exp) const { return field->Pow(a, exp); }

  void AddArray(U32Array1D dst, ConstU32Array1D a, ConstU32Array1D b) const {
    const std::size_t n = dst.shape(0);
    if (a.shape(0) != n || b.shape(0) != n) {
      throw std::invalid_argument("add_array: array lengths must match");
    }
    uint32_t* dst_ptr = dst.data();
    const uint32_t* a_ptr = a.data();
    const uint32_t* b_ptr = b.data();
    nb::gil_scoped_release release;
    field->AddArray(std::span<uint32_t>(dst_ptr, n),
                    std::span<const uint32_t>(a_ptr, n),
                    std::span<const uint32_t>(b_ptr, n));
  }

  void MulArray(U32Array1D dst, ConstU32Array1D a, ConstU32Array1D b) const {
    const std::size_t n = dst.shape(0);
    if (a.shape(0) != n || b.shape(0) != n) {
      throw std::invalid_argument("mul_array: array lengths must match");
    }
    uint32_t* dst_ptr = dst.data();
    const uint32_t* a_ptr = a.data();
    const uint32_t* b_ptr = b.data();
    nb::gil_scoped_release release;
    field->MulArray(std::span<uint32_t>(dst_ptr, n),
                    std::span<const uint32_t>(a_ptr, n),
                    std::span<const uint32_t>(b_ptr, n));
  }

  void FmaArray(U32Array1D dst, ConstU32Array1D a, uint32_t scalar) const {
    const std::size_t n = dst.shape(0);
    if (a.shape(0) != n) {
      throw std::invalid_argument("fma_array: array lengths must match");
    }
    uint32_t* dst_ptr = dst.data();
    const uint32_t* a_ptr = a.data();
    nb::gil_scoped_release release;
    field->FmaArray(std::span<uint32_t>(dst_ptr, n),
                    std::span<const uint32_t>(a_ptr, n), scalar);
  }
};

struct PyGF2X {
  std::shared_ptr<BinaryExtField> field;

  explicit PyGF2X(uint8_t m, const std::string& modulus = "",
                  bool check_irreducible = true) {
    if (m < 2 || m > 16) {
      throw std::invalid_argument("GF2X degree m must satisfy 2 <= m <= 16");
    }
    field = std::make_shared<BinaryExtField>(m, "int", modulus, "α",
                                             check_irreducible);
  }

  PyGF2X(uint8_t m, const std::vector<uint32_t>& modulus_coeffs,
         bool check_irreducible = true) {
    if (m < 2 || m > 16) {
      throw std::invalid_argument("GF2X degree m must satisfy 2 <= m <= 16");
    }
    const auto [mask, poly_str] = BinaryCoeffsToMaskAndPoly(modulus_coeffs, m);
    field = std::make_shared<BinaryExtField>(m, "int", poly_str, "α",
                                             check_irreducible);
    if (field->Modulus() != mask) {
      throw std::runtime_error("GF2X modulus bitmask mismatch");
    }
  }

  uint64_t Order() const { return field->Order(); }
  uint64_t Characteristic() const { return field->Characteristic(); }
  uint32_t Degree() const { return field->Degree(); }
  uint64_t Modulus() const { return field->Modulus(); }

  uint32_t Add(uint32_t a, uint32_t b) const { return field->Add(a, b); }
  uint32_t Sub(uint32_t a, uint32_t b) const { return field->Sub(a, b); }
  uint32_t Mul(uint32_t a, uint32_t b) const { return field->Mul(a, b); }
  uint32_t Div(uint32_t a, uint32_t b) const { return field->Div(a, b); }
  uint32_t Inv(uint32_t a) const { return field->Inv(a); }
  uint32_t Pow(uint32_t a, uint64_t exp) const { return field->Pow(a, exp); }

  void AddArray(U32Array1D dst, ConstU32Array1D a, ConstU32Array1D b) const {
    const std::size_t n = dst.shape(0);
    if (a.shape(0) != n || b.shape(0) != n) {
      throw std::invalid_argument("add_array: array lengths must match");
    }
    uint32_t* dst_ptr = dst.data();
    const uint32_t* a_ptr = a.data();
    const uint32_t* b_ptr = b.data();
    nb::gil_scoped_release release;
    field->AddArray(std::span<uint32_t>(dst_ptr, n),
                    std::span<const uint32_t>(a_ptr, n),
                    std::span<const uint32_t>(b_ptr, n));
  }

  void MulArray(U32Array1D dst, ConstU32Array1D a, ConstU32Array1D b) const {
    const std::size_t n = dst.shape(0);
    if (a.shape(0) != n || b.shape(0) != n) {
      throw std::invalid_argument("mul_array: array lengths must match");
    }
    uint32_t* dst_ptr = dst.data();
    const uint32_t* a_ptr = a.data();
    const uint32_t* b_ptr = b.data();
    nb::gil_scoped_release release;
    field->MulArray(std::span<uint32_t>(dst_ptr, n),
                    std::span<const uint32_t>(a_ptr, n),
                    std::span<const uint32_t>(b_ptr, n));
  }

  void FmaArray(U32Array1D dst, ConstU32Array1D a, uint32_t scalar) const {
    const std::size_t n = dst.shape(0);
    if (a.shape(0) != n) {
      throw std::invalid_argument("fma_array: array lengths must match");
    }
    uint32_t* dst_ptr = dst.data();
    const uint32_t* a_ptr = a.data();
    nb::gil_scoped_release release;
    field->FmaArray(std::span<uint32_t>(dst_ptr, n),
                    std::span<const uint32_t>(a_ptr, n), scalar);
  }
};

// ============================================================================
// Tamo-Barg & Availability LRC Wrapper
// ============================================================================

struct PyTamoBargLRC {
  using LrcPrime = coding::LocallyRecoverableCode<PrimeField>;
  using LrcBinary = coding::LocallyRecoverableCode<BinaryExtField>;

  std::variant<LrcPrime, LrcBinary> code;

  explicit PyTamoBargLRC(LrcPrime c) : code(std::move(c)) {}
  explicit PyTamoBargLRC(LrcBinary c) : code(std::move(c)) {}

  std::size_t Length() const {
    return std::visit([](const auto& c) { return c.Length(); }, code);
  }
  std::size_t Dimension() const {
    return std::visit([](const auto& c) { return c.Dimension(); }, code);
  }
  std::size_t Locality() const {
    return std::visit([](const auto& c) { return c.Locality(0); }, code);
  }
  std::size_t LocalDistance() const {
    return std::visit([](const auto& c) { return c.LocalDistance(0); }, code);
  }
  std::size_t Availability() const {
    return std::visit([](const auto& c) { return c.Availability(); }, code);
  }
  std::size_t MinimumDistance() const {
    return std::visit([](const auto& c) { return c.MinimumDistance(); }, code);
  }
  std::size_t DesignedDistance() const {
    return std::visit([](const auto& c) { return c.DesignedDistance(); }, code);
  }

  std::vector<std::vector<std::size_t>> RepairSets(std::size_t coord) const {
    return std::visit([coord](const auto& c) { return c.RepairSets(coord); },
                      code);
  }

  std::vector<std::vector<std::size_t>> FiberPartition(
      std::size_t repair_set_idx = 0) const {
    return std::visit(
        [repair_set_idx](const auto& c) {
          return c.FiberPartition(repair_set_idx);
        },
        code);
  }

  std::vector<uint32_t> EvaluationPoints() const {
    return std::visit(
        [](const auto& c) {
          return GFArrayToVec(c.Field(), c.EvaluationPoints());
        },
        code);
  }

  std::vector<uint32_t> GoodPoly() const {
    return std::visit(
        [](const auto& c) { return PolyToVec(c.Field(), c.GoodPolynomial()); },
        code);
  }

  std::vector<uint32_t> Encode(const std::vector<uint32_t>& msg) const {
    return std::visit(
        [&msg](const auto& c) {
          auto gf = c.Field();
          auto arr = MakeGFArray(gf, msg);
          auto cw = c.Encode(arr);
          return GFArrayToVec(gf, cw);
        },
        code);
  }

  std::vector<uint32_t> Unencode(const std::vector<uint32_t>& cw) const {
    return std::visit(
        [&cw](const auto& c) {
          auto gf = c.Field();
          auto arr = MakeGFArray(gf, cw);
          auto msg = c.Unencode(arr);
          return GFArrayToVec(gf, msg);
        },
        code);
  }

  nb::dict LocalRepair(const std::vector<uint32_t>& received,
                       std::size_t erased_index,
                       std::size_t repair_set_idx = 0) const {
    return std::visit(
        [&received, erased_index, repair_set_idx](const auto& c) {
          auto gf = c.Field();
          auto recv = MakeGFArray(gf, received);
          auto res = [&]() {
            nb::gil_scoped_release release;
            return c.LocalRepair(recv, erased_index, repair_set_idx);
          }();
          nb::dict out;
          out["success"] = res.Ok();
          out["recovered_symbol"] =
              res.Ok() ? static_cast<uint32_t>(gf->GetElementValue(
                             res.recovered_symbol.RawValue()))
                       : uint32_t{0};
          out["erased_index"] = res.erased_index;
          out["repair_set_index"] = res.repair_set_index;
          out["symbols_read"] = res.symbols_read;
          out["error_message"] = res.error_message;
          return out;
        },
        code);
  }

  nb::dict LocalRepairMultiple(const std::vector<uint32_t>& received,
                               const std::vector<std::size_t>& erased_positions,
                               bool allow_global_fallback = true) const {
    return std::visit(
        [&received, &erased_positions, allow_global_fallback](const auto& c) {
          auto gf = c.Field();
          auto recv = MakeGFArray(gf, received);
          auto res = [&]() {
            nb::gil_scoped_release release;
            return c.LocalRepairMultiple(recv, erased_positions,
                                         allow_global_fallback);
          }();
          nb::dict out;
          out["success"] = res.Ok();
          out["codeword"] = res.Ok() ? GFArrayToVec(gf, res.codeword)
                                     : std::vector<uint32_t>{};
          out["message"] = res.Ok() ? GFArrayToVec(gf, res.message)
                                    : std::vector<uint32_t>{};
          out["num_local_repairs"] = res.num_local_repairs;
          out["num_global_repairs"] = res.num_global_repairs;
          out["total_symbols_read"] = res.total_symbols_read;
          out["used_global_fallback"] = res.used_global_fallback;
          out["error_message"] = res.error_message;
          return out;
        },
        code);
  }
};

PyTamoBargLRC CreateTamoBargLrcPrime(
    const PyGFPrime& gf, std::size_t n, std::size_t k, std::size_t r,
    std::size_t delta = 2,
    const std::optional<std::string>& partition_type = std::nullopt) {
  return PyTamoBargLRC(
      coding::LocallyRecoverableCode<PrimeField>::CreateTamoBarg(
          gf.field, n, k, r, delta, ParsePartitionType(partition_type)));
}

PyTamoBargLRC CreateTamoBargLrcBinary(
    const PyGF2X& gf, std::size_t n, std::size_t k, std::size_t r,
    std::size_t delta = 2,
    const std::optional<std::string>& partition_type = std::nullopt) {
  return PyTamoBargLRC(
      coding::LocallyRecoverableCode<BinaryExtField>::CreateTamoBarg(
          gf.field, n, k, r, delta, ParsePartitionType(partition_type)));
}

PyTamoBargLRC CreateAvailabilityLrcPrime(
    const PyGFPrime& gf, std::size_t n, std::size_t k,
    const std::vector<std::size_t>& localities,
    const std::optional<std::vector<std::size_t>>& local_distances =
        std::nullopt) {
  const std::vector<std::size_t> deltas =
      local_distances.value_or(std::vector<std::size_t>{});
  return PyTamoBargLRC(
      coding::LocallyRecoverableCode<
          PrimeField>::CreateMultiplicativeAvailability(gf.field, n, k,
                                                        localities, deltas));
}

PyTamoBargLRC CreateAvailabilityLrcBinary(
    const PyGF2X& gf, std::size_t n, std::size_t k,
    const std::vector<std::size_t>& localities,
    const std::optional<std::vector<std::size_t>>& local_distances =
        std::nullopt) {
  const std::vector<std::size_t> deltas =
      local_distances.value_or(std::vector<std::size_t>{});
  return PyTamoBargLRC(
      coding::LocallyRecoverableCode<
          BinaryExtField>::CreateMultiplicativeAvailability(gf.field, n, k,
                                                            localities,
                                                            deltas));
}

// ============================================================================
// Interleaved GRS Wrapper
// ============================================================================

struct PyInterleavedRS {
  using IrsPrime = coding::InterleavedGRSCode<PrimeField>;
  using IrsBinary = coding::InterleavedGRSCode<BinaryExtField>;

  std::variant<IrsPrime, IrsBinary> code;

  explicit PyInterleavedRS(IrsPrime c) : code(std::move(c)) {}
  explicit PyInterleavedRS(IrsBinary c) : code(std::move(c)) {}

  template <typename GF>
  static coding::InterleavedGRSCode<GF> BuildWithExplicitPoints(
      const std::shared_ptr<GF>& gf, const std::vector<uint32_t>& eval_pts,
      const std::vector<std::size_t>& dims,
      const std::optional<std::vector<uint32_t>>& col_mults) {
    auto pts = MakeGFArray(gf, eval_pts);
    auto mults = col_mults.has_value()
                     ? MakeGFArray(gf, *col_mults)
                     : GFArray<GF>::Ones(gf, {eval_pts.size()});
    return coding::InterleavedGRSCode<GF>(gf, pts, mults, dims);
  }

  static PyInterleavedRS FromPrimeHomogeneous(
      const PyGFPrime& gf, const std::vector<uint32_t>& eval_pts, std::size_t k,
      std::size_t ell,
      const std::optional<std::vector<uint32_t>>& col_mults = std::nullopt) {
    return PyInterleavedRS(BuildWithExplicitPoints(
        gf.field, eval_pts, std::vector<std::size_t>(ell, k), col_mults));
  }

  static PyInterleavedRS FromPrimeHeterogeneous(
      const PyGFPrime& gf, const std::vector<uint32_t>& eval_pts,
      const std::vector<std::size_t>& dims,
      const std::optional<std::vector<uint32_t>>& col_mults = std::nullopt) {
    return PyInterleavedRS(
        BuildWithExplicitPoints(gf.field, eval_pts, dims, col_mults));
  }

  static PyInterleavedRS FromBinaryHomogeneous(
      const PyGF2X& gf, const std::vector<uint32_t>& eval_pts, std::size_t k,
      std::size_t ell,
      const std::optional<std::vector<uint32_t>>& col_mults = std::nullopt) {
    return PyInterleavedRS(BuildWithExplicitPoints(
        gf.field, eval_pts, std::vector<std::size_t>(ell, k), col_mults));
  }

  static PyInterleavedRS FromBinaryHeterogeneous(
      const PyGF2X& gf, const std::vector<uint32_t>& eval_pts,
      const std::vector<std::size_t>& dims,
      const std::optional<std::vector<uint32_t>>& col_mults = std::nullopt) {
    return PyInterleavedRS(
        BuildWithExplicitPoints(gf.field, eval_pts, dims, col_mults));
  }

  std::size_t Length() const {
    return std::visit([](const auto& c) { return c.Length(); }, code);
  }
  std::size_t InterleavingOrder() const {
    return std::visit([](const auto& c) { return c.InterleavingOrder(); },
                      code);
  }
  std::size_t CollaborativeRadius() const {
    return std::visit(
        [](const auto& c) { return c.CollaborativeDecodingRadius(); }, code);
  }
  std::size_t UniqueRadius() const {
    return std::visit([](const auto& c) { return c.UniqueDecodingRadius(); },
                      code);
  }
  std::vector<std::size_t> Dimensions() const {
    return std::visit([](const auto& c) { return c.Dimensions(); }, code);
  }

  std::vector<std::vector<uint32_t>> Encode(
      const std::vector<std::vector<uint32_t>>& messages) const {
    return std::visit(
        [&messages](const auto& c) {
          using GF = typename std::decay_t<decltype(c)>::field_type;
          auto gf = c.Field();
          std::vector<GFArray<GF>> msg_vecs;
          msg_vecs.reserve(messages.size());
          for (const auto& m : messages) {
            msg_vecs.push_back(MakeGFArray(gf, m));
          }
          auto cw = c.Encode(msg_vecs);
          return GFMatrixToVec(gf, cw);
        },
        code);
  }

  nb::dict Decode(
      const std::vector<std::vector<uint32_t>>& received_matrix,
      const std::string& mode = "ssb", bool fallback = true,
      const std::vector<std::size_t>& erasure_positions = {}) const {
    const coding::IRSDecoderMode parsed_mode = ParseIrsMode(mode);
    return std::visit(
        [&received_matrix, parsed_mode, fallback,
         &erasure_positions](const auto& c) {
          auto gf = c.Field();
          auto recv = MakeGFMatrix(gf, received_matrix);
          auto res = [&]() {
            nb::gil_scoped_release release;
            return erasure_positions.empty()
                       ? c.Decode(recv, parsed_mode, fallback)
                       : c.DecodeWithErasures(recv, erasure_positions,
                                              parsed_mode, fallback);
          }();
          nb::dict out;
          out["success"] = res.Ok();
          out["collaborative_radius"] = c.CollaborativeDecodingRadius();
          out["unique_radius"] = c.UniqueDecodingRadius();
          out["used_fallback"] = res.used_fallback;
          out["burst_error_positions"] = res.burst_error_positions;
          if (res.Ok()) {
            out["codeword_matrix"] = GFMatrixToVec(gf, res.codeword_matrix);
            std::vector<std::vector<uint32_t>> msgs;
            msgs.reserve(res.messages.size());
            for (const auto& m : res.messages) {
              msgs.push_back(GFArrayToVec(gf, m));
            }
            out["messages"] = std::move(msgs);
            out["error_message"] = std::string{};
          } else {
            out["codeword_matrix"] = std::vector<std::vector<uint32_t>>{};
            out["messages"] = std::vector<std::vector<uint32_t>>{};
            out["error_message"] = res.error_message;
          }
          return out;
        },
        code);
  }
};

// ============================================================================
// Generalized Reed-Solomon Wrapper
// ============================================================================

struct PyGRSCode {
  using GrsPrime = coding::GeneralizedReedSolomonCode<PrimeField>;
  using GrsBinary = coding::GeneralizedReedSolomonCode<BinaryExtField>;

  std::variant<GrsPrime, GrsBinary> code;

  explicit PyGRSCode(GrsPrime c) : code(std::move(c)) {}
  explicit PyGRSCode(GrsBinary c) : code(std::move(c)) {}

  template <typename GF>
  static coding::GeneralizedReedSolomonCode<GF> Build(
      const std::shared_ptr<GF>& gf, const std::vector<uint32_t>& eval_pts,
      std::size_t k, const std::optional<std::vector<uint32_t>>& col_mults) {
    auto pts = MakeGFArray(gf, eval_pts);
    auto mults = col_mults.has_value()
                     ? MakeGFArray(gf, *col_mults)
                     : GFArray<GF>::Ones(gf, {eval_pts.size()});
    return coding::GeneralizedReedSolomonCode<GF>(gf, pts, mults, k);
  }

  static PyGRSCode FromPrime(
      const PyGFPrime& gf, const std::vector<uint32_t>& eval_pts, std::size_t k,
      const std::optional<std::vector<uint32_t>>& col_mults = std::nullopt) {
    return PyGRSCode(Build(gf.field, eval_pts, k, col_mults));
  }

  static PyGRSCode FromBinary(
      const PyGF2X& gf, const std::vector<uint32_t>& eval_pts, std::size_t k,
      const std::optional<std::vector<uint32_t>>& col_mults = std::nullopt) {
    return PyGRSCode(Build(gf.field, eval_pts, k, col_mults));
  }

  std::size_t Length() const {
    return std::visit([](const auto& c) { return c.Length(); }, code);
  }
  std::size_t Dimension() const {
    return std::visit([](const auto& c) { return c.Dimension(); }, code);
  }
  std::size_t MinimumDistance() const {
    return std::visit([](const auto& c) { return c.MinimumDistance(); }, code);
  }

  std::vector<uint32_t> Encode(const std::vector<uint32_t>& msg) const {
    return std::visit(
        [&msg](const auto& c) {
          auto gf = c.Field();
          auto arr = MakeGFArray(gf, msg);
          auto cw = c.Encode(arr);
          return GFArrayToVec(gf, cw);
        },
        code);
  }

  std::vector<uint32_t> Unencode(const std::vector<uint32_t>& cw) const {
    return std::visit(
        [&cw](const auto& c) {
          auto gf = c.Field();
          auto arr = MakeGFArray(gf, cw);
          auto msg = c.Unencode(arr);
          return GFArrayToVec(gf, msg);
        },
        code);
  }

  nb::dict DecodeBM(const std::vector<uint32_t>& received) const {
    return std::visit(
        [&received](const auto& c) {
          auto gf = c.Field();
          auto recv = MakeGFArray(gf, received);
          auto res = [&]() {
            nb::gil_scoped_release release;
            return c.DecodeWithResult(recv, "BerlekampMassey");
          }();
          nb::dict out;
          out["success"] = res.Ok();
          if (res.Ok()) {
            out["codeword"] = GFArrayToVec(gf, res.codeword);
            out["message"] = GFArrayToVec(gf, res.message);
            out["num_errors"] = res.num_errors;
            out["error_message"] = std::string{};
          } else {
            out["codeword"] = std::vector<uint32_t>{};
            out["message"] = std::vector<uint32_t>{};
            out["num_errors"] = std::size_t{0};
            out["error_message"] = res.error_message;
          }
          return out;
        },
        code);
  }

  nb::dict DecodePower(
      const std::vector<uint32_t>& received,
      std::optional<std::size_t> powering_degree = std::nullopt,
      const std::string& mode = "ssb", bool fallback_to_bm = true,
      const std::vector<std::size_t>& erasure_positions = {}) const {
    const coding::IRSDecoderMode parsed_mode = ParseIrsMode(mode);
    return std::visit(
        [&received, powering_degree, parsed_mode, fallback_to_bm,
         &erasure_positions](const auto& c) {
          using GF = typename std::decay_t<decltype(c)>::field_type;
          auto gf = c.Field();
          auto recv = MakeGFArray(gf, received);
          coding::PowerDecodingDecoder<GF> decoder(c, powering_degree,
                                                   parsed_mode, fallback_to_bm);
          auto res = [&]() {
            nb::gil_scoped_release release;
            return erasure_positions.empty()
                       ? decoder.DecodeWithResult(recv)
                       : decoder.DecodeWithErasures(recv, erasure_positions);
          }();
          nb::dict out;
          out["success"] = res.Ok();
          out["decoding_radius"] = decoder.DecodingRadius();
          out["powering_degree"] = decoder.PoweringDegree();
          if (res.Ok()) {
            out["codeword"] = GFArrayToVec(gf, res.codeword);
            out["message"] = GFArrayToVec(gf, res.message);
            out["num_errors"] = res.num_errors;
            out["error_message"] = std::string{};
          } else {
            out["codeword"] = std::vector<uint32_t>{};
            out["message"] = std::vector<uint32_t>{};
            out["num_errors"] = std::size_t{0};
            out["error_message"] = res.error_message;
          }
          return out;
        },
        code);
  }

  nb::dict DecodeGuruswamiSudan(const std::vector<uint32_t>& received,
                                std::size_t tau) const {
    return std::visit(
        [&received, tau](const auto& c) {
          using GF = typename std::decay_t<decltype(c)>::field_type;
          auto gf = c.Field();
          auto recv = MakeGFArray(gf, received);
          auto decoder =
              coding::GuruswamiSudanDecoder<GF>::FromTargetRadius(c, tau);
          auto cw_list = [&]() {
            nb::gil_scoped_release release;
            return decoder.DecodeToList(recv);
          }();
          std::vector<std::vector<uint32_t>> codewords;
          std::vector<std::vector<uint32_t>> messages;
          codewords.reserve(cw_list.size());
          messages.reserve(cw_list.size());
          for (const auto& cw : cw_list) {
            codewords.push_back(GFArrayToVec(gf, cw));
            messages.push_back(GFArrayToVec(gf, c.Unencode(cw)));
          }
          nb::dict out;
          out["codewords"] = std::move(codewords);
          out["messages"] = std::move(messages);
          out["decoding_radius"] = decoder.DecodingRadius();
          out["johnson_radius"] = decoder.JohnsonRadius();
          out["multiplicity"] = decoder.Multiplicity();
          out["list_size"] = decoder.ListSize();
          return out;
        },
        code);
  }
};

// ============================================================================
// Gabidulin Binary Code Wrapper
// ============================================================================

struct PyGabidulinBinaryCode {
  using GabCode = coding::GabidulinCode<GaloisFieldBinary, BinaryExtField>;

  std::shared_ptr<BinaryExtField> ext_field;
  GabCode code;

  PyGabidulinBinaryCode(const PyGF2X& gf, std::size_t n, std::size_t k)
      : ext_field(gf.field),
        code(GabCode::CreateStandard(std::make_shared<GaloisFieldBinary>(),
                                     gf.field, n, k)) {}

  std::size_t Length() const { return code.Length(); }
  std::size_t Dimension() const { return code.Dimension(); }
  std::size_t MinimumDistance() const { return code.MinimumDistance(); }

  std::size_t RankWeight(const std::vector<uint32_t>& word) const {
    auto arr = MakeGFArray(ext_field, word);
    return code.Weight(arr);
  }

  std::vector<uint32_t> Encode(const std::vector<uint32_t>& msg) const {
    auto arr = MakeGFArray(ext_field, msg);
    auto cw = code.Encode(arr);
    return GFArrayToVec(ext_field, cw);
  }

  nb::dict Decode(const std::vector<uint32_t>& received,
                  const std::string& decoder = "GabidulinSkewEEA") const {
    auto recv = MakeGFArray(ext_field, received);
    auto res = [&]() {
      nb::gil_scoped_release release;
      return code.DecodeWithResult(recv, decoder);
    }();
    nb::dict out;
    out["success"] = res.Ok();
    out["minimum_distance"] = code.MinimumDistance();
    if (res.Ok()) {
      out["codeword"] = GFArrayToVec(ext_field, res.codeword);
      out["message"] = GFArrayToVec(ext_field, res.message);
      out["error_rank"] = code.Weight(recv - res.codeword);
      out["error_message"] = std::string{};
    } else {
      out["codeword"] = std::vector<uint32_t>{};
      out["message"] = std::vector<uint32_t>{};
      out["error_rank"] = std::size_t{0};
      out["error_message"] = res.error_message;
    }
    return out;
  }
};

// ============================================================================
// Algebra & Rings Free Functions
// ============================================================================

std::string SymmetricGroupOrder(uint32_t n) {
  const algebra::SymmetricGroup sym(n);
  return sym.Order().ToDecimal();
}

std::string AlternatingGroupOrder(uint32_t n) {
  const algebra::AlternatingGroup alt(n);
  return alt.Order().ToDecimal();
}

std::string DihedralGroupOrder(uint64_t n) {
  const algebra::DihedralGroup dih(n);
  return dih.Order().ToDecimal();
}

std::vector<uint64_t> NegacyclicMul(uint64_t q, std::size_t n,
                                    const std::vector<uint64_t>& a,
                                    const std::vector<uint64_t>& b) {
  const ring::NegacyclicRing<ring::ModularRing<uint64_t>> nr(
      ring::ModularRing<uint64_t>(q), n);
  return nr.Mul(a, b);
}

std::optional<std::pair<std::string, std::string>> CrtInt(
    const std::vector<int64_t>& residues, const std::vector<int64_t>& moduli) {
  const ring::IntegerRing<BigInt> zz;
  std::vector<BigInt> res_big;
  std::vector<BigInt> mod_big;
  res_big.reserve(residues.size());
  mod_big.reserve(moduli.size());
  for (int64_t r : residues) {
    res_big.emplace_back(r);
  }
  for (int64_t m : moduli) {
    mod_big.emplace_back(m);
  }
  const auto ans = ring::ChineseRemainderTheorem(zz, res_big, mod_big);
  if (!ans.has_value()) {
    return std::nullopt;
  }
  return std::make_pair(ans->first.ToDecimal(), ans->second.ToDecimal());
}

}  // namespace
}  // namespace python
}  // namespace xg

NB_MODULE(_xgalois, m) {
  using namespace xg::python;

  m.doc() = "XGalois C++20 nanobind Python bindings";

  // GFPrime
  nb::class_<PyGFPrime>(m, "GFPrime")
      .def(nb::init<uint32_t, bool>(), nb::arg("prime"),
           nb::arg("prime_testing") = true)
      .def_prop_ro("order", &PyGFPrime::Order)
      .def_prop_ro("characteristic", &PyGFPrime::Characteristic)
      .def_prop_ro("degree", &PyGFPrime::Degree)
      .def_prop_ro("modulus", &PyGFPrime::Modulus)
      .def("add", &PyGFPrime::Add, nb::arg("a"), nb::arg("b"))
      .def("sub", &PyGFPrime::Sub, nb::arg("a"), nb::arg("b"))
      .def("mul", &PyGFPrime::Mul, nb::arg("a"), nb::arg("b"))
      .def("div", &PyGFPrime::Div, nb::arg("a"), nb::arg("b"))
      .def("inv", &PyGFPrime::Inv, nb::arg("a"))
      .def("pow", &PyGFPrime::Pow, nb::arg("a"), nb::arg("exp"))
      .def("add_array", &PyGFPrime::AddArray, nb::arg("dst"), nb::arg("a"),
           nb::arg("b"))
      .def("mul_array", &PyGFPrime::MulArray, nb::arg("dst"), nb::arg("a"),
           nb::arg("b"))
      .def("fma_array", &PyGFPrime::FmaArray, nb::arg("dst"), nb::arg("a"),
           nb::arg("scalar"));

  // GF2X
  nb::class_<PyGF2X>(m, "GF2X")
      .def(nb::init<uint8_t, const std::string&, bool>(), nb::arg("m"),
           nb::arg("modulus") = "", nb::arg("check_irreducible") = true)
      .def(nb::init<uint8_t, const std::vector<uint32_t>&, bool>(),
           nb::arg("m"), nb::arg("modulus"),
           nb::arg("check_irreducible") = true)
      .def_prop_ro("order", &PyGF2X::Order)
      .def_prop_ro("characteristic", &PyGF2X::Characteristic)
      .def_prop_ro("degree", &PyGF2X::Degree)
      .def_prop_ro("modulus", &PyGF2X::Modulus)
      .def("add", &PyGF2X::Add, nb::arg("a"), nb::arg("b"))
      .def("sub", &PyGF2X::Sub, nb::arg("a"), nb::arg("b"))
      .def("mul", &PyGF2X::Mul, nb::arg("a"), nb::arg("b"))
      .def("div", &PyGF2X::Div, nb::arg("a"), nb::arg("b"))
      .def("inv", &PyGF2X::Inv, nb::arg("a"))
      .def("pow", &PyGF2X::Pow, nb::arg("a"), nb::arg("exp"))
      .def("add_array", &PyGF2X::AddArray, nb::arg("dst"), nb::arg("a"),
           nb::arg("b"))
      .def("mul_array", &PyGF2X::MulArray, nb::arg("dst"), nb::arg("a"),
           nb::arg("b"))
      .def("fma_array", &PyGF2X::FmaArray, nb::arg("dst"), nb::arg("a"),
           nb::arg("scalar"));

  // TamoBargLRC
  nb::class_<PyTamoBargLRC>(m, "TamoBargLRC")
      .def_prop_ro("length", &PyTamoBargLRC::Length)
      .def_prop_ro("dimension", &PyTamoBargLRC::Dimension)
      .def_prop_ro("locality", &PyTamoBargLRC::Locality)
      .def_prop_ro("local_distance", &PyTamoBargLRC::LocalDistance)
      .def_prop_ro("availability", &PyTamoBargLRC::Availability)
      .def_prop_ro("minimum_distance", &PyTamoBargLRC::MinimumDistance)
      .def_prop_ro("designed_distance", &PyTamoBargLRC::DesignedDistance)
      .def("repair_sets", &PyTamoBargLRC::RepairSets, nb::arg("coord"))
      .def("fiber_partition", &PyTamoBargLRC::FiberPartition,
           nb::arg("repair_set_idx") = 0)
      .def("evaluation_points", &PyTamoBargLRC::EvaluationPoints)
      .def("good_poly", &PyTamoBargLRC::GoodPoly)
      .def("encode", &PyTamoBargLRC::Encode, nb::arg("msg"))
      .def("unencode", &PyTamoBargLRC::Unencode, nb::arg("codeword"))
      .def("local_repair", &PyTamoBargLRC::LocalRepair, nb::arg("received"),
           nb::arg("erased_index"), nb::arg("repair_set_idx") = 0)
      .def("local_repair_multiple", &PyTamoBargLRC::LocalRepairMultiple,
           nb::arg("received"), nb::arg("erased_positions"),
           nb::arg("allow_global_fallback") = true);

  m.def("create_tamo_barg_lrc_prime", &CreateTamoBargLrcPrime, nb::arg("field"),
        nb::arg("n"), nb::arg("k"), nb::arg("r"), nb::arg("delta") = 2,
        nb::arg("partition_type") = nb::none());
  m.def("create_tamo_barg_lrc_binary", &CreateTamoBargLrcBinary,
        nb::arg("field"), nb::arg("n"), nb::arg("k"), nb::arg("r"),
        nb::arg("delta") = 2, nb::arg("partition_type") = nb::none());
  m.def("create_availability_lrc_prime", &CreateAvailabilityLrcPrime,
        nb::arg("field"), nb::arg("n"), nb::arg("k"), nb::arg("localities"),
        nb::arg("local_distances") = nb::none());
  m.def("create_availability_lrc_binary", &CreateAvailabilityLrcBinary,
        nb::arg("field"), nb::arg("n"), nb::arg("k"), nb::arg("localities"),
        nb::arg("local_distances") = nb::none());

  // InterleavedRS
  nb::class_<PyInterleavedRS>(m, "InterleavedRS")
      .def(nb::new_(&PyInterleavedRS::FromPrimeHomogeneous), nb::arg("field"),
           nb::arg("evaluation_points"), nb::arg("k"), nb::arg("ell"),
           nb::arg("column_multipliers") = nb::none())
      .def(nb::new_(&PyInterleavedRS::FromPrimeHeterogeneous), nb::arg("field"),
           nb::arg("evaluation_points"), nb::arg("dims"),
           nb::arg("column_multipliers") = nb::none())
      .def(nb::new_(&PyInterleavedRS::FromBinaryHomogeneous), nb::arg("field"),
           nb::arg("evaluation_points"), nb::arg("k"), nb::arg("ell"),
           nb::arg("column_multipliers") = nb::none())
      .def(nb::new_(&PyInterleavedRS::FromBinaryHeterogeneous),
           nb::arg("field"), nb::arg("evaluation_points"), nb::arg("dims"),
           nb::arg("column_multipliers") = nb::none())
      .def_prop_ro("length", &PyInterleavedRS::Length)
      .def_prop_ro("interleaving_order", &PyInterleavedRS::InterleavingOrder)
      .def_prop_ro("collaborative_radius",
                   &PyInterleavedRS::CollaborativeRadius)
      .def_prop_ro("unique_radius", &PyInterleavedRS::UniqueRadius)
      .def_prop_ro("dimensions", &PyInterleavedRS::Dimensions)
      .def("encode", &PyInterleavedRS::Encode, nb::arg("messages"))
      .def("decode", &PyInterleavedRS::Decode, nb::arg("received_matrix"),
           nb::arg("mode") = "ssb", nb::arg("fallback") = true,
           nb::arg("erasure_positions") = std::vector<std::size_t>{});

  // GRSCode
  nb::class_<PyGRSCode>(m, "GRSCode")
      .def(nb::new_(&PyGRSCode::FromPrime), nb::arg("field"),
           nb::arg("evaluation_points"), nb::arg("k"),
           nb::arg("column_multipliers") = nb::none())
      .def(nb::new_(&PyGRSCode::FromBinary), nb::arg("field"),
           nb::arg("evaluation_points"), nb::arg("k"),
           nb::arg("column_multipliers") = nb::none())
      .def_prop_ro("length", &PyGRSCode::Length)
      .def_prop_ro("dimension", &PyGRSCode::Dimension)
      .def_prop_ro("minimum_distance", &PyGRSCode::MinimumDistance)
      .def("encode", &PyGRSCode::Encode, nb::arg("msg"))
      .def("unencode", &PyGRSCode::Unencode, nb::arg("codeword"))
      .def("decode_bm", &PyGRSCode::DecodeBM, nb::arg("received"))
      .def("decode_power", &PyGRSCode::DecodePower, nb::arg("received"),
           nb::arg("powering_degree") = nb::none(), nb::arg("mode") = "ssb",
           nb::arg("fallback_to_bm") = true,
           nb::arg("erasure_positions") = std::vector<std::size_t>{})
      .def("decode_guruswami_sudan", &PyGRSCode::DecodeGuruswamiSudan,
           nb::arg("received"), nb::arg("tau"));

  // GabidulinBinaryCode
  nb::class_<PyGabidulinBinaryCode>(m, "GabidulinBinaryCode")
      .def(nb::init<const PyGF2X&, std::size_t, std::size_t>(),
           nb::arg("ext_field"), nb::arg("n"), nb::arg("k"))
      .def_prop_ro("length", &PyGabidulinBinaryCode::Length)
      .def_prop_ro("dimension", &PyGabidulinBinaryCode::Dimension)
      .def_prop_ro("minimum_distance", &PyGabidulinBinaryCode::MinimumDistance)
      .def("rank_weight", &PyGabidulinBinaryCode::RankWeight, nb::arg("word"))
      .def("encode", &PyGabidulinBinaryCode::Encode, nb::arg("msg"))
      .def("decode", &PyGabidulinBinaryCode::Decode, nb::arg("received"),
           nb::arg("decoder") = "GabidulinSkewEEA");

  // Algebra & Rings
  m.def("symmetric_group_order", &SymmetricGroupOrder, nb::arg("n"));
  m.def("alternating_group_order", &AlternatingGroupOrder, nb::arg("n"));
  m.def("dihedral_group_order", &DihedralGroupOrder, nb::arg("n"));
  m.def("negacyclic_mul", &NegacyclicMul, nb::arg("q"), nb::arg("N"),
        nb::arg("a"), nb::arg("b"));
  m.def("crt_int", &CrtInt, nb::arg("residues"), nb::arg("moduli"));
}
