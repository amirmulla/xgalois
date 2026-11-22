#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "test/oracle/tiny_json.hpp"
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
#include "xgalois/ring/matrix_ring.hpp"
#include "xgalois/ring/modular_ring.hpp"
#include "xgalois/ring/quotient_ring.hpp"
#include "xgalois/utils/bigint.hpp"

namespace xg {
namespace oracle {
namespace {

constexpr const char* kOracleCliVersion = "0.4.0";

std::pair<uint64_t, std::string> ParseBinaryModulusMaskAndPoly(
    const JsonValue& modulus_arr, uint32_t degree) {
  const auto& coeffs = modulus_arr.AsArray();
  if (coeffs.size() != static_cast<std::size_t>(degree) + 1) {
    throw std::invalid_argument(
        "Binary extension field modulus must have exactly m + 1 ascending "
        "coefficients");
  }
  uint64_t mask = 0;
  std::string poly;
  for (std::size_t deg = 0; deg <= static_cast<std::size_t>(degree); ++deg) {
    const uint64_t c = coeffs[deg].AsUint64();
    if (c > 1ULL) {
      throw std::invalid_argument(
          "Binary extension field modulus coefficients must be 0 or 1");
    }
    if (c == 1ULL) {
      mask |= (uint64_t{1} << deg);
    }
  }
  if (((mask >> degree) & 1ULL) == 0ULL) {
    throw std::invalid_argument(
        "Binary extension field modulus must be monic of degree m");
  }
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

template <typename Fn>
JsonValue WithField(const JsonValue& field_json, Fn&& fn) {
  const uint64_t p = field_json.At("p").AsUint64();
  const uint64_t m =
      field_json.Contains("m") ? field_json.At("m").AsUint64() : 1ULL;

  if (m == 1ULL) {
    if (p < 2ULL || p > 0x7FFFFFFFULL) {
      throw std::invalid_argument("Unsupported prime field characteristic p");
    }
    auto gf =
        std::make_shared<GaloisFieldPrime<uint32_t>>(static_cast<uint32_t>(p));
    return fn(gf);
  }

  if (p == 2ULL && m >= 2ULL && m <= 16ULL) {
    if (!field_json.Contains("modulus")) {
      throw std::invalid_argument(
          "Binary extension field requires explicit ascending 'modulus' array");
    }
    const auto [mask, poly_str] = ParseBinaryModulusMaskAndPoly(
        field_json.At("modulus"), static_cast<uint32_t>(m));
    auto gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(
        static_cast<uint8_t>(m), "int", poly_str, "α",
        /*check_irreducible=*/true);
    if (gf->Modulus() != mask) {
      throw std::runtime_error(
          "Constructed binary extension field modulus bitmask mismatch");
    }
    return fn(gf);
  }

  throw std::invalid_argument(
      "Unsupported field specification: expected prime field (m == 1) or "
      "binary extension field (p == 2, 2 <= m <= 16)");
}

std::vector<std::size_t> ExtractSizeArray(const JsonValue& arr) {
  std::vector<std::size_t> out;
  out.reserve(arr.Size());
  for (const auto& elem : arr.AsArray()) {
    out.push_back(static_cast<std::size_t>(elem.AsUint64()));
  }
  return out;
}

std::vector<uint64_t> ExtractUint64Array(const JsonValue& arr) {
  std::vector<uint64_t> out;
  out.reserve(arr.Size());
  for (const auto& elem : arr.AsArray()) {
    out.push_back(elem.AsUint64());
  }
  return out;
}

JsonValue SizeVecToJson(const std::vector<std::size_t>& vec) {
  JsonValue::Array arr;
  arr.reserve(vec.size());
  for (std::size_t v : vec) {
    arr.emplace_back(static_cast<int64_t>(v));
  }
  return JsonValue(std::move(arr));
}

JsonValue SizeMatToJson(const std::vector<std::vector<std::size_t>>& mat) {
  JsonValue::Array arr;
  arr.reserve(mat.size());
  for (const auto& row : mat) {
    arr.push_back(SizeVecToJson(row));
  }
  return JsonValue(std::move(arr));
}

JsonValue Uint64VecToJson(const std::vector<uint64_t>& vec) {
  JsonValue::Array arr;
  arr.reserve(vec.size());
  for (uint64_t v : vec) {
    arr.emplace_back(static_cast<int64_t>(v));
  }
  return JsonValue(std::move(arr));
}

template <typename GaloisField>
GFArray<GaloisField> MakeGFArrayFromJson(
    const std::shared_ptr<GaloisField>& field, const JsonValue& arr) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFArray<GaloisField>::raw_array_type;
  using ShapeType = typename GFArray<GaloisField>::shape_type;
  const auto& elems = arr.AsArray();
  RawArray raw = RawArray::from_shape(ShapeType{elems.size()});
  for (std::size_t i = 0; i < elems.size(); ++i) {
    raw(i) =
        field->SetElementValue(static_cast<StorageType>(elems[i].AsUint64()));
  }
  return GFArray<GaloisField>::FromRaw(field, std::move(raw));
}

template <typename GaloisField>
GFMatrix<GaloisField> MakeGFMatrixFromJson(
    const std::shared_ptr<GaloisField>& field, const JsonValue& mat_json) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;
  const auto& rows = mat_json.AsArray();
  if (rows.empty()) {
    throw std::invalid_argument("Matrix JSON cannot be empty");
  }
  const std::size_t num_rows = rows.size();
  const std::size_t num_cols = rows[0].AsArray().size();
  RawArray raw = RawArray::from_shape(ShapeType{num_rows, num_cols});
  for (std::size_t r = 0; r < num_rows; ++r) {
    const auto& row = rows[r].AsArray();
    if (row.size() != num_cols) {
      throw std::invalid_argument("Ragged 2D matrix JSON");
    }
    for (std::size_t c = 0; c < num_cols; ++c) {
      raw(r, c) =
          field->SetElementValue(static_cast<StorageType>(row[c].AsUint64()));
    }
  }
  return GFMatrix<GaloisField>::FromRaw(field, std::move(raw));
}

template <typename GaloisField>
JsonValue GFArrayToJson(const std::shared_ptr<GaloisField>& field,
                        const GFArray<GaloisField>& vec) {
  JsonValue::Array arr;
  const auto span = vec.RawSpan();
  arr.reserve(span.size());
  for (const auto& raw_v : span) {
    arr.emplace_back(static_cast<int64_t>(field->GetElementValue(raw_v)));
  }
  return JsonValue(std::move(arr));
}

template <typename GaloisField>
JsonValue GFMatrixToJson(const std::shared_ptr<GaloisField>& field,
                         const GFMatrix<GaloisField>& mat) {
  const std::size_t rows = mat.Shape(0);
  const std::size_t cols = mat.Shape(1);
  JsonValue::Array out;
  out.reserve(rows);
  for (std::size_t r = 0; r < rows; ++r) {
    JsonValue::Array row_arr;
    row_arr.reserve(cols);
    for (std::size_t c = 0; c < cols; ++c) {
      row_arr.emplace_back(
          static_cast<int64_t>(field->GetElementValue(mat.Raw(r, c))));
    }
    out.emplace_back(std::move(row_arr));
  }
  return JsonValue(std::move(out));
}

template <typename GaloisField>
JsonValue PolyToJson(const std::shared_ptr<GaloisField>& field,
                     const PolynomialDense<GaloisField>& poly) {
  JsonValue::Array arr;
  if (poly.IsZero()) {
    return JsonValue(std::move(arr));
  }
  const int deg = poly.Degree();
  arr.reserve(static_cast<std::size_t>(deg) + 1);
  for (int i = 0; i <= deg; ++i) {
    arr.emplace_back(static_cast<int64_t>(
        field->GetElementValue(poly.RawCoeff(static_cast<std::size_t>(i)))));
  }
  return JsonValue(std::move(arr));
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
// RPC Handlers
// ============================================================================

JsonValue HandlePing() {
  JsonValue::Object res;
  res["version"] = std::string(kOracleCliVersion);
  res["status"] = std::string("ok");
  return JsonValue(std::move(res));
}

JsonValue HandleLrcEncodeAndRepair(const JsonValue& params) {
  return WithField(params.At("field"), [&](auto gf) -> JsonValue {
    using GF = typename decltype(gf)::element_type;

    const std::size_t n = static_cast<std::size_t>(params.At("n").AsUint64());
    const std::size_t k = static_cast<std::size_t>(params.At("k").AsUint64());
    const std::size_t r = static_cast<std::size_t>(params.At("r").AsUint64());
    const std::size_t delta =
        params.Contains("delta")
            ? static_cast<std::size_t>(params.At("delta").AsUint64())
            : 2U;

    std::optional<coding::LRCPartitionType> ptype = std::nullopt;
    if (params.Contains("partition_type") &&
        !params.At("partition_type").IsNull()) {
      const std::string& pt = params.At("partition_type").AsString();
      if (pt == "multiplicative") {
        ptype = coding::LRCPartitionType::kMultiplicative;
      } else if (pt == "additive") {
        ptype = coding::LRCPartitionType::kAdditive;
      } else {
        throw std::invalid_argument("Unsupported partition_type: '" + pt + "'");
      }
    }

    const auto code = coding::LocallyRecoverableCode<GF>::CreateTamoBarg(
        gf, n, k, r, delta, ptype);

    const auto msg = MakeGFArrayFromJson(gf, params.At("message"));
    const auto codeword = code.Encode(msg);
    const auto unencoded = code.Unencode(codeword);

    std::vector<std::size_t> erased_positions;
    if (params.Contains("erased_positions")) {
      erased_positions = ExtractSizeArray(params.At("erased_positions"));
    }
    const bool allow_global_fallback =
        params.Contains("allow_global_fallback")
            ? params.At("allow_global_fallback").AsBool()
            : true;

    // Perturb every erased symbol so repair cannot accidentally rely on the
    // original value at any erased index.
    auto corrupted = codeword;
    const auto one = gf->MultiplicativeIdentity();
    for (std::size_t pos : erased_positions) {
      if (pos < n) {
        corrupted.Raw(pos) = gf->Add(corrupted.Raw(pos), one);
      }
    }

    auto repaired = corrupted;
    bool repair_success = true;
    std::size_t symbols_read = 0;
    std::size_t num_local_repairs = 0;
    std::size_t num_global_repairs = 0;
    bool used_global_fallback = false;
    std::string repair_error;

    if (erased_positions.size() == 1) {
      const std::size_t pos = erased_positions[0];
      const auto res = code.LocalRepair(corrupted, pos, 0);
      repair_success = res.Ok();
      symbols_read = res.symbols_read;
      repair_error = res.error_message;
      if (res.Ok()) {
        repaired.Raw(pos) = res.recovered_symbol.RawValue();
        num_local_repairs = 1;
      }
    } else if (erased_positions.size() > 1) {
      const auto res = code.LocalRepairMultiple(corrupted, erased_positions,
                                                allow_global_fallback);
      repair_success = res.Ok();
      symbols_read = res.total_symbols_read;
      num_local_repairs = res.num_local_repairs;
      num_global_repairs = res.num_global_repairs;
      used_global_fallback = res.used_global_fallback;
      repair_error = res.error_message;
      if (res.Ok()) {
        repaired = res.codeword;
      }
    }

    JsonValue::Array repair_sets_json;
    repair_sets_json.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      repair_sets_json.push_back(SizeVecToJson(code.RepairSet(i, 0)));
    }

    JsonValue::Object out;
    out["codeword"] = GFArrayToJson(gf, codeword);
    out["repaired_codeword"] = GFArrayToJson(gf, repaired);
    out["minimum_distance"] = static_cast<int64_t>(code.MinimumDistance());
    out["designed_distance"] = static_cast<int64_t>(code.DesignedDistance());
    out["symbols_read"] = static_cast<int64_t>(symbols_read);
    out["evaluation_points"] = GFArrayToJson(gf, code.EvaluationPoints());
    out["repair_sets"] = JsonValue(std::move(repair_sets_json));
    out["fiber_partition"] = SizeMatToJson(code.FiberPartition(0));
    out["good_poly"] = PolyToJson(gf, code.GoodPolynomial());
    out["repair_success"] = repair_success;
    out["used_global_fallback"] = used_global_fallback;
    out["num_local_repairs"] = static_cast<int64_t>(num_local_repairs);
    out["num_global_repairs"] = static_cast<int64_t>(num_global_repairs);
    out["unencoded_message"] = GFArrayToJson(gf, unencoded);
    if (!repair_error.empty()) {
      out["repair_error"] = repair_error;
    }
    return JsonValue(std::move(out));
  });
}

JsonValue HandleLrcAvailabilityVerify(const JsonValue& params) {
  return WithField(params.At("field"), [&](auto gf) -> JsonValue {
    using GF = typename decltype(gf)::element_type;

    const std::string construction = params.Contains("construction")
                                         ? params.At("construction").AsString()
                                         : std::string("multiplicative");

    const std::vector<std::size_t> localities =
        ExtractSizeArray(params.At("localities"));
    std::vector<std::size_t> local_distances;
    if (params.Contains("local_distances") &&
        !params.At("local_distances").IsNull()) {
      local_distances = ExtractSizeArray(params.At("local_distances"));
    }

    auto build_code = [&]() {
      if (construction == "multiplicative") {
        const std::size_t n =
            static_cast<std::size_t>(params.At("n").AsUint64());
        const std::size_t k =
            static_cast<std::size_t>(params.At("k").AsUint64());
        return coding::LocallyRecoverableCode<
            GF>::CreateMultiplicativeAvailability(gf, n, k, localities,
                                                  local_distances);
      }
      if (construction == "additive") {
        const std::size_t num_superblocks =
            static_cast<std::size_t>(params.At("num_superblocks").AsUint64());
        const std::size_t message_superblocks = static_cast<std::size_t>(
            params.At("message_superblocks").AsUint64());
        return coding::LocallyRecoverableCode<GF>::CreateAdditiveAvailability(
            gf, num_superblocks, message_superblocks, localities,
            local_distances);
      }
      throw std::invalid_argument("Unsupported availability construction: '" +
                                  construction + "'");
    };

    const auto code = build_code();
    const std::size_t n = code.Length();
    const std::size_t k = code.Dimension();
    const std::size_t avail = code.Availability();

    const auto msg = MakeGFArrayFromJson(gf, params.At("message"));
    const auto codeword = code.Encode(msg);

    const std::size_t test_coord =
        params.Contains("test_coord")
            ? static_cast<std::size_t>(params.At("test_coord").AsUint64())
            : 0U;
    if (test_coord >= n) {
      throw std::invalid_argument("test_coord out of range");
    }

    JsonValue::Array all_repair_sets;
    all_repair_sets.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      all_repair_sets.push_back(SizeMatToJson(code.RepairSets(i)));
    }

    const auto one = gf->MultiplicativeIdentity();
    auto verify_coord_via_set = [&](std::size_t coord, std::size_t u) {
      const auto& rset = code.RepairSet(coord, u);
      std::vector<bool> in_rset(n, false);
      for (std::size_t idx : rset) {
        in_rset[idx] = true;
      }
      // Corrupt coord AND every coordinate outside repair set u so local repair
      // is forced to depend strictly on symbols in R_u(coord).
      auto corrupted = codeword;
      for (std::size_t j = 0; j < n; ++j) {
        if (!in_rset[j]) {
          corrupted.Raw(j) = gf->Add(corrupted.Raw(j), one);
        }
      }
      return code.LocalRepair(corrupted, coord, u);
    };

    JsonValue::Array repaired_via_each_set;
    repaired_via_each_set.reserve(avail);
    for (std::size_t u = 0; u < avail; ++u) {
      const auto res = verify_coord_via_set(test_coord, u);
      const uint64_t expected_val =
          gf->GetElementValue(codeword.Raw(test_coord));
      const uint64_t recovered_val =
          res.Ok() ? gf->GetElementValue(res.recovered_symbol.RawValue())
                   : 0ULL;
      JsonValue::Object entry;
      entry["repair_set_idx"] = static_cast<int64_t>(u);
      entry["success"] = res.Ok();
      entry["recovered_symbol"] = static_cast<int64_t>(recovered_val);
      entry["expected_symbol"] = static_cast<int64_t>(expected_val);
      entry["symbols_read"] = static_cast<int64_t>(res.symbols_read);
      entry["matches_original"] = res.Ok() && (recovered_val == expected_val);
      repaired_via_each_set.emplace_back(std::move(entry));
    }

    bool all_coords_all_sets_ok = true;
    for (std::size_t i = 0; i < n && all_coords_all_sets_ok; ++i) {
      for (std::size_t u = 0; u < avail; ++u) {
        const auto res = verify_coord_via_set(i, u);
        if (!res.Ok() || res.symbols_read != code.Locality(u) ||
            res.recovered_symbol.RawValue() != codeword.Raw(i)) {
          all_coords_all_sets_ok = false;
          break;
        }
      }
    }

    JsonValue::Object out;
    out["n"] = static_cast<int64_t>(n);
    out["k"] = static_cast<int64_t>(k);
    out["availability"] = static_cast<int64_t>(avail);
    out["designed_distance"] = static_cast<int64_t>(code.DesignedDistance());
    out["minimum_distance"] = static_cast<int64_t>(code.MinimumDistance());
    out["codeword"] = GFArrayToJson(gf, codeword);
    out["evaluation_points"] = GFArrayToJson(gf, code.EvaluationPoints());
    out["all_repair_sets"] = JsonValue(std::move(all_repair_sets));
    out["repaired_via_each_set"] = JsonValue(std::move(repaired_via_each_set));
    out["all_coords_all_sets_ok"] = all_coords_all_sets_ok;
    out["test_coord"] = static_cast<int64_t>(test_coord);
    return JsonValue(std::move(out));
  });
}

JsonValue HandleIrsDecode(const JsonValue& params) {
  return WithField(params.At("field"), [&](auto gf) -> JsonValue {
    using GF = typename decltype(gf)::element_type;

    const auto eval_pts =
        MakeGFArrayFromJson(gf, params.At("evaluation_points"));
    const std::size_t n = eval_pts.Size();
    const auto col_mults =
        params.Contains("column_multipliers") &&
                !params.At("column_multipliers").IsNull()
            ? MakeGFArrayFromJson(gf, params.At("column_multipliers"))
            : GFArray<GF>::Ones(gf, {n});

    std::vector<std::size_t> dims;
    if (params.Contains("dims")) {
      dims = ExtractSizeArray(params.At("dims"));
    } else {
      const std::size_t k = static_cast<std::size_t>(params.At("k").AsUint64());
      const std::size_t ell =
          static_cast<std::size_t>(params.At("ell").AsUint64());
      dims.assign(ell, k);
    }

    const coding::InterleavedGRSCode<GF> irs(gf, eval_pts, col_mults, dims);
    const auto received =
        MakeGFMatrixFromJson(gf, params.At("received_matrix"));

    const std::string mode_str =
        params.Contains("mode") ? params.At("mode").AsString() : "ssb";
    const coding::IRSDecoderMode mode = ParseIrsMode(mode_str);
    const bool fallback = params.Contains("fallback_to_independent")
                              ? params.At("fallback_to_independent").AsBool()
                              : true;

    std::vector<std::size_t> erasure_positions;
    if (params.Contains("erasure_positions")) {
      erasure_positions = ExtractSizeArray(params.At("erasure_positions"));
    }

    const auto res = erasure_positions.empty()
                         ? irs.Decode(received, mode, fallback)
                         : irs.DecodeWithErasures(received, erasure_positions,
                                                  mode, fallback);

    JsonValue::Object out;
    out["success"] = res.Ok();
    out["collaborative_radius"] =
        static_cast<int64_t>(irs.CollaborativeDecodingRadius());
    out["unique_radius"] = static_cast<int64_t>(irs.UniqueDecodingRadius());
    out["used_fallback"] = res.used_fallback;
    out["burst_error_positions"] = SizeVecToJson(res.burst_error_positions);
    if (res.Ok()) {
      out["codeword_matrix"] = GFMatrixToJson(gf, res.codeword_matrix);
      JsonValue::Array msgs;
      msgs.reserve(res.messages.size());
      for (std::size_t idx = 0; idx < res.messages.size(); ++idx) {
        msgs.push_back(GFArrayToJson(gf, res.messages[idx]));
      }
      out["messages"] = JsonValue(std::move(msgs));
    } else {
      out["codeword_matrix"] = JsonValue(JsonValue::Array{});
      out["messages"] = JsonValue(JsonValue::Array{});
      out["error_message"] = res.error_message;
    }
    return JsonValue(std::move(out));
  });
}

JsonValue HandleGrsPowerDecode(const JsonValue& params) {
  return WithField(params.At("field"), [&](auto gf) -> JsonValue {
    using GF = typename decltype(gf)::element_type;

    const auto eval_pts =
        MakeGFArrayFromJson(gf, params.At("evaluation_points"));
    const std::size_t n = eval_pts.Size();
    const auto col_mults =
        params.Contains("column_multipliers") &&
                !params.At("column_multipliers").IsNull()
            ? MakeGFArrayFromJson(gf, params.At("column_multipliers"))
            : GFArray<GF>::Ones(gf, {n});
    const std::size_t k = static_cast<std::size_t>(params.At("k").AsUint64());

    const coding::GeneralizedReedSolomonCode<GF> grs(gf, eval_pts, col_mults,
                                                     k);
    std::optional<std::size_t> powering_degree = std::nullopt;
    if (params.Contains("powering_degree") &&
        !params.At("powering_degree").IsNull()) {
      powering_degree =
          static_cast<std::size_t>(params.At("powering_degree").AsUint64());
    }
    const std::string mode_str =
        params.Contains("mode") ? params.At("mode").AsString() : "ssb";
    const coding::IRSDecoderMode mode = ParseIrsMode(mode_str);
    const bool fallback_to_bm = params.Contains("fallback_to_bm")
                                    ? params.At("fallback_to_bm").AsBool()
                                    : true;

    const coding::PowerDecodingDecoder<GF> decoder(grs, powering_degree, mode,
                                                   fallback_to_bm);
    const auto received = MakeGFArrayFromJson(gf, params.At("received"));

    std::vector<std::size_t> erasure_positions;
    if (params.Contains("erasure_positions")) {
      erasure_positions = ExtractSizeArray(params.At("erasure_positions"));
    }

    const auto res =
        erasure_positions.empty()
            ? decoder.DecodeWithResult(received)
            : decoder.DecodeWithErasures(received, erasure_positions);

    JsonValue::Object out;
    out["success"] = res.Ok();
    out["decoding_radius"] = static_cast<int64_t>(decoder.DecodingRadius());
    out["powering_degree"] = static_cast<int64_t>(decoder.PoweringDegree());
    if (res.Ok()) {
      out["codeword"] = GFArrayToJson(gf, res.codeword);
      out["message"] = GFArrayToJson(gf, res.message);
      out["num_errors"] = static_cast<int64_t>(res.num_errors);
    } else {
      out["codeword"] = JsonValue(JsonValue::Array{});
      out["message"] = JsonValue(JsonValue::Array{});
      out["error_message"] = res.error_message;
    }
    return JsonValue(std::move(out));
  });
}

JsonValue HandleGrsGuruswamiSudan(const JsonValue& params) {
  return WithField(params.At("field"), [&](auto gf) -> JsonValue {
    using GF = typename decltype(gf)::element_type;

    const auto eval_pts =
        MakeGFArrayFromJson(gf, params.At("evaluation_points"));
    const std::size_t n = eval_pts.Size();
    const auto col_mults =
        params.Contains("column_multipliers") &&
                !params.At("column_multipliers").IsNull()
            ? MakeGFArrayFromJson(gf, params.At("column_multipliers"))
            : GFArray<GF>::Ones(gf, {n});
    const std::size_t k = static_cast<std::size_t>(params.At("k").AsUint64());

    const coding::GeneralizedReedSolomonCode<GF> grs(gf, eval_pts, col_mults,
                                                     k);

    auto build_decoder = [&]() {
      if (params.Contains("tau") && !params.At("tau").IsNull()) {
        const std::size_t tau =
            static_cast<std::size_t>(params.At("tau").AsUint64());
        return coding::GuruswamiSudanDecoder<GF>::FromTargetRadius(grs, tau);
      }
      const std::size_t mult =
          params.Contains("multiplicity") && !params.At("multiplicity").IsNull()
              ? static_cast<std::size_t>(params.At("multiplicity").AsUint64())
              : 1U;
      if (params.Contains("list_size") && !params.At("list_size").IsNull()) {
        const std::size_t list_size =
            static_cast<std::size_t>(params.At("list_size").AsUint64());
        return coding::GuruswamiSudanDecoder<GF>(grs, mult, list_size);
      }
      return coding::GuruswamiSudanDecoder<GF>(grs, mult);
    };

    const auto decoder = build_decoder();
    const auto received = MakeGFArrayFromJson(gf, params.At("received"));
    const auto cw_list = decoder.DecodeToList(received);

    JsonValue::Array codewords_json;
    JsonValue::Array messages_json;
    codewords_json.reserve(cw_list.size());
    messages_json.reserve(cw_list.size());
    for (std::size_t idx = 0; idx < cw_list.size(); ++idx) {
      codewords_json.push_back(GFArrayToJson(gf, cw_list[idx]));
      messages_json.push_back(GFArrayToJson(gf, grs.Unencode(cw_list[idx])));
    }

    JsonValue::Object out;
    out["codewords"] = JsonValue(std::move(codewords_json));
    out["messages"] = JsonValue(std::move(messages_json));
    out["decoding_radius"] = static_cast<int64_t>(decoder.DecodingRadius());
    out["johnson_radius"] = static_cast<int64_t>(decoder.JohnsonRadius());
    out["multiplicity"] = static_cast<int64_t>(decoder.Multiplicity());
    out["list_size"] = static_cast<int64_t>(decoder.ListSize());
    return JsonValue(std::move(out));
  });
}

JsonValue HandleGroupsAndRingsCheck(const JsonValue& params) {
  const std::string& sub_op = params.At("sub_op").AsString();

  if (sub_op == "symmetric_group") {
    const uint32_t deg = static_cast<uint32_t>(params.At("degree").AsUint64());
    const algebra::SymmetricGroup sym(deg);
    JsonValue::Object out;
    out["order"] = sym.Order().ToDecimal();
    out["bsgs_order"] = sym.Bsgs().Order().ToDecimal();
    if (deg <= 7U) {
      out["is_solvable"] = IsSolvable(sym);
    }
    return JsonValue(std::move(out));
  }

  if (sub_op == "alternating_group") {
    const uint32_t deg = static_cast<uint32_t>(params.At("degree").AsUint64());
    const algebra::AlternatingGroup alt(deg);
    JsonValue::Object out;
    out["order"] = alt.Order().ToDecimal();
    out["bsgs_order"] = alt.Bsgs().Order().ToDecimal();
    if (deg <= 7U) {
      out["is_solvable"] = IsSolvable(alt);
    }
    return JsonValue(std::move(out));
  }

  if (sub_op == "dihedral_group") {
    const uint64_t n = params.At("n").AsUint64();
    const algebra::DihedralGroup dih(n);
    JsonValue::Object out;
    out["order"] = dih.Order().ToDecimal();
    out["center_order"] = static_cast<int64_t>(Center(dih).Order());
    return JsonValue(std::move(out));
  }

  if (sub_op == "negacyclic_mul") {
    const uint64_t q = params.At("modulus_q").AsUint64();
    const std::size_t deg =
        static_cast<std::size_t>(params.At("degree").AsUint64());
    const ring::NegacyclicRing<ring::ModularRing<uint64_t>> nr(
        ring::ModularRing<uint64_t>(q), deg);
    const auto a = ExtractUint64Array(params.At("a"));
    const auto b = ExtractUint64Array(params.At("b"));
    const auto prod = nr.Mul(a, b);
    JsonValue::Object out;
    out["product"] = Uint64VecToJson(prod);
    return JsonValue(std::move(out));
  }

  if (sub_op == "crt_int") {
    const ring::IntegerRing<BigInt> zz;
    std::vector<BigInt> residues;
    std::vector<BigInt> moduli;
    for (const auto& r : params.At("residues").AsArray()) {
      residues.emplace_back(r.AsInt64());
    }
    for (const auto& m : params.At("moduli").AsArray()) {
      moduli.emplace_back(m.AsInt64());
    }
    const auto res = ring::ChineseRemainderTheorem(zz, residues, moduli);
    JsonValue::Object out;
    out["consistent"] = res.has_value();
    if (res.has_value()) {
      out["solution"] = res->first.ToDecimal();
      out["lcm"] = res->second.ToDecimal();
    }
    return JsonValue(std::move(out));
  }

  if (sub_op == "matrix_berkowitz") {
    const uint64_t mod = params.At("modulus").AsUint64();
    const std::size_t n = static_cast<std::size_t>(params.At("n").AsUint64());
    const ring::MatrixRing<ring::ModularRing<uint64_t>> mr(
        ring::ModularRing<uint64_t>(mod), n);
    const auto mat = ExtractUint64Array(params.At("matrix"));
    JsonValue::Object out;
    out["det"] = static_cast<int64_t>(mr.Det(mat));
    out["char_poly_desc"] = Uint64VecToJson(mr.CharacteristicPolynomial(mat));
    return JsonValue(std::move(out));
  }

  throw std::invalid_argument("Unknown groups_and_rings_check sub_op: '" +
                              sub_op + "'");
}

JsonValue HandleGabidulinDecode(const JsonValue& params) {
  const auto& fld = params.At("field");
  const uint64_t p = fld.At("p").AsUint64();
  const uint32_t m = static_cast<uint32_t>(fld.At("m").AsUint64());
  if (p != 2ULL || m < 2U || m > 16U) {
    throw std::invalid_argument(
        "gabidulin_decode currently supports binary extensions GF(2^m) with "
        "2 <= m <= 16");
  }
  const auto [mask, poly_str] =
      ParseBinaryModulusMaskAndPoly(fld.At("modulus"), m);
  auto base_gf = std::make_shared<GaloisFieldBinary>();
  auto ext_gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(
      static_cast<uint8_t>(m), "int", poly_str, "α",
      /*check_irreducible=*/true);
  if (static_cast<uint64_t>(ext_gf->Modulus()) != mask) {
    throw std::invalid_argument("Binary extension field modulus mismatch");
  }

  const std::size_t n = static_cast<std::size_t>(params.At("n").AsUint64());
  const std::size_t k = static_cast<std::size_t>(params.At("k").AsUint64());
  const auto code = coding::GabidulinCode<
      GaloisFieldBinary,
      GaloisFieldBinaryExtension<uint32_t>>::CreateStandard(base_gf, ext_gf, n,
                                                            k);
  const auto received = MakeGFArrayFromJson(ext_gf, params.At("received"));
  const std::string decoder_name = params.Contains("decoder")
                                       ? params.At("decoder").AsString()
                                       : "GabidulinSkewEEA";
  const auto res = code.DecodeWithResult(received, decoder_name);

  JsonValue::Object out;
  out["success"] = res.Ok();
  out["minimum_distance"] = static_cast<int64_t>(code.MinimumDistance());
  if (res.Ok()) {
    out["codeword"] = GFArrayToJson(ext_gf, res.codeword);
    out["message"] = GFArrayToJson(ext_gf, res.message);
    out["error_rank"] =
        static_cast<int64_t>(code.Weight(received - res.codeword));
  } else {
    out["codeword"] = JsonValue(JsonValue::Array{});
    out["message"] = JsonValue(JsonValue::Array{});
    out["error_message"] = res.error_message;
  }
  return JsonValue(std::move(out));
}

JsonValue DispatchMethod(const std::string& method, const JsonValue& params) {
  if (method == "ping") {
    return HandlePing();
  }
  if (method == "lrc_encode_and_repair") {
    return HandleLrcEncodeAndRepair(params);
  }
  if (method == "lrc_availability_verify") {
    return HandleLrcAvailabilityVerify(params);
  }
  if (method == "irs_decode") {
    return HandleIrsDecode(params);
  }
  if (method == "grs_power_decode") {
    return HandleGrsPowerDecode(params);
  }
  if (method == "grs_guruswami_sudan") {
    return HandleGrsGuruswamiSudan(params);
  }
  if (method == "gabidulin_decode") {
    return HandleGabidulinDecode(params);
  }
  if (method == "groups_and_rings_check") {
    return HandleGroupsAndRingsCheck(params);
  }
  throw std::invalid_argument("Unknown RPC method: '" + method + "'");
}

bool IsWhitespaceOnly(std::string_view sv) noexcept {
  for (char c : sv) {
    if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
      return false;
    }
  }
  return true;
}

}  // namespace
}  // namespace oracle
}  // namespace xg

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::string line;
  while (std::getline(std::cin, line)) {
    if (xg::oracle::IsWhitespaceOnly(line)) {
      continue;
    }

    xg::oracle::JsonValue req_id(nullptr);
    xg::oracle::JsonValue::Object resp;

    try {
      const xg::oracle::JsonValue req = xg::oracle::ParseJson(line);
      if (req.Contains("id")) {
        req_id = req.At("id");
      }
      const std::string& method = req.At("method").AsString();
      const xg::oracle::JsonValue params =
          req.Contains("params")
              ? req.At("params")
              : xg::oracle::JsonValue(xg::oracle::JsonValue::Object{});
      xg::oracle::JsonValue result = xg::oracle::DispatchMethod(method, params);

      resp["id"] = req_id;
      resp["ok"] = true;
      resp["result"] = std::move(result);
    } catch (const std::exception& ex) {
      resp["id"] = req_id;
      resp["ok"] = false;
      resp["error"] = std::string(ex.what());
    }

    std::cout << xg::oracle::DumpJson(xg::oracle::JsonValue(std::move(resp)), 0)
              << '\n'
              << std::flush;
  }

  return 0;
}
