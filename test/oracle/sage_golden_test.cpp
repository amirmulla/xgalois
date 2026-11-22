#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "test/oracle/tiny_json.hpp"
#include "xgalois/algebra/group_algorithms.hpp"
#include "xgalois/algebra/groups.hpp"
#include "xgalois/algebra/snf.hpp"
#include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
#include "xgalois/coding/grs.hpp"
#include "xgalois/coding/interleaved_rs.hpp"
#include "xgalois/coding/lrc.hpp"
#include "xgalois/field/gf_bigint.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_binary_large.hpp"
#include "xgalois/field/gf_galois_group.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/field/gf_tower.hpp"
#include "xgalois/ring/ideal.hpp"
#include "xgalois/ring/integer_ring.hpp"
#include "xgalois/ring/matrix_ring.hpp"
#include "xgalois/ring/modular_ring.hpp"
#include "xgalois/ring/quadratic_integers.hpp"
#include "xgalois/ring/quotient_ring.hpp"
#include "xgalois/ring/skew_polynomial_ring.hpp"
#include "xgalois/utils/bigint.hpp"

#if !defined(XGALOIS_GOLDEN_DIR) && defined(__has_include)
#if __has_include("rules_cc/cc/runfiles/runfiles.h")
#define XGALOIS_HAS_BAZEL_RUNFILES 1
#include "rules_cc/cc/runfiles/runfiles.h"
#endif
#endif

namespace xg {
namespace oracle {
namespace {

constexpr const char* kGoldenFilename = "v0.3_algebra_fields_poly.json";

std::string ResolveGoldenPath(const std::string& filename) {
  namespace fs = std::filesystem;
  std::error_code ec;

  // 1. Explicit environment variable override.
  if (const char* env_dir = std::getenv("XGALOIS_GOLDEN_DIR");
      env_dir != nullptr && env_dir[0] != '\0') {
    const fs::path candidate = fs::path(env_dir) / filename;
    if (fs::exists(candidate, ec) && !ec) {
      return candidate.string();
    }
  }

  // 2. Compile-time macro path (CMake build/test tree).
#ifdef XGALOIS_GOLDEN_DIR
  {
    ec.clear();
    const fs::path candidate = fs::path(XGALOIS_GOLDEN_DIR) / filename;
    if (fs::exists(candidate, ec) && !ec) {
      return candidate.string();
    }
  }
#endif

  // 3. Bazel rules_cc Runfiles lookup when XGALOIS_GOLDEN_DIR is not defined.
#if defined(XGALOIS_HAS_BAZEL_RUNFILES)
  {
    std::string error;
#ifdef BAZEL_CURRENT_REPOSITORY
    std::unique_ptr<rules_cc::cc::runfiles::Runfiles> runfiles(
        rules_cc::cc::runfiles::Runfiles::CreateForTest(
            BAZEL_CURRENT_REPOSITORY, &error));
#else
    std::unique_ptr<rules_cc::cc::runfiles::Runfiles> runfiles(
        rules_cc::cc::runfiles::Runfiles::CreateForTest(&error));
#endif
    if (runfiles != nullptr) {
      const std::string rel_candidates[] = {
          "_main/test/oracle/golden/" + filename,
          "__main__/test/oracle/golden/" + filename,
          "test/oracle/golden/" + filename,
      };
      for (const auto& rel : rel_candidates) {
        const std::string resolved = runfiles->Rlocation(rel);
        ec.clear();
        if (!resolved.empty() && fs::exists(resolved, ec) && !ec) {
          return resolved;
        }
      }
    }
  }
#endif

  // 4. Bazel environment directory fallback (TEST_SRCDIR / RUNFILES_DIR).
  for (const char* env_name : {"TEST_SRCDIR", "RUNFILES_DIR"}) {
    if (const char* base = std::getenv(env_name);
        base != nullptr && base[0] != '\0') {
      const fs::path base_path(base);
      const fs::path candidates[] = {
          base_path / "_main" / "test" / "oracle" / "golden" / filename,
          base_path / "__main__" / "test" / "oracle" / "golden" / filename,
          base_path / "test" / "oracle" / "golden" / filename,
      };
      for (const auto& candidate : candidates) {
        ec.clear();
        if (fs::exists(candidate, ec) && !ec) {
          return candidate.string();
        }
      }
    }
  }

  // 5. Workspace-relative fallbacks.
  const fs::path rel_fallbacks[] = {
      fs::path("test") / "oracle" / "golden" / filename,
      fs::path("..") / "test" / "oracle" / "golden" / filename,
      fs::path("..") / ".." / "test" / "oracle" / "golden" / filename,
  };
  for (const auto& candidate : rel_fallbacks) {
    ec.clear();
    if (fs::exists(candidate, ec) && !ec) {
      return candidate.string();
    }
  }

  throw std::runtime_error("Failed to resolve golden oracle JSON file: " +
                           filename);
}

std::string BinaryModulusPolyString(const JsonValue& modulus_arr) {
  const auto& coeffs = modulus_arr.AsArray();
  std::string poly;
  for (std::size_t deg = coeffs.size(); deg-- > 0;) {
    if ((coeffs[deg].AsUint64() & 1ULL) == 0ULL) {
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
  return poly.empty() ? "0" : poly;
}

std::vector<uint64_t> ExtractUint64Array(const JsonValue& arr) {
  std::vector<uint64_t> out;
  out.reserve(arr.Size());
  for (const auto& elem : arr.AsArray()) {
    out.push_back(elem.AsUint64());
  }
  return out;
}

std::vector<uint32_t> ExtractUint32Array(const JsonValue& arr) {
  std::vector<uint32_t> out;
  out.reserve(arr.Size());
  for (const auto& elem : arr.AsArray()) {
    out.push_back(static_cast<uint32_t>(elem.AsUint64()));
  }
  return out;
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

std::vector<std::size_t> ExtractSizeArray(const JsonValue& arr) {
  std::vector<std::size_t> out;
  out.reserve(arr.Size());
  for (const auto& elem : arr.AsArray()) {
    out.push_back(static_cast<std::size_t>(elem.AsUint64()));
  }
  return out;
}

template <typename GaloisField>
GFMatrix<GaloisField> MakeGFMatrixFromJson(
    const std::shared_ptr<GaloisField>& field, const JsonValue& mat_json) {
  using StorageType = typename GaloisField::element_type;
  using RawArray = typename GFMatrix<GaloisField>::raw_array_type;
  using ShapeType = typename GFMatrix<GaloisField>::shape_type;
  const auto& rows = mat_json.AsArray();
  const std::size_t num_rows = rows.size();
  const std::size_t num_cols = rows.empty() ? 0U : rows[0].AsArray().size();
  RawArray raw = RawArray::from_shape(ShapeType{num_rows, num_cols});
  for (std::size_t r = 0; r < num_rows; ++r) {
    const auto& row = rows[r].AsArray();
    for (std::size_t c = 0; c < num_cols; ++c) {
      raw(r, c) =
          field->SetElementValue(static_cast<StorageType>(row[c].AsUint64()));
    }
  }
  return GFMatrix<GaloisField>::FromRaw(field, std::move(raw));
}

// ============================================================================
// Vector Executors (one per canonical operation)
// ============================================================================

void VerifyGroupVector(const JsonValue& vec) {
  const std::string& op = vec.At("operation").AsString();

  if (op == "symmetric_group_order") {
    const uint32_t deg = static_cast<uint32_t>(vec.At("degree").AsUint64());
    const BigInt expected =
        BigInt::FromDecimal(vec.At("expected_order").AsString());
    const algebra::SymmetricGroup sym(deg);
    EXPECT_EQ(sym.Order(), expected) << vec.At("id").AsString();
    EXPECT_EQ(sym.Bsgs().Order(), expected) << vec.At("id").AsString();
  } else if (op == "alternating_group_properties") {
    const uint32_t deg = static_cast<uint32_t>(vec.At("degree").AsUint64());
    const BigInt expected =
        BigInt::FromDecimal(vec.At("expected_order").AsString());
    const algebra::AlternatingGroup alt(deg);
    EXPECT_EQ(alt.Order(), expected) << vec.At("id").AsString();
    EXPECT_EQ(alt.Bsgs().Order(), expected) << vec.At("id").AsString();
    EXPECT_EQ(IsSolvable(alt), vec.At("expected_is_solvable").AsBool())
        << vec.At("id").AsString();
  } else if (op == "dihedral_group_properties") {
    const uint64_t n = vec.At("n").AsUint64();
    const BigInt expected =
        BigInt::FromDecimal(vec.At("expected_order").AsString());
    const algebra::DihedralGroup dih(n);
    EXPECT_EQ(dih.Order(), expected) << vec.At("id").AsString();
    EXPECT_EQ(Center(dih).Order(), vec.At("expected_center_order").AsUint64())
        << vec.At("id").AsString();
  } else if (op == "symmetric_group_structure") {
    const uint32_t deg = static_cast<uint32_t>(vec.At("degree").AsUint64());
    const algebra::SymmetricGroup sym(deg);
    EXPECT_EQ(sym.Order(),
              BigInt::FromDecimal(vec.At("expected_order").AsString()));
    EXPECT_EQ(CommutatorSubgroup(sym).Order(),
              vec.At("expected_commutator_order").AsUint64());
    EXPECT_EQ(SylowSubgroup(sym, 2).Order(),
              vec.At("expected_sylow_2_order").AsUint64());
    EXPECT_EQ(SylowSubgroup(sym, 3).Order(),
              vec.At("expected_sylow_3_order").AsUint64());
    EXPECT_EQ(IsSolvable(sym), vec.At("expected_is_solvable").AsBool());
  } else if (op == "snf_abelian_decomposition") {
    std::vector<std::vector<BigInt>> relations;
    for (const auto& row : vec.At("relations").AsArray()) {
      std::vector<BigInt> r;
      for (const auto& entry : row.AsArray()) {
        r.push_back(BigInt(entry.AsInt64()));
      }
      relations.push_back(std::move(r));
    }
    const auto decomp =
        DecomposeFinitelyGeneratedAbelianGroup<BigInt>(relations);
    EXPECT_EQ(decomp.free_rank, vec.At("expected_free_rank").AsUint64());
    const auto& exp_factors = vec.At("expected_invariant_factors").AsArray();
    ASSERT_EQ(decomp.invariant_factors.size(), exp_factors.size());
    for (std::size_t i = 0; i < exp_factors.size(); ++i) {
      EXPECT_EQ(decomp.invariant_factors[i],
                BigInt::FromDecimal(exp_factors[i].AsString()));
    }
  } else {
    FAIL() << "Unhandled required group operation: " << op;
  }
}

void VerifyRingVector(const JsonValue& vec) {
  const std::string& op = vec.At("operation").AsString();

  if (op == "modular_ring_units") {
    const uint64_t mod = vec.At("modulus").AsUint64();
    const ring::ModularRing<uint64_t> z_mod(mod);
    const auto units = z_mod.UnitGroupElements();
    EXPECT_EQ(units, ExtractUint64Array(vec.At("expected_units")));
    EXPECT_EQ(units.size(), vec.At("expected_unit_count").AsUint64());
    EXPECT_EQ(z_mod.EulerTotient(),
              vec.At("expected_euler_totient").AsUint64());
    for (uint64_t u : units) {
      const auto inv = z_mod.TryInv(u);
      ASSERT_TRUE(inv.has_value());
      EXPECT_EQ(z_mod.Mul(u, *inv), 1ULL);
    }
  } else if (op == "gaussian_integer_gcd") {
    const ring::GaussianIntegerRing zi;
    for (const auto& tc : vec.At("cases").AsArray()) {
      const ring::GaussianInt a{BigInt(tc.At("a").At(0).AsInt64()),
                                BigInt(tc.At("a").At(1).AsInt64())};
      const ring::GaussianInt b{BigInt(tc.At("b").At(0).AsInt64()),
                                BigInt(tc.At("b").At(1).AsInt64())};
      const auto [g, s, t] = zi.ExtendedGcd(a, b);
      const ring::GaussianInt expected_g{
          BigInt(tc.At("expected_gcd_canonical").At(0).AsInt64()),
          BigInt(tc.At("expected_gcd_canonical").At(1).AsInt64())};
      EXPECT_EQ(g, expected_g);
      EXPECT_EQ(zi.Norm(g), BigInt(tc.At("expected_gcd_norm").AsInt64()));
      EXPECT_EQ(zi.Add(zi.Mul(s, a), zi.Mul(t, b)), g);
      EXPECT_TRUE(zi.IsZero(zi.DivRem(a, g).second));
      EXPECT_TRUE(zi.IsZero(zi.DivRem(b, g).second));
    }
  } else if (op == "eisenstein_integer_gcd") {
    const ring::EisensteinIntegerRing zw;
    for (const auto& tc : vec.At("cases").AsArray()) {
      const ring::EisensteinInt a{BigInt(tc.At("a").At(0).AsInt64()),
                                  BigInt(tc.At("a").At(1).AsInt64())};
      const ring::EisensteinInt b{BigInt(tc.At("b").At(0).AsInt64()),
                                  BigInt(tc.At("b").At(1).AsInt64())};
      const auto [g, s, t] = zw.ExtendedGcd(a, b);
      const ring::EisensteinInt expected_g{
          BigInt(tc.At("expected_gcd_canonical").At(0).AsInt64()),
          BigInt(tc.At("expected_gcd_canonical").At(1).AsInt64())};
      EXPECT_EQ(g, expected_g);
      EXPECT_EQ(zw.Norm(g), BigInt(tc.At("expected_gcd_norm").AsInt64()));
      EXPECT_EQ(zw.Add(zw.Mul(s, a), zw.Mul(t, b)), g);
      EXPECT_TRUE(zw.IsZero(zw.DivRem(a, g).second));
      EXPECT_TRUE(zw.IsZero(zw.DivRem(b, g).second));
    }
  } else if (op == "negacyclic_ring_mul") {
    const uint64_t q = vec.At("modulus_q").AsUint64();
    const std::size_t deg =
        static_cast<std::size_t>(vec.At("degree").AsUint64());
    const ring::NegacyclicRing<ring::ModularRing<uint64_t>> nr(
        ring::ModularRing<uint64_t>(q), deg);
    const auto a = ExtractUint64Array(vec.At("a"));
    const auto b = ExtractUint64Array(vec.At("b"));
    const auto expected = ExtractUint64Array(vec.At("expected_product"));
    EXPECT_EQ(nr.Mul(a, b), expected);
  } else if (op == "matrix_ring_berkowitz_det") {
    const uint64_t mod = vec.At("modulus").AsUint64();
    const std::size_t n = static_cast<std::size_t>(vec.At("n").AsUint64());
    const ring::MatrixRing<ring::ModularRing<uint64_t>> mr(
        ring::ModularRing<uint64_t>(mod), n);
    const auto mat = ExtractUint64Array(vec.At("matrix"));
    EXPECT_EQ(mr.Det(mat), vec.At("expected_det").AsUint64());
    EXPECT_EQ(mr.CharacteristicPolynomial(mat),
              ExtractUint64Array(vec.At("expected_char_poly_desc")));
  } else if (op == "generalized_crt_int") {
    const ring::IntegerRing<BigInt> zz;
    for (const auto& tc : vec.At("cases").AsArray()) {
      std::vector<BigInt> residues;
      std::vector<BigInt> moduli;
      for (const auto& r : tc.At("residues").AsArray()) {
        residues.push_back(BigInt(r.AsInt64()));
      }
      for (const auto& m : tc.At("moduli").AsArray()) {
        moduli.push_back(BigInt(m.AsInt64()));
      }
      const auto res = ring::ChineseRemainderTheorem(zz, residues, moduli);
      const bool expected_consistent = tc.At("expected_consistent").AsBool();
      ASSERT_EQ(res.has_value(), expected_consistent)
          << tc.At("name").AsString();
      if (expected_consistent) {
        EXPECT_EQ(res->first, BigInt(tc.At("expected_solution").AsInt64()));
        EXPECT_EQ(res->second, BigInt(tc.At("expected_lcm").AsInt64()));
      }
    }
  } else if (op == "frobenius_skew_ring") {
    const auto& fld_json = vec.At("field");
    const uint8_t m = static_cast<uint8_t>(fld_json.At("m").AsUint64());
    const std::string poly_str =
        BinaryModulusPolyString(fld_json.At("modulus"));
    const GaloisFieldBinaryExtension<uint32_t> gf(m, "int", poly_str);
    const uint32_t step =
        static_cast<uint32_t>(vec.At("frobenius_step").AsUint64());
    const ring::FrobeniusSkewRing<GaloisFieldBinaryExtension<uint32_t>> skew(
        gf, step);

    const auto eval_poly = ExtractUint32Array(vec.At("eval_poly"));
    const uint32_t eval_pt =
        static_cast<uint32_t>(vec.At("eval_point").AsUint64());
    EXPECT_EQ(skew.EvalOperator(eval_poly, eval_pt),
              static_cast<uint32_t>(vec.At("expected_eval").AsUint64()));

    for (const auto& mc : vec.At("minimal_poly_cases").AsArray()) {
      const auto pts = ExtractUint32Array(mc.At("points"));
      const auto min_poly = skew.MinimalPolynomial(pts);
      EXPECT_EQ(skew.Degree(min_poly), mc.At("expected_degree").AsInt64());
      EXPECT_EQ(min_poly, ExtractUint32Array(mc.At("expected_coeffs")));
      for (uint32_t pt : pts) {
        EXPECT_EQ(skew.EvalOperator(min_poly, pt), 0U);
      }
    }
  } else {
    FAIL() << "Unhandled required ring operation: " << op;
  }
}

void VerifyFieldVector(const JsonValue& vec) {
  const std::string& op = vec.At("operation").AsString();

  if (op == "gf2_extension_arith") {
    const auto& fld_json = vec.At("field");
    const uint8_t m = static_cast<uint8_t>(fld_json.At("m").AsUint64());
    const std::string poly_str =
        BinaryModulusPolyString(fld_json.At("modulus"));
    const GaloisFieldBinaryExtension<uint32_t> gf(m, "int", poly_str);
    for (const auto& tc : vec.At("cases").AsArray()) {
      const uint32_t a = static_cast<uint32_t>(tc.At("a").AsUint64());
      const uint32_t b = static_cast<uint32_t>(tc.At("b").AsUint64());
      const uint64_t exp = tc.At("exp").AsUint64();
      EXPECT_EQ(gf.Mul(a, b),
                static_cast<uint32_t>(tc.At("expected_mul").AsUint64()));
      EXPECT_EQ(gf.Inv(a),
                static_cast<uint32_t>(tc.At("expected_inv").AsUint64()));
      EXPECT_EQ(gf.Pow(a, exp),
                static_cast<uint32_t>(tc.At("expected_pow").AsUint64()));
    }
  } else if (op == "field_tower_relative_ops") {
    const auto& fld_json = vec.At("field");
    const uint32_t p = static_cast<uint32_t>(fld_json.At("p").AsUint64());
    const GaloisFieldPrime<uint32_t> base_fp(p);
    const auto mod_coeffs = ExtractUint32Array(fld_json.At("modulus"));
    const FieldTower<GaloisFieldPrime<uint32_t>, 2> tower(base_fp, mod_coeffs);

    const auto a = tower.FromCoeffs(ExtractUint32Array(vec.At("a")));
    const auto b = tower.FromCoeffs(ExtractUint32Array(vec.At("b")));
    const auto exp_mul =
        tower.FromCoeffs(ExtractUint32Array(vec.At("expected_mul")));
    const auto exp_inv =
        tower.FromCoeffs(ExtractUint32Array(vec.At("expected_inv")));

    EXPECT_TRUE(tower.Equals(tower.Mul(a, b), exp_mul));
    EXPECT_TRUE(tower.Equals(tower.Inv(a), exp_inv));
    EXPECT_EQ(
        tower.RelativeTrace(a),
        static_cast<uint32_t>(vec.At("expected_relative_trace").AsUint64()));
    EXPECT_EQ(
        tower.RelativeNorm(a),
        static_cast<uint32_t>(vec.At("expected_relative_norm").AsUint64()));
  } else if (op == "cyclotomic_cosets") {
    const uint64_t n = vec.At("n").AsUint64();
    const uint64_t q = vec.At("q").AsUint64();
    const auto actual = CyclotomicCosets(n, q);
    const auto& exp_arr = vec.At("expected_cosets").AsArray();
    ASSERT_EQ(actual.size(), exp_arr.size());
    for (std::size_t i = 0; i < exp_arr.size(); ++i) {
      EXPECT_EQ(actual[i], ExtractUint64Array(exp_arr[i]));
    }
  } else if (op == "prime_big_256_arith") {
    const auto fp =
        GaloisFieldPrimeBig<256>::FromHex(vec.At("p_hex").AsString(), false);
    const auto a = BigUInt<256>::FromHex(vec.At("a_hex").AsString());
    const auto b = BigUInt<256>::FromHex(vec.At("b_hex").AsString());
    const auto exp_mul =
        BigUInt<256>::FromHex(vec.At("expected_mul_hex").AsString());
    const auto exp_inv =
        BigUInt<256>::FromHex(vec.At("expected_inv_hex").AsString());

    EXPECT_EQ(fp.Mul(a, b), exp_mul) << vec.At("curve").AsString();
    EXPECT_EQ(fp.Inv(a), exp_inv) << vec.At("curve").AsString();
    EXPECT_EQ(fp.Mul(a, exp_inv), BigUInt<256>::One());

    const auto qr = BigUInt<256>::FromHex(vec.At("qr_hex").AsString());
    const auto exp_sqrt =
        BigUInt<256>::FromHex(vec.At("expected_sqrt_hex").AsString());
    const auto sqrt_opt = fp.Sqrt(qr);
    ASSERT_TRUE(sqrt_opt.has_value()) << vec.At("curve").AsString();
    EXPECT_EQ(*sqrt_opt, exp_sqrt);
    EXPECT_EQ(fp.Mul(*sqrt_opt, *sqrt_opt), qr);

    const auto qnr = BigUInt<256>::FromHex(vec.At("qnr_hex").AsString());
    EXPECT_FALSE(fp.HasSqrt(qnr));
    EXPECT_FALSE(fp.Sqrt(qnr).has_value());
  } else if (op == "gf2_binary_large_arith") {
    const auto& fld_json = vec.At("field");
    ASSERT_EQ(fld_json.At("m").AsUint64(), 163ULL);
    std::vector<std::size_t> exponents;
    const auto& mod_arr = fld_json.At("modulus").AsArray();
    for (std::size_t deg = 0; deg < mod_arr.size(); ++deg) {
      if ((mod_arr[deg].AsUint64() & 1ULL) != 0ULL) {
        exponents.push_back(deg);
      }
    }
    const GaloisFieldBinaryLarge<163> gf163(std::move(exponents));
    using Elem = GaloisFieldBinaryLarge<163>::element_type;

    const Elem a = Elem::FromHex(vec.At("a_hex").AsString());
    const Elem b = Elem::FromHex(vec.At("b_hex").AsString());
    const Elem exp_mul = Elem::FromHex(vec.At("expected_mul_hex").AsString());
    const Elem exp_inv = Elem::FromHex(vec.At("expected_inv_hex").AsString());
    EXPECT_EQ(gf163.Mul(a, b), exp_mul);
    EXPECT_EQ(gf163.Inv(a), exp_inv);
    EXPECT_EQ(gf163.Trace(a),
              static_cast<uint8_t>(vec.At("expected_trace_a").AsUint64()));

    const Elem tr0 = Elem::FromHex(vec.At("trace0_elem_hex").AsString());
    const Elem tr1 = Elem::FromHex(vec.At("trace1_elem_hex").AsString());
    const Elem exp_ht =
        Elem::FromHex(vec.At("expected_half_trace_hex").AsString());
    EXPECT_EQ(gf163.Trace(tr0), 0U);
    EXPECT_EQ(gf163.Trace(tr1), 1U);
    const Elem ht = gf163.HalfTrace(tr0);
    EXPECT_EQ(ht, exp_ht);
    EXPECT_EQ(gf163.Add(gf163.Sqr(ht), ht), tr0);
  } else {
    FAIL() << "Unhandled required field operation: " << op;
  }
}

template <typename GaloisField>
void VerifyTamoBargLrcImpl(const std::shared_ptr<GaloisField>& gf,
                           const JsonValue& vec) {
  const std::string& id = vec.At("id").AsString();
  const std::size_t n = static_cast<std::size_t>(vec.At("n").AsUint64());
  const std::size_t k = static_cast<std::size_t>(vec.At("k").AsUint64());
  const std::size_t r =
      static_cast<std::size_t>(vec.At("locality_r").AsUint64());
  const std::size_t delta =
      vec.Contains("local_distance")
          ? static_cast<std::size_t>(vec.At("local_distance").AsUint64())
          : 2U;
  const std::string& pt_str = vec.At("partition_type").AsString();
  const coding::LRCPartitionType ptype =
      (pt_str == "additive") ? coding::LRCPartitionType::kAdditive
                             : coding::LRCPartitionType::kMultiplicative;

  const auto code = coding::LocallyRecoverableCode<GaloisField>::CreateTamoBarg(
      gf, n, k, r, delta, ptype);

  EXPECT_EQ(code.MinimumDistance(), vec.At("expected_min_distance").AsUint64())
      << id;
  EXPECT_EQ(code.GoodPolynomial().Degree(),
            vec.At("good_poly_degree").AsInt64())
      << id;
  EXPECT_EQ(code.EvaluationPoints(),
            MakeGFArrayFromJson(gf, vec.At("eval_points")))
      << id;

  const auto msg = MakeGFArrayFromJson(gf, vec.At("message"));
  const auto expected_cw = MakeGFArrayFromJson(gf, vec.At("expected_codeword"));
  const auto encoded = code.Encode(msg);
  EXPECT_EQ(encoded, expected_cw) << id;
  EXPECT_TRUE(code.Contains(encoded)) << id;
  EXPECT_EQ(code.Unencode(encoded), msg) << id;

  // Single-erasure strict local repair
  const std::size_t erased_pos =
      static_cast<std::size_t>(vec.At("erased_position").AsUint64());
  auto corrupted = encoded;
  corrupted.Raw(erased_pos) =
      gf->Add(corrupted.Raw(erased_pos), gf->MultiplicativeIdentity());
  const auto repair = code.LocalRepair(corrupted, erased_pos, 0);
  ASSERT_TRUE(repair.Ok()) << id << ": " << repair.error_message;
  EXPECT_EQ(repair.symbols_read, vec.At("expected_symbols_read").AsUint64())
      << id;
  EXPECT_EQ(repair.recovered_symbol.RawValue(), expected_cw.Raw(erased_pos))
      << id;

  // Multi-erasure local-only repair for (r, delta) with delta >= 3
  if (vec.Contains("multi_erased_positions")) {
    const auto multi_pos = ExtractSizeArray(vec.At("multi_erased_positions"));
    auto multi_corrupted = encoded;
    for (std::size_t pos : multi_pos) {
      multi_corrupted.Raw(pos) =
          gf->Add(multi_corrupted.Raw(pos), gf->MultiplicativeIdentity());
    }
    const auto mres = code.LocalRepairMultiple(multi_corrupted, multi_pos,
                                               /*allow_global_fallback=*/false);
    ASSERT_TRUE(mres.Ok()) << id << ": " << mres.error_message;
    EXPECT_FALSE(mres.used_global_fallback) << id;
    EXPECT_EQ(mres.total_symbols_read,
              vec.At("expected_symbols_read").AsUint64())
        << id;
    EXPECT_EQ(mres.codeword, expected_cw) << id;
    EXPECT_EQ(mres.message, msg) << id;
  }
}

template <typename GaloisField>
void VerifyLrcAvailabilityImpl(const std::shared_ptr<GaloisField>& gf,
                               const JsonValue& vec) {
  const std::string& id = vec.At("id").AsString();
  const std::size_t n = static_cast<std::size_t>(vec.At("n").AsUint64());
  const std::size_t k = static_cast<std::size_t>(vec.At("k").AsUint64());
  const auto localities = ExtractSizeArray(vec.At("localities"));
  const auto local_distances = ExtractSizeArray(vec.At("local_distances"));

  const auto code = coding::LocallyRecoverableCode<
      GaloisField>::CreateMultiplicativeAvailability(gf, n, k, localities,
                                                     local_distances);

  EXPECT_EQ(code.Availability(), vec.At("expected_availability").AsUint64())
      << id;
  EXPECT_EQ(code.DesignedDistance(),
            vec.At("expected_designed_distance").AsUint64())
      << id;
  EXPECT_EQ(code.EvaluationPoints(),
            MakeGFArrayFromJson(gf, vec.At("eval_points")))
      << id;

  const auto msg = MakeGFArrayFromJson(gf, vec.At("message"));
  const auto expected_cw = MakeGFArrayFromJson(gf, vec.At("expected_codeword"));
  const auto encoded = code.Encode(msg);
  EXPECT_EQ(encoded, expected_cw) << id;
  EXPECT_TRUE(code.Contains(encoded)) << id;
  EXPECT_EQ(code.Unencode(encoded), msg) << id;

  // Verify pairwise disjointness of repair sets at every coordinate
  const std::size_t avail = code.Availability();
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<bool> used(n, false);
    for (std::size_t u = 0; u < avail; ++u) {
      for (std::size_t idx : code.RepairSet(i, u)) {
        EXPECT_NE(idx, i) << id;
        EXPECT_FALSE(used[idx]) << id << " overlapping repair set at " << i;
        used[idx] = true;
      }
    }
  }

  // Verify independent local repair from EVERY disjoint repair set
  const std::size_t test_coord =
      static_cast<std::size_t>(vec.At("test_coord").AsUint64());
  for (std::size_t u = 0; u < avail; ++u) {
    const auto& rset = code.RepairSet(test_coord, u);
    std::vector<bool> in_rset(n, false);
    for (std::size_t idx : rset) {
      in_rset[idx] = true;
    }
    auto corrupted = encoded;
    for (std::size_t j = 0; j < n; ++j) {
      if (!in_rset[j]) {
        corrupted.Raw(j) =
            gf->Add(corrupted.Raw(j), gf->MultiplicativeIdentity());
      }
    }
    const auto res = code.LocalRepair(corrupted, test_coord, u);
    ASSERT_TRUE(res.Ok()) << id << " set " << u << ": " << res.error_message;
    EXPECT_EQ(res.symbols_read, localities[u]) << id;
    EXPECT_EQ(res.recovered_symbol.RawValue(), expected_cw.Raw(test_coord))
        << id;
  }
}

void VerifyCodingVector(const JsonValue& vec) {
  const std::string& op = vec.At("operation").AsString();

  if (op == "rs_encode_bm_decode") {
    const auto& fld_json = vec.At("field");
    const uint8_t m = static_cast<uint8_t>(fld_json.At("m").AsUint64());
    const std::string poly_str =
        BinaryModulusPolyString(fld_json.At("modulus"));
    auto gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m, "int",
                                                                     poly_str);

    const std::size_t k = static_cast<std::size_t>(vec.At("k").AsUint64());
    const auto eval_pts = MakeGFArrayFromJson(gf, vec.At("eval_points"));
    const auto rs = coding::GeneralizedReedSolomonCode<
        GaloisFieldBinaryExtension<uint32_t>>::StandardReedSolomon(gf, eval_pts,
                                                                   k);

    const auto msg = MakeGFArrayFromJson(gf, vec.At("message"));
    const auto expected_cw =
        MakeGFArrayFromJson(gf, vec.At("expected_codeword"));
    const auto received = MakeGFArrayFromJson(gf, vec.At("received"));

    const auto encoded = rs.Encode(msg);
    EXPECT_EQ(encoded, expected_cw);
    EXPECT_TRUE(rs.Contains(encoded));

    const auto decoded_cw = rs.DecodeToCode(received, "BerlekampMassey");
    EXPECT_EQ(decoded_cw, expected_cw);
    EXPECT_EQ(rs.Unencode(decoded_cw), msg);
  } else if (op == "grs_guruswami_sudan_list_decode") {
    const auto& fld_json = vec.At("field");
    const uint8_t m = static_cast<uint8_t>(fld_json.At("m").AsUint64());
    const std::string poly_str =
        BinaryModulusPolyString(fld_json.At("modulus"));
    auto gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m, "int",
                                                                     poly_str);

    const std::size_t k = static_cast<std::size_t>(vec.At("k").AsUint64());
    const auto eval_pts = MakeGFArrayFromJson(gf, vec.At("eval_points"));
    const auto col_mults =
        MakeGFArrayFromJson(gf, vec.At("column_multipliers"));
    const coding::GeneralizedReedSolomonCode<
        GaloisFieldBinaryExtension<uint32_t>>
        grs(gf, eval_pts, col_mults, k);

    const auto msg = MakeGFArrayFromJson(gf, vec.At("message"));
    const auto expected_cw =
        MakeGFArrayFromJson(gf, vec.At("expected_codeword"));
    const auto received = MakeGFArrayFromJson(gf, vec.At("received"));

    EXPECT_EQ(grs.Encode(msg), expected_cw);

    const std::size_t target_errors =
        static_cast<std::size_t>(vec.At("target_errors").AsUint64());
    const auto gs_decoder = coding::GuruswamiSudanDecoder<
        GaloisFieldBinaryExtension<uint32_t>>::FromTargetRadius(grs,
                                                                target_errors);
    EXPECT_EQ(gs_decoder.JohnsonRadius(),
              vec.At("expected_johnson_radius").AsUint64());
    EXPECT_GE(gs_decoder.DecodingRadius(), target_errors);

    const auto cw_list = gs_decoder.DecodeToList(received);
    ASSERT_FALSE(cw_list.empty());
    EXPECT_NE(std::find(cw_list.begin(), cw_list.end(), expected_cw),
              cw_list.end());
    for (const auto& cand : cw_list) {
      EXPECT_TRUE(grs.Contains(cand));
      EXPECT_LE(grs.Distance(cand, received), gs_decoder.DecodingRadius());
    }
  } else if (op == "tamo_barg_lrc") {
    const auto& fld_json = vec.At("field");
    const uint64_t p = fld_json.At("p").AsUint64();
    const uint64_t m = fld_json.At("m").AsUint64();
    if (m == 1ULL) {
      auto gf = std::make_shared<GaloisFieldPrime<uint32_t>>(
          static_cast<uint32_t>(p));
      VerifyTamoBargLrcImpl(gf, vec);
    } else {
      const std::string poly_str =
          BinaryModulusPolyString(fld_json.At("modulus"));
      auto gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(
          static_cast<uint8_t>(m), "int", poly_str);
      VerifyTamoBargLrcImpl(gf, vec);
    }
  } else if (op == "lrc_availability") {
    const auto& fld_json = vec.At("field");
    const uint64_t p = fld_json.At("p").AsUint64();
    const uint64_t m = fld_json.At("m").AsUint64();
    if (m == 1ULL) {
      auto gf = std::make_shared<GaloisFieldPrime<uint32_t>>(
          static_cast<uint32_t>(p));
      VerifyLrcAvailabilityImpl(gf, vec);
    } else {
      const std::string poly_str =
          BinaryModulusPolyString(fld_json.At("modulus"));
      auto gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(
          static_cast<uint8_t>(m), "int", poly_str);
      VerifyLrcAvailabilityImpl(gf, vec);
    }
  } else if (op == "interleaved_rs_collaborative_decode") {
    const auto& fld_json = vec.At("field");
    const uint8_t m = static_cast<uint8_t>(fld_json.At("m").AsUint64());
    const std::string poly_str =
        BinaryModulusPolyString(fld_json.At("modulus"));
    auto gf = std::make_shared<GaloisFieldBinaryExtension<uint32_t>>(m, "int",
                                                                     poly_str);

    const std::size_t k = static_cast<std::size_t>(vec.At("k").AsUint64());
    const std::size_t ell =
        static_cast<std::size_t>(vec.At("interleaving_order").AsUint64());
    const auto eval_pts = MakeGFArrayFromJson(gf, vec.At("eval_points"));

    const auto irs = coding::InterleavedGRSCode<
        GaloisFieldBinaryExtension<uint32_t>>::CreateHomogeneous(gf, eval_pts,
                                                                 k, ell);
    EXPECT_EQ(irs.UniqueDecodingRadius(), vec.At("unique_radius").AsUint64());
    EXPECT_EQ(irs.CollaborativeDecodingRadius(),
              vec.At("collaborative_radius").AsUint64());

    std::vector<GFArray<GaloisFieldBinaryExtension<uint32_t>>> msgs;
    for (const auto& m_arr : vec.At("messages").AsArray()) {
      msgs.push_back(MakeGFArrayFromJson(gf, m_arr));
    }
    const auto expected_cw =
        MakeGFMatrixFromJson(gf, vec.At("expected_codeword_matrix"));
    const auto received = MakeGFMatrixFromJson(gf, vec.At("received_matrix"));
    const auto expected_burst =
        ExtractSizeArray(vec.At("burst_error_positions"));

    const auto encoded = irs.Encode(msgs);
    EXPECT_EQ(encoded, expected_cw);
    EXPECT_TRUE(irs.Contains(encoded));

    for (const auto mode :
         {coding::IRSDecoderMode::kSsb, coding::IRSDecoderMode::kSsbExact,
          coding::IRSDecoderMode::kBkyHankel}) {
      const auto res =
          irs.Decode(received, mode, /*fallback_to_independent=*/false);
      ASSERT_TRUE(res.Ok()) << res.error_message;
      EXPECT_FALSE(res.used_fallback);
      EXPECT_EQ(res.codeword_matrix, expected_cw);
      EXPECT_EQ(res.burst_error_positions, expected_burst);
      EXPECT_EQ(res.messages, msgs);
    }
  } else {
    FAIL() << "Unhandled required coding operation: " << op;
  }
}

// ============================================================================
// Tests for tiny_json.hpp
// ============================================================================

TEST(TinyJsonTest, RoundTripPrimitivesContainersAndEscapes) {
  const std::string raw = R"({
    "z_null": null,
    "a_bool_t": true,
    "b_bool_f": false,
    "c_int_neg": -9223372036854775808,
    "d_int_pos": 9223372036854775807,
    "e_sci": -1.25e-6,
    "f_str": "hello \"world\" \\ \n\t\r\b\f \u0041 \u03b1 \uD83D\uDE00",
    "g_arr": [1, "two", [3, 4], {"nested": true}]
  })";

  const JsonValue parsed = ParseJson(raw);
  ASSERT_TRUE(parsed.IsObject());
  EXPECT_EQ(parsed.Size(), 8U);
  EXPECT_TRUE(parsed.Contains("z_null"));
  EXPECT_FALSE(parsed.Contains("missing_key"));

  EXPECT_TRUE(parsed.At("z_null").IsNull());
  EXPECT_TRUE(parsed.At("a_bool_t").AsBool());
  EXPECT_FALSE(parsed.At("b_bool_f").AsBool());
  EXPECT_EQ(parsed.At("c_int_neg").AsInt64(),
            std::numeric_limits<int64_t>::min());
  EXPECT_EQ(parsed.At("d_int_pos").AsInt64(),
            std::numeric_limits<int64_t>::max());
  EXPECT_EQ(parsed.At("d_int_pos").AsUint64(),
            static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
  EXPECT_NEAR(parsed.At("e_sci").AsDouble(), -1.25e-6, 1e-18);

  const std::string& s = parsed.At("f_str").AsString();
  EXPECT_NE(s.find("\"world\""), std::string::npos);
  EXPECT_NE(s.find("\n\t\r\b\f"), std::string::npos);
  EXPECT_NE(s.find('A'), std::string::npos);
  // UTF-8 encoding of U+03B1 (Greek alpha): 0xCE 0xB1
  EXPECT_NE(s.find("\xCE\xB1"), std::string::npos);
  // UTF-8 encoding of U+1F600 (grinning face): 0xF0 0x9F 0x98 0x80
  EXPECT_NE(s.find("\xF0\x9F\x98\x80"), std::string::npos);

  EXPECT_EQ(parsed.At("g_arr").Size(), 4U);
  EXPECT_TRUE(parsed.At("g_arr").At(3).At("nested").AsBool());

  // Deterministic sorted-key serialization & round-trip equality
  const std::string compact = DumpJson(parsed, 0);
  const std::string pretty = DumpJson(parsed, 2);
  EXPECT_LT(compact.find("\"a_bool_t\""), compact.find("\"b_bool_f\""));
  EXPECT_LT(compact.find("\"g_arr\""), compact.find("\"z_null\""));
  EXPECT_EQ(ParseJson(compact), parsed);
  EXPECT_EQ(ParseJson(pretty), parsed);
}

TEST(TinyJsonTest, CheckedAccessorsAndParseErrors) {
  const JsonValue neg = ParseJson("-5");
  EXPECT_TRUE(neg.IsInt());
  EXPECT_TRUE(neg.IsNumber());
  EXPECT_THROW((void)neg.AsUint64(), std::runtime_error);
  EXPECT_THROW((void)neg.AsBool(), std::runtime_error);
  EXPECT_THROW((void)neg.AsString(), std::runtime_error);
  EXPECT_THROW((void)neg.Size(), std::runtime_error);

  const JsonValue obj = ParseJson(R"({"k": 1})");
  EXPECT_THROW((void)obj.At("nope"), std::runtime_error);

  const char* invalid_inputs[] = {
      "",
      "{",
      "[1, 2",
      "{\"a\": 1} trailing",
      "01",
      "\"unterminated",
      "\"bad\\xescape\"",
      "\"\\u123\"",
      "\"\\uD800\"",
      "tru",
      "[1,]",
  };
  for (const char* bad : invalid_inputs) {
    try {
      (void)ParseJson(bad);
      FAIL() << "Expected parse error for input: " << bad;
    } catch (const std::runtime_error& err) {
      EXPECT_NE(std::string(err.what()).find("offset"), std::string::npos)
          << "Error missing offset for input: " << bad;
    }
  }

  EXPECT_THROW((void)ParseJsonFile("/nonexistent/path/to/golden.json"),
               std::runtime_error);
}

// ============================================================================
// Tier-1 SageMath / Reference Golden Oracle Suite
// ============================================================================

TEST(SageGoldenOracleTest, ExecuteAllRequiredAndSkipPlannedVectors) {
  const std::string golden_path = ResolveGoldenPath(kGoldenFilename);
  const JsonValue root = ParseJsonFile(golden_path);

  ASSERT_TRUE(root.IsObject());
  EXPECT_EQ(root.At("schema_version").AsString(), "1.0");
  EXPECT_EQ(root.At("seed").AsInt64(), 2026);

  const auto& summary = root.At("summary");
  const auto& vectors = root.At("vectors").AsArray();
  EXPECT_EQ(vectors.size(), summary.At("total_vectors").AsUint64());

  std::size_t required_executed = 0;
  std::size_t planned_skipped = 0;

  for (const auto& vec : vectors) {
    const std::string& id = vec.At("id").AsString();
    const std::string& status = vec.At("status").AsString();
    const std::string& category = vec.At("category").AsString();

    if (status == "planned") {
      EXPECT_TRUE(vec.Contains("target_version")) << id;
      ++planned_skipped;
      continue;
    }

    ASSERT_EQ(status, "required") << "Unexpected status for vector " << id;
    if (category == "groups") {
      VerifyGroupVector(vec);
    } else if (category == "rings") {
      VerifyRingVector(vec);
    } else if (category == "fields") {
      VerifyFieldVector(vec);
    } else if (category == "coding") {
      VerifyCodingVector(vec);
    } else {
      FAIL() << "Unknown category '" << category << "' for vector " << id;
    }
    ++required_executed;
  }

  EXPECT_GE(required_executed, 27U);
  EXPECT_EQ(required_executed, summary.At("required_vectors").AsUint64());
  EXPECT_EQ(planned_skipped, summary.At("planned_vectors").AsUint64());
  EXPECT_EQ(planned_skipped, 1U);
}

}  // namespace
}  // namespace oracle
}  // namespace xg
