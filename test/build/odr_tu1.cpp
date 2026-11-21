#include <gtest/gtest.h>

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Include EVERY public header in xgalois/ to verify zero ODR violations
// across multiple translation units.
#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
#include "xgalois/databases/embedded_conway_polys.hpp"
#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_base.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_concepts.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_factory.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/linalg/linalg.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/field.hpp"
#include "xgalois/utils/math.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace test_odr {

// Defined in odr_tu2.cpp.
bool OdrCheckTU2();

bool OdrCheckTU1() {
  static_assert(FiniteField<GF2>);
  static_assert(FiniteField<GFP<uint8_t>>);
  static_assert(FiniteField<GFPLOG<uint8_t>>);
  static_assert(FiniteField<GF2X<uint8_t>>);
  static_assert(FiniteField<GF2XLOG<uint8_t>>);
  static_assert(FiniteField<GF2XZECH<uint8_t>>);
  static_assert(FiniteField<GFPX<uint8_t>>);

  // utils/math.hpp
  if (!utils::IsPrime(65521) || utils::ModInv(3, 7) != 5 ||
      utils::ModPow(3, 6, 7) != 1) {
    return false;
  }
  const auto [eg, st] = utils::ExtendedGcd<int64_t>(21, 15);
  if (eg != 3 || 21 * st.first + 15 * st.second != 3) {
    return false;
  }

  // databases/embedded_conway_polys.hpp & databases/interface.hpp
  const auto embedded = databases::detail::LookupEmbeddedConwayPoly(2, 8);
  if (!embedded.has_value() || embedded->empty()) {
    return false;
  }
  databases::PrimeFactorsDatabase pf_db;
  const auto pf = pf_db.fetch(121);
  if (pf.factors.size() != 1 || pf.factors[0] != 11 ||
      pf.multiplicities[0] != 2) {
    return false;
  }
  const std::string conway_str = databases::GetConwayPolynomial(2, 4);
  if (conway_str.empty()) {
    return false;
  }

  // field/* & utils/field.hpp
  if (utils::ConvertRepresentation("poly") != FieldRepresentation::POLY) {
    return false;
  }
  auto gf2 = std::make_shared<GF2>();
  auto gf7 = std::make_shared<GFP<uint8_t>>(7);
  auto gf7_tbl = std::make_shared<GFPLOG<uint8_t>>(7);
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  auto gf16_log = std::make_shared<GF2XLOG<uint8_t>>(4);
  auto gf16_zech = std::make_shared<GF2XZECH<uint8_t>>(4);
  auto gf9 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL));

  const auto parsed_pow = utils::ParsePowerString("g^2", *gf7);
  if (parsed_pow != gf7->Pow(gf7->MultiplicativeGenerator(), 2)) {
    return false;
  }
  const auto variant_field = GF(7);
  (void)FetchElement(variant_field, 3);

  GaloisFieldElement<GFP<uint8_t>> e3(3, gf7);
  GaloisFieldElement<GFP<uint8_t>> e5(5, gf7);
  if ((e3 * e5).Value() != 1) {
    return false;
  }
  if (gf7_tbl->Mul(3, 5) != 1 || gf16->Mul(2, gf16->Inv(2)) != 1 ||
      !gf16_zech->IsOne(
          gf16_zech->Mul(gf16_zech->SetElementValue(3),
                         gf16_zech->Inv(gf16_zech->SetElementValue(3)))) ||
      !gf9->IsOne(gf9->Mul(gf9->FromInt(4), gf9->Inv(gf9->FromInt(4))))) {
    return false;
  }

  // poly/poly_dense.hpp & utils/poly.hpp
  const auto p = utils::ParsePolynomial(gf7, "x^2 + 1");
  if (!utils::IsIrreducible(p)) {
    return false;
  }
  const auto [lc, factors] = utils::Factor(p * p);
  if (!lc.IsOne() || factors.size() != 1 || factors[0].second != 2) {
    return false;
  }

  // linalg/linalg.hpp
  const auto eye3 = linalg::eye(gf7, 3);
  if (!linalg::det(eye3).IsOne() || linalg::inv(eye3) != eye3) {
    return false;
  }
  const auto cp = linalg::characteristic_poly(eye3);
  if (cp.Degree() != 3) {
    return false;
  }

  // coding/coding.hpp & channel/channel.hpp
  const auto ham = coding::HammingCode(gf2, 3);
  if (!coding::SatisfiesHammingBound(ham) ||
      coding::SingletonBound(ham.Length(), ham.Dimension()) != 4) {
    return false;
  }
  const auto rs = coding::CreateReedSolomonCode(gf16_log, 15, 11);
  const auto msg = GFArray<GF2XLOG<uint8_t>>::Ones(gf16_log, {11});
  const auto cw = rs.Encode(msg);
  channel::ErrorErasureChannel<GF2XLOG<uint8_t>> eec(gf16_log, 15, 1, 2, 123);
  const auto tx = eec.transmit_with_erasures(cw);
  const auto dec = rs.DecodeWithErasures(tx.received, tx.erasure_positions);
  if (!dec.Ok() || dec.codeword != cw) {
    return false;
  }

  return true;
}

}  // namespace test_odr
}  // namespace xg

TEST(OdrMultiTuTest, BothTranslationUnitsLinkAndExecuteCleanly) {
  EXPECT_TRUE(xg::test_odr::OdrCheckTU1());
  EXPECT_TRUE(xg::test_odr::OdrCheckTU2());
}
