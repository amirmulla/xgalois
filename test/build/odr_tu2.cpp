#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Include EVERY public header in xgalois/ in a second translation unit to
// detect any missing `inline` / ODR multiple-definition linker errors.
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

bool OdrCheckTU2() {
  static_assert(FiniteField<GF2>);
  static_assert(FiniteField<GFP<uint8_t>>);
  static_assert(FiniteField<GFPLOG<uint8_t>>);
  static_assert(FiniteField<GF2X<uint8_t>>);
  static_assert(FiniteField<GF2XLOG<uint8_t>>);
  static_assert(FiniteField<GF2XZECH<uint8_t>>);
  static_assert(FiniteField<GFPX<uint8_t>>);

  // utils/math.hpp
  if (!utils::IsPrime(104729) || utils::ModMul(5, 6, 7) != 2 ||
      utils::Gcd<uint64_t>(48, 18) != 6) {
    return false;
  }
  const auto factors = utils::PrimeFactorizeWithMultiplicity(823543);
  if (factors.size() != 1 || factors[0].first != 7 || factors[0].second != 7) {
    return false;
  }

  // databases/embedded_conway_polys.hpp & databases/interface.hpp
  const auto embedded = databases::detail::LookupEmbeddedConwayPoly(3, 2);
  if (!embedded.has_value() || embedded->empty()) {
    return false;
  }
  databases::ConwayPolyDatabase cdb;
  const auto cres = cdb.fetch(2, 4);
  const std::string pstr = databases::ConvertPolyResultToString(
      cres.nonzero_degrees, cres.nonzero_coeffs, 2);
  if (pstr.empty()) {
    return false;
  }

  // field/* & utils/field.hpp
  if (utils::ConvertRepresentation("hex") != FieldRepresentation::HEX) {
    return false;
  }
  auto gf2 = std::make_shared<GF2>();
  auto gf11 = std::make_shared<GFP<uint8_t>>(11);
  auto gf11_tbl = std::make_shared<GFPLOG<uint8_t>>(11);
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  auto gf16_log = std::make_shared<GF2XLOG<uint8_t>>(4);
  auto gf16_zech = std::make_shared<GF2XZECH<uint8_t>>(4);
  auto gf9 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL));

  const auto parsed = utils::ParsePowerString("g^3", *gf11);
  if (parsed != gf11->Pow(gf11->MultiplicativeGenerator(), 3)) {
    return false;
  }
  const auto vfield = GF({2, 4});
  (void)FetchElement(vfield, 5);

  GaloisFieldElement<GFP<uint8_t>> a(4, gf11);
  GaloisFieldElement<GFP<uint8_t>> b(3, gf11);
  if ((a * b).Value() != 1) {
    return false;
  }
  if (gf11_tbl->Mul(4, 3) != 1 || gf16->Add(5, 5) != 0 ||
      !gf16_zech->IsZero(gf16_zech->Add(gf16_zech->SetElementValue(7),
                                        gf16_zech->SetElementValue(7))) ||
      !gf9->IsZero(gf9->Sub(gf9->FromInt(5), gf9->FromInt(5)))) {
    return false;
  }

  // poly/poly_dense.hpp & utils/poly.hpp
  const auto p = utils::ParsePolynomial(gf2, "x^4 + x + 1");
  if (!utils::IsPrimitive(p)) {
    return false;
  }
  const auto g = utils::PolynomialGcd(p, p);
  if (g != p) {
    return false;
  }

  // linalg/linalg.hpp
  const auto I4 = linalg::eye(gf2, 4);
  const auto plu_res = linalg::plu(I4);
  if (linalg::dot(plu_res.P, I4) != linalg::dot(plu_res.L, plu_res.U)) {
    return false;
  }

  // coding/coding.hpp & channel/channel.hpp
  const auto golay = coding::BinaryGolayCode();
  const auto tgolay = coding::TernaryGolayCode();
  const auto rm = coding::ReedMullerCode(1, 3);
  if (golay.Length() != 23 || tgolay.Length() != 11 || rm.Length() != 8) {
    return false;
  }
  const auto bch =
      coding::BCHCode<GF2, GF2XLOG<uint8_t>>::FromErrorCorrectionCapability(
          gf2, gf16_log, 2);
  channel::BinarySymmetricChannel bsc(gf2, bch.Length(), 0.0, 456);
  const auto stats = channel::SimulateCodePerformance(
      bch, bsc, channel::SimulationConfig{5, 0, 456});
  if (stats.frame_errors != 0) {
    return false;
  }

  return true;
}

}  // namespace test_odr
}  // namespace xg
