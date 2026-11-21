#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Include EVERY public header in xgalois/ in a second translation unit to
// detect any missing `inline` / ODR multiple-definition linker errors
// (all 21 public headers).
#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
#include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
#include "xgalois/coding/ldpc.hpp"
#include "xgalois/coding/polar.hpp"
#include "xgalois/databases/embedded_conway_polys.hpp"
#include "xgalois/databases/interface.hpp"
#include "xgalois/field/gf128.hpp"
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
  static_assert(FiniteField<GF128>);

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
  const auto cres = cdb.Fetch(2, 4);
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
  auto gf128 = std::make_shared<GF128>();
  const UInt128 a128 = UInt128::FromHex("0x66e94bd4ef8a2c3b884cfa59ca342b2e");
  const UInt128 b128 = UInt128::FromHex("0x0388dace60b6a392f328c2b971b2fe78");
  if (gf128->GhashMul(a128, b128) != gf128->GhashMul(b128, a128)) {
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
  const auto I4 = linalg::Eye(gf2, 4);
  const auto plu_res = linalg::Plu(I4);
  if (linalg::Dot(plu_res.p, I4) != linalg::Dot(plu_res.l, plu_res.u)) {
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
  if (stats.frame_errors != 0 || stats.Fer() != 0.0) {
    return false;
  }
  auto serc =
      channel::StaticErrorRateChannel<GF2>::Fixed(gf2, bch.Length(), 1, 789ULL);
  const auto bch_rx =
      serc.Transmit(bch.Encode(GFArray<GF2>::Zeros(gf2, {bch.Dimension()})));
  if (!bch.DecodeWithResult(bch_rx).Ok()) {
    return false;
  }

  // M6 coding headers: ldpc.hpp, polar.hpp, guruswami_sudan_decoder.hpp
  const auto ldpc = coding::LDPCCode::FromQuasiCyclic(gf2, {{0, 1, 2}}, 4);
  coding::BeliefPropagationDecoder bp(ldpc);
  const auto ldpc_cw =
      ldpc.Encode(GFArray<GF2>::Zeros(gf2, {ldpc.Dimension()}));
  if (!bp.DecodeWithResult(ldpc_cw).Ok()) {
    return false;
  }
  const auto polar = coding::PolarCode::FromBhattacharyya(gf2, 8, 4, 0.5);
  coding::SuccessiveCancellationDecoder sc(polar);
  const auto polar_cw = polar.Encode(GFArray<GF2>::Zeros(gf2, {4}));
  if (!sc.DecodeWithResult(polar_cw).Ok()) {
    return false;
  }
  const auto rs_gs = coding::CreateReedSolomonCode(gf16_log, 15, 5);
  const auto gs_dec =
      coding::GuruswamiSudanDecoder<GF2XLOG<uint8_t>>::FromTargetRadius(rs_gs,
                                                                        6);
  const auto gs_cw =
      rs_gs.Encode(GFArray<GF2XLOG<uint8_t>>::Zeros(gf16_log, {5}));
  if (!gs_dec.DecodeWithResult(gs_cw).Ok()) {
    return false;
  }

  return true;
}

}  // namespace test_odr
}  // namespace xg
