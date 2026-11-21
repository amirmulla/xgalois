#include <gtest/gtest.h>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/poly/poly_dense.hpp"
#include "xgalois/utils/poly.hpp"

namespace xg {
namespace {

// Generates deterministic adversarial ASCII, UTF-8, malformed, and near-valid
// strings for parser fuzzing while avoiding multi-gigabyte polynomial degree
// allocations (medium valid-looking exponent runs are clamped, whereas >20
// digit overflow runs are preserved to exercise overflow rejection).
std::vector<std::string> GenerateFuzzStrings(size_t count, uint64_t seed) {
  std::vector<std::string> corpus = {
      "",
      " ",
      "\t\n\r",
      "+",
      "-",
      "*",
      "^",
      "()",
      "(",
      ")",
      "((x))",
      "x^",
      "x^-1",
      "x^+2",
      "x^0",
      "x^1",
      "x^2 + x + 1",
      "-x^3 + 2*x - 5",
      "(1+α)*x^2 + (α)*x + 1",
      "0x",
      "0xFF",
      "0xFFFFFFFFFFFFFFFFFFFFF",
      "0x1G",
      "g^",
      "g^-",
      "g^0",
      "g^-1",
      "g^99999999999999999999999",
      "x^99999999999999999999999",
      "α^2 + α + 1",
      "α^-3",
      "β^2 + 1",
      "++x",
      "--x",
      "x**2",
      "2*x^2+",
      "+*x",
      "(2+3)x^2",
      "()x",
      "0",
      "1",
      "-0",
      "999999999999999999999999999",
      "\xFF\xFE\x80\xC0",
      "x^\xC3\xA9",
      "α\x80+1",
      std::string(1, '\0') + "x+1",
      "x+1" + std::string(1, '\0'),
  };

  const std::vector<std::string> tokens = {
      "0",   "1",    "2",    "3",        "7",   "9",
      "x",   "X",    "g",    "α",        "β",   "+",
      "-",   "*",    "^",    "(",        ")",   "0x",
      "0xA", " ",    "\t",   "!",        "?",   "_",
      "abc", "\x80", "\xFF", "\xCE\xB1", "255", "9999999999999999999999",
  };

  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<size_t> mode_dist(0, 3);
  std::uniform_int_distribution<size_t> len_dist(1, 14);
  std::uniform_int_distribution<size_t> tok_dist(0, tokens.size() - 1);
  std::uniform_int_distribution<int> byte_dist(0, 255);

  while (corpus.size() < count) {
    const size_t mode = mode_dist(rng);
    std::string s;
    if (mode == 0) {
      // Structured-ish token soup.
      const size_t ntok = len_dist(rng);
      for (size_t i = 0; i < ntok; ++i) {
        const std::string &tok = tokens[tok_dist(rng)];
        // Prevent accidental 4..10 digit positive exponents right after '^'.
        if (!s.empty() && s.back() == '^' && !tok.empty() &&
            std::isdigit(static_cast<unsigned char>(tok[0])) != 0 &&
            tok.size() < 20) {
          s.push_back(tok[0]);
        } else {
          s += tok;
        }
      }
    } else if (mode == 1) {
      // Raw random bytes including invalid UTF-8 and control characters.
      const size_t nbytes = len_dist(rng);
      s.reserve(nbytes);
      for (size_t i = 0; i < nbytes; ++i) {
        char ch = static_cast<char>(byte_dist(rng));
        if (!s.empty() && s.back() == '^' &&
            std::isdigit(static_cast<unsigned char>(ch)) != 0) {
          ch = '3';
        }
        s.push_back(ch);
      }
    } else if (mode == 2) {
      // Near-valid polynomial / element strings with random mutations.
      const size_t deg = rng() % 5;
      for (size_t d = 0; d <= deg; ++d) {
        if (d > 0) {
          s += (rng() & 1U) ? " + " : " - ";
        }
        s += std::to_string(rng() % 17);
        if (rng() & 1U) {
          s += "*x^" + std::to_string(d);
        } else if (rng() & 1U) {
          s += "*α^" + std::to_string(d);
        }
      }
      if (!s.empty() && (rng() % 3 != 0)) {
        const size_t pos = rng() % s.size();
        s[pos] = static_cast<char>(33 + (rng() % 90));
      }
    } else {
      // Hex / power / parenthesized mutations.
      const char *prefixes[] = {"0x", "g^", "g^-", "(", "α^", "x^"};
      s = prefixes[rng() % 6];
      const size_t ntok = 1 + (rng() % 5);
      for (size_t i = 0; i < ntok; ++i) {
        s += tokens[tok_dist(rng)];
      }
    }

    // Clamp any digit run immediately following '^' of length 4..18 to 2 digits
    // so ParsePolynomial never attempts a huge legitimate vector allocation.
    for (size_t i = 0; i < s.size(); ++i) {
      if (s[i] == '^') {
        size_t j = i + 1;
        while (j < s.size() &&
               std::isdigit(static_cast<unsigned char>(s[j])) != 0) {
          ++j;
        }
        const size_t run = j - (i + 1);
        if (run >= 4 && run <= 18) {
          s.erase(i + 3, run - 2);
        }
      }
    }

    corpus.push_back(std::move(s));
  }

  return corpus;
}

template <typename FieldType>
PolynomialDense<FieldType> RandomPoly(const std::shared_ptr<FieldType> &field,
                                      size_t max_deg, std::mt19937_64 &rng) {
  using Raw = typename FieldType::element_type;
  std::uniform_int_distribution<size_t> deg_dist(0, max_deg);
  const size_t deg = deg_dist(rng);
  std::vector<Raw> coeffs(deg + 1);
  for (size_t i = 0; i <= deg; ++i) {
    coeffs[i] = field->Random(rng);
  }
  return PolynomialDense<FieldType>::FromRaw(std::move(coeffs), field);
}

// ============================================================================
// 1. Fuzz ParsePolynomial (1,000+ strings across multiple field families)
// ============================================================================

TEST(FuzzSmokeTest, FuzzParsePolynomial) {
  auto gf2 = std::make_shared<GF2>();
  auto gf7 = std::make_shared<GFP<uint8_t>>(7);
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
  auto gf9 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL));

  const auto strings = GenerateFuzzStrings(1000, 0xF00D0001ULL);
  ASSERT_EQ(strings.size(), 1000U);

  size_t parsed_ok = 0;
  size_t rejected_cleanly = 0;

  auto check_parse = [&](const auto &field, const std::string &s) {
    try {
      const auto poly = utils::ParsePolynomial(field, s);
      // Returned polynomial must always be validly trimmed.
      EXPECT_TRUE(poly.IsZero() || !field->IsZero(poly.LeadingRawCoeff()))
          << "Untrimmed polynomial parsed from: " << s;
      ++parsed_ok;
    } catch (const std::invalid_argument &) {
      ++rejected_cleanly;
    } catch (const std::out_of_range &) {
      ++rejected_cleanly;
    } catch (const std::exception &e) {
      FAIL() << "Unexpected exception type for input '" << s
             << "': " << e.what();
    }
  };

  for (size_t i = 0; i < strings.size(); ++i) {
    const std::string &s = strings[i];
    switch (i % 4) {
      case 0:
        check_parse(gf2, s);
        break;
      case 1:
        check_parse(gf7, s);
        break;
      case 2:
        check_parse(gf256, s);
        break;
      case 3:
        check_parse(gf9, s);
        break;
    }
  }

  EXPECT_GT(parsed_ok, 0U);
  EXPECT_GT(rejected_cleanly, 0U);
  EXPECT_EQ(parsed_ok + rejected_cleanly, strings.size());
}

// ============================================================================
// 2. Fuzz SetElementValue(std::string) across GF(7), GF2X(8), GF2XZECH(8),
//    GFPX(3, 2) (500 strings each)
// ============================================================================

TEST(FuzzSmokeTest, FuzzSetElementValueString) {
  auto gf7 = std::make_shared<GFP<uint8_t>>(7);
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
  auto gf256_zech = std::make_shared<GF2XZECH<uint8_t>>(8);
  auto gf9 = std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL));

  const auto strings = GenerateFuzzStrings(500, 0xF00D0002ULL);
  ASSERT_EQ(strings.size(), 500U);

  auto fuzz_field = [&](const auto &field) {
    size_t ok = 0;
    size_t err = 0;
    for (const std::string &s : strings) {
      try {
        const auto elem = field->SetElementValue(s);
        // Verify returned element is in canonical range and printable.
        (void)field->ToString(elem);
        ++ok;
      } catch (const std::invalid_argument &) {
        ++err;
      } catch (const std::out_of_range &) {
        ++err;
      } catch (const std::exception &e) {
        FAIL() << "Unexpected exception for SetElementValue('" << s
               << "'): " << e.what();
      }
    }
    EXPECT_GT(ok, 0U);
    EXPECT_GT(err, 0U);
    EXPECT_EQ(ok + err, strings.size());
  };

  fuzz_field(gf7);
  fuzz_field(gf256);
  fuzz_field(gf256_zech);
  fuzz_field(gf9);
}

// ============================================================================
// 3. Fuzz PolynomialDense::DivRem & utils::Factor (200 pairs across 4 fields)
// ============================================================================

template <typename FieldType>
void FuzzPolyDivRemAndFactor(const std::shared_ptr<FieldType> &field,
                             size_t num_pairs, uint64_t seed) {
  using Poly = PolynomialDense<FieldType>;
  std::mt19937_64 rng(seed);

  for (size_t i = 0; i < num_pairs; ++i) {
    const Poly a = RandomPoly(field, 8, rng);
    const Poly b = RandomPoly(field, 5, rng);

    if (b.IsZero()) {
      EXPECT_THROW((void)a.DivRem(b), std::invalid_argument);
    } else {
      const auto [q, r] = a.DivRem(b);
      EXPECT_EQ(q * b + r, a);
      EXPECT_LT(r.Degree(), b.Degree());
      EXPECT_TRUE(q.IsZero() || !field->IsZero(q.LeadingRawCoeff()));
      EXPECT_TRUE(r.IsZero() || !field->IsZero(r.LeadingRawCoeff()));
    }

    // Factor a moderate-degree nonzero polynomial (deg <= 5 for fast ASan).
    const Poly f = RandomPoly(field, 5, rng);
    if (f.IsZero()) {
      EXPECT_THROW((void)utils::Factor(f), std::invalid_argument);
    } else {
      const auto [lc, factors] = utils::Factor(f, seed + i);
      Poly reconstructed = Poly::Monomial(field, 0, lc.RawValue());
      for (const auto &[irred, mult] : factors) {
        EXPECT_TRUE(irred.IsMonic());
        EXPECT_TRUE(utils::IsIrreducible(irred));
        EXPECT_GT(mult, 0U);
        reconstructed *= irred.Pow(mult);
      }
      EXPECT_EQ(reconstructed, f);
    }
  }
}

TEST(FuzzSmokeTest, FuzzPolyDivRemAndFactor) {
  FuzzPolyDivRemAndFactor(std::make_shared<GF2>(), 50, 0xF00D0101ULL);
  FuzzPolyDivRemAndFactor(std::make_shared<GFP<uint8_t>>(5), 50, 0xF00D0102ULL);
  FuzzPolyDivRemAndFactor(std::make_shared<GF2X<uint8_t>>(4), 50,
                          0xF00D0103ULL);
  FuzzPolyDivRemAndFactor(
      std::make_shared<GFPX<uint8_t>>(std::make_pair(3ULL, 2ULL)), 50,
      0xF00D0104ULL);
}

// ============================================================================
// 4. Fuzz BerlekampMasseyDecoder, SugiyamaDecoder, and GRSDecoder (PGZ)
//    on 300 random words with random error weights 0..n and random erasures
// ============================================================================

template <typename FieldType>
void FuzzGRSDecodersOnCode(
    const coding::GeneralizedReedSolomonCode<FieldType> &code, size_t trials,
    uint64_t seed) {
  const auto field = code.Field();
  const size_t n = code.Length();
  const size_t k = code.Dimension();
  const size_t redundancy = n - k;

  coding::BerlekampMasseyDecoder<FieldType> bm_dec(code);
  coding::SugiyamaDecoder<FieldType> sug_dec(code);
  coding::GRSDecoder<FieldType> pgz_dec(code);

  std::mt19937_64 rng(seed);
  std::vector<size_t> all_pos(n);
  std::iota(all_pos.begin(), all_pos.end(), size_t{0});

  for (size_t trial = 0; trial < trials; ++trial) {
    const auto msg = GFArray<FieldType>::Random(field, {k}, rng);
    const auto cw = code.Encode(msg);

    // First half of trials bias heavily toward within-capacity (2e + s <= n -
    // k); second half span arbitrary error/erasure weights across 0..n.
    size_t num_erasures = 0;
    size_t num_errors = 0;
    if (trial % 2 == 0) {
      num_erasures = rng() % (redundancy + 1);
      const size_t max_e = (redundancy - num_erasures) / 2;
      num_errors = rng() % (max_e + 1);
    } else {
      num_erasures = rng() % (n + 1);
      num_errors = rng() % (n - num_erasures + 1);
    }

    std::shuffle(all_pos.begin(), all_pos.end(), rng);
    std::vector<size_t> erasures(
        all_pos.begin(),
        all_pos.begin() + static_cast<std::ptrdiff_t>(num_erasures));
    std::sort(erasures.begin(), erasures.end());

    auto received = cw;
    auto rx_span = received.RawSpan();
    const auto zero = field->AdditiveIdentity();

    for (size_t pos : erasures) {
      rx_span[pos] = zero;
    }
    for (size_t i = 0; i < num_errors; ++i) {
      const size_t pos = all_pos[num_erasures + i];
      auto delta = field->Random(rng);
      while (field->IsZero(delta)) {
        delta = field->Random(rng);
      }
      rx_span[pos] = field->Add(rx_span[pos], delta);
    }

    const bool guaranteed_correctable =
        (2 * num_errors + num_erasures <= redundancy);

    auto verify_decoder = [&](const auto &decoder, const char *name) {
      const auto res = decoder.DecodeWithErasures(received, erasures);
      if (guaranteed_correctable) {
        ASSERT_TRUE(res.Ok())
            << name << " failed within capacity (e=" << num_errors
            << ", s=" << num_erasures << ", r=" << redundancy
            << "): " << res.error_message;
        EXPECT_EQ(res.codeword, cw) << name << " wrong codeword";
        EXPECT_EQ(res.message, msg) << name << " wrong message";
      } else {
        // Beyond capacity: never crash, and if it reports Ok(), the returned
        // codeword MUST be a valid codeword whose unencoding matches message!
        if (res.Ok()) {
          EXPECT_TRUE(code.Contains(res.codeword))
              << name << " returned Ok() with non-codeword";
          EXPECT_EQ(code.Unencode(res.codeword), res.message)
              << name << " returned inconsistent message/codeword";
        }
      }
    };

    verify_decoder(bm_dec, "BerlekampMassey");
    verify_decoder(sug_dec, "Sugiyama");
    verify_decoder(pgz_dec, "PGZ");
  }
}

TEST(FuzzSmokeTest, FuzzAllGRSDecoders300Words) {
  auto gf16 = std::make_shared<GF2XLOG<uint8_t>>(4);
  auto rs = coding::CreateReedSolomonCode(gf16, 15, 9);
  auto ext_rs = coding::CreateExtendedReedSolomonCode(gf16, 16, 10);

  // 200 words on standard RS(15, 9) + 100 words on extended RS(16, 10) with
  // alpha_0 = 0 => 300 fuzzed words across all 3 decoders.
  FuzzGRSDecodersOnCode(rs, 200, 0xF00D0201ULL);
  FuzzGRSDecodersOnCode(ext_rs, 100, 0xF00D0202ULL);
}

// ============================================================================
// 5. Fuzz ErrorErasureChannel (500 random valid configs)
// ============================================================================

TEST(FuzzSmokeTest, FuzzErrorErasureChannel500Configs) {
  auto gf32 = std::make_shared<GF2XLOG<uint8_t>>(5);
  std::mt19937_64 rng(0xF00D0301ULL);

  for (size_t trial = 0; trial < 500; ++trial) {
    const size_t n = 4 + (rng() % 61);  // block length in [4, 64]

    if (trial % 2 == 0) {
      // Exact (num_errors, num_erasures) configuration.
      const size_t num_erasures = rng() % (n + 1);
      const size_t num_errors = rng() % (n - num_erasures + 1);

      channel::ErrorErasureChannel<GF2XLOG<uint8_t>> ch(
          gf32, n, num_errors, num_erasures, 0xA00000ULL + trial);
      const auto cw = GFArray<GF2XLOG<uint8_t>>::Random(gf32, {n}, rng);
      const auto tx = ch.transmit_with_erasures(cw);

      ASSERT_EQ(tx.received.Size(), n);
      ASSERT_EQ(tx.erasure_mask.size(), n);
      ASSERT_EQ(tx.erasure_positions.size(), num_erasures);

      std::vector<bool> seen_erased(n, false);
      for (size_t i = 0; i < tx.erasure_positions.size(); ++i) {
        const size_t pos = tx.erasure_positions[i];
        ASSERT_LT(pos, n);
        if (i > 0) {
          EXPECT_LT(tx.erasure_positions[i - 1], pos);
        }
        seen_erased[pos] = true;
        EXPECT_TRUE(tx.erasure_mask[pos]);
        EXPECT_TRUE(gf32->IsZero(tx.received.Raw(pos)));
      }

      size_t non_erased_errors = 0;
      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(tx.erasure_mask[i], seen_erased[i]);
        if (!seen_erased[i] && tx.received.Raw(i) != cw.Raw(i)) {
          ++non_erased_errors;
        }
      }
      // Exact error count outside erasures guarantees error/erasure
      // disjointness!
      EXPECT_EQ(non_erased_errors, num_errors);
    } else {
      // Range ([min_e, max_e], [min_s, max_s]) configuration with max_e +
      // max_s <= n.
      const size_t max_erasures = rng() % (n + 1);
      const size_t min_erasures = rng() % (max_erasures + 1);
      const size_t max_errors = rng() % (n - max_erasures + 1);
      const size_t min_errors = rng() % (max_errors + 1);

      channel::ErrorErasureChannel<GF2XLOG<uint8_t>> ch(
          gf32, n, min_errors, max_errors, min_erasures, max_erasures,
          0xB00000ULL + trial);
      const auto cw = GFArray<GF2XLOG<uint8_t>>::Random(gf32, {n}, rng);
      const auto tx = ch.transmit_with_erasures(cw);

      ASSERT_EQ(tx.received.Size(), n);
      ASSERT_EQ(tx.erasure_mask.size(), n);
      EXPECT_GE(tx.erasure_positions.size(), min_erasures);
      EXPECT_LE(tx.erasure_positions.size(), max_erasures);

      std::vector<bool> seen_erased(n, false);
      for (size_t i = 0; i < tx.erasure_positions.size(); ++i) {
        const size_t pos = tx.erasure_positions[i];
        ASSERT_LT(pos, n);
        if (i > 0) {
          EXPECT_LT(tx.erasure_positions[i - 1], pos);
        }
        seen_erased[pos] = true;
        EXPECT_TRUE(tx.erasure_mask[pos]);
        EXPECT_TRUE(gf32->IsZero(tx.received.Raw(pos)));
      }

      size_t non_erased_errors = 0;
      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(tx.erasure_mask[i], seen_erased[i]);
        if (!seen_erased[i] && tx.received.Raw(i) != cw.Raw(i)) {
          ++non_erased_errors;
        }
      }
      EXPECT_GE(non_erased_errors, min_errors);
      EXPECT_LE(non_erased_errors, max_errors);
    }
  }
}

}  // namespace
}  // namespace xg
