#include "xgalois/channel/channel.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

namespace {

template <typename Field>
size_t CountArrayDiffs(const GFArray<Field> &a, const GFArray<Field> &b) {
  EXPECT_EQ(a.Size(), b.Size());
  const auto sa = a.RawSpan();
  const auto sb = b.RawSpan();
  size_t diffs = 0;
  for (size_t i = 0; i < sa.size(); ++i) {
    if (sa[i] != sb[i]) {
      ++diffs;
    }
  }
  return diffs;
}

template <typename Field>
struct MockDecodingResult {
  GFArray<Field> codeword;
  GFArray<Field> message;
  bool success = false;
  bool Ok() const noexcept { return success; }
};

// Minimal duck-typed code modeling a (3, 1) repetition code over GF(2),
// exercising the SimulateCodePerformance(code, channel) contract.
struct MockRepetition3Code {
  using field_type = GaloisFieldBinary;
  using codeword_type = GFArray<GaloisFieldBinary>;
  using message_type = GFArray<GaloisFieldBinary>;

  explicit MockRepetition3Code(std::shared_ptr<GaloisFieldBinary> f)
      : field(std::move(f)) {}

  std::shared_ptr<GaloisFieldBinary> Field() const { return field; }
  size_t Dimension() const { return 1; }

  codeword_type Encode(const message_type &msg) const {
    auto cw = codeword_type::Zeros(field, {3});
    const uint8_t bit = msg.Raw(0);
    for (size_t i = 0; i < 3; ++i) {
      cw.Raw(i) = bit;
    }
    return cw;
  }

  message_type Unencode(const codeword_type &cw) const {
    auto msg = message_type::Zeros(field, {1});
    msg.Raw(0) = cw.Raw(0);
    return msg;
  }

  MockDecodingResult<GaloisFieldBinary> DecodeWithResult(
      const codeword_type &rx) const {
    size_t ones = 0;
    for (size_t i = 0; i < 3; ++i) {
      if (rx.Raw(i) != 0) {
        ++ones;
      }
    }
    const uint8_t dec_bit = (ones >= 2) ? uint8_t{1} : uint8_t{0};
    auto cw = codeword_type::Zeros(field, {3});
    for (size_t i = 0; i < 3; ++i) {
      cw.Raw(i) = dec_bit;
    }
    auto msg = message_type::Zeros(field, {1});
    msg.Raw(0) = dec_bit;
    return MockDecodingResult<GaloisFieldBinary>{std::move(cw), std::move(msg),
                                                 true};
  }

  MockDecodingResult<GaloisFieldBinary> DecodeWithErasures(
      const codeword_type &rx,
      const std::vector<size_t> &erasure_positions) const {
    if (erasure_positions.size() >= 3) {
      return MockDecodingResult<GaloisFieldBinary>{};
    }
    std::vector<bool> erased(3, false);
    for (size_t pos : erasure_positions) {
      if (pos < 3) {
        erased[pos] = true;
      }
    }
    size_t ones = 0;
    size_t zeros = 0;
    for (size_t i = 0; i < 3; ++i) {
      if (!erased[i]) {
        if (rx.Raw(i) != 0) {
          ++ones;
        } else {
          ++zeros;
        }
      }
    }
    if (ones == zeros) {
      return MockDecodingResult<GaloisFieldBinary>{};
    }
    const uint8_t dec_bit = (ones > zeros) ? uint8_t{1} : uint8_t{0};
    auto cw = codeword_type::Zeros(field, {3});
    for (size_t i = 0; i < 3; ++i) {
      cw.Raw(i) = dec_bit;
    }
    auto msg = message_type::Zeros(field, {1});
    msg.Raw(0) = dec_bit;
    return MockDecodingResult<GaloisFieldBinary>{std::move(cw), std::move(msg),
                                                 true};
  }

  std::shared_ptr<GaloisFieldBinary> field;
};

}  // namespace

class ChannelTest : public ::testing::Test {
 protected:
  void SetUp() override {
    gf16 = std::make_shared<GF2X<uint8_t>>(4);
    gf17 = std::make_shared<GaloisFieldPrime<uint8_t>>(17);

    message_length = 8;
    message = xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>>::from_shape(
        {message_length});
    for (size_t i = 0; i < message_length; ++i) {
      message[i] = GaloisFieldElementBase<GF2X<uint8_t>>(i + 1, gf16);
    }

    prime_message =
        xt::xarray<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>::
            from_shape({message_length});
    for (size_t i = 0; i < message_length; ++i) {
      prime_message[i] = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>(
          (i * 2 + 1) % 17, gf17);
    }
  }

  std::shared_ptr<GF2X<uint8_t>> gf16;
  std::shared_ptr<GaloisFieldPrime<uint8_t>> gf17;
  xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>> message;
  xt::xarray<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> prime_message;
  size_t message_length;
};

// =============================================================================
// Existing Backwards-Compatible Legacy API Tests
// =============================================================================

TEST_F(ChannelTest, StaticErrorRateChannelFixed) {
  auto channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 2);

  EXPECT_EQ(channel->input_size(), message_length);
  EXPECT_EQ(channel->output_size(), message_length);
  EXPECT_EQ(channel->number_errors().first, 2);
  EXPECT_EQ(channel->number_errors().second, 2);

  auto transmitted = channel->transmit(message);
  EXPECT_EQ(transmitted.size(), message.size());

  size_t error_count = 0;
  for (size_t i = 0; i < message_length; ++i) {
    if (message[i] != transmitted[i]) error_count++;
  }
  EXPECT_EQ(error_count, 2);

  auto no_error_channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 0);
  auto no_error_transmitted = no_error_channel->transmit(message);
  for (size_t i = 0; i < message_length; ++i) {
    EXPECT_EQ(message[i], no_error_transmitted[i]);
  }
}

TEST_F(ChannelTest, StaticErrorRateChannelVariable) {
  auto channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 1, 3, 12345ULL);

  EXPECT_EQ(channel->number_errors().first, 1);
  EXPECT_EQ(channel->number_errors().second, 3);

  std::vector<size_t> error_counts;
  for (int i = 0; i < 100; ++i) {
    auto transmitted = channel->transmit(message);
    size_t errors = 0;
    for (size_t j = 0; j < message_length; ++j) {
      if (message[j] != transmitted[j]) errors++;
    }
    error_counts.push_back(errors);
  }

  for (size_t count : error_counts) {
    EXPECT_GE(count, 1);
    EXPECT_LE(count, 3);
  }

  auto min_errors = *std::min_element(error_counts.begin(), error_counts.end());
  auto max_errors = *std::max_element(error_counts.begin(), error_counts.end());
  EXPECT_LT(min_errors, max_errors);
}

TEST_F(ChannelTest, StaticErrorRateChannelInvalidInput) {
  auto channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 2);

  xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>> wrong_size_message =
      xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>>::from_shape({3});
  for (size_t i = 0; i < 3; ++i) {
    wrong_size_message[i] = GaloisFieldElementBase<GF2X<uint8_t>>(i + 1, gf16);
  }
  EXPECT_THROW(channel->transmit(wrong_size_message), std::invalid_argument);

  EXPECT_THROW(
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 3, 2),
      std::invalid_argument);

  // Constructor eagerly validates max_errors <= input_size_.
  EXPECT_THROW(
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, message_length + 1),
      std::invalid_argument);
  EXPECT_THROW(
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 1, message_length + 1),
      std::invalid_argument);
}

TEST_F(ChannelTest, ErrorErasureChannelFixed) {
  auto channel = std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(
      gf16, message_length, 2, 1);

  EXPECT_EQ(channel->number_errors().first, 2);
  EXPECT_EQ(channel->number_errors().second, 2);
  EXPECT_EQ(channel->number_erasures().first, 1);
  EXPECT_EQ(channel->number_erasures().second, 1);

  auto [transmitted, erasures] =
      channel->transmit_unsafe_with_erasures(message);
  EXPECT_EQ(transmitted.size(), message.size());
  EXPECT_EQ(erasures.size(), message.size());

  size_t error_count = 0;
  size_t erasure_count = 0;
  for (size_t i = 0; i < message_length; ++i) {
    if (erasures[i]) {
      erasure_count++;
      GaloisFieldElementBase<GF2X<uint8_t>> zero_element(
          gf16->AdditiveIdentity(), gf16);
      EXPECT_EQ(transmitted[i], zero_element);
    } else if (message[i] != transmitted[i]) {
      error_count++;
    }
  }
  EXPECT_EQ(error_count, 2);
  EXPECT_EQ(erasure_count, 1);

  auto regular_transmitted = channel->transmit(message);
  EXPECT_EQ(regular_transmitted.size(), message.size());
}

TEST_F(ChannelTest, ErrorErasureChannelVariable) {
  auto channel = std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(
      gf16, message_length, 1, 2, 0, 2, 777ULL);

  for (int i = 0; i < 50; ++i) {
    auto [transmitted, erasures] =
        channel->transmit_unsafe_with_erasures(message);

    size_t error_count = 0;
    size_t erasure_count = 0;
    for (size_t j = 0; j < message_length; ++j) {
      if (erasures[j]) {
        erasure_count++;
      } else if (message[j] != transmitted[j]) {
        error_count++;
      }
    }

    EXPECT_GE(error_count, 1);
    EXPECT_LE(error_count, 2);
    EXPECT_GE(erasure_count, 0);
    EXPECT_LE(erasure_count, 2);
  }
}

TEST_F(ChannelTest, QarySymmetricChannelBasic) {
  const double epsilon = 0.1;
  auto channel =
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, epsilon);

  EXPECT_EQ(channel->input_size(), message_length);
  EXPECT_EQ(channel->output_size(), message_length);
  EXPECT_DOUBLE_EQ(channel->error_probability(), epsilon);

  auto transmitted = channel->transmit(message);
  EXPECT_EQ(transmitted.size(), message.size());

  auto perfect_channel =
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, 0.0);
  auto perfect_transmitted = perfect_channel->transmit(message);
  for (size_t i = 0; i < message_length; ++i) {
    EXPECT_EQ(message[i], perfect_transmitted[i]);
  }
}

TEST_F(ChannelTest, QarySymmetricChannelProbabilities) {
  const double epsilon = 0.1;
  auto channel =
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, epsilon);

  double prob_0 = channel->probability_of_exactly_t_errors(0);
  double prob_1 = channel->probability_of_exactly_t_errors(1);
  double prob_all = channel->probability_of_exactly_t_errors(message_length);

  EXPECT_GT(prob_0, 0.0);
  EXPECT_GT(prob_1, 0.0);
  EXPECT_GT(prob_all, 0.0);

  double cum_prob_0 = channel->probability_of_at_most_t_errors(0);
  double cum_prob_1 = channel->probability_of_at_most_t_errors(1);

  EXPECT_DOUBLE_EQ(cum_prob_0, prob_0);
  EXPECT_DOUBLE_EQ(cum_prob_1, prob_0 + prob_1);

  double total_prob = 0.0;
  for (size_t t = 0; t <= message_length; ++t) {
    total_prob += channel->probability_of_exactly_t_errors(t);
  }
  EXPECT_NEAR(total_prob, 1.0, 1e-10);
}

TEST_F(ChannelTest, QarySymmetricChannelStatistics) {
  const double epsilon = 0.1;
  const int num_trials = 1000;
  auto channel =
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, epsilon, 42ULL);

  std::vector<size_t> error_counts(message_length + 1, 0);

  for (int trial = 0; trial < num_trials; ++trial) {
    auto transmitted = channel->transmit(message);
    size_t errors = 0;
    for (size_t i = 0; i < message_length; ++i) {
      if (message[i] != transmitted[i]) errors++;
    }
    error_counts[errors]++;
  }

  for (size_t t = 0; t <= message_length; ++t) {
    double empirical_prob = static_cast<double>(error_counts[t]) / num_trials;
    double theoretical_prob = channel->probability_of_exactly_t_errors(t);

    if (theoretical_prob > 0.01) {
      EXPECT_NEAR(empirical_prob, theoretical_prob, 0.05);
    }
  }
}

TEST_F(ChannelTest, QarySymmetricChannelInvalidInput) {
  EXPECT_THROW(std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
                   gf16, message_length, -0.1),
               std::invalid_argument);
  EXPECT_THROW(std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
                   gf16, message_length, 1.1),
               std::invalid_argument);
  EXPECT_THROW(
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, std::numeric_limits<double>::quiet_NaN()),
      std::invalid_argument);
  EXPECT_THROW(
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, std::numeric_limits<double>::infinity()),
      std::invalid_argument);

  EXPECT_NO_THROW(
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, 0.0));
  EXPECT_NO_THROW(
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, 1.0));
}

TEST_F(ChannelTest, ChannelWithPrimeField) {
  auto channel = std::make_shared<
      channels::StaticErrorRateChannel<GaloisFieldPrime<uint8_t>>>(
      gf17, message_length, 1);

  auto transmitted = channel->transmit(prime_message);
  EXPECT_EQ(transmitted.size(), prime_message.size());

  size_t error_count = 0;
  for (size_t i = 0; i < message_length; ++i) {
    if (prime_message[i] != transmitted[i]) error_count++;
  }
  EXPECT_EQ(error_count, 1);
}

TEST_F(ChannelTest, ChannelStringRepresentation) {
  auto static_channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 2);
  auto var_channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 1, 3);
  auto ee_channel =
      std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(
          gf16, message_length, 2, 1);
  auto qary_channel =
      std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(
          gf16, message_length, 0.1);

  std::string static_str = static_channel->to_string();
  std::string var_str = var_channel->to_string();
  std::string ee_str = ee_channel->to_string();
  std::string qary_str = qary_channel->to_string();

  EXPECT_FALSE(static_str.empty());
  EXPECT_FALSE(var_str.empty());
  EXPECT_FALSE(ee_str.empty());
  EXPECT_FALSE(qary_str.empty());

  EXPECT_NE(static_str.find("Static error rate"), std::string::npos);
  EXPECT_NE(var_str.find("between"), std::string::npos);
  EXPECT_NE(ee_str.find("Error-and-erasure"), std::string::npos);
  EXPECT_NE(qary_str.find("q-ary symmetric"), std::string::npos);
}

TEST_F(ChannelTest, ChannelOperatorCall) {
  auto channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 1);

  auto transmitted1 = channel->transmit(message);
  auto transmitted2 = (*channel)(message);

  size_t errors1 = 0, errors2 = 0;
  for (size_t i = 0; i < message_length; ++i) {
    if (message[i] != transmitted1[i]) errors1++;
    if (message[i] != transmitted2[i]) errors2++;
  }

  EXPECT_EQ(errors1, 1);
  EXPECT_EQ(errors2, 1);
}

TEST_F(ChannelTest, FactoryFunctions) {
  auto static_channel =
      channels::CreateStaticErrorRateChannel(gf16, message_length, 2);
  EXPECT_EQ(static_channel->number_errors().first, 2);
  EXPECT_EQ(static_channel->number_errors().second, 2);

  auto var_channel =
      channels::CreateStaticErrorRateChannel(gf16, message_length, 1, 3);
  EXPECT_EQ(var_channel->number_errors().first, 1);
  EXPECT_EQ(var_channel->number_errors().second, 3);

  auto ee_channel =
      channels::CreateErrorErasureChannel(gf16, message_length, 2, 1);
  EXPECT_EQ(ee_channel->number_errors().first, 2);
  EXPECT_EQ(ee_channel->number_errors().second, 2);
  EXPECT_EQ(ee_channel->number_erasures().first, 1);
  EXPECT_EQ(ee_channel->number_erasures().second, 1);

  auto var_ee_channel =
      channels::CreateErrorErasureChannel(gf16, message_length, 1, 2, 0, 1);
  EXPECT_EQ(var_ee_channel->number_errors().first, 1);
  EXPECT_EQ(var_ee_channel->number_errors().second, 2);
  EXPECT_EQ(var_ee_channel->number_erasures().first, 0);
  EXPECT_EQ(var_ee_channel->number_erasures().second, 1);

  auto qary_channel =
      channels::CreateQarySymmetricChannel(gf16, message_length, 0.1);
  EXPECT_DOUBLE_EQ(qary_channel->error_probability(), 0.1);
}

TEST_F(ChannelTest, ChannelFieldCompatibility) {
  auto channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, 1);

  EXPECT_EQ(channel->field(), gf16);

  auto transmitted = channel->transmit(message);
  EXPECT_EQ(transmitted.size(), message.size());

  for (size_t i = 0; i < transmitted.size(); ++i) {
    EXPECT_EQ(transmitted[i].Field(), gf16);
  }
}

TEST_F(ChannelTest, ErrorHandling) {
  EXPECT_THROW(std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(
                   gf16, message_length, 3, 2, 1, 1),
               std::invalid_argument);
  EXPECT_THROW(std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(
                   gf16, message_length, 1, 1, 3, 2),
               std::invalid_argument);

  auto channel =
      std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(
          gf16, message_length, message_length);
  EXPECT_NO_THROW(channel->transmit(message));
}

// =============================================================================
// Wave 3C: GFArray Native API Across GF(2), GF(7), GF(2^8), Zech, GF(3^2)
// =============================================================================

TEST(ChannelGFArrayTest, NativeTransmissionAcrossAllFieldFamilies) {
  // 1. GF(2)
  {
    auto gf2 = std::make_shared<GaloisFieldBinary>();
    std::mt19937_64 rng(101);
    auto msg = GFArray<GaloisFieldBinary>::Random(gf2, {16}, rng);
    channel::StaticErrorRateChannel<GaloisFieldBinary> ch(
        gf2, 16, 5, std::optional<uint64_t>{42ULL});
    auto rx = ch.transmit(msg);
    EXPECT_EQ(rx.Size(), 16U);
    EXPECT_EQ(rx.Field(), gf2);
    EXPECT_EQ(CountArrayDiffs(msg, rx), 5U);
  }

  // 2. GF(7)
  {
    auto gf7 = std::make_shared<GaloisFieldPrime<uint8_t>>(7);
    std::mt19937_64 rng(102);
    auto msg = GFArray<GaloisFieldPrime<uint8_t>>::Random(gf7, {20}, rng);
    channel::StaticErrorRateChannel<GaloisFieldPrime<uint8_t>> ch(
        gf7, 20, 6, std::optional<uint64_t>{42ULL});
    auto rx = ch.transmit(msg);
    EXPECT_EQ(CountArrayDiffs(msg, rx), 6U);
  }

  // 3. GF(2^8) standard binary extension
  {
    auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
    std::mt19937_64 rng(103);
    auto msg = GFArray<GF2X<uint8_t>>::Random(gf256, {32}, rng);
    channel::StaticErrorRateChannel<GF2X<uint8_t>> ch(
        gf256, 32, 7, std::optional<uint64_t>{42ULL});
    auto rx = ch.transmit(msg);
    EXPECT_EQ(CountArrayDiffs(msg, rx), 7U);
  }

  // 4. GF(2^8) Zech-log internal representation
  {
    auto gf256_zech = std::make_shared<GFBEZechLogTables<uint8_t>>(8);
    std::mt19937_64 rng(104);
    auto msg =
        GFArray<GFBEZechLogTables<uint8_t>>::Random(gf256_zech, {32}, rng);
    channel::StaticErrorRateChannel<GFBEZechLogTables<uint8_t>> ch(
        gf256_zech, 32, 8, std::optional<uint64_t>{42ULL});
    auto rx = ch.transmit(msg);
    EXPECT_EQ(CountArrayDiffs(msg, rx), 8U);
  }

  // 5. GF(3^2) general GaloisFieldExtension
  {
    auto gf9 =
        std::make_shared<GaloisFieldExtension<uint32_t>>(std::pair{3ULL, 2ULL});
    std::mt19937_64 rng(105);
    auto msg = GFArray<GaloisFieldExtension<uint32_t>>::Random(gf9, {24}, rng);
    channel::StaticErrorRateChannel<GaloisFieldExtension<uint32_t>> ch(
        gf9, 24, 9, std::optional<uint64_t>{42ULL});
    auto rx = ch.transmit(msg);
    EXPECT_EQ(CountArrayDiffs(msg, rx), 9U);

    // Also test QarySymmetricChannel with epsilon == 1.0 over GF(3^2).
    channel::QarySymmetricChannel<GaloisFieldExtension<uint32_t>> qsc(
        gf9, 24, 1.0, 42ULL);
    auto all_flipped = qsc.transmit(msg);
    EXPECT_EQ(CountArrayDiffs(msg, all_flipped), 24U);
  }
}

TEST(ChannelGFArrayTest, LargeDimensionFloydSampling) {
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
  constexpr size_t kLargeN = 50000;
  constexpr size_t kErrors = 5;
  auto msg = GFArray<GF2X<uint8_t>>::Zeros(gf256, {kLargeN});
  channel::StaticErrorRateChannel<GF2X<uint8_t>> ch(
      gf256, kLargeN, kErrors, std::optional<uint64_t>{99ULL});
  auto rx = ch.transmit(msg);
  EXPECT_EQ(CountArrayDiffs(msg, rx), kErrors);
}

// =============================================================================
// Bug B-14 Regression: ErrorErasureChannel Validation, Disjointness, Zero UB
// =============================================================================

TEST(ChannelBugFixTest, ErrorErasureConstructorBoundsAndOverflow) {
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  constexpr size_t kN = 10;

  // max_errors > n
  EXPECT_THROW(channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, kN + 1, 0),
               std::invalid_argument);
  EXPECT_THROW(
      channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, 0, kN + 1, 0, 0),
      std::invalid_argument);

  // max_erasures > n
  EXPECT_THROW(channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, 0, kN + 1),
               std::invalid_argument);
  EXPECT_THROW(
      channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, 0, 0, 0, kN + 1),
      std::invalid_argument);

  // max_errors + max_erasures > n
  EXPECT_THROW(channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, 6, 5),
               std::invalid_argument);
  EXPECT_THROW(
      channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, 1, 6, 1, 5),
      std::invalid_argument);

  // Overflow-safe check (size_t::max() + 1 would wrap around if added naively).
  constexpr size_t kMax = std::numeric_limits<size_t>::max();
  EXPECT_THROW(channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, kMax, 1),
               std::invalid_argument);

  // Exact boundary max_errors + max_erasures == n is valid.
  EXPECT_NO_THROW(channel::ErrorErasureChannel<GF2X<uint8_t>>(gf16, kN, 6, 4));
}

TEST(ChannelBugFixTest, ErrorErasureSafeTransmitAndDisjointness) {
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
  constexpr size_t kN = 16;
  channel::ErrorErasureChannel<GF2X<uint8_t>> ch(gf256, kN, 5, 4, 2025ULL);

  // Size mismatch throws on safe transmit_with_erasures.
  auto bad_size = GFArray<GF2X<uint8_t>>::Ones(gf256, {kN - 1});
  EXPECT_THROW(ch.transmit_with_erasures(bad_size), std::invalid_argument);

  // Field mismatch throws on safe transmit_with_erasures.
  auto other_gf256 = std::make_shared<GF2X<uint8_t>>(8);
  auto wrong_field = GFArray<GF2X<uint8_t>>::Ones(other_gf256, {kN});
  EXPECT_THROW(ch.transmit_with_erasures(wrong_field), std::invalid_argument);

  // Verify exact disjointness, sorted erasure_positions, mask consistency, and
  // zeroed erased positions across 50 trials.
  auto msg = GFArray<GF2X<uint8_t>>::Ones(gf256, {kN});
  for (int trial = 0; trial < 50; ++trial) {
    auto [received, erasure_positions, erasure_mask] =
        ch.transmit_with_erasures(msg);
    ASSERT_EQ(received.Size(), kN);
    ASSERT_EQ(erasure_mask.size(), kN);
    ASSERT_EQ(erasure_positions.size(), 4U);
    EXPECT_TRUE(
        std::is_sorted(erasure_positions.begin(), erasure_positions.end()));

    // Check positions have no duplicates and match mask exactly.
    for (size_t i = 0; i < erasure_positions.size(); ++i) {
      EXPECT_LT(erasure_positions[i], kN);
      EXPECT_TRUE(erasure_mask[erasure_positions[i]]);
      EXPECT_EQ(received.Raw(erasure_positions[i]), gf256->AdditiveIdentity());
      if (i > 0) {
        EXPECT_LT(erasure_positions[i - 1], erasure_positions[i]);
      }
    }

    size_t mask_erasure_count = 0;
    size_t non_erased_error_count = 0;
    for (size_t i = 0; i < kN; ++i) {
      if (erasure_mask[i]) {
        ++mask_erasure_count;
      } else if (received.Raw(i) != msg.Raw(i)) {
        ++non_erased_error_count;
      }
    }
    EXPECT_EQ(mask_erasure_count, 4U);
    EXPECT_EQ(non_erased_error_count, 5U);
  }

  // Exact full partition: num_errors + num_erasures == n.
  channel::ErrorErasureChannel<GF2X<uint8_t>> full_ch(gf256, kN, 10, 6, 7ULL);
  auto [rx, pos, mask] = full_ch.transmit_with_erasures(msg);
  EXPECT_EQ(pos.size(), 6U);
  size_t errs = 0;
  for (size_t i = 0; i < kN; ++i) {
    if (!mask[i]) {
      EXPECT_NE(rx.Raw(i), msg.Raw(i));
      ++errs;
    } else {
      EXPECT_EQ(rx.Raw(i), gf256->AdditiveIdentity());
    }
  }
  EXPECT_EQ(errs, 10U);
}

// =============================================================================
// Bug B-26 Regression: QarySymmetricChannel Numerical Stability & Edge Cases
// =============================================================================

TEST(ChannelBugFixTest, QarySymmetricChannelLogSpaceAndEpsilonOne) {
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);

  // Large n = 2000, t = 1000 must not overflow to inf*0 = NaN.
  channel::QarySymmetricChannel<GF2X<uint8_t>> large_qsc(gf256, 2000, 0.5);
  const double p_mid = large_qsc.probability_of_exactly_t_errors(1000);
  EXPECT_TRUE(std::isfinite(p_mid));
  EXPECT_GT(p_mid, 0.0);
  // By Stirling's approximation, C(2000, 1000) * 2^-2000 ~= 1 / sqrt(pi * 1000)
  // ~= 0.01784124.
  EXPECT_NEAR(p_mid, 1.0 / std::sqrt(std::numbers::pi * 1000.0), 1e-4);

  // Even n = 10000, t = 5000 stays finite and accurate.
  channel::QarySymmetricChannel<GF2X<uint8_t>> huge_qsc(gf256, 10000, 0.5);
  const double p_huge = huge_qsc.probability_of_exactly_t_errors(5000);
  EXPECT_TRUE(std::isfinite(p_huge));
  EXPECT_NEAR(p_huge, 1.0 / std::sqrt(std::numbers::pi * 5000.0), 1e-4);

  // epsilon == 0.0 exact probabilities.
  channel::QarySymmetricChannel<GF2X<uint8_t>> zero_qsc(gf256, 16, 0.0);
  EXPECT_DOUBLE_EQ(zero_qsc.probability_of_exactly_t_errors(0), 1.0);
  EXPECT_DOUBLE_EQ(zero_qsc.probability_of_exactly_t_errors(1), 0.0);
  EXPECT_DOUBLE_EQ(zero_qsc.probability_of_at_most_t_errors(0), 1.0);

  // epsilon == 1.0 flips 100% of symbols and has exact probabilities.
  channel::QarySymmetricChannel<GF2X<uint8_t>> one_qsc(gf256, 64, 1.0, 42ULL);
  EXPECT_DOUBLE_EQ(one_qsc.probability_of_exactly_t_errors(64), 1.0);
  EXPECT_DOUBLE_EQ(one_qsc.probability_of_exactly_t_errors(63), 0.0);
  EXPECT_DOUBLE_EQ(one_qsc.probability_of_at_most_t_errors(63), 0.0);
  EXPECT_DOUBLE_EQ(one_qsc.probability_of_at_most_t_errors(64), 1.0);

  std::mt19937_64 rng(55);
  auto msg = GFArray<GF2X<uint8_t>>::Random(gf256, {64}, rng);
  auto rx = one_qsc.transmit(msg);
  EXPECT_EQ(CountArrayDiffs(msg, rx), 64U);
}

// =============================================================================
// Deterministic Reproducibility via SetSeed Across All Channel Types
// =============================================================================

TEST(ChannelReproducibilityTest, SetSeedProducesIdenticalOutputs) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
  std::mt19937_64 rng(999);
  auto msg_q = GFArray<GF2X<uint8_t>>::Random(gf256, {32}, rng);
  auto msg_b = GFArray<GaloisFieldBinary>::Random(gf2, {32}, rng);

  // StaticErrorRateChannel
  {
    channel::StaticErrorRateChannel<GF2X<uint8_t>> c1(gf256, 32, 2, 8);
    channel::StaticErrorRateChannel<GF2X<uint8_t>> c2(gf256, 32, 2, 8);
    c1.SetSeed(42);
    c2.set_seed(42);
    for (int i = 0; i < 10; ++i) {
      EXPECT_EQ(c1.transmit(msg_q), c2.transmit(msg_q));
    }
  }

  // ErrorErasureChannel
  {
    channel::ErrorErasureChannel<GF2X<uint8_t>> c1(gf256, 32, 1, 4, 1, 4);
    channel::ErrorErasureChannel<GF2X<uint8_t>> c2(gf256, 32, 1, 4, 1, 4);
    c1.SetSeed(42);
    c2.SetSeed(42);
    for (int i = 0; i < 10; ++i) {
      auto r1 = c1.transmit_with_erasures(msg_q);
      auto r2 = c2.transmit_with_erasures(msg_q);
      EXPECT_EQ(r1.received, r2.received);
      EXPECT_EQ(r1.erasure_positions, r2.erasure_positions);
      EXPECT_EQ(r1.erasure_mask, r2.erasure_mask);
    }
  }

  // QarySymmetricChannel
  {
    channel::QarySymmetricChannel<GF2X<uint8_t>> c1(gf256, 32, 0.25);
    channel::QarySymmetricChannel<GF2X<uint8_t>> c2(gf256, 32, 0.25);
    c1.SetSeed(42);
    c2.SetSeed(42);
    for (int i = 0; i < 10; ++i) {
      EXPECT_EQ(c1.transmit(msg_q), c2.transmit(msg_q));
    }
  }

  // BinarySymmetricChannel
  {
    channel::BinarySymmetricChannel c1(gf2, 32, 0.2);
    channel::BinarySymmetricChannel c2(gf2, 32, 0.2);
    c1.SetSeed(42);
    c2.SetSeed(42);
    for (int i = 0; i < 10; ++i) {
      EXPECT_EQ(c1.transmit(msg_b), c2.transmit(msg_b));
    }
  }

  // QaryErasureChannel
  {
    channel::QaryErasureChannel<GF2X<uint8_t>> c1(gf256, 32, 0.3);
    channel::QaryErasureChannel<GF2X<uint8_t>> c2(gf256, 32, 0.3);
    c1.SetSeed(42);
    c2.SetSeed(42);
    for (int i = 0; i < 10; ++i) {
      auto r1 = c1.transmit_with_erasures(msg_q);
      auto r2 = c2.transmit_with_erasures(msg_q);
      EXPECT_EQ(r1.received, r2.received);
      EXPECT_EQ(r1.erasure_positions, r2.erasure_positions);
    }
  }

  // AWGNChannel
  {
    channel::AWGNChannel c1(gf2, 32, 2.0);
    channel::AWGNChannel c2(gf2, 32, 2.0);
    c1.SetSeed(42);
    c2.SetSeed(42);
    for (int i = 0; i < 10; ++i) {
      auto r1 = c1.transmit_detailed(msg_b);
      auto r2 = c2.transmit_detailed(msg_b);
      EXPECT_EQ(r1.soft_received, r2.soft_received);
      EXPECT_EQ(r1.llr, r2.llr);
      EXPECT_EQ(r1.hard_decision, r2.hard_decision);
    }
  }
}

// =============================================================================
// New Channel Models: BinarySymmetricChannel, QaryErasureChannel, AWGNChannel
// =============================================================================

TEST(BinarySymmetricChannelTest, BehaviorCapacityAndStatistics) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // Invalid crossover probabilities
  EXPECT_THROW(channel::BinarySymmetricChannel(gf2, 16, -0.01),
               std::invalid_argument);
  EXPECT_THROW(channel::BinarySymmetricChannel(gf2, 16, 1.01),
               std::invalid_argument);
  EXPECT_THROW(channel::BinarySymmetricChannel(
                   gf2, 16, std::numeric_limits<double>::quiet_NaN()),
               std::invalid_argument);

  // Capacity checks
  channel::BinarySymmetricChannel bsc0(gf2, 64, 0.0, 1ULL);
  channel::BinarySymmetricChannel bsc05(gf2, 64, 0.5, 1ULL);
  channel::BinarySymmetricChannel bsc1(gf2, 64, 1.0, 1ULL);
  EXPECT_DOUBLE_EQ(bsc0.capacity(), 1.0);
  EXPECT_NEAR(bsc05.capacity(), 0.0, 1e-12);
  EXPECT_DOUBLE_EQ(bsc1.capacity(), 1.0);

  auto msg = GFArray<GaloisFieldBinary>::Zeros(gf2, {64});
  EXPECT_EQ(CountArrayDiffs(msg, bsc0.transmit(msg)), 0U);
  EXPECT_EQ(CountArrayDiffs(msg, bsc1.transmit(msg)), 64U);

  // Empirical bit-flip rate check with fixed seed using default-field factory.
  constexpr size_t kBits = 10000;
  constexpr double kP = 0.15;
  auto bsc = channel::CreateBinarySymmetricChannel(kBits, kP, 2026ULL);
  auto bits = GFArray<GaloisFieldBinary>::Zeros(bsc->field(), {kBits});
  auto rx = bsc->transmit(bits);
  const double empirical_p =
      static_cast<double>(CountArrayDiffs(bits, rx)) / kBits;
  EXPECT_NEAR(empirical_p, kP, 0.02);
  EXPECT_NE(bsc->to_string().find("Binary symmetric"), std::string::npos);
}

TEST(QaryErasureChannelTest, ErasureRateAndConsistency) {
  auto gf256 = std::make_shared<GF2X<uint8_t>>(8);
  auto gf9 =
      std::make_shared<GaloisFieldExtension<uint32_t>>(std::pair{3ULL, 2ULL});

  // Invalid epsilon
  EXPECT_THROW(channel::QaryErasureChannel<GF2X<uint8_t>>(gf256, 16, -0.1),
               std::invalid_argument);
  EXPECT_THROW(channel::QaryErasureChannel<GF2X<uint8_t>>(gf256, 16, 1.1),
               std::invalid_argument);
  EXPECT_THROW(channel::QaryErasureChannel<GF2X<uint8_t>>(
                   gf256, 16, std::numeric_limits<double>::quiet_NaN()),
               std::invalid_argument);

  // Edge cases epsilon == 0.0 and epsilon == 1.0
  auto msg = GFArray<GF2X<uint8_t>>::Ones(gf256, {32});
  channel::QaryErasureChannel<GF2X<uint8_t>> qec0(gf256, 32, 0.0, 1ULL);
  auto tx0 = qec0.transmit_with_erasures(msg);
  EXPECT_EQ(tx0.received, msg);
  EXPECT_TRUE(tx0.erasure_positions.empty());

  channel::QaryErasureChannel<GF2X<uint8_t>> qec1(gf256, 32, 1.0, 1ULL);
  auto tx1 = qec1.transmit_with_erasures(msg);
  EXPECT_EQ(tx1.erasure_positions.size(), 32U);
  for (size_t i = 0; i < 32; ++i) {
    EXPECT_EQ(tx1.erasure_positions[i], i);
    EXPECT_TRUE(tx1.erasure_mask[i]);
    EXPECT_EQ(tx1.received.Raw(i), gf256->AdditiveIdentity());
  }

  // Size mismatch rejected
  auto bad_msg = GFArray<GF2X<uint8_t>>::Ones(gf256, {10});
  EXPECT_THROW(qec1.transmit_with_erasures(bad_msg), std::invalid_argument);

  // Statistical erasure rate check with fixed seed over GF(3^2)
  constexpr size_t kN = 4000;
  constexpr double kEps = 0.25;
  auto qec = channel::CreateQaryErasureChannel(gf9, kN, kEps, 31415ULL);
  auto ext_msg = GFArray<GaloisFieldExtension<uint32_t>>::Ones(gf9, {kN});
  auto [rx, positions, mask] = qec->transmit_with_erasures(ext_msg);
  EXPECT_EQ(rx.Size(), kN);
  EXPECT_EQ(mask.size(), kN);
  EXPECT_TRUE(std::is_sorted(positions.begin(), positions.end()));

  size_t mask_count = 0;
  for (size_t i = 0; i < kN; ++i) {
    if (mask[i]) {
      ++mask_count;
      EXPECT_EQ(rx.Raw(i), gf9->AdditiveIdentity());
    } else {
      EXPECT_EQ(rx.Raw(i), ext_msg.Raw(i));
    }
  }
  EXPECT_EQ(mask_count, positions.size());
  const double empirical_eps = static_cast<double>(positions.size()) / kN;
  EXPECT_NEAR(empirical_eps, kEps, 0.025);
}

TEST(AWGNChannelTest, BPSKModulationLLRAndEmpiricalBER) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // Invalid parameters
  EXPECT_THROW(
      channel::AWGNChannel(gf2, 16, std::numeric_limits<double>::quiet_NaN()),
      std::invalid_argument);
  EXPECT_THROW(
      channel::AWGNChannel(gf2, 16, std::numeric_limits<double>::infinity()),
      std::invalid_argument);
  EXPECT_THROW(channel::AWGNChannel(gf2, 16, 3.0, 0.0), std::invalid_argument);
  EXPECT_THROW(channel::AWGNChannel(gf2, 16, 3.0, 1.5), std::invalid_argument);
  EXPECT_THROW(channel::AWGNChannel::FromSigma(gf2, 16, -0.5),
               std::invalid_argument);
  EXPECT_THROW(channel::AWGNChannel::FromSigma(
                   gf2, 16, std::numeric_limits<double>::quiet_NaN()),
               std::invalid_argument);

  // BPSK mapping: 0 -> +1.0, 1 -> -1.0
  auto bits = GFArray<GaloisFieldBinary>::Zeros(gf2, {4});
  bits.Raw(0) = 0;
  bits.Raw(1) = 1;
  bits.Raw(2) = 0;
  bits.Raw(3) = 1;
  const auto modulated = channel::AWGNChannel::ModulateBPSK(bits);
  ASSERT_EQ(modulated.size(), 4U);
  EXPECT_DOUBLE_EQ(modulated[0], 1.0);
  EXPECT_DOUBLE_EQ(modulated[1], -1.0);
  EXPECT_DOUBLE_EQ(modulated[2], 1.0);
  EXPECT_DOUBLE_EQ(modulated[3], -1.0);

  // Noiseless sigma == 0.0 produces exact BPSK symbols and +/-inf LLRs
  auto noiseless = channel::AWGNChannel::FromSigma(gf2, 4, 0.0, 1ULL);
  auto det0 = noiseless.transmit_detailed(bits);
  EXPECT_EQ(det0.soft_received, modulated);
  EXPECT_EQ(det0.hard_decision, bits);
  EXPECT_TRUE(std::isinf(det0.llr[0]) && det0.llr[0] > 0.0);
  EXPECT_TRUE(std::isinf(det0.llr[1]) && det0.llr[1] < 0.0);
  EXPECT_DOUBLE_EQ(noiseless.theoretical_ber(), 0.0);

  // Exact LLR formula: LLR = 2 * y / sigma^2, and hard decision matches sign.
  auto ch_sigma = channel::AWGNChannel::FromSigma(gf2, 4, 0.5, 42ULL);
  EXPECT_DOUBLE_EQ(ch_sigma.sigma(), 0.5);
  EXPECT_DOUBLE_EQ(ch_sigma.noise_variance(), 0.25);
  auto det = ch_sigma.transmit_detailed(bits);
  for (size_t i = 0; i < 4; ++i) {
    const double expected_llr = 2.0 * det.soft_received[i] / 0.25;
    EXPECT_DOUBLE_EQ(det.llr[i], expected_llr);
    const uint8_t expected_hard = (det.llr[i] >= 0.0) ? uint8_t{0} : uint8_t{1};
    EXPECT_EQ(det.hard_decision.Raw(i), expected_hard);
  }

  // Eb/N0 with code rate 0.5 has 2x noise variance vs Es/N0 at same dB.
  auto ch_r1 = channel::AWGNChannel::FromEsN0(gf2, 16, 3.0);
  auto ch_r05 = channel::AWGNChannel::FromEbN0(gf2, 16, 3.0, 0.5);
  EXPECT_NEAR(ch_r05.noise_variance(), 2.0 * ch_r1.noise_variance(), 1e-12);

  // Low SNR (0 dB) vs High SNR (8 dB) empirical BER checks against theoretical.
  constexpr size_t kN = 20000;
  auto zeros = GFArray<GaloisFieldBinary>::Zeros(gf2, {kN});

  auto low_snr = channel::AWGNChannel::FromEsN0(gf2, kN, 0.0, 1234ULL);
  auto rx_low = low_snr.transmit(zeros);
  const double ber_low = static_cast<double>(CountArrayDiffs(zeros, rx_low)) /
                         static_cast<double>(kN);
  // At Es/N0 = 0 dB, theoretical BPSK BER = Q(sqrt(2)) ~= 0.07865.
  EXPECT_NEAR(ber_low, low_snr.theoretical_ber(), 0.01);

  auto high_snr = channel::AWGNChannel::FromEsN0(gf2, kN, 8.0, 1234ULL);
  auto rx_high = high_snr.transmit(zeros);
  const double ber_high = static_cast<double>(CountArrayDiffs(zeros, rx_high)) /
                          static_cast<double>(kN);
  EXPECT_LT(ber_high, 0.002);
  EXPECT_LT(ber_high, ber_low);
}

// =============================================================================
// Monte Carlo BER/SER/FER Simulation Framework
// =============================================================================

TEST(SimulationFrameworkTest, DeterministicMonteCarloAndEarlyStopping) {
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  // 1. Deterministic functional simulation: (3, 1) repetition code over BSC.
  auto run_sim = [&](uint64_t seed) {
    channel::BinarySymmetricChannel bsc(gf2, 3, 0.1, seed ^ 0x55AAULL);
    auto encode = [&](const GFArray<GaloisFieldBinary> &msg) {
      auto cw = GFArray<GaloisFieldBinary>::Zeros(gf2, {3});
      for (size_t i = 0; i < 3; ++i) {
        cw.Raw(i) = msg.Raw(0);
      }
      return cw;
    };
    auto ch_fn = [&](const GFArray<GaloisFieldBinary> &cw) {
      return bsc.transmit(cw);
    };
    auto decode = [&](const GFArray<GaloisFieldBinary> &rx) {
      size_t ones = (rx.Raw(0) != 0) + (rx.Raw(1) != 0) + (rx.Raw(2) != 0);
      auto out = GFArray<GaloisFieldBinary>::Zeros(gf2, {1});
      out.Raw(0) = (ones >= 2) ? uint8_t{1} : uint8_t{0};
      return out;
    };
    return channel::SimulateCodePerformance(
        gf2, 1, encode, ch_fn, decode,
        channel::SimulationConfig{
            .max_frames = 5000, .target_frame_errors = 0, .seed = seed});
  };

  const auto stats1 = run_sim(42);
  const auto stats2 = run_sim(42);
  EXPECT_EQ(stats1.frames_sent, 5000U);
  EXPECT_EQ(stats1.frame_errors, stats2.frame_errors);
  EXPECT_EQ(stats1.symbol_errors, stats2.symbol_errors);
  EXPECT_EQ(stats1.bit_errors, stats2.bit_errors);
  EXPECT_EQ(stats1.decoding_failures, 0U);
  // Exact theoretical FER for (3,1) rep code over BSC(p=0.1) is
  // 3*p^2*(1-p) + p^3 = 0.028.
  EXPECT_NEAR(stats1.FER(), 0.028, 0.01);
  EXPECT_DOUBLE_EQ(stats1.FER(), stats1.BER());
  EXPECT_DOUBLE_EQ(stats1.SER(), stats1.BER());

  // 2. Early stopping when target_frame_errors is reached.
  channel::BinarySymmetricChannel noisy_bsc(gf2, 4, 0.4, 77ULL);
  auto identity = [](const GFArray<GaloisFieldBinary> &x) { return x; };
  auto noisy_ch = [&](const GFArray<GaloisFieldBinary> &cw) {
    return noisy_bsc.transmit(cw);
  };
  const auto early_stats = channel::SimulateCodePerformance(
      gf2, 4, identity, noisy_ch, identity, 10000, 25, 42ULL);
  EXPECT_EQ(early_stats.frame_errors, 25U);
  EXPECT_LT(early_stats.frames_sent, 10000U);
  EXPECT_GT(early_stats.BER(), 0.0);
  EXPECT_LE(early_stats.BER(), 1.0);

  // 3. Decoding failures counted for DecodingResult-like and std::optional.
  auto always_fail_opt = [](const GFArray<GaloisFieldBinary> &)
      -> std::optional<GFArray<GaloisFieldBinary>> { return std::nullopt; };
  auto perfect_ch = [](const GFArray<GaloisFieldBinary> &cw) { return cw; };
  const auto fail_stats = channel::SimulateCodePerformance(
      gf2, 4, identity, perfect_ch, always_fail_opt, 10, 0, 7ULL);
  EXPECT_EQ(fail_stats.frames_sent, 10U);
  EXPECT_EQ(fail_stats.decoding_failures, 10U);
  EXPECT_EQ(fail_stats.frame_errors, 10U);
  EXPECT_DOUBLE_EQ(fail_stats.FER(), 1.0);
  EXPECT_DOUBLE_EQ(fail_stats.SER(), 1.0);
  EXPECT_DOUBLE_EQ(fail_stats.BER(), 1.0);

  // 4. Direct (code, channel) overload with BSC and QaryErasureChannel.
  MockRepetition3Code rep_code(gf2);
  channel::BinarySymmetricChannel bsc_code(gf2, 3, 0.1);
  const auto code_stats_a = channel::SimulateCodePerformance(
      rep_code, bsc_code, channel::SimulationConfig{2000, 0, 123ULL});
  const auto code_stats_b = channel::SimulateCodePerformance(
      rep_code, bsc_code, channel::SimulationConfig{2000, 0, 123ULL});
  EXPECT_EQ(code_stats_a.frame_errors, code_stats_b.frame_errors);
  EXPECT_NEAR(code_stats_a.FER(), 0.028, 0.015);

  // With BEC(eps=0.5), rep-3 fails iff all 3 bits are erased (prob = 0.125).
  channel::BinaryErasureChannel bec_code(gf2, 3, 0.5);
  const auto bec_stats = channel::SimulateCodePerformance(
      rep_code, bec_code, channel::SimulationConfig{2000, 0, 456ULL});
  EXPECT_EQ(bec_stats.decoding_failures, bec_stats.frame_errors);
  EXPECT_NEAR(bec_stats.FER(), 0.125, 0.025);

  // 5. Edge cases: max_frames == 0 and null field.
  const auto zero_stats = channel::SimulateCodePerformance(
      gf2, 4, identity, perfect_ch, identity, 0, 0, 1ULL);
  EXPECT_EQ(zero_stats.frames_sent, 0U);
  EXPECT_DOUBLE_EQ(zero_stats.FER(), 0.0);
  EXPECT_DOUBLE_EQ(zero_stats.SER(), 0.0);
  EXPECT_DOUBLE_EQ(zero_stats.BER(), 0.0);

  std::shared_ptr<GaloisFieldBinary> null_gf2;
  EXPECT_THROW(channel::SimulateCodePerformance(null_gf2, 4, identity,
                                                perfect_ch, identity),
               std::invalid_argument);
}

// =============================================================================
// Null Field & Zero Input Size Edge Cases Across All Channel Types
// =============================================================================

TEST(ChannelEdgeCasesTest, NullFieldRejected) {
  std::shared_ptr<GF2X<uint8_t>> null_gf16;
  std::shared_ptr<GaloisFieldBinary> null_gf2;

  EXPECT_THROW(channel::StaticErrorRateChannel<GF2X<uint8_t>>(null_gf16, 8, 1),
               std::invalid_argument);
  EXPECT_THROW(channel::ErrorErasureChannel<GF2X<uint8_t>>(null_gf16, 8, 1, 1),
               std::invalid_argument);
  EXPECT_THROW(channel::QarySymmetricChannel<GF2X<uint8_t>>(null_gf16, 8, 0.1),
               std::invalid_argument);
  EXPECT_THROW(channel::BinarySymmetricChannel(null_gf2, 8, 0.1),
               std::invalid_argument);
  EXPECT_THROW(channel::QaryErasureChannel<GF2X<uint8_t>>(null_gf16, 8, 0.1),
               std::invalid_argument);
  EXPECT_THROW(channel::AWGNChannel(null_gf2, 8, 3.0), std::invalid_argument);
  EXPECT_THROW(channel::AWGNChannel::FromSigma(null_gf2, 8, 0.5),
               std::invalid_argument);
}

TEST(ChannelEdgeCasesTest, ZeroInputSizeHandledCleanly) {
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  auto gf2 = std::make_shared<GaloisFieldBinary>();

  auto empty_q = GFArray<GF2X<uint8_t>>::Zeros(gf16, {0});
  auto empty_b = GFArray<GaloisFieldBinary>::Zeros(gf2, {0});
  auto empty_legacy =
      xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>>::from_shape({0});

  // StaticErrorRateChannel
  channel::StaticErrorRateChannel<GF2X<uint8_t>> s_ch(gf16, 0, 0);
  EXPECT_EQ(s_ch.transmit(empty_q).Size(), 0U);
  EXPECT_EQ(s_ch.transmit(empty_legacy).size(), 0U);
  EXPECT_THROW(channel::StaticErrorRateChannel<GF2X<uint8_t>>(gf16, 0, 1),
               std::invalid_argument);

  // ErrorErasureChannel
  channel::ErrorErasureChannel<GF2X<uint8_t>> ee_ch(gf16, 0, 0, 0);
  auto ee_res = ee_ch.transmit_with_erasures(empty_q);
  EXPECT_EQ(ee_res.received.Size(), 0U);
  EXPECT_TRUE(ee_res.erasure_positions.empty());
  EXPECT_EQ(ee_res.erasure_mask.size(), 0U);
  auto ee_legacy = ee_ch.transmit_with_erasures(empty_legacy);
  EXPECT_EQ(ee_legacy.first.size(), 0U);
  EXPECT_EQ(ee_legacy.second.size(), 0U);

  // QarySymmetricChannel
  channel::QarySymmetricChannel<GF2X<uint8_t>> qsc(gf16, 0, 0.5);
  EXPECT_EQ(qsc.transmit(empty_q).Size(), 0U);
  EXPECT_DOUBLE_EQ(qsc.probability_of_exactly_t_errors(0), 1.0);
  EXPECT_DOUBLE_EQ(qsc.probability_of_exactly_t_errors(1), 0.0);

  // BinarySymmetricChannel
  channel::BinarySymmetricChannel bsc(gf2, 0, 0.5);
  EXPECT_EQ(bsc.transmit(empty_b).Size(), 0U);

  // QaryErasureChannel
  channel::QaryErasureChannel<GF2X<uint8_t>> qec(gf16, 0, 0.5);
  auto qec_res = qec.transmit_with_erasures(empty_q);
  EXPECT_EQ(qec_res.received.Size(), 0U);
  EXPECT_TRUE(qec_res.erasure_positions.empty());
  EXPECT_EQ(qec_res.erasure_mask.size(), 0U);

  // AWGNChannel
  channel::AWGNChannel awgn(gf2, 0, 3.0);
  auto awgn_res = awgn.transmit_detailed(empty_b);
  EXPECT_TRUE(awgn_res.soft_received.empty());
  EXPECT_TRUE(awgn_res.llr.empty());
  EXPECT_EQ(awgn_res.hard_decision.Size(), 0U);
}
