/**
 * @file channel_test.cpp
 * @brief Comprehensive tests for channel implementations
 */

#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>

#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_prime.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/channel/channel.hpp"

using namespace xg;

class ChannelTest : public ::testing::Test {
protected:
    void SetUp() override {
        gf16 = std::make_shared<GF2X<uint8_t>>(4);
        gf17 = std::make_shared<GaloisFieldPrime<uint8_t>>(17);

        // Create test message using GaloisFieldElementBase
        message_length = 8;
        message = xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>>::from_shape({message_length});
        for (size_t i = 0; i < message_length; ++i) {
            message[i] = GaloisFieldElementBase<GF2X<uint8_t>>(i + 1, gf16);
        }

        // Create message for prime field
        prime_message = xt::xarray<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>>::from_shape({message_length});
        for (size_t i = 0; i < message_length; ++i) {
            prime_message[i] = GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>((i * 2 + 1) % 17, gf17);
        }
    }

    std::shared_ptr<GF2X<uint8_t>> gf16;
    std::shared_ptr<GaloisFieldPrime<uint8_t>> gf17;
    xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>> message;
    xt::xarray<GaloisFieldElementBase<GaloisFieldPrime<uint8_t>>> prime_message;
    size_t message_length;
};

// Test StaticErrorRateChannel
TEST_F(ChannelTest, StaticErrorRateChannelFixed) {
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 2);

    // Check channel properties
    EXPECT_EQ(channel->input_size(), message_length);
    EXPECT_EQ(channel->output_size(), message_length);
    EXPECT_EQ(channel->number_errors().first, 2);
    EXPECT_EQ(channel->number_errors().second, 2);

    // Test transmission
    auto transmitted = channel->transmit(message);
    EXPECT_EQ(transmitted.size(), message.size());

    // Count errors
    size_t error_count = 0;
    for (size_t i = 0; i < message_length; ++i) {
        if (message[i] != transmitted[i]) error_count++;
    }
    EXPECT_EQ(error_count, 2);

    // Test with zero errors
    auto no_error_channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 0);
    auto no_error_transmitted = no_error_channel->transmit(message);
    for (size_t i = 0; i < message_length; ++i) {
        EXPECT_EQ(message[i], no_error_transmitted[i]);
    }
}

TEST_F(ChannelTest, StaticErrorRateChannelVariable) {
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 1, 3);

    // Check channel properties
    EXPECT_EQ(channel->number_errors().first, 1);
    EXPECT_EQ(channel->number_errors().second, 3);

    // Test multiple transmissions
    std::vector<size_t> error_counts;
    for (int i = 0; i < 100; ++i) {
        auto transmitted = channel->transmit(message);
        size_t errors = 0;
        for (size_t j = 0; j < message_length; ++j) {
            if (message[j] != transmitted[j]) errors++;
        }
        error_counts.push_back(errors);
    }

    // Check that all error counts are within range
    for (size_t count : error_counts) {
        EXPECT_GE(count, 1);
        EXPECT_LE(count, 3);
    }

    // Check that we get some variation
    auto min_errors = *std::min_element(error_counts.begin(), error_counts.end());
    auto max_errors = *std::max_element(error_counts.begin(), error_counts.end());
    EXPECT_LT(min_errors, max_errors);
}

TEST_F(ChannelTest, StaticErrorRateChannelInvalidInput) {
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 2);

    // Test wrong message size
    xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>> wrong_size_message =
        xt::xarray<GaloisFieldElementBase<GF2X<uint8_t>>>::from_shape({3});
    for (size_t i = 0; i < 3; ++i) {
        wrong_size_message[i] = GaloisFieldElementBase<GF2X<uint8_t>>(i + 1, gf16);
    }
    EXPECT_THROW(channel->transmit(wrong_size_message), std::invalid_argument);

    // Test invalid error range
    EXPECT_THROW(std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 3, 2), std::invalid_argument);

    // Test too many errors - this should throw when transmitting, not during construction
    auto too_many_errors_channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, message_length + 1);
    EXPECT_THROW(too_many_errors_channel->transmit(message), std::invalid_argument);
}

// Test ErrorErasureChannel
TEST_F(ChannelTest, ErrorErasureChannelFixed) {
    auto channel = std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(gf16, message_length, 2, 1);

    // Check channel properties
    EXPECT_EQ(channel->number_errors().first, 2);
    EXPECT_EQ(channel->number_errors().second, 2);
    EXPECT_EQ(channel->number_erasures().first, 1);
    EXPECT_EQ(channel->number_erasures().second, 1);

    // Test transmission with erasures
    auto [transmitted, erasures] = channel->transmit_unsafe_with_erasures(message);
    EXPECT_EQ(transmitted.size(), message.size());
    EXPECT_EQ(erasures.size(), message.size());

    // Count errors and erasures
    size_t error_count = 0;
    size_t erasure_count = 0;
    for (size_t i = 0; i < message_length; ++i) {
        if (erasures[i]) {
            erasure_count++;
            GaloisFieldElementBase<GF2X<uint8_t>> zero_element(gf16->AdditiveIdentity(), gf16);
            EXPECT_EQ(transmitted[i], zero_element);
        } else if (message[i] != transmitted[i]) {
            error_count++;
        }
    }
    EXPECT_EQ(error_count, 2);
    EXPECT_EQ(erasure_count, 1);

    // Test regular transmission (without erasure information)
    auto regular_transmitted = channel->transmit(message);
    EXPECT_EQ(regular_transmitted.size(), message.size());
}

TEST_F(ChannelTest, ErrorErasureChannelVariable) {
    auto channel = std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(gf16, message_length, 1, 2, 0, 2);

    // Test multiple transmissions
    for (int i = 0; i < 50; ++i) {
        auto [transmitted, erasures] = channel->transmit_unsafe_with_erasures(message);

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

// Test QarySymmetricChannel
TEST_F(ChannelTest, QarySymmetricChannelBasic) {
    const double epsilon = 0.1;
    auto channel = std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, epsilon);

    // Check channel properties
    EXPECT_EQ(channel->input_size(), message_length);
    EXPECT_EQ(channel->output_size(), message_length);
    EXPECT_DOUBLE_EQ(channel->error_probability(), epsilon);

    // Test transmission
    auto transmitted = channel->transmit(message);
    EXPECT_EQ(transmitted.size(), message.size());

    // Test with zero error probability
    auto perfect_channel = std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, 0.0);
    auto perfect_transmitted = perfect_channel->transmit(message);
    for (size_t i = 0; i < message_length; ++i) {
        EXPECT_EQ(message[i], perfect_transmitted[i]);
    }
}

TEST_F(ChannelTest, QarySymmetricChannelProbabilities) {
    const double epsilon = 0.1;
    auto channel = std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, epsilon);

    // Test probability calculations
    double prob_0 = channel->probability_of_exactly_t_errors(0);
    double prob_1 = channel->probability_of_exactly_t_errors(1);
    double prob_all = channel->probability_of_exactly_t_errors(message_length);

    EXPECT_GT(prob_0, 0.0);
    EXPECT_GT(prob_1, 0.0);
    EXPECT_GT(prob_all, 0.0);

    // Test cumulative probabilities
    double cum_prob_0 = channel->probability_of_at_most_t_errors(0);
    double cum_prob_1 = channel->probability_of_at_most_t_errors(1);

    EXPECT_DOUBLE_EQ(cum_prob_0, prob_0);
    EXPECT_DOUBLE_EQ(cum_prob_1, prob_0 + prob_1);

    // Test that total probability is approximately 1
    double total_prob = 0.0;
    for (size_t t = 0; t <= message_length; ++t) {
        total_prob += channel->probability_of_exactly_t_errors(t);
    }
    EXPECT_NEAR(total_prob, 1.0, 1e-10);
}

TEST_F(ChannelTest, QarySymmetricChannelStatistics) {
    const double epsilon = 0.1;
    const int num_trials = 1000;
    auto channel = std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, epsilon);

    std::vector<size_t> error_counts(message_length + 1, 0);

    for (int trial = 0; trial < num_trials; ++trial) {
        auto transmitted = channel->transmit(message);
        size_t errors = 0;
        for (size_t i = 0; i < message_length; ++i) {
            if (message[i] != transmitted[i]) errors++;
        }
        error_counts[errors]++;
    }

    // Check that empirical probabilities are reasonably close to theoretical
    for (size_t t = 0; t <= message_length; ++t) {
        double empirical_prob = static_cast<double>(error_counts[t]) / num_trials;
        double theoretical_prob = channel->probability_of_exactly_t_errors(t);

        // Allow for some statistical variation
        if (theoretical_prob > 0.01) {
            EXPECT_NEAR(empirical_prob, theoretical_prob, 0.05);
        }
    }
}

TEST_F(ChannelTest, QarySymmetricChannelInvalidInput) {
    // Test invalid epsilon values
    EXPECT_THROW(std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, -0.1), std::invalid_argument);
    EXPECT_THROW(std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, 1.1), std::invalid_argument);

    // Valid boundary values
    EXPECT_NO_THROW(std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, 0.0));
    EXPECT_NO_THROW(std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, 1.0));
}

// Test with different field types
TEST_F(ChannelTest, ChannelWithPrimeField) {
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GaloisFieldPrime<uint8_t>>>(gf17, message_length, 1);

    auto transmitted = channel->transmit(prime_message);
    EXPECT_EQ(transmitted.size(), prime_message.size());

    // Count errors
    size_t error_count = 0;
    for (size_t i = 0; i < message_length; ++i) {
        if (prime_message[i] != transmitted[i]) error_count++;
    }
    EXPECT_EQ(error_count, 1);
}

// Test channel string representation
TEST_F(ChannelTest, ChannelStringRepresentation) {
    auto static_channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 2);
    auto var_channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 1, 3);
    auto ee_channel = std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(gf16, message_length, 2, 1);
    auto qary_channel = std::make_shared<channels::QarySymmetricChannel<GF2X<uint8_t>>>(gf16, message_length, 0.1);

    std::string static_str = static_channel->to_string();
    std::string var_str = var_channel->to_string();
    std::string ee_str = ee_channel->to_string();
    std::string qary_str = qary_channel->to_string();

    EXPECT_FALSE(static_str.empty());
    EXPECT_FALSE(var_str.empty());
    EXPECT_FALSE(ee_str.empty());
    EXPECT_FALSE(qary_str.empty());

    // Check that descriptions contain expected keywords
    EXPECT_NE(static_str.find("Static error rate"), std::string::npos);
    EXPECT_NE(var_str.find("between"), std::string::npos);
    EXPECT_NE(ee_str.find("Error-and-erasure"), std::string::npos);
    EXPECT_NE(qary_str.find("q-ary symmetric"), std::string::npos);
}

// Test channel operator() functionality
TEST_F(ChannelTest, ChannelOperatorCall) {
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 1);

    // Test that operator() works the same as transmit()
    auto transmitted1 = channel->transmit(message);
    auto transmitted2 = (*channel)(message);

    // Both should introduce the same number of errors (1)
    size_t errors1 = 0, errors2 = 0;
    for (size_t i = 0; i < message_length; ++i) {
        if (message[i] != transmitted1[i]) errors1++;
        if (message[i] != transmitted2[i]) errors2++;
    }

    EXPECT_EQ(errors1, 1);
    EXPECT_EQ(errors2, 1);
}

// Test factory functions
TEST_F(ChannelTest, FactoryFunctions) {
    // Test CreateStaticErrorRateChannel
    auto static_channel = channels::CreateStaticErrorRateChannel(gf16, message_length, 2);
    EXPECT_EQ(static_channel->number_errors().first, 2);
    EXPECT_EQ(static_channel->number_errors().second, 2);

    // Test CreateStaticErrorRateChannel with range
    auto var_channel = channels::CreateStaticErrorRateChannel(gf16, message_length, 1, 3);
    EXPECT_EQ(var_channel->number_errors().first, 1);
    EXPECT_EQ(var_channel->number_errors().second, 3);

    // Test CreateErrorErasureChannel
    auto ee_channel = channels::CreateErrorErasureChannel(gf16, message_length, 2, 1);
    EXPECT_EQ(ee_channel->number_errors().first, 2);
    EXPECT_EQ(ee_channel->number_errors().second, 2);
    EXPECT_EQ(ee_channel->number_erasures().first, 1);
    EXPECT_EQ(ee_channel->number_erasures().second, 1);

    // Test CreateErrorErasureChannel with range
    auto var_ee_channel = channels::CreateErrorErasureChannel(gf16, message_length, 1, 2, 0, 1);
    EXPECT_EQ(var_ee_channel->number_errors().first, 1);
    EXPECT_EQ(var_ee_channel->number_errors().second, 2);
    EXPECT_EQ(var_ee_channel->number_erasures().first, 0);
    EXPECT_EQ(var_ee_channel->number_erasures().second, 1);

    // Test CreateQarySymmetricChannel
    auto qary_channel = channels::CreateQarySymmetricChannel(gf16, message_length, 0.1);
    EXPECT_DOUBLE_EQ(qary_channel->error_probability(), 0.1);
}

// Test channel field compatibility
TEST_F(ChannelTest, ChannelFieldCompatibility) {
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, 1);

    // Test that channel uses the correct field
    EXPECT_EQ(channel->field(), gf16);

    // Test that channel works with field elements
    auto transmitted = channel->transmit(message);
    EXPECT_EQ(transmitted.size(), message.size());

    // Verify all elements are valid field elements
    for (size_t i = 0; i < transmitted.size(); ++i) {
        EXPECT_EQ(transmitted[i].Field(), gf16);
    }
}

// Test error handling with invalid parameters
TEST_F(ChannelTest, ErrorHandling) {
    // Test ErrorErasureChannel with invalid ranges
    EXPECT_THROW(std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(gf16, message_length, 3, 2, 1, 1), std::invalid_argument);
    EXPECT_THROW(std::make_shared<channels::ErrorErasureChannel<GF2X<uint8_t>>>(gf16, message_length, 1, 1, 3, 2), std::invalid_argument);

    // Test that too many errors throws exception during position generation
    // This would be caught in the generate_error_positions method
    auto channel = std::make_shared<channels::StaticErrorRateChannel<GF2X<uint8_t>>>(gf16, message_length, message_length);
    EXPECT_NO_THROW(channel->transmit(message)); // Should work with errors = input_size
}
