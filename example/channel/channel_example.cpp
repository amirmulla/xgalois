/**
 * @file channel_example.cpp
 * @brief Examples demonstrating channel usage with XGalois fields
 */

#include <bitset>
#include <iomanip>
#include <iostream>
#include <memory>

#include "xgalois/channel/channel.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/field/gf_element.hpp"
#include "xgalois/field/gf_prime.hpp"

using namespace xg;

int main() {
  std::cout << "XGalois Channel Examples" << std::endl;
  std::cout << "=======================" << std::endl;

  // Example 1: Binary field channels
  std::cout << "\n=== Binary Field Channel Examples ===" << std::endl;

  // Create GF(2^4) field
  auto gf16 = std::make_shared<GF2X<uint8_t>>(4);
  std::cout << "Created ";
  gf16->Print(std::cout);
  std::cout << std::endl;

  // Create a test message using field elements
  const size_t message_length = 8;
  using ElementType = GaloisFieldElementBase<GF2X<uint8_t>>;
  xt::xarray<ElementType> message = xt::empty<ElementType>({message_length});

  // Initialize message with field elements
  for (size_t i = 0; i < message_length; ++i) {
    message[i] = ElementType(static_cast<uint8_t>(i + 1), gf16);
  }

  std::cout << "\nOriginal message (in GF(2^4)): [";
  for (size_t i = 0; i < message.size(); ++i) {
    std::cout << static_cast<int>(message[i].Value());
    if (i < message.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Example 1a: Static Error Rate Channel
  std::cout << "\n--- Static Error Rate Channel ---" << std::endl;
  auto static_channel =
      channels::CreateStaticErrorRateChannel(gf16, message_length, 2);
  std::cout << "Channel: " << *static_channel << std::endl;

  auto transmitted1 = static_channel->transmit(message);
  std::cout << "transmitted: [";
  for (size_t i = 0; i < transmitted1.size(); ++i) {
    std::cout << static_cast<int>(transmitted1[i].Value());
    if (i < transmitted1.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Count errors
  size_t error_count = 0;
  for (size_t i = 0; i < message_length; ++i) {
    if (message[i] != transmitted1[i]) error_count++;
  }
  std::cout << "Number of errors introduced: " << error_count << std::endl;

  // Example 1b: Variable Error Rate Channel
  std::cout << "\n--- Variable Error Rate Channel ---" << std::endl;
  auto var_channel =
      channels::CreateStaticErrorRateChannel(gf16, message_length, 1, 3);
  std::cout << "Channel: " << *var_channel << std::endl;

  for (int trial = 0; trial < 3; ++trial) {
    auto transmitted = var_channel->transmit(message);
    size_t errors = 0;
    for (size_t i = 0; i < message_length; ++i) {
      if (message[i] != transmitted[i]) errors++;
    }
    std::cout << "Trial " << (trial + 1) << " - Errors: " << errors
              << std::endl;
  }

  // Example 1c: Error-Erasure Channel
  std::cout << "\n--- Error-Erasure Channel ---" << std::endl;
  auto ee_channel =
      channels::CreateErrorErasureChannel(gf16, message_length, 2, 1);
  std::cout << "Channel: " << *ee_channel << std::endl;

  auto transmitted_pair = ee_channel->transmit_unsafe_with_erasures(message);
  auto& transmitted2 = transmitted_pair.first;
  auto& erasures = transmitted_pair.second;

  std::cout << "transmitted: [";
  for (size_t i = 0; i < transmitted2.size(); ++i) {
    std::cout << static_cast<int>(transmitted2[i].Value());
    if (i < transmitted2.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "erasures: [";
  for (size_t i = 0; i < erasures.size(); ++i) {
    std::cout << (erasures[i] ? "1" : "0");
    if (i < erasures.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Example 2: Prime field channels
  std::cout << "\n=== Prime Field Channel Examples ===" << std::endl;

  // Create GF(31) field
  auto gf31 = std::make_shared<GFP<uint32_t>>(31);
  std::cout << "Created ";
  gf31->Print(std::cout);
  std::cout << std::endl;

  // Create a test message using field elements
  using PrimeElementType = GaloisFieldElementBase<GFP<uint32_t>>;
  xt::xarray<PrimeElementType> prime_message =
      xt::empty<PrimeElementType>({message_length});

  // Initialize message with field elements
  for (size_t i = 0; i < message_length; ++i) {
    prime_message[i] = PrimeElementType(static_cast<uint32_t>(i + 1), gf31);
  }

  std::cout << "\nOriginal message (in GF(31)): [";
  for (size_t i = 0; i < prime_message.size(); ++i) {
    std::cout << prime_message[i].Value();
    if (i < prime_message.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Example 2a: q-ary Symmetric Channel
  std::cout << "\n--- q-ary Symmetric Channel ---" << std::endl;
  double error_prob = 0.1;
  auto qsc =
      channels::CreateQarySymmetricChannel(gf31, message_length, error_prob);
  std::cout << "Channel: " << *qsc << std::endl;

  auto transmitted3 = qsc->transmit(prime_message);
  std::cout << "transmitted: [";
  for (size_t i = 0; i < transmitted3.size(); ++i) {
    std::cout << transmitted3[i].Value();
    if (i < transmitted3.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;

  // Count errors
  error_count = 0;
  for (size_t i = 0; i < message_length; ++i) {
    if (prime_message[i] != transmitted3[i]) error_count++;
  }
  std::cout << "Number of errors introduced: " << error_count << std::endl;

  return 0;
}
