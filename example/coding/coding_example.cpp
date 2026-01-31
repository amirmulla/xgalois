#include "xgalois/coding/coding.hpp"
#include "xgalois/field/gf_extension.hpp"
#include "xgalois/field/gf_prime.hpp"
#include <iostream>
#include <vector>
#include <xtensor/containers/xadapt.hpp>
#include <xtensor/containers/xarray.hpp>

using namespace xg;
using namespace xg::coding;

// Example 1: Reed-Solomon code over GF(16)
void example_reed_solomon() {
    std::cout << "=== Reed-Solomon Code Example ===" << std::endl;

    // Create GF(16) = GF(2^4)
    auto field = std::make_shared<GaloisFieldBinary<uint8_t>>(4);

    // Create RS(15, 11) code - can correct 2 errors
    size_t length = 15;
    size_t dimension = 11;

    auto rs_code = CreateReedSolomonCode(field, length, dimension);

    std::cout << "Code: " << rs_code->ToString() << std::endl;
    std::cout << "Minimum distance: " << rs_code->MinimumDistance() << std::endl;
    std::cout << "Error correction capability: " << ErrorCorrectionCapability(*rs_code) << std::endl;
    std::cout << "Code rate: " << CodeRate(*rs_code) << std::endl;
    std::cout << "Is MDS: " << (IsMDS(*rs_code) ? "Yes" : "No") << std::endl;

    // Create a message using xtensor
    xt::xarray<uint8_t> message = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    // Encode the message
    auto codeword = rs_code->Encode(message);

    std::cout << "Original message: ";
    for (size_t i = 0; i < message.size(); ++i) {
        std::cout << static_cast<int>(message(i)) << " ";
    }
    std::cout << std::endl;

    std::cout << "Encoded codeword: ";
    for (size_t i = 0; i < codeword.size(); ++i) {
        std::cout << static_cast<int>(codeword(i)) << " ";
    }
    std::cout << std::endl;

    // Introduce errors
    auto received = codeword;
    received(3) = 15; // Error at position 3
    received(7) = 14; // Error at position 7

    std::cout << "Received word (with 2 errors): ";
    for (size_t i = 0; i < received.size(); ++i) {
        std::cout << static_cast<int>(received(i)) << " ";
    }
    std::cout << std::endl;

    // Decode
    try {
        auto decoded_message = rs_code->DecodeToMessage(received);

        std::cout << "Decoded message: ";
        for (size_t i = 0; i < decoded_message.size(); ++i) {
            std::cout << static_cast<int>(decoded_message(i)) << " ";
        }
        std::cout << std::endl;

        // Check if decoding was successful (compare xtensor arrays)
        bool success = true;
        if (decoded_message.size() == message.size()) {
            for (size_t i = 0; i < message.size(); ++i) {
                if (decoded_message(i) != message(i)) {
                    success = false;
                    break;
                }
            }
        } else {
            success = false;
        }
        std::cout << "Decoding " << (success ? "successful" : "failed") << std::endl;

    } catch (const std::exception& e) {
        std::cout << "Decoding failed: " << e.what() << std::endl;
    }

    std::cout << std::endl;
}

// Example 2: Generalized Reed-Solomon code
void example_generalized_reed_solomon() {
    std::cout << "=== Generalized Reed-Solomon Code Example ===" << std::endl;

    // Create GF(8) = GF(2^3)
    auto field = std::make_shared<GaloisFieldBinary<uint8_t>>(3);

    // Define evaluation points and column multipliers using xtensor
    xt::xarray<uint8_t> eval_points = {1, 2, 3, 4, 5, 6, 7}; // 7 points
    xt::xarray<uint8_t> multipliers = {1, 1, 2, 2, 3, 3, 4}; // Non-uniform multipliers
    size_t dimension = 4;

    auto grs_code = CreateGeneralizedReedSolomonCode(field, eval_points, multipliers, dimension);

    std::cout << "Code: " << grs_code->ToString() << std::endl;
    std::cout << "Minimum distance: " << grs_code->MinimumDistance() << std::endl;
    std::cout << "Error correction capability: " << ErrorCorrectionCapability(*grs_code) << std::endl;
    std::cout << "Is MDS: " << (IsMDS(*grs_code) ? "Yes" : "No") << std::endl;

    // Test encoding and decoding
    xt::xarray<uint8_t> message = {1, 3, 5, 7};
    auto codeword = grs_code->Encode(message);
    auto decoded = grs_code->DecodeToMessage(codeword);

    // Compare xtensor arrays
    bool success = true;
    if (decoded.size() == message.size()) {
        for (size_t i = 0; i < message.size(); ++i) {
            if (decoded(i) != message(i)) {
                success = false;
                break;
            }
        }
    } else {
        success = false;
    }

    std::cout << "Encoding/decoding test: " <<
                 (success ? "Passed" : "Failed") << std::endl;

    std::cout << std::endl;
}

// Example 3: Cyclic code
void example_cyclic_code() {
    std::cout << "=== Cyclic Code Example ===" << std::endl;

    // Create GF(2)
    auto field = std::make_shared<GaloisFieldPrime<uint8_t>>(2);

    // Create generator polynomial g(x) = x^3 + x + 1
    std::vector<uint8_t> gen_coeffs = {1, 1, 0, 1}; // 1 + x + x^3
    auto gen_poly = xg::poly::PolyDense<uint8_t>(gen_coeffs, field);

    // Create (7, 4) cyclic code
    size_t length = 7;
    auto cyclic_code = CreateCyclicCode(field, length, gen_poly);

    std::cout << "Code: " << cyclic_code->ToString() << std::endl;
    std::cout << "Generator polynomial degree: " << gen_poly.Degree() << std::endl;
    std::cout << "Minimum distance: " << cyclic_code->MinimumDistance() << std::endl;

    // Test cyclic property
    xt::xarray<uint8_t> message = {1, 0, 1, 1};
    auto codeword = cyclic_code->Encode(message);
    auto shifted = cyclic_code->CyclicShift(codeword, 1);

    std::cout << "Original codeword: ";
    for (size_t i = 0; i < codeword.size(); ++i) {
        std::cout << static_cast<int>(codeword(i)) << " ";
    }
    std::cout << std::endl;

    std::cout << "Cyclically shifted: ";
    for (size_t i = 0; i < shifted.size(); ++i) {
        std::cout << static_cast<int>(shifted(i)) << " ";
    }
    std::cout << std::endl;

    // Check if shifted version is still a codeword
    bool is_codeword = cyclic_code->Contains(shifted);
    std::cout << "Shifted version is codeword: " << (is_codeword ? "Yes" : "No") << std::endl;

    std::cout << std::endl;
}

// Example 4: Code bounds and properties
void example_code_bounds() {
    std::cout << "=== Code Bounds and Properties ===" << std::endl;

    // Create various codes and check bounds
    auto field = std::make_shared<GaloisFieldBinary<uint8_t>>(4);

    // Reed-Solomon code
    auto rs_code = CreateReedSolomonCode(field, 15, 11);

    std::cout << "Reed-Solomon (15, 11) code:" << std::endl;
    std::cout << "  Satisfies Singleton bound: " <<
                 (SatisfiesSingletonBound(*rs_code) ? "Yes" : "No") << std::endl;
    std::cout << "  Satisfies Hamming bound: " <<
                 (SatisfiesHammingBound(*rs_code) ? "Yes" : "No") << std::endl;
    std::cout << "  Is perfect: " <<
                 (IsPerfectCode(*rs_code) ? "Yes" : "No") << std::endl;
    std::cout << "  Is MDS: " <<
                 (IsMDS(*rs_code) ? "Yes" : "No") << std::endl;

    std::cout << std::endl;
}

// Example 5: Encoder/Decoder comparison
void example_encoder_decoder_comparison() {
    std::cout << "=== Encoder/Decoder Comparison ===" << std::endl;

    auto field = std::make_shared<GaloisFieldBinary<uint8_t>>(3);
    auto rs_code = CreateReedSolomonCode(field, 7, 4);

    std::cout << "Available encoders: ";
    for (const auto& name : rs_code->AvailableEncoders()) {
        std::cout << name << " ";
    }
    std::cout << std::endl;

    std::cout << "Available decoders: ";
    for (const auto& name : rs_code->AvailableDecoders()) {
        std::cout << name << " ";
    }
    std::cout << std::endl;

    // Test different encoders
    xt::xarray<uint8_t> message = {1, 2, 3, 4};

    try {
        auto encoder1 = rs_code->GetEncoder("GRSEncoder");
        auto codeword1 = encoder1->Encode(message);

        std::cout << "GRS Encoder result: ";
        for (size_t i = 0; i < codeword1.size(); ++i) {
            std::cout << static_cast<int>(codeword1(i)) << " ";
        }
        std::cout << std::endl;

    } catch (const std::exception& e) {
        std::cout << "GRS Encoder error: " << e.what() << std::endl;
    }

    std::cout << std::endl;
}

int main() {
    std::cout << "XGalois Coding Theory Examples" << std::endl;
    std::cout << "==============================" << std::endl << std::endl;

    try {
        example_reed_solomon();
        example_generalized_reed_solomon();
        example_cyclic_code();
        example_code_bounds();
        example_encoder_decoder_comparison();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "All examples completed successfully!" << std::endl;
    return 0;
}
