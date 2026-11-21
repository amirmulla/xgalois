# xgalois

**xgalois** is a modern, header-only C++20 library for high-performance finite-field arithmetic, polynomial algebra, finite-field linear algebra, communication channel modeling, and algebraic & modern error-correcting codes over arbitrary Galois fields $\mathrm{GF}(p^m)$.

## Key Capabilities

- **Finite Fields (`xg::GaloisFieldBase`, `xg::GF`)**
  - **Binary Field $\mathrm{GF}(2)$** (`GF2`): specialized bitwise arithmetic and bulk array kernels.
  - **Prime Fields $\mathrm{GF}(p)$** (`GFP`, `GFPLOG`): division-free Barrett reduction up to $p \le 2^{63} - 1$, Tonelli–Shanks square roots, Baby-Step Giant-Step (BSGS) and Pohlig–Hellman discrete logarithms, and log-table acceleration (`GFPLOG`).
  - **Binary Extension Fields $\mathrm{GF}(2^m)$** (`GF2X`, `GF2XLOG`, `GF2XZECH`, `GF128`): hardware-accelerated carry-less multiplication (`PCLMULQDQ` / `SSSE3` with portable fallback) up to $\mathrm{GF}(2^{64})$ and dedicated $\mathrm{GF}(2^{128})$ (`xgalois/field/gf128.hpp`, `UInt128`, Karatsuba `PCLMULQDQ`, AES-GCM NIST SP 800-38D `GhashMul` modulo $x^{128} + x^7 + x^2 + x + 1$, and general carry-less Barrett reduction); optimized log/antilog tables, Zech logarithm tables, and built-in embedded Conway and irreducible polynomial databases.
  - **General Prime Extension Fields $\mathrm{GF}(p^m)$** (`GFPX`): compact `ExtensionElement` polynomial representation with Frobenius, Trace, Norm, Minimal Polynomial, BSGS discrete log, and Tonelli–Shanks square root.
  - **Unified Runtime Factory** (`xg::GF(...)`) returning a type-safe `std::variant` (`GaloisFieldVariant`), plus zero-overhead direct constructors and C++20 concepts (`FiniteField`, `UnsignedFieldWord`).
- **Contiguous Field Arrays (`xg::GFArray`, `xg::GFMatrix`)**
  - Stores a single `std::shared_ptr<GaloisField>` alongside a packed `xt::xarray<storage_type>` buffer, avoiding per-element reference-count overhead and unlocking bulk span/SIMD kernels (`AddArray`, `SubArray`, `MulArray`, `ScalarMulArray`, `FmaArray`, `DotProduct`, Montgomery `BatchInv`).
- **Polynomial Algebra (`xg::PolynomialDense`, `xg::utils`)**
  - Dense polynomials with Karatsuba multiplication, Euclidean `DivRem`, `PowMod`, `Derivative`, Horner/batch evaluation (`EvalBatch`), and composition.
  - Exact algorithms over arbitrary $\mathrm{GF}(q)$: `PolynomialGcd`, `PolynomialLcm`, `PolynomialDenseExtendedGcd`, Rabin irreducibility (`IsIrreducible`) and primitivity (`IsPrimitive`) tests, square-free / distinct-degree / Cantor–Zassenhaus complete factorization (`Factor`), `FindRoots`, `ChienSearch`, and $O(k^2)$ Newton interpolation (`Interpolate`).
- **Finite-Field Linear Algebra (`xg::linalg`)**
  - Bit-packed `uint64_t` **M4RI** engine for $\mathrm{GF}(2)$ matrices alongside cache-friendly FMA kernels for arbitrary fields.
  - Strict Google C++ Style `PascalCase` API: `Dot`, `Vdot`, `Outer`, `Cross`, `Kron`, `MatrixPower`, `RowEchelon`, `Rref`, `MatrixRank`, `Trace`, `Det`, `Inv`, `Solve`, `Nullspace`, `LeftNullspace`, `RowSpace`, `ColumnSpace`, `SystematicForm`, `Plu`, and division-free Berkowitz `CharacteristicPoly`.
- **Coding Theory (`xg::coding`)**
  - `LinearCode` (generator/parity-check duality, puncturing, shortening) with classic factories: `HammingCode`, `ExtendedHammingCode`, `BinaryGolayCode`, `ExtendedBinaryGolayCode`, `TernaryGolayCode`, and `ReedMullerCode`.
  - `CyclicCode` and `BCHCode` (binary and non-binary, narrow- and general-sense, with algebraic Berlekamp–Massey / Chien / Forney decoding).
  - `GeneralizedReedSolomonCode` (`ReedSolomon`, `ExtendedReedSolomon`, general GRS) with `BerlekampMassey`, `Sugiyama`, and `PGZ` decoders supporting full **Errors-and-Erasures** decoding up to $2e + s \le n - k$, plus **Guruswami–Sudan list decoding** (`xgalois/coding/decoder/guruswami_sudan_decoder.hpp`, `GuruswamiSudanDecoder`, `FromTargetRadius`, `DecodeToList`, `JohnsonRadius`) decoding beyond $\lfloor (d - 1) / 2 \rfloor$ all the way up to the Johnson bound $n - 1 - \lfloor \sqrt{n(k - 1)} \rfloor$.
  - **LDPC Codes** (`xgalois/coding/ldpc.hpp`): `SparseBinaryMatrix` (CSR/CSC Tanner graph), `LDPCCode` (`FromQuasiCyclic`, `MakeRegularLDPC`), hard-decision `BitFlippingDecoder`, and soft-decision `BeliefPropagationDecoder` supporting Sum-Product, Normalized Min-Sum (`BpAlgorithm::kNormalizedMinSum`), and Offset Min-Sum (`BpAlgorithm::kOffsetMinSum`) over channel LLRs (`DecodeLlr`).
  - **Polar Codes** (`xgalois/coding/polar.hpp`): `PolarCode` (`FromDesignSnr` Gaussian Approximation and `FromBhattacharyya` construction, systematic & non-systematic, optional CRC-aided CA-Polar), $O(N \log N)$ `PolarEncoder`, `SuccessiveCancellationDecoder` (SC), and `SuccessiveCancellationListDecoder` (SCL / CA-SCL).
  - Exact arbitrary-precision Singleton/MDS and Hamming/perfection bound checks (`IsMDS`, `SatisfiesHammingBound`, `IsPerfectCode`).
- **Channels & Simulation (`xg::channel`)**
  - `BinarySymmetricChannel` (BSC), `QarySymmetricChannel` (QSC), `QaryErasureChannel` / `BinaryErasureChannel` (QEC/BEC), `ErrorErasureChannel`, `StaticErrorRateChannel`, and BPSK `AWGNChannel` with exact LLR computation (`Transmit`, `TransmitWithErasures`, `TransmitSoft`, `TransmitLlr`, `TransmitDetailed`).
  - Deterministic, seeded Monte Carlo BER/SER/FER evaluation via `SimulateCodePerformance` (`SimulationStats::Fer()`, `Ser()`, `Ber()`).

## Requirements

- **Compiler**: C++20-capable compiler (`GCC 11+` or `Clang 14+`)
- **Build System**:
  - [Bazel](https://bazel.build/) 8 or 9 (pure Bzlmod via `MODULE.bazel`), **or**
  - [CMake](https://cmake.org/) `>= 3.20` (`FetchContent` or system `xtl` / `xtensor`)
- **Documentation (optional)**: Doxygen and Python 3.10+ (`sphinx`, `furo`, `breathe`, `exhale`)

## Quickstart

```cpp
#include <iostream>
#include <memory>
#include <random>

#include "xgalois/channel/channel.hpp"
#include "xgalois/coding/coding.hpp"
#include "xgalois/field/gf_array.hpp"
#include "xgalois/field/gf_binary.hpp"
#include "xgalois/linalg/linalg.hpp"

int main() {
  using Field = xg::GF2XLOG<uint8_t>;  // GF(2^8) with optimized log tables
  auto gf256 = std::make_shared<Field>(8);

  // Classic CCSDS Reed-Solomon RS(255, 223) code (d = 33, corrects up to 16 errors)
  auto rs = xg::coding::GeneralizedReedSolomonCode<Field>::ReedSolomon(
      gf256, /*length=*/255, /*dimension=*/223);

  // Encode a random 223-symbol message and verify parity checks H * c^T == 0
  std::mt19937_64 rng(42);
  auto msg = xg::GFArray<Field>::Random(gf256, {223}, rng);
  auto codeword = rs.Encode(msg);
  auto syndrome = xg::linalg::Dot(rs.ParityCheckMatrix(), codeword);

  // Transmit through an Error-and-Erasure channel (10 errors + 12 erasures; 2*10 + 12 = 32 <= 32)
  xg::channel::ErrorErasureChannel<Field> ch(gf256, 255, 10, 12, /*seed=*/7);
  auto tx = ch.TransmitWithErasures(codeword);

  // Decode with Berlekamp-Massey + Chien search + Forney (errors & erasures)
  auto res = rs.DecodeWithErasures(tx.received, tx.erasure_positions);
  std::cout << "Valid codeword: " << std::boolalpha
            << (syndrome == xg::GFArray<Field>::Zeros(gf256, {32}))
            << ", Decoded OK: " << res.Ok()
            << " (errors=" << res.num_errors
            << ", erasures=" << res.num_erasures
            << ", matches=" << (res.message == msg) << ")\n";
  return 0;
}
```

## Build, Test, Install, Benchmarks & Docs

### Building with Bazel (Pure Bzlmod, Bazel 8/9)

```bash
# Build library and all examples
bazel build //xgalois/... //example/...

# Run all unit tests
bazel test //test/...

# Run the aggregated test suite with error output
bazel test //test:galois_field_tests --test_output=errors

# Run an example binary
bazel run //example:coding_example
```

### Building & Installing with CMake (`>= 3.20`)

CMake automatically uses system-installed `xtl` / `xtensor` if found, or fetches pinned releases via `FetchContent`:

```bash
# Configure and build (add -DXGALOIS_ENABLE_NATIVE=ON for -march=native / PCLMULQDQ)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run unit tests via CTest
ctest --test-dir build --output-on-failure

# Run an example binary
./build/example/coding_example

# Install headers, databases, and CMake package config
cmake --install build --prefix /usr/local
```

**Consuming an installed `xgalois` via `find_package`:**

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_app LANGUAGES CXX)

find_package(xgalois CONFIG REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE xgalois::xgalois)
```

**Consuming directly via `FetchContent`:**

```cmake
include(FetchContent)
FetchContent_Declare(
  xgalois
  GIT_REPOSITORY https://github.com/amirmulla/xgalois.git
  GIT_TAG        main
)
FetchContent_MakeAvailable(xgalois)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE xgalois::xgalois)
```

### Sanitizers & Coverage (Bazel)

Configured in `.bazelrc`:

```bash
# AddressSanitizer
bazel test --config=asan //test:galois_field_tests

# UndefinedBehaviorSanitizer
bazel test --config=ubsan //test:galois_field_tests

# ThreadSanitizer
bazel test --config=tsan //test:galois_field_tests

# Code coverage (LCOV)
bazel coverage --config=coverage //test:galois_field_tests
```

### Native Benchmarks

Use `--config=native` (`-c opt -O3 -march=native`) to enable `PCLMULQDQ` / `SSSE3` / native SIMD tuning:

```bash
bazel run --config=native //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark
bazel run --config=native //benchmark/gf_prime/src:gf_prime_benchmark
```

### Building Documentation

Requires Doxygen and Python 3.10+:

```bash
cd doc
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
./build.sh
```

The generated HTML documentation is written to `doc/_build/html/index.html`.

## Architecture

```text
xgalois/
├── xgalois/
│   ├── field/       # gf_base, gf_concepts, gf_element, gf_array, gf_prime, gf_binary, gf128, gf_extension, gf_factory
│   ├── poly/        # PolynomialDense (Karatsuba, DivRem, PowMod, Derivative, EvalBatch)
│   ├── utils/       # math (primality, factorization), field helpers, poly (GCD, Rabin, Cantor-Zassenhaus, Interpolate)
│   ├── linalg/      # Contiguous GFArray linear algebra + bit-packed GF(2) M4RI engine (PascalCase Google Style)
│   ├── channel/     # BSC, QSC, QEC/BEC, ErrorErasure, StaticErrorRate, AWGN-LLR, SimulateCodePerformance
│   ├── coding/      # LinearCode, Hamming, Golay, ReedMuller, CyclicCode, BCHCode, GRS/RS, Guruswami-Sudan, LDPC, Polar
│   ├── databases/   # Embedded Conway and irreducible polynomial databases
│   └── BUILD.bazel  # Granular targets (:field_core, :gf_array, :field, :poly, :linalg, :channel, :coding, :xgalois)
├── cmake/           # CMake package configuration template (xgaloisConfig.cmake.in)
├── test/            # Unit, differential, algebraic-property, ODR, DB-integrity, and fuzz tests (GoogleTest)
├── benchmark/       # Google Benchmark suites for prime, binary extension, and Zech implementations
├── example/         # End-to-end runnable examples for every module (8 binaries)
├── doc/             # Sphinx + Breathe + Exhale API documentation
├── CMakeLists.txt   # CMake >= 3.20 header-only build, CTest, FetchContent, and install config
├── MODULE.bazel     # Pure Bzlmod dependency manifest (Bazel 8/9)
└── .bazelrc         # C++20, strict warnings, opt/native/dbg, ASan/UBSan/TSan, and coverage configs
```

## License

Licensed under the [MIT License](LICENSE).
