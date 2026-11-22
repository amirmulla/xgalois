# XGalois (v1.0.0)

**XGalois** is a modern, header-only C++20 computational algebra and coding-theory library for high-performance finite-field arithmetic, abstract algebra, rings, elliptic curves, fast polynomial & linear algebra, communication channel modeling, and classical & modern error-correcting codes over arbitrary Galois fields $\mathrm{GF}(p^m)$.

## Feature Matrix

| Domain | Modules | Highlights |
| --- | --- | --- |
| **Algebraic Structures & Groups** | `xgalois/algebra/` | C++20 algebraic concepts (`Magma`..`Field`, `Module`, `VectorSpace`, `Algebra`), `Permutation`, deterministic & randomized Schreier–Sims `BsgsChain`, `CyclicGroup`, `SymmetricGroup`, `AlternatingGroup`, `DihedralGroup`, `MatrixGroup` ($\mathrm{GL}_n, \mathrm{SL}_n$), `DirectProductGroup`, `SemidirectProductGroup`, `Subgroup`, `QuotientGroup`, `SylowSubgroup`, `GroupHomomorphism`, `SmithNormalForm` Abelian decomposition. |
| **Rings, Ideals & Skew Polynomials** | `xgalois/ring/` | `IntegerRing`, `ModularRing`, `GaussianIntegerRing`, `EisensteinIntegerRing`, `PolynomialRing`, `QuotientRing`, `NegacyclicRing` ($\mathbb{Z}_q[x]/(x^N+1)$), `CyclicRing`, `MatrixRing` (division-free Berkowitz det), `SkewPolynomialRing` / `FrobeniusSkewRing`, `PrincipalIdeal`, generalized CRT. |
| **Finite Fields, Towers & Elliptic Curves** | `xgalois/field/`, `xgalois/ec/` | Prime (`GFP`, `GFPLOG`), Binary (`GF2`, `GF2X`, `GF2XLOG`, `GF2XZECH`), dedicated `GF128` ($\mathrm{GF}(2^{128})$ / AES-GCM GHASH), Extension (`GFPX`), Recursive `FieldTower`, Canright `CompositeGF256Tower`, `GaloisGroup` / `FieldIsomorphism`, ONB Type I & II `GaloisFieldNormalBinary`, `GaloisFieldDualBasis`, `GaloisFieldPrimeBig<Bits>` / `GaloisFieldPrimeDyn` (`xg::primes`), `GaloisFieldBinaryLarge<M>`, `GFArray`/`GFMatrix`, `EllipticCurve` (Weierstrass, Binary, Montgomery, Edwards), `WeilPairing`, `TatePairing`, `EllipticRiemannRochMatrix`. |
| **Fast Polynomials, Linear Algebra & SIMD** | `xgalois/poly/`, `xgalois/linalg/`, `xgalois/simd/` | Shoup `NttPlan`, 3-prime `CrtNttPlan`, `BluesteinDft`, Lin–Chung–Han `AdditiveFftPlan`, Möller `HalfGcd`, `SubproductTree`, `PolynomialDense`, `PolynomialSparse`, PID `HermiteNormalForm` / `SmithNormalForm`, `FrobeniusNormalForm`, `BlockWiedemannNullspace`, `BlockLanczosNullspaceBinary`, `M4riMatMul` / `M4riRref`, `StrassenMatMul`, runtime-dispatched AVX2/PCLMULQDQ, AVX-512 GFNI/VPCLMULQDQ, ARM64 NEON/PMULL, and `std::jthread` `ThreadPool`. |
| **Classical & Modern Coding Theory** | `xgalois/coding/` | `LinearCode` (Hamming, Golay, Reed–Muller), `CyclicCode`, `BCHCode`, GRS/RS (`BerlekampMassey`, `Sugiyama` + Half-GCD, **Guruswami–Sudan** list decoding, **Power Decoding**), **Tamo–Barg LRC** (single-set, $(r,\delta)$, multi-set availability), **Interleaved RS** (SSB, Kötter–Nielsen, BKY Hankel), **Gabidulin** & **Kötter–Kschischang Subspace**, **LDPC** (Bit-Flipping & BP Sum-Product/Min-Sum), **Polar** (SC/SCL/CA-SCL), **Convolutional** (Viterbi/BCJR), **3GPP LTE Turbo**, **Product** (Hard & Chase–Pyndiah TPC), **Concatenated** (Forney GMD), **Hermitian AG** (Feng–Rao + BMS & Sudan), **LT & RaptorQ Fountain** (inactivation), **Code Modifiers** (`PlotkinCode`, puncturing, shortening, extending, expurgating), and `BigInt` **MacWilliams** weight enumerator. |
| **Channels, Simulation, Python & Oracle** | `xgalois/channel/`, `python/`, `test/oracle/` | AWGN, BSC, QSC, BEC/QEC, Error–Erasure, Gilbert–Elliott Burst, Gray BPSK/QPSK/8-PSK/16-QAM/64-QAM, Rayleigh/Rician Fading, BAC/Z-Channel, IDS (Insertion/Deletion/Substitution), deterministic parallel BER/SER/FER sweeps, `nanobind` Python/NumPy/SageMath interop, and 2-tier SageMath 10.5 differential verification. |

## Requirements

- **Compiler**: C++20-capable compiler (`GCC 11+` or `Clang 14+`)
- **Build System**:
  - [Bazel](https://bazel.build/) 8 or 9 (pure Bzlmod via `MODULE.bazel`), **or**
  - [CMake](https://cmake.org/) `>= 3.20` (`FetchContent` or system `xtl` / `xtensor`)
- **Python bindings (optional)**: Python 3.10+, `scikit-build-core`, `nanobind`, `numpy`
- **Live Differential Oracle (optional)**: SageMath 10.5 (`sagemath/sagemath:10.5`)
- **Documentation (optional)**: Doxygen and Python 3.10+ (`sphinx`, `furo`, `breathe`, `exhale`)

## C++20 Quickstart

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

## Quick-Start Builds & Workflows

### Bazel 9 (Pure Bzlmod)

```bash
# Build the library and all examples
bazel build --cxxopt=-std=c++20 //xgalois/... //example/...

# Run all unit tests (including Tier-1 SageMath golden JSON tests)
bazel test --cxxopt=-std=c++20 //test:galois_field_tests --test_output=errors

# Run with ASan / UBSan / TSan or native SIMD benchmarks
bazel test --config=asan //test:galois_field_tests
bazel run --config=native //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark
```

### CMake (`>= 3.20`) & `find_package(xgalois CONFIG REQUIRED)`

```bash
# Configure, build, test, and install (add -DXGALOIS_ENABLE_NATIVE=ON for -march=native)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure --parallel "$(nproc)"
cmake --install build --prefix /usr/local
```

Consume an installed `xgalois` from any downstream CMake project:

```cmake
cmake_minimum_required(VERSION 3.20)
project(my_app LANGUAGES CXX)

find_package(xgalois CONFIG REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE xgalois::xgalois)
```

Or pull directly via `FetchContent`:

```cmake
include(FetchContent)
FetchContent_Declare(
  xgalois
  GIT_REPOSITORY https://github.com/amirmulla/xgalois.git
  GIT_TAG        v1.0.0
)
FetchContent_MakeAvailable(xgalois)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE xgalois::xgalois)
```

### Python (`nanobind` / NumPy / SageMath Interop)

```bash
# Install the Python package and compiled C++20 extension from source
pip install .
```

### SageMath 10.5 Live Differential Oracle (2-Tier Verification)

XGalois verifies mathematical correctness against **SageMath 10.5** via two tiers:
- **Tier 1 (Offline Golden JSON)**: Deterministic golden vectors under `test/oracle/golden/` executed in standard C++ unit tests (`sage_golden_test`).
- **Tier 2 (Live Differential Oracle)**: Builds `xgalois_oracle_cli` and runs randomized differential checks directly against SageMath 10.5 (also executed in CI via `sagemath/sagemath:10.5`).

```bash
# 1. Build the oracle CLI (ABI-safe standalone target)
cmake -S . -B build-oracle \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_STANDARD=20 \
  -DXGALOIS_BUILD_TESTS=OFF \
  -DXGALOIS_BUILD_EXAMPLES=OFF \
  -DXGALOIS_BUILD_ORACLE=ON
cmake --build build-oracle --target xgalois_oracle_cli --parallel 2

# 2. Verify golden JSON freshness (zero drift)
sage -python test/oracle/sage/generate_golden.py \
  --seed=2026 --check-dir=test/oracle/golden

# 3. Run the live differential oracle suite
sage -python test/oracle/sage/live_oracle.py \
  --cli=build-oracle/xgalois_oracle_cli \
  --seed=2026 --trials=250
```

### Building Documentation

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
│   ├── algebra/     # Algebraic concepts, Permutation, BsgsChain, Cyclic/Symmetric/Alternating/Dihedral/Matrix/Product/Quotient/Sylow groups
│   ├── ring/        # Integer, Modular, Gaussian/Eisenstein, Polynomial, Quotient, Negacyclic, Cyclic, Matrix, Skew rings & ideals
│   ├── field/       # Prime, Binary, GF(2^128), Extension, FieldTower, CompositeGF256Tower, ONB/Dual basis, Big/Dyn primes, GFArray
│   ├── ec/          # EllipticCurve (Weierstrass/Binary/Montgomery/Edwards), WeilPairing, TatePairing, EllipticRiemannRochMatrix
│   ├── poly/        # Dense & Sparse polynomials, Shoup NTT, 3-prime CRT-NTT, Bluestein DFT, Additive FFT, Half-GCD, SubproductTree
│   ├── linalg/      # GFArray linear algebra, M4RI, Strassen, PID HNF/SNF, Frobenius normal form, Block Wiedemann & Block Lanczos
│   ├── simd/        # Runtime-dispatched AVX2/PCLMULQDQ, AVX-512 GFNI/VPCLMULQDQ, ARM64 NEON/PMULL, std::jthread ThreadPool
│   ├── coding/      # Linear, Cyclic, BCH, GRS/RS (BM/Sugiyama/GS/Power), LRC, IRS, Gabidulin/Subspace, LDPC, Polar, Conv, Turbo, TPC, GMD, AG, Fountain, Modifiers, MacWilliams
│   ├── channel/     # AWGN, BSC, QSC, BEC, Gilbert-Elliott, Gray PSK/QAM, Rayleigh/Rician, BAC/Z-Channel, IDS, parallel simulation
│   ├── utils/       # Number theory, BigInt, polynomial factorization, interpolation, field & coding helpers
│   ├── databases/   # Embedded Conway and irreducible polynomial databases
│   └── BUILD.bazel  # Granular Bazel targets (:field_core, :gf_array, :field, :poly, :linalg, :channel, :coding, :xgalois, ...)
├── python/          # nanobind C++20 Python bindings with NumPy & SageMath interop
├── test/            # GoogleTest unit, algebraic-property, ODR, fuzz, and 2-tier SageMath 10.5 oracle suites (test/oracle/)
├── benchmark/       # Google Benchmark suites
├── example/         # End-to-end runnable C++20 examples
├── cmake/           # CMake package configuration template (xgaloisConfig.cmake.in)
├── doc/             # Sphinx + Breathe + Exhale API documentation
├── CMakeLists.txt   # CMake >= 3.20 header-only build, CTest, FetchContent, oracle CLI, and install config
├── MODULE.bazel     # Pure Bzlmod dependency manifest (Bazel 8/9)
└── .bazelrc         # C++20, strict warnings, opt/native/dbg, ASan/UBSan/TSan, and coverage configs
```

## License

Licensed under the [MIT License](LICENSE).
