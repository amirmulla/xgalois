.. _dev_guide:

Developer's Guide
=================

Thank you for your interest in contributing to **xgalois**! This guide covers
the repository architecture, C++20 and strict Google C++ Style conventions,
Bazel and CMake target structure, formatting/static-analysis tools, sanitizers,
and contribution workflow.

Project Structure
-----------------

.. code-block:: text

   xgalois/
   ├── MODULE.bazel                     # Pure Bzlmod dependency manifest (Bazel 8/9)
   ├── CMakeLists.txt                   # CMake >= 3.20 header-only build, CTest, and install
   ├── cmake/
   │   └── xgaloisConfig.cmake.in       # Relocatable find_package(xgalois CONFIG) template
   ├── .bazelrc                         # C++20, warnings, opt/native/dbg, sanitizers, coverage
   ├── .clang-format                    # Google-style formatting rules
   ├── .clang-tidy                      # Static analysis check configuration
   ├── .github/workflows/ci.yml         # CI: GCC/Clang builds, CMake, ASan/UBSan, format, docs
   ├── xgalois/                         # Header-only C++20 library source
   │   ├── BUILD.bazel                  # Granular cc_library targets + umbrella :xgalois
   │   ├── field/
   │   │   ├── gf_base.hpp              # Abstract GaloisFieldBase<T> + bulk span kernels
   │   │   ├── gf_concepts.hpp          # C++20 concepts (UnsignedFieldWord, FiniteField)
   │   │   ├── gf_element.hpp           # Value-type GaloisFieldElement<Field> wrapper
   │   │   ├── gf_array.hpp             # Contiguous GFArray<Field> / GFMatrix<Field>
   │   │   ├── gf_prime.hpp             # GF(p): Barrett GFP, table-backed GFPLOG
   │   │   ├── gf_binary.hpp            # GF(2) & GF(2^m): PCLMULQDQ/SSSE3/Log/Zech
   │   │   ├── gf128.hpp                # GF(2^128): UInt128, Karatsuba PCLMULQDQ, AES-GCM GHASH
   │   │   ├── gf_extension.hpp         # GF(p^m): ExtensionElement & GFPX
   │   │   └── gf_factory.hpp           # Runtime xg::GF(...) variant factory
   │   ├── poly/
   │   │   └── poly_dense.hpp           # PolynomialDense<Field> (Karatsuba, DivRem, PowMod)
   │   ├── utils/
   │   │   ├── math.hpp                 # Primality, factorization, safe integer powers
   │   │   ├── field.hpp                # Representation parsing and field helpers
   │   │   └── poly.hpp                 # GCD/ExtGCD, Rabin, Cantor–Zassenhaus, Interpolate
   │   ├── linalg/
   │   │   └── linalg.hpp               # Contiguous GFArray linalg + GF(2) M4RI + Berkowitz
   │   ├── channel/
   │   │   └── channel.hpp              # BSC, QSC, QEC, ErrorErasure, AWGN-LLR, Simulation
   │   ├── coding/
   │   │   ├── coding.hpp               # Umbrella header + factories & Singleton/Hamming bounds
   │   │   ├── abstract_code.hpp        # AbstractCode base + encoder/decoder registry
   │   │   ├── abstract_linear_code.hpp # AbstractLinearCode base (G, H, syndromes, weight)
   │   │   ├── decoding_result.hpp      # Structured DecodingResult<Field>
   │   │   ├── linear_code.hpp          # LinearCode + Hamming, Golay, Reed–Muller factories
   │   │   ├── cyclic_code.hpp          # CyclicCode
   │   │   ├── bch.hpp                  # BCHCode & BCHDecoder (binary and non-binary)
   │   │   ├── grs.hpp                  # GeneralizedReedSolomonCode (GRS / RS / Extended RS)
   │   │   ├── ldpc.hpp                 # SparseBinaryMatrix, LDPCCode, BitFlipping & BP decoders
   │   │   ├── polar.hpp                # PolarCode, PolarEncoder, SC & SCL/CA-SCL decoders
   │   │   ├── encoder/                 # Encoder, GeneratorMatrix, Cyclic, GRS, GRSSystematic
   │   │   └── decoder/                 # Decoder, Syndrome, CyclicSyndrome, BM, Sugiyama, PGZ,
   │   │                                # guruswami_sudan_decoder.hpp (Koetter + Roth–Ruckenstein)
   │   └── databases/
   │       ├── interface.hpp            # Conway & irreducible polynomial database interface
   │       ├── embedded_conway_polys.hpp
   │       ├── conway_polys.db
   │       └── irreducible_polys.db
   ├── test/                            # GoogleTest unit, differential, property, ODR & fuzz suites
   ├── benchmark/                       # Google Benchmark suites
   ├── example/                         # End-to-end runnable examples (8 binaries)
   └── doc/                             # Sphinx + Breathe + Exhale documentation

Granular Bazel Targets (``xgalois/BUILD.bazel``)
------------------------------------------------

Headers are partitioned into fine-grained ``cc_library`` targets to keep build
graphs minimal and avoid circular dependencies:

- ``//xgalois:utils_math`` — standalone number-theory helpers
- ``//xgalois:databases`` — Conway & irreducible polynomial databases
- ``//xgalois:field_core`` — ``gf_base.hpp``, ``gf_concepts.hpp``, ``gf_element.hpp``
- ``//xgalois:gf_array`` — contiguous ``GFArray`` / ``GFMatrix``
- ``//xgalois:utils_field`` — representation and field utility helpers
- ``//xgalois:poly`` — ``PolynomialDense``
- ``//xgalois:utils_poly`` — GCD, Rabin, Cantor–Zassenhaus, root finding, interpolation
- ``//xgalois:utils`` — aggregates ``:utils_math``, ``:utils_field``, ``:utils_poly``
- ``//xgalois:field`` — all concrete field implementations (including ``gf128.hpp``) and ``gf_factory.hpp``
- ``//xgalois:linalg`` — finite-field linear algebra and ``GF(2)`` M4RI engine
- ``//xgalois:channel`` — channel models and ``SimulateCodePerformance``
- ``//xgalois:coding`` — linear, cyclic, BCH, GRS/RS, Guruswami–Sudan, LDPC, and Polar codes, encoders, and decoders
- ``//xgalois:xgalois`` — umbrella target depending on all modules

CMake exports a single header-only target ``xgalois::xgalois`` (declared in
``CMakeLists.txt``) with ``cxx_std_20`` and automatic ``xtl``/``xtensor``
integration.

Coding Conventions
------------------

- **Language standard**: **C++20** (``-std=c++20`` configured in ``.bazelrc`` and ``cxx_std_20`` in ``CMakeLists.txt``).
- **Warnings**: All files under ``xgalois/`` are compiled with
  ``-Wall -Wextra -Wpedantic``; code must compile warning-free on both GCC 11+
  and Clang 14+.
- **Strict Google C++ Style Naming across ALL modules**:

  - Classes / structs / type aliases: ``PascalCase`` (e.g., ``GaloisFieldPrime``, ``GaloisFieldBinary128``, ``GFArray``, ``DecodingResult``, ``PluResult``)
  - Methods & free functions across **every** module (``field``, ``poly``, ``utils``, ``linalg``, ``channel``, ``coding``, ``databases``): ``PascalCase`` — e.g., ``Add()``, ``DivRem()``, ``DecodeWithErasures()``, ``linalg::Dot()``, ``linalg::Solve()``, ``linalg::Rref()``, ``linalg::Plu()``, ``linalg::CharacteristicPoly()``, ``channel.TransmitWithErasures()``, ``channel.TransmitLlr()``, ``stats.Fer()``, ``db.Fetch()``
  - Variables, struct fields, & parameters: ``snake_case`` (e.g., ``evaluation_points``, ``erasure_positions``, ``plu_res.p``, ``plu_res.l``, ``plu_res.u``); private data members use a trailing underscore (``field_``)
  - Constants & enumerators: ``kPascalCase`` (e.g., ``kLogZero``, ``kKaratsubaThreshold``, ``BpAlgorithm::kNormalizedMinSum``)
  - Namespaces: ``lowercase`` (``xg``, ``xg::utils``, ``xg::linalg``, ``xg::channel``, ``xg::coding``, ``xg::databases``)

- **Headers**: Use ``.hpp`` extension for all C++ headers.
- **Include guards**: Follow the trailing-underscore Google convention
  ``XGALOIS_<MODULE>_<FILE>_HPP_`` (for example,
  ``#ifndef XGALOIS_FIELD_GF_ARRAY_HPP_``).

Formatting & Static Analysis
----------------------------

Code formatting is enforced in CI via ``.clang-format`` (Google style):

.. code-block:: bash

   # Check formatting across library, tests, and examples
   find xgalois test example -type f \( -name '*.hpp' -o -name '*.cpp' \) \
     -exec clang-format --dry-run --Werror {} +

   # Auto-format modified files in place
   clang-format -i path/to/file.hpp path/to/file.cpp

Static analysis checks are configured in ``.clang-tidy`` (covering
``bugprone-*``, ``clang-analyzer-*``, ``google-*``, ``modernize-*``,
``performance-*``, and ``readability-*`` checks on ``xgalois/`` headers).

Testing & Sanitizers
--------------------

.. code-block:: bash

   # Run every test target (Bazel)
   bazel test //test/...

   # Run the full aggregated test suite
   bazel test //test:galois_field_tests --test_output=errors

   # Run an individual test target (e.g. prime fields, GF128, GRS, Guruswami-Sudan, LDPC/Polar)
   bazel test //test:gf_prime_test //test:gf128_test //test:grs_test //test:guruswami_sudan_test //test:ldpc_polar_test

   # Run with AddressSanitizer, UndefinedBehaviorSanitizer, or ThreadSanitizer
   bazel test --config=asan //test:galois_field_tests --test_output=errors
   bazel test --config=ubsan //test:galois_field_tests --test_output=errors
   bazel test --config=tsan //test:galois_field_tests --test_output=errors

   # Generate LCOV coverage report
   bazel coverage --config=coverage //test:galois_field_tests

   # Build and run all tests via CMake / CTest
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ctest --test-dir build --output-on-failure

Contribution Workflow
---------------------

1. Fork or branch from ``main``.
2. Add or update headers under ``xgalois/<module>/`` and register any new
   targets or dependencies in ``xgalois/BUILD.bazel`` (and ``CMakeLists.txt`` /
   ``test/CMakeLists.txt`` if adding tests or examples).
3. Add comprehensive unit and property tests under ``test/`` (updating both
   ``test/BUILD.bazel`` and ``test/CMakeLists.txt``).
4. Format all touched C++ files with ``clang-format -i``.
5. Verify both Bazel and CMake builds, tests, examples, and sanitizers pass cleanly:

   .. code-block:: bash

      bazel build //xgalois/... //example/...
      bazel test //test/...
      bazel test --config=asan //test:galois_field_tests
      bazel test --config=ubsan //test:galois_field_tests
      cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && ctest --test-dir build --output-on-failure

6. If public APIs or guides changed, rebuild documentation with
   ``cd doc && ./build.sh`` (which runs ``sphinx-build -W``) and ensure zero
   warnings.
7. Open a pull request describing the motivation, mathematical/algorithmic
   details, and test coverage.
