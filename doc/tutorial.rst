.. _tutorial:

Tutorial
========

This tutorial walks through the core workflows in **xgalois v0.2**. Every code
listing below is a complete, self-contained **C++20** program targeting the
strict Google-style ``PascalCase`` v0.2 API.

By the end of this guide, you will know how to:

- Construct prime, binary extension (up to :math:`\mathrm{GF}(2^{128})`), and general extension fields (directly or via ``xg::GF``)
- Perform contiguous array and linear algebra operations with ``xg::GFArray`` and ``xg::linalg``
- Manipulate, factor, find roots of, and interpolate polynomials with ``xg::PolynomialDense`` and ``xg::utils``
- Encode and decode Reed–Solomon codes with errors, erasures, and Guruswami–Sudan list decoding
- Encode and soft-decode LDPC and Polar codes over AWGN channels
- Run deterministic Monte Carlo BER/SER/FER simulations with ``SimulateCodePerformance``

Prerequisites
-------------

- A **C++20**-capable compiler (``GCC 11+`` or ``Clang 14+``)
- `Bazel <https://bazel.build/>`_ 8 or 9 (pure Bzlmod) **or** `CMake <https://cmake.org/>`_ ``>= 3.20``
- A local clone of the xgalois repository

Creating Finite Fields
----------------------

xgalois supports five field families:

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Family
     - Direct Class / Alias
     - Example
   * - Binary field :math:`\mathrm{GF}(2)`
     - ``GaloisFieldBinary`` / ``GF2``
     - :math:`\mathrm{GF}(2)`
   * - Prime field :math:`\mathrm{GF}(p)`
     - ``GaloisFieldPrime<T>`` / ``GFP<T>``, ``GFPLOG<T>``
     - :math:`\mathrm{GF}(7)`, :math:`\mathrm{GF}(2^{31} - 1)`
   * - Binary extension :math:`\mathrm{GF}(2^m)`
     - ``GF2X<T>``, ``GF2XLOG<T>``, ``GF2XZECH<T>``
     - :math:`\mathrm{GF}(2^8)`, :math:`\mathrm{GF}(2^{64})`
   * - 128-bit binary extension :math:`\mathrm{GF}(2^{128})`
     - ``GaloisFieldBinary128`` / ``GF128``
     - AES-GCM / GHASH (:math:`x^{128} + x^7 + x^2 + x + 1`)
   * - General extension :math:`\mathrm{GF}(p^m)`
     - ``GaloisFieldExtension<T>`` / ``GFPX<T>``
     - :math:`\mathrm{GF}(3^2)`, :math:`\mathrm{GF}(5^3)`

**Direct construction (zero-overhead static types):**

.. code-block:: cpp

   #include <iostream>
   #include <memory>

   #include "xgalois/field/gf_binary.hpp"
   #include "xgalois/field/gf_element.hpp"
   #include "xgalois/field/gf_prime.hpp"

   int main() {
     using PrimeField = xg::GaloisFieldPrime<uint8_t>;
     auto gf7 = std::make_shared<PrimeField>(7);

     xg::GaloisFieldElement<PrimeField> a(3, gf7);
     xg::GaloisFieldElement<PrimeField> b(5, gf7);
     auto c = a + b;                                                    // (3 + 5) mod 7 = 1
     auto inv_a = a.Inv();                                              // 3^-1 mod 7 = 5
     auto sqrt_2 = xg::GaloisFieldElement<PrimeField>(2, gf7).Sqrt();   // 3 (since 3^2 == 2 mod 7)

     std::cout << "3 + 5 = " << c << ", 3^-1 = " << inv_a
               << ", sqrt(2) = " << sqrt_2 << '\n';

     // Binary extension field GF(2^8) with optimized log/antilog tables
     auto gf256 = std::make_shared<xg::GF2XLOG<uint8_t>>(8);
     uint8_t prod = gf256->Mul(0x57, 0x83);
     std::cout << "GF(2^8) order = " << gf256->Order()
               << ", 0x57 * 0x83 = " << static_cast<int>(prod) << '\n';
     return 0;
   }

**Runtime construction with the ``xg::GF(...)`` factory:**

When field parameters are chosen at runtime, use ``xg::GF(...)`` to obtain a
``xg::GaloisFieldVariant`` and dispatch with ``std::visit``:

.. code-block:: cpp

   #include <iostream>
   #include <utility>
   #include <variant>

   #include "xgalois/field/gf_factory.hpp"

   int main() {
     xg::GaloisFieldVariant gf256 =
         xg::GF(std::pair<uint64_t, uint64_t>{2, 8}, "int", "", "α", "log_opt");
     xg::GaloisFieldVariant gf9 = xg::GF(9, "poly");

     std::visit(
         [](const auto& field) {
           std::cout << "Order=" << field->Order()
                     << ", Characteristic=" << field->Characteristic()
                     << ", Degree=" << field->Degree() << '\n';
         },
         gf256);

     auto elem = xg::FetchElement(gf9, 5);
     std::visit(
         [](const auto& e) { std::cout << "GF(9) element(5) = " << e << '\n'; },
         elem);
     return 0;
   }

Contiguous Arrays & Linear Algebra
----------------------------------

``xg::GFArray<Field>`` (aliased as ``xg::GFMatrix<Field>``) packs raw field
words contiguously alongside a single shared field pointer, enabling fast bulk
array kernels and Google-style ``PascalCase`` linear algebra in ``xg::linalg``:

.. code-block:: cpp

   #include <iostream>
   #include <memory>

   #include "xgalois/field/gf_array.hpp"
   #include "xgalois/field/gf_prime.hpp"
   #include "xgalois/linalg/linalg.hpp"

   int main() {
     using Field = xg::GaloisFieldPrime<uint8_t>;
     auto gf7 = std::make_shared<Field>(7);

     // Invertible 3x3 system A * x = b over GF(7)
     xg::GFArray<Field> A(gf7, xt::xarray<uint8_t>{{1, 2, 3},
                                                   {2, 5, 1},
                                                   {3, 1, 4}});
     xg::GFArray<Field> b(gf7, xt::xarray<uint8_t>{4, 6, 2});

     auto x = xg::linalg::Solve(A, b);
     auto check = xg::linalg::Dot(A, x);  // A * x == b

     std::cout << "Det(A) = " << xg::linalg::Det(A)
               << ", A * x == b: " << std::boolalpha << (check == b) << '\n';

     // Nullspace of a 2x4 matrix M over GF(7)
     xg::GFArray<Field> M(gf7, xt::xarray<uint8_t>{{1, 2, 3, 4},
                                                   {2, 4, 1, 3}});
     auto N = xg::linalg::Nullspace(M);  // shape: (2, 4), satisfies M * N^T == 0
     auto zero_check = xg::linalg::Dot(M, xg::linalg::Transpose(N));
     std::cout << "Nullspace dimension: " << N.Shape(0)
               << ", M * N^T == 0: "
               << (zero_check == xg::GFArray<Field>::Zeros(gf7, {2, 2})) << '\n';
     return 0;
   }

Polynomial Arithmetic & Algorithms
----------------------------------

``xg::PolynomialDense<Field>`` stores dense coefficients in ascending degree
order and pairs with exact algebraic algorithms in ``xgalois/utils/poly.hpp``:

.. code-block:: cpp

   #include <iostream>
   #include <memory>
   #include <vector>

   #include "xgalois/field/gf_element.hpp"
   #include "xgalois/field/gf_prime.hpp"
   #include "xgalois/poly/poly_dense.hpp"
   #include "xgalois/utils/poly.hpp"

   int main() {
     using Field = xg::GaloisFieldPrime<uint8_t>;
     using Elem = xg::GaloisFieldElement<Field>;
     using Poly = xg::PolynomialDense<Field>;

     auto gf7 = std::make_shared<Field>(7);

     // f(x) = 3x^2 + 2x + 1, g(x) = x + 1 over GF(7)
     Poly f({Elem(1, gf7), Elem(2, gf7), Elem(3, gf7)});
     Poly g({Elem(1, gf7), Elem(1, gf7)});

     auto product = f * g;
     auto [q, r] = f.DivRem(g);
     auto df = f.Derivative();

     std::cout << "f(x) * g(x) = " << product.ToString() << '\n';
     std::cout << "f / g = " << q.ToString() << ", rem = " << r.ToString() << '\n';
     std::cout << "f'(x) = " << df.ToString() << '\n';

     // GCD, irreducibility, complete factorization, and root finding
     auto gcd_fg = xg::utils::PolynomialGcd(product, g);
     auto [lc, factors] = xg::utils::Factor(product);
     auto roots = xg::utils::FindRoots(product);

     std::cout << "gcd(f*g, g) = " << gcd_fg.ToString()
               << ", f irreducible? " << std::boolalpha
               << xg::utils::IsIrreducible(f) << '\n';
     std::cout << "Leading coeff = " << lc << ", #factors = " << factors.size()
               << ", #roots = " << roots.size() << '\n';

     // Exact O(k^2) Newton interpolation from 3 evaluation points
     std::vector<Elem> xs = {Elem(0, gf7), Elem(1, gf7), Elem(2, gf7)};
     std::vector<Elem> ys = {f(xs[0]), f(xs[1]), f(xs[2])};
     auto recovered = xg::utils::Interpolate(xs, ys);
     std::cout << "Interpolated == f: " << (recovered == f) << '\n';
     return 0;
   }

Reed–Solomon Codes, Errors-and-Erasures & Simulation
----------------------------------------------------

Construct a primitive Reed–Solomon code over :math:`\mathrm{GF}(2^4)`, transmit
through a seeded ``ErrorErasureChannel``, decode with ``BerlekampMasseyDecoder``,
and benchmark code performance with ``SimulateCodePerformance``:

.. code-block:: cpp

   #include <iostream>
   #include <memory>
   #include <random>

   #include "xgalois/channel/channel.hpp"
   #include "xgalois/coding/coding.hpp"
   #include "xgalois/field/gf_array.hpp"
   #include "xgalois/field/gf_binary.hpp"

   int main() {
     using Field = xg::GF2XLOG<uint8_t>;
     auto gf16 = std::make_shared<Field>(4);  // GF(2^4)

     // Primitive RS(15, 11): d = 5, corrects any (e, s) with 2e + s <= 4
     auto rs = xg::coding::GeneralizedReedSolomonCode<Field>::ReedSolomon(
         gf16, /*length=*/15, /*dimension=*/11);

     std::mt19937_64 rng(42);
     auto message = xg::GFArray<Field>::Random(gf16, {11}, rng);
     auto codeword = rs.Encode(message);

     // Introduce 1 symbol error + 2 erasures (2*1 + 2 = 4 <= 4)
     xg::channel::ErrorErasureChannel<Field> ee_channel(
         gf16, /*input_size=*/15, /*num_errors=*/1, /*num_erasures=*/2,
         /*seed=*/12345);
     auto tx = ee_channel.TransmitWithErasures(codeword);

     // Decode errors and erasures with BerlekampMasseyDecoder
     xg::coding::BerlekampMasseyDecoder<Field> bm_decoder(rs);
     auto res = bm_decoder.DecodeWithErasures(tx.received, tx.erasure_positions);

     std::cout << "Decode OK: " << std::boolalpha << res.Ok()
               << " (errors=" << res.num_errors
               << ", erasures=" << res.num_erasures
               << ", recovered=" << (res.message == message) << ")\n";

     // Deterministic Monte Carlo FER/SER/BER simulation over a Q-ary Symmetric Channel
     xg::channel::QarySymmetricChannel<Field> qsc(gf16, /*input_size=*/15,
                                                  /*epsilon=*/0.05, /*seed=*/99);
     xg::channel::SimulationConfig cfg{
         .max_frames = 500, .target_frame_errors = 50, .seed = 2026};
     auto stats = xg::channel::SimulateCodePerformance(rs, qsc, cfg);

     std::cout << "Simulated " << stats.frames_sent << " frames: FER="
               << stats.Fer() << ", SER=" << stats.Ser()
               << ", BER=" << stats.Ber() << '\n';
     return 0;
   }

GF(2^128), Guruswami–Sudan List Decoding, LDPC & Polar Codes
------------------------------------------------------------

xgalois v0.2 also includes cryptographic :math:`\mathrm{GF}(2^{128})` arithmetic
with AES-GCM GHASH multiplication, Guruswami–Sudan list decoding beyond the
half-minimum-distance bound, iterative LDPC Belief Propagation, and Arikan Polar
SC/SCL/CA-SCL codes:

.. code-block:: cpp

   #include <algorithm>
   #include <iostream>
   #include <memory>
   #include <random>
   #include <vector>

   #include "xgalois/channel/channel.hpp"
   #include "xgalois/coding/coding.hpp"
   #include "xgalois/coding/decoder/guruswami_sudan_decoder.hpp"
   #include "xgalois/coding/ldpc.hpp"
   #include "xgalois/coding/polar.hpp"
   #include "xgalois/field/gf128.hpp"
   #include "xgalois/field/gf_array.hpp"
   #include "xgalois/field/gf_binary.hpp"

   int main() {
     // 1. GF(2^128) & NIST SP 800-38D AES-GCM GHASH multiplication
     xg::GF128 gf128;
     const auto h = xg::UInt128::FromHex("0x66e94bd4ef8a2c3b884cfa59ca342b2e");
     const auto c1 = xg::UInt128::FromHex("0x0388dace60b6a392f328c2b971b2fe78");
     const auto x1 = gf128.GhashMul(c1, h);
     std::cout << "GHASH(C1, H) = " << x1.ToHex(/*zero_pad=*/true) << '\n';

     // 2. Guruswami-Sudan list decoding on RS(15, 3): unique radius t = 6, Johnson = 9
     using Field16 = xg::GF2XLOG<uint8_t>;
     auto gf16 = std::make_shared<Field16>(4);
     auto rs = xg::coding::GeneralizedReedSolomonCode<Field16>::ReedSolomon(
         gf16, /*length=*/15, /*dimension=*/3);
     auto gs = xg::coding::GuruswamiSudanDecoder<Field16>::FromTargetRadius(rs, 8);
     std::mt19937_64 rng(42);
     auto rs_msg = xg::GFArray<Field16>::Random(gf16, {3}, rng);
     auto rs_cw = rs.Encode(rs_msg);
     auto rs_ch = xg::channel::StaticErrorRateChannel<Field16>::Fixed(
         gf16, 15, /*num_errors=*/8, /*seed=*/17);
     auto candidates = gs.DecodeToList(rs_ch.Transmit(rs_cw));
     const bool found =
         std::find(candidates.begin(), candidates.end(), rs_cw) != candidates.end();
     std::cout << "Guruswami-Sudan (e=8 > t=6, Johnson=" << gs.JohnsonRadius()
               << "): list_size=" << candidates.size()
               << ", contains_tx=" << std::boolalpha << found << '\n';

     // 3. Quasi-Cyclic LDPC with Normalized Min-Sum Belief Propagation over AWGN
     auto gf2 = std::make_shared<xg::GF2>();
     const std::vector<std::vector<int>> qc_base = {
         {0, 0, 0, 0, 0, 0},
         {0, 1, 2, 3, 4, 5},
         {0, 2, 4, 6, 8, 10},
     };
     auto ldpc = xg::coding::LDPCCode::FromQuasiCyclic(gf2, qc_base, /*z=*/17);
     xg::coding::BeliefPropagationDecoder bp(
         ldpc, xg::coding::BpAlgorithm::kNormalizedMinSum, /*max_iter=*/50, 0.75);
     auto ldpc_msg = xg::GFArray<xg::GF2>::Random(gf2, {ldpc.Dimension()}, rng);
     auto ldpc_awgn = xg::channel::AWGNChannel::FromEbN0(
         gf2, ldpc.Length(), /*eb_n0_db=*/5.0, ldpc.Rate(), /*seed=*/101);
     auto ldpc_res = bp.DecodeLlr(ldpc_awgn.TransmitLlr(ldpc.Encode(ldpc_msg)));
     std::cout << "LDPC QC(102, 53) Min-Sum OK=" << ldpc_res.Ok()
               << ", recovered=" << (ldpc_res.message == ldpc_msg) << '\n';

     // 4. Polar(64, 32) with CRC-8 and CA-SCL(L=8) over AWGN
     auto polar = xg::coding::PolarCode::FromDesignSnr(
         gf2, /*n=*/64, /*k=*/32, /*design_snr_db=*/2.5, /*systematic=*/false,
         /*crc_bits=*/8);
     xg::coding::SuccessiveCancellationListDecoder scl(polar, /*list_size=*/8);
     auto polar_msg = xg::GFArray<xg::GF2>::Random(gf2, {32}, rng);
     auto polar_awgn = xg::channel::AWGNChannel::FromEbN0(
         gf2, 64, /*eb_n0_db=*/4.5, 0.5, /*seed=*/202);
     auto polar_res = scl.DecodeLlr(polar_awgn.TransmitLlr(polar.Encode(polar_msg)));
     std::cout << "Polar(64, 32) CA-SCL(L=8) OK=" << polar_res.Ok()
               << ", recovered=" << (polar_res.message == polar_msg) << '\n';
     return 0;
   }

Next Steps
----------

- Run the end-to-end examples in ``example/`` via ``bazel run //example:<name>``
  or ``./build/example/<name>`` (CMake).
- Browse the :doc:`C++ API Reference <api/library_root>` for complete class documentation.
- Read the :doc:`Developer's Guide <dev_guide>` to learn how to contribute.
