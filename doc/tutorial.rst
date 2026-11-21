.. _tutorial:

Tutorial
========

This tutorial walks through the core workflows in **xgalois v0.2**. Every code
listing below is a complete, self-contained C++23 program targeting the actual
v0.2 API.

By the end of this guide, you will know how to:

- Construct prime, binary extension, and general extension fields (directly or via ``xg::GF``)
- Perform contiguous array and linear algebra operations with ``xg::GFArray`` and ``xg::linalg``
- Manipulate, factor, find roots of, and interpolate polynomials with ``xg::PolynomialDense`` and ``xg::utils``
- Encode and decode Reed–Solomon codes with errors and erasures using ``BerlekampMasseyDecoder``
- Run deterministic Monte Carlo BER/SER/FER simulations with ``SimulateCodePerformance``

Prerequisites
-------------

- A C++23-capable compiler (``g++ 13+`` or ``clang++ 17+``)
- `Bazel <https://bazel.build/>`_ 8 or 9 (pure Bzlmod)
- A local clone of the xgalois repository

Creating Finite Fields
----------------------

xgalois supports four field families:

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
array kernels and linear algebra in ``xg::linalg``:

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

     auto x = xg::linalg::solve(A, b);
     auto check = xg::linalg::dot(A, x);  // A * x == b

     std::cout << "det(A) = " << xg::linalg::det(A)
               << ", A * x == b: " << std::boolalpha << (check == b) << '\n';

     // Nullspace of a 2x4 matrix M over GF(7)
     xg::GFArray<Field> M(gf7, xt::xarray<uint8_t>{{1, 2, 3, 4},
                                                   {2, 4, 1, 3}});
     auto N = xg::linalg::nullspace(M);  // shape: (2, 4), satisfies M * N^T == 0
     auto zero_check = xg::linalg::dot(M, xg::linalg::transpose(N));
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
     auto tx = ee_channel.transmit_with_erasures(codeword);

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
               << stats.FER() << ", SER=" << stats.SER()
               << ", BER=" << stats.BER() << '\n';
     return 0;
   }

Next Steps
----------

- Run the end-to-end examples in ``example/`` via ``bazel run //example:<name>``.
- Browse the :doc:`C++ API Reference <api/library_root>` for complete class documentation.
- Read the :doc:`Developer's Guide <dev_guide>` to learn how to contribute.
