.. _tutorial:

Tutorial
========

This tutorial walks you through the fundamental capabilities of the
**xgalois** library. By the end, you will know how to:

- Create and manipulate finite fields
- Perform polynomial arithmetic over those fields
- Construct and use error-correcting codes
- Simulate noisy communication channels

Prerequisites
-------------

- A C++17 compatible compiler (``clang++`` or ``g++``)
- `Bazel <https://bazel.build/>`_ build system installed
- The xgalois source code (``git clone``)


Creating Finite Fields
-----------------------

xgalois supports three families of finite fields:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Type
     - Class
     - Example
   * - Prime field :math:`\mathrm{GF}(p)`
     - ``GaloisFieldPrime``
     - :math:`\mathrm{GF}(7)`, :math:`\mathrm{GF}(31)`
   * - Binary extension :math:`\mathrm{GF}(2^n)`
     - ``GaloisFieldBinaryExtension``
     - :math:`\mathrm{GF}(2^8)`, :math:`\mathrm{GF}(2^{16})`
   * - General extension :math:`\mathrm{GF}(p^n)`
     - ``GaloisFieldExtension``
     - :math:`\mathrm{GF}(3^4)`, :math:`\mathrm{GF}(5^3)`

**Example — Creating GF(7):**

.. code-block:: cpp

   #include "xgalois/field/gf_factory.hpp"

   auto field = xg::CreateGaloisField(7);       // GF(7)
   auto a = field->SetElementValue(3);
   auto b = field->SetElementValue(5);
   auto c = field->Add(a, b);                   // (3 + 5) mod 7 = 1
   std::cout << field->ToString(c) << std::endl; // prints: 1


**Example — Creating GF(2^8):**

.. code-block:: cpp

   auto gf256 = xg::CreateGaloisField(2, 8);  // GF(2^8) = GF(256)


Polynomial Arithmetic
----------------------

Polynomials over Galois fields are first-class citizens in xgalois.

.. code-block:: cpp

   #include "xgalois/poly/poly_dense.hpp"

   auto field = xg::CreateGaloisField(7);

   // Create f(x) = 3x^2 + 2x + 1  over GF(7)
   auto f = xg::PolynomialDense<...>({one, two, three});
   auto g = xg::PolynomialDense<...>({one, one});           // g(x) = x + 1

   auto product = f * g;         // polynomial multiplication
   auto [q, r] = f.DivRem(g);    // quotient and remainder
   auto df = f.Derivative();     // formal derivative


Reed–Solomon Codes
-------------------

Construct a standard Reed–Solomon code over a prime field:

.. code-block:: cpp

   #include "xgalois/coding/grs.hpp"

   auto field = xg::CreateGaloisField(7);
   auto rs = xg::coding::GeneralizedReedSolomonCode<uint32_t>::ReedSolomon(
       field, /*length=*/6, /*dimension=*/3, /*primitive_element=*/3);

   // Encode
   auto msg = ...;
   auto codeword = rs->Encode(msg);

   // Decode
   auto decoded = rs->DecodeToMessage(received_word);


Simulating a Noisy Channel
----------------------------

.. code-block:: cpp

   #include "xgalois/channel/channel.hpp"

   auto channel = xg::channels::StaticErrorRateChannel(field, /*n=*/6, /*errors=*/2);
   auto received = channel.transmit(codeword);

   // The received word now contains up to 2 symbol errors.


Next Steps
----------

- Browse the :doc:`C++ API Reference <api/library_root>` for complete class documentation.
- Read the :doc:`Developer's Guide <dev_guide>` to learn how to contribute.
