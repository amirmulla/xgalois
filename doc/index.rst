.. _xgalois-docs:

xgalois Documentation
=====================

**xgalois** is a high-performance C++20 library for computations over
`Galois Fields <https://en.wikipedia.org/wiki/Finite_field>`_ (finite fields),
providing a comprehensive toolkit for polynomial algebra, finite-field linear
algebra, communication channel modeling, and algebraic & modern coding theory —
all natively over :math:`\mathrm{GF}(p^m)`.

.. note::

   This is the documentation for **xgalois v0.2.0**.
   The library targets **C++20** and supports both **Bazel 8/9 (pure Bzlmod)**
   and **CMake >= 3.20** (``FetchContent`` and ``find_package(xgalois CONFIG REQUIRED)``).

----

Tutorials & Getting Started
----------------------------

.. topic:: Tutorial

   A hands-on introduction to the xgalois library. Learn how to construct
   finite fields (including :math:`\mathrm{GF}(2^{128})` with AES-GCM GHASH),
   work with contiguous ``GFArray`` buffers, factor polynomials, solve linear
   systems, decode Reed–Solomon codes with errors, erasures, and Guruswami–Sudan
   list decoding, decode LDPC and Polar codes over AWGN, and run Monte Carlo
   channel simulations.

   :doc:`Go to Tutorial → <tutorial>`

.. topic:: Installation

   Build, test, install, sanitizer, benchmark, and documentation instructions
   using Bazel 8/9 (pure Bzlmod) or CMake >= 3.20 on Linux and macOS.

   :doc:`Go to Installation Guide → <installation>`

----

Comprehensive Reference Manual
-------------------------------

The reference manual contains automatically generated documentation
for every class, function, and namespace in the xgalois C++20 library.

.. topic:: Finite Fields & Contiguous Arrays — ``xg::GaloisFieldBase``, ``xg::GFArray``

   Prime fields :math:`\mathrm{GF}(p)` (division-free Barrett reduction,
   Tonelli–Shanks, BSGS/Pohlig–Hellman, log tables), binary and binary
   extension fields :math:`\mathrm{GF}(2^m)` up to :math:`\mathrm{GF}(2^{64})`
   (PCLMULQDQ/SSSE3, optimized log/antilog tables, Zech logarithm tables),
   dedicated :math:`\mathrm{GF}(2^{128})` (``GF128``, ``UInt128``, Karatsuba
   PCLMULQDQ, NIST SP 800-38D ``GhashMul``, general carry-less Barrett
   reduction), general extension fields :math:`\mathrm{GF}(p^m)`, runtime
   factory ``xg::GF(...)``, and packed ``xg::GFArray`` / ``xg::GFMatrix``
   containers.

.. topic:: Polynomials — ``xg::PolynomialDense``, ``xg::utils``

   Dense polynomials over arbitrary Galois fields with Karatsuba
   multiplication, Euclidean ``DivRem``, ``PowMod``, formal derivatives,
   batch evaluation, GCD/Extended GCD, Rabin irreducibility/primitivity tests,
   Cantor–Zassenhaus factorization, ``FindRoots``, ``ChienSearch``, and
   :math:`O(k^2)` Newton ``Interpolate``.

.. topic:: Linear Algebra — ``xg::linalg``

   Contiguous ``GFArray`` matrix operations and a bit-packed ``uint64_t``
   M4RI engine for :math:`\mathrm{GF}(2)` using strict Google-style
   ``PascalCase`` naming: ``Dot``, ``Vdot``, ``Outer``, ``Cross``, ``Kron``,
   ``MatrixPower``, ``RowEchelon``, ``Rref``, ``MatrixRank``, ``Trace``,
   ``Det``, ``Inv``, ``Solve``, ``Nullspace``, ``LeftNullspace``, ``RowSpace``,
   ``ColumnSpace``, ``SystematicForm``, ``Plu``, and division-free Berkowitz
   ``CharacteristicPoly``.

.. topic:: Coding Theory — ``xg::coding``

   ``LinearCode`` (with Hamming, Extended Hamming, Binary/Extended Binary/Ternary
   Golay, and Reed–Muller factories), ``CyclicCode``, ``BCHCode`` (binary and
   non-binary with algebraic BM/Chien/Forney decoding),
   ``GeneralizedReedSolomonCode`` (GRS/RS/Extended RS) with Berlekamp–Massey,
   Sugiyama, PGZ (errors-and-erasures), and **Guruswami–Sudan list decoding**
   (``GuruswamiSudanDecoder`` up to the Johnson bound), plus **LDPC codes**
   (``LDPCCode``, ``BitFlippingDecoder``, ``BeliefPropagationDecoder`` with
   Sum-Product and Normalized/Offset Min-Sum) and **Polar codes**
   (``PolarCode``, ``SuccessiveCancellationDecoder``,
   ``SuccessiveCancellationListDecoder`` for SC/SCL/CA-SCL).

.. topic:: Communication Channels & Simulation — ``xg::channel``

   Seeded channel models (``BinarySymmetricChannel``, ``QarySymmetricChannel``,
   ``QaryErasureChannel``, ``BinaryErasureChannel``, ``ErrorErasureChannel``,
   ``StaticErrorRateChannel``, BPSK ``AWGNChannel`` with exact LLR) and
   deterministic Monte Carlo BER/SER/FER simulation via
   ``SimulateCodePerformance``.

.. topic:: Databases — ``xg::databases``

   Embedded and file-backed lookup tables for Conway and irreducible
   polynomials used for fast canonical field construction (``Fetch``,
   ``FetchAll``, ``FileExists``, ``GetDbPath``).

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference Manual

   api/library_root

----

Developer's Guide
------------------

.. topic:: Contributing to xgalois

   Guidelines for repository layout, C++20 and strict Google C++ Style
   conventions, granular Bazel and CMake targets, formatting, sanitizers,
   testing, and submitting contributions.

   :doc:`Go to Developer's Guide → <dev_guide>`

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Guides

   tutorial
   installation
   dev_guide

----

Indices and Tables
==================

* :ref:`genindex`
* :ref:`search`
