.. _xgalois-docs:

xgalois Documentation
=====================

**xgalois** is a high-performance C++23 library for computations over
`Galois Fields <https://en.wikipedia.org/wiki/Finite_field>`_ (finite fields),
providing a comprehensive toolkit for polynomial algebra, finite-field linear
algebra, communication channel modeling, and algebraic coding theory — all
natively over :math:`\mathrm{GF}(p^m)`.

.. note::

   This is the documentation for **xgalois v0.2.0**.
   The library targets C++23 and pure Bzlmod (Bazel 8/9).

----

Tutorials & Getting Started
----------------------------

.. topic:: Tutorial

   A hands-on introduction to the xgalois library. Learn how to construct
   finite fields, work with contiguous ``GFArray`` buffers, factor polynomials,
   solve linear systems, decode Reed–Solomon codes with errors and erasures,
   and run Monte Carlo channel simulations.

   :doc:`Go to Tutorial → <tutorial>`

.. topic:: Installation

   Build, test, sanitizer, benchmark, and documentation instructions using
   Bazel 8/9 (pure Bzlmod) on Linux and macOS.

   :doc:`Go to Installation Guide → <installation>`

----

Comprehensive Reference Manual
-------------------------------

The reference manual contains automatically generated documentation
for every class, function, and namespace in the xgalois C++23 library.

.. topic:: Finite Fields & Contiguous Arrays — ``xg::GaloisFieldBase``, ``xg::GFArray``

   Prime fields :math:`\mathrm{GF}(p)` (division-free Barrett reduction,
   Tonelli–Shanks, BSGS/Pohlig–Hellman, log tables), binary and binary
   extension fields :math:`\mathrm{GF}(2^m)` up to :math:`\mathrm{GF}(2^{64})`
   (PCLMULQDQ/SSSE3, optimized log/antilog tables, Zech logarithm tables),
   general extension fields :math:`\mathrm{GF}(p^m)`, runtime factory
   ``xg::GF(...)``, and packed ``xg::GFArray`` / ``xg::GFMatrix`` containers.

.. topic:: Polynomials — ``xg::PolynomialDense``, ``xg::utils``

   Dense polynomials over arbitrary Galois fields with Karatsuba
   multiplication, Euclidean ``DivRem``, ``PowMod``, formal derivatives,
   batch evaluation, GCD/Extended GCD, Rabin irreducibility/primitivity tests,
   Cantor–Zassenhaus factorization, ``FindRoots``, ``ChienSearch``, and
   :math:`O(k^2)` Newton ``Interpolate``.

.. topic:: Linear Algebra — ``xg::linalg``

   Contiguous ``GFArray`` matrix operations and a bit-packed ``uint64_t``
   M4RI engine for :math:`\mathrm{GF}(2)`: ``dot``, ``vdot``, ``outer``,
   ``cross``, ``kron``, ``matrix_power``, ``row_echelon``, ``rref``,
   ``matrix_rank``, ``trace``, ``det``, ``inv``, ``solve``, ``nullspace``,
   ``left_nullspace``, ``row_space``, ``column_space``, ``systematic_form``,
   ``plu``, and division-free Berkowitz ``characteristic_poly``.

.. topic:: Coding Theory — ``xg::coding``

   ``LinearCode`` (with Hamming, Extended Hamming, Binary/Extended Binary/Ternary
   Golay, and Reed–Muller factories), ``CyclicCode``, ``BCHCode`` (binary and
   non-binary with algebraic BM/Chien/Forney decoding), and
   ``GeneralizedReedSolomonCode`` (GRS/RS/Extended RS) with Berlekamp–Massey,
   Sugiyama, and PGZ decoders supporting full errors-and-erasures decoding.

.. topic:: Communication Channels & Simulation — ``xg::channel``

   Seeded channel models (``BinarySymmetricChannel``, ``QarySymmetricChannel``,
   ``QaryErasureChannel``, ``BinaryErasureChannel``, ``ErrorErasureChannel``,
   ``StaticErrorRateChannel``, BPSK ``AWGNChannel`` with exact LLR) and
   deterministic Monte Carlo BER/SER/FER simulation via
   ``SimulateCodePerformance``.

.. topic:: Databases — ``xg::databases``

   Embedded and file-backed lookup tables for Conway and irreducible
   polynomials used for fast canonical field construction.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference Manual

   api/library_root

----

Developer's Guide
------------------

.. topic:: Contributing to xgalois

   Guidelines for repository layout, C++23 coding conventions, granular Bazel
   targets, formatting, sanitizers, testing, and submitting contributions.

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
