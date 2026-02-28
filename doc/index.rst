.. _xgalois-docs:

xgalois Documentation
=====================

**xgalois** is a high-performance C++ library for computations over
`Galois Fields <https://en.wikipedia.org/wiki/Finite_field>`_ (finite fields),
providing a comprehensive toolkit for algebraic coding theory, linear algebra,
and polynomial arithmetic — all natively over :math:`\mathrm{GF}(p^n)`.

.. note::

   This is the documentation for **xgalois v0.1.0**.
   The library is under active development; the API may evolve.

----

Tutorials & Getting Started
----------------------------

.. topic:: Tutorial

   A hands-on introduction to the xgalois library. Learn how to create
   finite fields, perform arithmetic, construct codes, and simulate
   channels — all in a few lines of C++.

   :doc:`Go to Tutorial → <tutorial>`

.. topic:: Installation

   Build instructions using Bazel for macOS and Linux.

   :doc:`Go to Installation Guide → <installation>`

----

Comprehensive Reference Manual
-------------------------------

The reference manual contains automatically generated documentation
for every class, function, and namespace in the xgalois C++ library.

.. topic:: Finite Fields — ``xg::GaloisFieldBase``

   Representations and arithmetic for prime fields :math:`\mathrm{GF}(p)`,
   binary fields :math:`\mathrm{GF}(2^n)`, and general extension fields
   :math:`\mathrm{GF}(p^n)`. Includes element representation modes
   (integer, polynomial, logarithmic, power).

.. topic:: Coding Theory — ``xg::coding``

   Abstract linear codes, Generalized Reed–Solomon (GRS) codes,
   cyclic codes, encoders and decoders. Mirrors the structure found in
   `SageMath's Coding Theory module <https://doc.sagemath.org/html/en/reference/coding/index.html>`_.

.. topic:: Communication Channels — ``xg::channels``

   Simulated noisy channels for error-correcting code evaluation.
   Includes static error-rate channels, erasure channels, and
   q-ary symmetric channels.

.. topic:: Linear Algebra — ``xg::linalg``

   Matrix operations over finite fields: dot products, Kronecker products,
   row echelon form, RREF, determinants, matrix rank, inversion,
   and system solving.

.. topic:: Polynomials — ``xg::PolynomialDense``

   Dense polynomial representation over arbitrary Galois fields with
   full arithmetic (``+``, ``-``, ``*``, ``/``, ``%``, ``^``), Horner
   evaluation, formal derivatives, and division with remainder.

.. topic:: Databases — ``xg::databases``

   Built-in lookup tables for prime factorizations, irreducible
   polynomials, and Conway polynomials — used internally for
   fast field construction.

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference Manual

   api/library_root

----

Developer's Guide
------------------

.. topic:: Contributing to xgalois

   Guidelines for development setup, coding conventions, testing,
   and submitting contributions.

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
