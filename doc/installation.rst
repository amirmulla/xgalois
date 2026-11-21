.. _installation:

Installation Guide
==================

This guide covers building, testing, benchmarking, and generating documentation
for **xgalois** from source on Linux and macOS.

Requirements
------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Dependency
     - Version / Notes
   * - C++ compiler
     - C++23-capable compiler (``g++ 13+`` or ``clang++ 17+``)
   * - Bazel
     - Bazel 8 or 9 (pure Bzlmod via ``MODULE.bazel``)
   * - Doxygen (docs only)
     - 1.9+ (required by Breathe/Exhale to generate XML API metadata)
   * - Python (docs only)
     - 3.10+ (for Sphinx, Furo, Breathe, and Exhale)

Building with Bazel
-------------------

Clone the repository and build the library and all examples:

.. code-block:: bash

   git clone https://github.com/amirmulla/xgalois.git
   cd xgalois
   bazel build //xgalois/... //example/...

Running Tests
-------------

Run the full unit test suite (all targets under ``//test/...``):

.. code-block:: bash

   bazel test //test/...

Run the aggregated suite with verbose error output:

.. code-block:: bash

   bazel test //test:galois_field_tests --test_output=errors

Run an individual runnable example:

.. code-block:: bash

   bazel run //example:coding_example

Build Configurations (``.bazelrc``)
-----------------------------------

xgalois ships pre-configured Bazel profiles in ``.bazelrc`` (C++23 and strict
``-Wall -Wextra -Wpedantic`` warnings on ``xgalois/`` are enabled by default):

.. code-block:: bash

   # Optimized release build (-c opt -O3)
   bazel build --config=opt //xgalois/...

   # Native ISA tuning (-c opt -O3 -march=native, enables PCLMULQDQ / SSSE3)
   bazel run --config=native //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark
   bazel run --config=native //benchmark/gf_prime/src:gf_prime_benchmark

   # Debug build (-c dbg -O0 -g)
   bazel build --config=dbg //xgalois/...

   # Sanitizers (AddressSanitizer, UndefinedBehaviorSanitizer, ThreadSanitizer)
   bazel test --config=asan //test:galois_field_tests
   bazel test --config=ubsan //test:galois_field_tests
   bazel test --config=tsan //test:galois_field_tests

   # Code coverage (LCOV combined report for //xgalois/...)
   bazel coverage --config=coverage //test:galois_field_tests

Building the Documentation
--------------------------

Install ``doxygen`` via your system package manager (e.g.,
``sudo apt-get install -y doxygen`` on Ubuntu/Debian or ``brew install doxygen``
on macOS), then build the Sphinx site:

.. code-block:: bash

   cd doc
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ./build.sh

Open ``doc/_build/html/index.html`` in your browser to view the generated site.
``./build.sh`` runs ``sphinx-build -W -b html . _build/html`` so documentation
warnings are caught as errors.

Dependencies
------------

All C++ build and test dependencies are declared in ``MODULE.bazel`` and
fetched automatically by Bazel (pure Bzlmod):

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Dependency
     - Version
     - Purpose
   * - `rules_cc <https://github.com/bazelbuild/rules_cc>`_
     - 0.2.17
     - C++ Bazel build rules and runfiles library
   * - `platforms <https://github.com/bazelbuild/platforms>`_
     - 1.0.0
     - Bazel target platform constraints
   * - `googletest <https://github.com/google/googletest>`_
     - 1.17.0.bcr.2
     - Unit testing framework
   * - `google_benchmark <https://github.com/google/benchmark>`_
     - 1.9.5
     - Microbenchmarking framework
   * - `xtl <https://github.com/xtensor-stack/xtl>`_
     - 0.8.0
     - Foundational C++ template library for xtensor
   * - `xtensor <https://github.com/xtensor-stack/xtensor>`_
     - 0.26.0
     - Multi-dimensional array containers backing ``GFArray``
