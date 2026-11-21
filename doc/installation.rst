.. _installation:

Installation Guide
==================

This guide covers building, testing, installing, benchmarking, and generating
documentation for **xgalois** from source on Linux and macOS using either
**Bazel** or **CMake**.

Requirements
------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Dependency
     - Version / Notes
   * - C++ compiler
     - C++20-capable compiler (``GCC 11+`` or ``Clang 14+``)
   * - Bazel (option 1)
     - Bazel 8 or 9 (pure Bzlmod via ``MODULE.bazel``)
   * - CMake (option 2)
     - CMake ``>= 3.20`` (auto-detects system ``xtl``/``xtensor`` or fetches via ``FetchContent``)
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

Running Tests (Bazel)
^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

xgalois ships pre-configured Bazel profiles in ``.bazelrc`` (**C++20** and
strict ``-Wall -Wextra -Wpedantic`` warnings on ``xgalois/`` are enabled by
default):

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

Building & Installing with CMake
--------------------------------

xgalois provides a first-class CMake ``>= 3.20`` build exporting the
header-only interface target ``xgalois::xgalois`` with ``cxx_std_20``.
Available CMake options:

- ``XGALOIS_BUILD_TESTS`` (default ``ON`` when top-level, ``OFF`` when included as a subproject)
- ``XGALOIS_BUILD_EXAMPLES`` (default ``ON`` when top-level, ``OFF`` when included as a subproject)
- ``XGALOIS_ENABLE_NATIVE`` (default ``OFF``; enables ``-march=native`` for SIMD/``PCLMULQDQ`` acceleration)

Configure, Build, and Test (CMake)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/amirmulla/xgalois.git
   cd xgalois

   # Configure (add -DXGALOIS_ENABLE_NATIVE=ON for -march=native)
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

   # Build tests and examples
   cmake --build build -j

   # Run unit tests via CTest
   ctest --test-dir build --output-on-failure

   # Run an example binary
   ./build/example/coding_example

Installing and Using ``find_package``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the headers, polynomial databases, and relocatable CMake package config:

.. code-block:: bash

   cmake --install build --prefix /usr/local

Then link against ``xgalois::xgalois`` in any downstream CMake project:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.20)
   project(my_app LANGUAGES CXX)

   find_package(xgalois CONFIG REQUIRED)

   add_executable(my_app main.cpp)
   target_link_libraries(my_app PRIVATE xgalois::xgalois)

If you installed to a custom prefix (for example ``$HOME/.local``), pass
``-DCMAKE_PREFIX_PATH=$HOME/.local`` when configuring your downstream project.

Using ``FetchContent``
^^^^^^^^^^^^^^^^^^^^^^

To embed xgalois directly in a CMake project without a separate install step:

.. code-block:: cmake

   include(FetchContent)
   FetchContent_Declare(
     xgalois
     GIT_REPOSITORY https://github.com/amirmulla/xgalois.git
     GIT_TAG        main
   )
   FetchContent_MakeAvailable(xgalois)

   add_executable(my_app main.cpp)
   target_link_libraries(my_app PRIVATE xgalois::xgalois)

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

All C++ build and test dependencies are declared in ``MODULE.bazel`` (for
Bazel) and ``CMakeLists.txt`` (for CMake, using ``find_package`` with
``FetchContent`` fallback):

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Dependency
     - Version
     - Purpose
   * - `rules_cc <https://github.com/bazelbuild/rules_cc>`_
     - 0.2.17
     - C++ Bazel build rules and runfiles library (Bazel only)
   * - `platforms <https://github.com/bazelbuild/platforms>`_
     - 1.0.0
     - Bazel target platform constraints (Bazel only)
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
