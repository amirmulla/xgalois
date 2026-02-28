.. _installation:

Installation Guide
==================

This guide covers building xgalois from source on macOS and Linux.

Requirements
------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Dependency
     - Version
   * - C++ compiler
     - C++17 or later (``clang++ 14+`` / ``g++ 11+``)
   * - Bazel
     - 7.0+
   * - Python (docs only)
     - 3.10+

Building with Bazel
--------------------

Clone the repository and build:

.. code-block:: bash

   git clone https://github.com/amirmulla/xgalois.git
   cd xgalois
   bazel build //...

Running all tests:

.. code-block:: bash

   bazel test //test/...

Building the Documentation
---------------------------

.. code-block:: bash

   cd doc
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ./build.sh

Open ``doc/_build/html/index.html`` in your browser to view the generated site.

Dependencies
------------

xgalois depends on the following libraries, all fetched automatically by Bazel:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Library
     - Version
     - Purpose
   * - `Google Test <https://github.com/google/googletest>`_
     - 1.17.0
     - Unit testing framework
   * - `xtensor <https://github.com/xtensor-stack/xtensor>`_
     - 0.26.0
     - Multi-dimensional array containers
   * - `xtl <https://github.com/xtensor-stack/xtl>`_
     - 0.8.0
     - Foundational library for xtensor
   * - `rules_foreign_cc <https://github.com/bazel-contrib/rules_foreign_cc>`_
     - 0.15.0
     - Build CMake projects with Bazel
