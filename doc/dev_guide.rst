.. _dev_guide:

Developer's Guide
=================

Thank you for your interest in contributing to **xgalois**! This guide
will help you set up your development environment and understand the
project's conventions.


Project Structure
-----------------

.. code-block:: text

   xgalois/
   ├── xgalois/                 # C++ library source code
   │   ├── field/               # Galois field implementations
   │   │   ├── gf_base.hpp      # Abstract base template
   │   │   ├── gf_prime.hpp     # GF(p) — prime fields
   │   │   ├── gf_binary.hpp    # GF(2^n) — binary extension fields
   │   │   ├── gf_extension.hpp # GF(p^n) — general extensions
   │   │   ├── gf_element.hpp   # Wrapper for field elements
   │   │   └── gf_factory.hpp   # Factory for creating fields
   │   ├── coding/              # Error-correcting codes
   │   │   ├── abstract_code.hpp
   │   │   ├── abstract_linear_code.hpp
   │   │   ├── grs.hpp          # Generalized Reed–Solomon
   │   │   ├── cyclic_code.hpp
   │   │   ├── encoder/         # Encoding strategies
   │   │   └── decoder/         # Decoding strategies
   │   ├── channel/             # Communication channel models
   │   ├── linalg/              # Linear algebra over finite fields
   │   ├── poly/                # Polynomial representations
   │   ├── databases/           # Lookup tables (Conway, irreducible polys)
   │   └── utils/               # Helper utilities
   ├── test/                    # Unit tests (Google Test)
   ├── benchmark/               # Performance benchmarks
   ├── example/                 # Usage examples
   ├── doc/                     # Sphinx documentation (you are here)
   └── WORKSPACE.bazel          # Bazel workspace configuration


Coding Conventions
-------------------

- **Language**: C++17.
- **Naming**:

  - Classes: ``PascalCase`` (e.g., ``GaloisFieldPrime``)
  - Functions / methods: ``PascalCase`` (e.g., ``Add()``, ``Encode()``)
  - Variables: ``snake_case`` (e.g., ``evaluation_points``)
  - Constants: ``kPascalCase`` (e.g., ``kGfZero``)
  - Namespaces: ``lowercase`` (``xg``, ``xg::coding``, ``xg::linalg``)

- **Headers**: Use ``.hpp`` extension for C++ headers.
- **Include guards**: ``#ifndef XGALOIS_MODULE_FILE_HPP``


Running Tests
--------------

.. code-block:: bash

   # Run all tests
   bazel test //test/...

   # Run a specific test
   bazel test //test:field_test

   # Run with verbose output
   bazel test //test:field_test --test_output=all


Adding a New Module
---------------------

1. Create your header(s) under ``xgalois/<module_name>/``.
2. Add the ``cc_library`` target in the corresponding ``BUILD.bazel``.
3. Write unit tests under ``test/``.
4. Run ``bazel test //test/...`` to verify.
5. Rebuild the documentation with ``cd doc && ./build.sh``.
   Your new classes will automatically appear in the API Reference.
