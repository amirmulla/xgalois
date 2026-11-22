find_package(Threads REQUIRED)
find_package(GTest CONFIG QUIET)

if(NOT GTest_FOUND)
  include(FetchContent)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.tar.gz
  )
  FetchContent_MakeAvailable(googletest)
endif()

function(xgalois_add_test target_name)
  add_executable(${target_name} ${ARGN})
  target_link_libraries(${target_name}
    PRIVATE
      xgalois::xgalois
      GTest::gtest
      GTest::gtest_main
      Threads::Threads
  )
  target_compile_options(${target_name}
    PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
      $<$<CXX_COMPILER_ID:MSVC>:/W4>
  )
  add_test(
    NAME ${target_name}
    COMMAND ${target_name}
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )
endfunction()

# Field tests
xgalois_add_test(gf_prime_test field/gf_prime_test.cpp)
xgalois_add_test(gf_binary_test field/gf_binary_test.cpp)
xgalois_add_test(gf_tables_test field/gf_tables_test.cpp)
xgalois_add_test(gf_extension_test field/gf_extension_test.cpp)
xgalois_add_test(gf_factory_test field/gf_factory_test.cpp)
xgalois_add_test(gf_element_test field/gf_element_test.cpp)
xgalois_add_test(gf128_test field/gf128_test.cpp)

# Polynomial tests
xgalois_add_test(poly_dense_test poly/poly_dense_test.cpp)
xgalois_add_test(poly_algo_test poly/poly_algo_test.cpp)

# Database tests
xgalois_add_test(interface_test databases/interface_test.cpp)
xgalois_add_test(db_integrity_test databases/db_integrity_test.cpp)

# Linear algebra tests
xgalois_add_test(linalg_prime_test linalg/linalg_prime_test.cpp)
xgalois_add_test(linalg_binary_test linalg/linalg_binary_test.cpp)
xgalois_add_test(linalg_gf2_test linalg/linalg_gf2_test.cpp)
xgalois_add_test(linalg_extension_test linalg/linalg_extension_test.cpp)

# Utility tests
xgalois_add_test(math_test utils/math_test.cpp)
xgalois_add_test(field_test utils/field_test.cpp)
xgalois_add_test(poly_test utils/poly_test.cpp)

# Channel tests
xgalois_add_test(channel_test channel/channel_test.cpp)

# Coding theory tests
xgalois_add_test(coding_test coding/coding_test.cpp)
xgalois_add_test(grs_test coding/grs_test.cpp)
xgalois_add_test(cyclic_bch_test coding/cyclic_bch_test.cpp)
xgalois_add_test(syndrome_linear_test coding/syndrome_linear_test.cpp)
xgalois_add_test(ldpc_polar_test coding/ldpc_polar_test.cpp)
xgalois_add_test(guruswami_sudan_test coding/guruswami_sudan_test.cpp)

# Property, fuzz, and multi-TU ODR tests
xgalois_add_test(algebraic_property_test property/algebraic_property_test.cpp)
xgalois_add_test(fuzz_smoke_test fuzz/fuzz_smoke_test.cpp)
xgalois_add_test(odr_multi_tu_test build/odr_tu1.cpp build/odr_tu2.cpp)

# v0.3 Wave A tests (BigInt, Algebra Concepts/SNF, SIMD Dispatch/ThreadPool)
xgalois_add_test(bigint_test utils/bigint_test.cpp)
xgalois_add_test(concepts_snf_test algebra/concepts_snf_test.cpp)
xgalois_add_test(simd_dispatch_test simd/simd_dispatch_test.cpp)

# v0.3 Wave B tests (Groups/BSGS, Rings/Skew, Field Towers/ONB/BigInt Fields)
xgalois_add_test(groups_test algebra/groups_test.cpp)
xgalois_add_test(rings_test ring/rings_test.cpp)
xgalois_add_test(advanced_fields_test field/advanced_fields_test.cpp)

# v0.3 Wave C tests (Fast Polynomials / NTT / Additive FFT / Half-GCD & Tier-1 Golden Oracle)
xgalois_add_test(fast_poly_test poly/fast_poly_test.cpp)
xgalois_add_test(sage_golden_test oracle/sage_golden_test.cpp)
target_compile_definitions(sage_golden_test
  PRIVATE
    XGALOIS_GOLDEN_DIR="${CMAKE_CURRENT_SOURCE_DIR}/oracle/golden"
)

# v0.4 Wave A tests (Tamo-Barg LRC, Interleaved RS/Power Decoding, Advanced Channels/Sim)
xgalois_add_test(lrc_test coding/lrc_test.cpp)
xgalois_add_test(interleaved_rs_test coding/interleaved_rs_test.cpp)
xgalois_add_test(advanced_channels_test channel/advanced_channels_test.cpp)

# v0.4 Wave B tests (Gabidulin Rank-Metric & Kotter-Kschischang Subspace Codes)
xgalois_add_test(gabidulin_subspace_test coding/gabidulin_subspace_test.cpp)

if(XGALOIS_BUILD_ORACLE)
  add_executable(xgalois_oracle_cli oracle/xgalois_oracle_cli.cpp)
  target_link_libraries(xgalois_oracle_cli PRIVATE xgalois::xgalois Threads::Threads)
  target_compile_options(xgalois_oracle_cli
    PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
      $<$<CXX_COMPILER_ID:MSVC>:/W4>
  )
endif()
