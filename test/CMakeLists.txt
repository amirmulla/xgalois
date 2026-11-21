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
