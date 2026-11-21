cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project(xgalois
  VERSION 0.2.0
  DESCRIPTION "High-performance C++20 Galois field, polynomial, linear algebra, and coding theory library"
  LANGUAGES CXX
)

# Enforce C++20 standard project-wide
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FetchContent)

string(COMPARE EQUAL "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}" _xgalois_is_top_level)

option(XGALOIS_BUILD_TESTS "Build XGalois unit tests" ${_xgalois_is_top_level})
option(XGALOIS_BUILD_EXAMPLES "Build XGalois example binaries" ${_xgalois_is_top_level})
option(XGALOIS_ENABLE_NATIVE "Enable -march=native SIMD/PCLMULQDQ optimizations" OFF)

# ------------------------------------------------------------------------------
# Dependencies: xtl (0.8.0) and xtensor (0.26.0)
# ------------------------------------------------------------------------------
find_package(xtl CONFIG QUIET)
if(xtl_FOUND)
  set(XGALOIS_USE_EXTERNAL_XTL ON)
else()
  set(XGALOIS_USE_EXTERNAL_XTL OFF)
  FetchContent_Declare(
    xtl
    URL https://github.com/xtensor-stack/xtl/archive/refs/tags/0.8.0.tar.gz
    URL_HASH SHA256=ee38153b7dd0ec84cee3361f5488a4e7e6ddd26392612ac8821cbc76e740273a
    SOURCE_SUBDIR include
  )
  FetchContent_MakeAvailable(xtl)
endif()

find_package(xtensor CONFIG QUIET)
if(xtensor_FOUND)
  set(XGALOIS_USE_EXTERNAL_XTENSOR ON)
else()
  set(XGALOIS_USE_EXTERNAL_XTENSOR OFF)
  FetchContent_Declare(
    xtensor
    URL https://github.com/xtensor-stack/xtensor/archive/refs/tags/0.26.0.tar.gz
    URL_HASH SHA256=f5f42267d850f781d71097b50567a480a82cd6875a5ec3e6238555e0ef987dc6
    SOURCE_SUBDIR include
  )
  FetchContent_MakeAvailable(xtensor)
endif()

# ------------------------------------------------------------------------------
# Header-only library target: xgalois::xgalois
# ------------------------------------------------------------------------------
add_library(xgalois INTERFACE)
add_library(xgalois::xgalois ALIAS xgalois)

target_compile_features(xgalois INTERFACE cxx_std_20)

target_include_directories(xgalois
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

if(XGALOIS_USE_EXTERNAL_XTL)
  target_link_libraries(xgalois INTERFACE xtl)
else()
  target_include_directories(xgalois
    SYSTEM INTERFACE
      "$<BUILD_INTERFACE:${xtl_SOURCE_DIR}/include>"
  )
endif()

if(XGALOIS_USE_EXTERNAL_XTENSOR)
  target_link_libraries(xgalois INTERFACE xtensor)
else()
  target_include_directories(xgalois
    SYSTEM INTERFACE
      "$<BUILD_INTERFACE:${xtensor_SOURCE_DIR}/include>"
  )
endif()

if(IS_ABSOLUTE "${CMAKE_INSTALL_DATADIR}")
  set(_xgalois_install_db_dir "${CMAKE_INSTALL_DATADIR}/xgalois/databases")
else()
  set(_xgalois_install_db_dir "$<INSTALL_PREFIX>/${CMAKE_INSTALL_DATADIR}/xgalois/databases")
endif()

target_compile_definitions(xgalois
  INTERFACE
    "$<BUILD_INTERFACE:XGALOIS_DEFAULT_DB_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/xgalois/databases\">"
    "$<INSTALL_INTERFACE:XGALOIS_DEFAULT_DB_DIR=\"${_xgalois_install_db_dir}\">"
)

if(XGALOIS_ENABLE_NATIVE)
  target_compile_options(xgalois
    INTERFACE
      "$<BUILD_INTERFACE:$<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-march=native>>"
  )
endif()

# ------------------------------------------------------------------------------
# Installation & CMake package config
# ------------------------------------------------------------------------------
set(_xgalois_cmake_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/xgalois")

install(
  TARGETS xgalois
  EXPORT xgaloisTargets
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/xgalois/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/xgalois"
  FILES_MATCHING PATTERN "*.hpp"
)

install(
  FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/xgalois/databases/conway_polys.db"
    "${CMAKE_CURRENT_SOURCE_DIR}/xgalois/databases/irreducible_polys.db"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/xgalois/databases"
)

if(NOT XGALOIS_USE_EXTERNAL_XTL)
  install(
    DIRECTORY "${xtl_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
endif()

if(NOT XGALOIS_USE_EXTERNAL_XTENSOR)
  install(
    DIRECTORY "${xtensor_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
endif()

install(
  EXPORT xgaloisTargets
  FILE xgaloisTargets.cmake
  NAMESPACE xgalois::
  DESTINATION "${_xgalois_cmake_install_dir}"
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/xgaloisConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/xgaloisConfig.cmake"
  INSTALL_DESTINATION "${_xgalois_cmake_install_dir}"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/xgaloisConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
  ARCH_INDEPENDENT
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/xgaloisConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/xgaloisConfigVersion.cmake"
  DESTINATION "${_xgalois_cmake_install_dir}"
)

# ------------------------------------------------------------------------------
# Tests and examples
# ------------------------------------------------------------------------------
if(XGALOIS_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(XGALOIS_BUILD_EXAMPLES)
  add_subdirectory(example)
endif()
