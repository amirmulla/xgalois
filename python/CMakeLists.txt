cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(Threads REQUIRED)

# Attempt to locate nanobind installed in the active Python environment first.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_VARIABLE _xgalois_nanobind_cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  RESULT_VARIABLE _xgalois_nanobind_probe_rc
  ERROR_QUIET
)

if(_xgalois_nanobind_probe_rc EQUAL 0 AND EXISTS "${_xgalois_nanobind_cmake_dir}")
  list(PREPEND CMAKE_PREFIX_PATH "${_xgalois_nanobind_cmake_dir}")
endif()

find_package(nanobind CONFIG QUIET)
if(NOT nanobind_FOUND)
  include(FetchContent)
  message(STATUS "nanobind not found via find_package; fetching wjakob/nanobind v2.4.0")
  FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG v2.4.0
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(nanobind)
endif()

nanobind_add_module(
  _xgalois
  STABLE_ABI
  NB_STATIC
  src/bindings.cpp
)

target_compile_features(_xgalois PRIVATE cxx_std_20)
target_link_libraries(_xgalois PRIVATE xgalois::xgalois Threads::Threads)
target_compile_options(_xgalois
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

install(TARGETS _xgalois LIBRARY DESTINATION xgalois)
