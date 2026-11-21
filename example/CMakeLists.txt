function(xgalois_add_example target_name)
  add_executable(${target_name} ${ARGN})
  target_link_libraries(${target_name} PRIVATE xgalois::xgalois)
  target_compile_options(${target_name}
    PRIVATE
      $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
      $<$<CXX_COMPILER_ID:MSVC>:/W4>
  )
endfunction()

# Field examples
xgalois_add_example(gf_factory_example field/gf_factory_example.cpp)
xgalois_add_example(gf_prime_example field/gf_prime_example.cpp)
xgalois_add_example(gf_binary_example field/gf_binary_example.cpp)
xgalois_add_example(gf_extension_example field/gf_extension_example.cpp)

# Polynomial examples
xgalois_add_example(poly_example poly/poly_example.cpp)

# Linear algebra examples
xgalois_add_example(linalg_example linalg/linalg_example.cpp)

# Channel and coding examples
xgalois_add_example(channel_example channel/channel_example.cpp)
xgalois_add_example(coding_example coding/coding_example.cpp)
