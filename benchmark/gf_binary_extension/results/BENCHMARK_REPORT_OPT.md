# GF2X Binary Extension Field Benchmark Report (Optimized)

## Executive Summary

This report presents comprehensive benchmark results for different implementations of binary extension field (GF2X) operations. The benchmarks compare four different implementation approaches:

- **GF2X**: Standard polynomial representation
- **GF2XLOG**: Logarithmic representation
- **GF2XLOGOPT**: Optimized logarithmic representation
- **GF2XZECH**: Zech logarithm representation

The results are categorized by field size: Small (2²-2⁸), Medium (2⁹-2¹⁶), and Large (2¹⁷-2²⁰).

## Test Environment

- **Date**: June 24, 2025
- **Time Unit**: Nanoseconds (ns)
- **Memory Usage**: Tracked in KB
- **Metrics**: Real time, CPU time, iterations, and peak memory usage

## Performance Overview

### Small Fields (2² to 2⁸)

#### Addition Performance
Updated results based on new benchmarks:
- **GF2X**: ~1.36 ns (consistent performance)
- **GF2XLOG**: ~1.36 ns (comparable to GF2X)
- **GF2XLOGOPT**: ~1.36 ns (comparable to GF2X)
- **GF2XZECH**: ~1.38 ns (slightly slower due to logarithm operations)

#### Multiplication Performance
Updated results based on new benchmarks:
- **GF2X**: 8.85-31.02 ns (increases with field size)
- **GF2XLOG**: ~1.90 ns (constant, excellent for small fields)
- **GF2XLOGOPT**: ~1.56 ns (best performance)
- **GF2XZECH**: ~1.93 ns (competitive)

#### Division Performance
Updated results based on new benchmarks:
- **GF2X**: 21.29-106.36 ns (significantly slower)
- **GF2XLOG**: ~2.08 ns (excellent)
- **GF2XLOGOPT**: ~1.92 ns (best performance)
- **GF2XZECH**: ~1.92 ns (competitive)

#### Inversion Performance
Updated results based on new benchmarks:
- **GF2X**: 16.17-75.90 ns (moderate performance)
- **GF2XLOG**: ~1.63 ns (excellent)
- **GF2XLOGOPT**: ~1.55 ns (best performance)
- **GF2XZECH**: ~1.64 ns (competitive)

#### Exponentiation Performance
Updated results based on new benchmarks:
- **GF2X**: 62.78-265.62 ns (poor scaling)
- **GF2XLOG**: ~1.43 ns (excellent)
- **GF2XLOGOPT**: ~1.49 ns (excellent)
- **GF2XZECH**: ~1.46 ns (best performance)

### Medium Fields (2⁹ to 2¹⁶)

#### Key Observations
- **GF2X**: Addition remains fast (~1.36 ns), but multiplication/division scale poorly
- **GF2XLOG**: Maintains excellent performance across all operations
- **GF2XLOGOPT**: Slightly better than GF2XLOG in most cases
- **GF2XZECH**: Consistent performance, particularly strong in inversion/exponentiation

### Large Fields (2¹⁷ to 2²⁰)

#### Performance Trends
- **GF2X**: Severe performance degradation for complex operations
- **Logarithmic approaches**: Maintain relatively stable performance
- **Memory usage**: Increases significantly for larger fields, especially with logarithmic tables

## Detailed Results

### Small Fields Performance Summary

| Operation | GF2X | GF2XLOG | GF2XLOGOPT | GF2XZECH |
|-----------|------|---------|------------|----------|
| Addition | 1.36 ns | 1.36 ns | 1.36 ns | 1.38 ns |
| Multiplication | 8.85-31.02 ns | 1.90 ns | 1.56 ns | 1.93 ns |
| Division | 21.29-106.36 ns | 2.08 ns | 1.92 ns | 1.92 ns |
| Inversion | 16.17-75.90 ns | 1.63 ns | 1.55 ns | 1.64 ns |
| Exponentiation | 62.78-265.62 ns | 1.43 ns | 1.49 ns | 1.46 ns |

### Memory Usage Analysis
Updated memory usage based on new benchmarks:
- **Small fields**: 6.1-6.4 MB baseline
- **Medium fields**: 6.4-12.6 MB (logarithmic tables increase memory)
- **Large fields**: Up to 107.8 MB (significant memory overhead for large logarithmic tables)

### Scaling Analysis

#### GF2X (Polynomial)
- **Strengths**: Minimal memory overhead, excellent addition performance
- **Weaknesses**: Poor scaling for multiplication, division, and exponentiation

#### GF2XLOG (Logarithmic)
- **Strengths**: Consistent performance across all operations
- **Weaknesses**: High memory usage for large fields

#### GF2XLOGOPT (Optimized Logarithmic)
- **Strengths**: Best overall performance for most operations
- **Weaknesses**: Similar memory overhead to GF2XLOG

#### GF2XZECH (Zech Logarithms)
- **Strengths**: Excellent inversion and exponentiation performance
- **Weaknesses**: Slightly slower addition compared to other methods

## Recommendations

### For Small Fields (2² to 2⁸)
- **General purpose**: Use **GF2XLOGOPT** for best overall performance
- **Memory constrained**: Use **GF2X** if only addition is needed
- **Inversion/Exponentiation heavy**: Use **GF2XZECH**

### For Medium Fields (2⁹ to 2¹⁶)
- **Recommended**: **GF2XLOGOPT** provides the best balance of performance and memory usage
- **Alternative**: **GF2XZECH** for inversion-heavy workloads

### For Large Fields (2¹⁷ to 2²⁰)
- **Performance critical**: **GF2XLOGOPT** or **GF2XZECH**
- **Memory critical**: **GF2X** only for addition-heavy workloads
- **Balanced**: Consider memory vs. performance trade-offs carefully

## Conclusion

The benchmark results demonstrate that logarithmic representations (GF2XLOG, GF2XLOGOPT, GF2XZECH) significantly outperform polynomial representation (GF2X) for complex operations like multiplication, division, inversion, and exponentiation. The optimized logarithmic implementation (GF2XLOGOPT) provides the best overall performance across most operations, while Zech logarithms (GF2XZECH) excel in inversion and exponentiation operations.

The trade-off between performance and memory usage becomes more pronounced for larger fields, where logarithmic tables can consume significant memory but provide substantial performance benefits for complex operations.

---

*Benchmark conducted on June 24, 2025 using Google Benchmark framework.*
