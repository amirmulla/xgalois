# GF(2^n) Binary Extension Field Arithmetic Benchmark Report

## Executive Summary

This report analyzes the performance of four different implementations of finite field arithmetic over GF(2^n) binary extension fields. The benchmark covers operations across small (n=2-8), medium (n=9-16), and large (n=17-20) field sizes, testing five fundamental operations: addition, multiplication, division, inversion, and exponentiation.

**Key Findings:**
- **GF2XZECH** (Zech logarithms) provides the best overall performance for multiplicative operations
- **GF2X** baseline is competitive only for addition operations
- **GF2XLOGOPT** offers an excellent balance of performance and implementation complexity
- Performance differences become more pronounced as field size increases

## Implementation Overview

### 1. GF2X (Baseline)
- Standard polynomial arithmetic implementation
- Direct computation without logarithmic optimizations
- Simple implementation but scales poorly for complex operations

### 2. GF2XLOG (Logarithmic)
- Uses discrete logarithm tables for multiplication and division
- Converts multiplicative operations to additive operations
- Significant performance improvement over baseline

### 3. GF2XLOGOPT (Optimized Logarithmic)
- Enhanced version of GF2XLOG with additional optimizations
- Better memory management and lookup table optimization
- Balanced performance across all operations

### 4. GF2XZECH (Zech Logarithms)
- Uses Zech logarithms for efficient field arithmetic
- Particularly optimized for inversion and exponentiation
- Highest memory usage but best performance for complex operations

## Performance Analysis by Field Size

### Small Fields (n=2-8, Field Orders: 4-256)

| Implementation | Addition | Multiplication | Division | Inversion | Exponentiation |
|:---------------|:--------:|:-------------:|:--------:|:---------:|:--------------:|
| GF2X           | 8.68     | 81.41         | 580.28   | 495.98    | 683.99         |
| GF2XLOG        | 8.48     | 21.72         | 23.71    | 12.10     | 21.70          |
| GF2XLOGOPT     | 8.46     | 13.29         | 14.16    | 10.69     | 21.19          |
| GF2XZECH       | 21.24    | 10.90         | 11.12    | 8.18      | 8.75           |

*Times in nanoseconds (ns) for n=8 (Field Order 256)*

### Medium Fields (n=9-16, Field Orders: 512-65,536)

| Implementation | Addition | Multiplication | Division | Inversion | Exponentiation |
|:---------------|:--------:|:-------------:|:--------:|:---------:|:--------------:|
| GF2X           | 8.48     | 158.62        | 2089.00  | 1938.97   | 2033.79        |
| GF2XLOG        | 8.48     | 21.56         | 24.44    | 12.19     | 17.12          |
| GF2XLOGOPT     | 8.50     | 13.46         | 14.22    | 10.83     | 17.07          |
| GF2XZECH       | 21.17    | 10.97         | 11.20    | 11.29     | 9.04           |

*Times in nanoseconds (ns) for n=16 (Field Order 65,536)*

### Large Fields (n=17-20, Field Orders: 131,072-1,048,576)

| Implementation | Addition | Multiplication | Division | Inversion | Exponentiation |
|:---------------|:--------:|:-------------:|:--------:|:---------:|:--------------:|
| GF2X           | 8.60     | 161.73        | 3049.14  | 2761.97   | 3177.38        |
| GF2XLOG        | 8.47     | 22.13         | 23.56    | 12.89     | 17.69          |
| GF2XLOGOPT     | 8.47     | 13.94         | 14.62    | 11.16     | 17.87          |
| GF2XZECH       | 21.34    | 10.97         | 11.20    | 8.29      | 9.04           |

*Times in nanoseconds (ns) for n=20 (Field Order ≈1,048,576)*

## Detailed Performance Comparison

### 1. Addition Performance
- **Best:** GF2X, GF2XLOG, GF2XLOGOPT (~8.4-8.6 ns) - All perform similarly
- **Worst:** GF2XZECH (~21.7 ns) - 2.6x slower due to Zech logarithm conversion overhead
- **Insight:** Addition is inherently fast in GF(2^n) and doesn't benefit from logarithmic optimizations

### 2. Multiplication Performance
- **Best:** GF2XZECH (10.90 ns) - Fastest across all field sizes
- **Good:** GF2XLOGOPT (13.94 ns) - 28% slower than Zech but still excellent
- **Decent:** GF2XLOG (22.13 ns) - Good improvement over baseline
- **Worst:** GF2X (161.73 ns) - 14.8x slower than best implementation

### 3. Division Performance
- **Best:** GF2XZECH (11.20 ns) - Leverages logarithmic properties effectively
- **Good:** GF2XLOGOPT (14.62 ns) - Actually fastest for division operations
- **Decent:** GF2XLOG (23.56 ns) - Significant improvement over baseline
- **Worst:** GF2X (3049.14 ns) - 207x slower, completely impractical for large fields

### 4. Inversion Performance
- **Best:** GF2XZECH (8.29 ns) - Zech logarithms excel at inversion
- **Good:** GF2XLOGOPT (11.16 ns) - 34% slower but still very fast
- **Decent:** GF2XLOG (12.89 ns) - Good logarithmic optimization
- **Worst:** GF2X (2761.97 ns) - 333x slower than best

### 5. Exponentiation Performance
- **Best:** GF2XZECH (9.04 ns) - Outstanding performance
- **Good:** GF2XLOG (17.69 ns) - Solid logarithmic approach
- **Decent:** GF2XLOGOPT (17.87 ns) - Similar to basic LOG
- **Worst:** GF2X (3177.38 ns) - 350x slower than best

## Memory Usage Analysis

Based on the MemoryPeak_KB data:

| Implementation | Small Fields | Medium Fields | Large Fields |
|:---------------|:------------:|:-------------:|:------------:|
| GF2X           | 6.2 MB       | 6.3 MB        | 10.1 MB      |
| GF2XLOG        | 6.3 MB       | 8.6-9.5 MB    | 69.8-94.4 MB |
| GF2XLOGOPT     | 6.3 MB       | 9.7-10.0 MB   | 94.4-103.3 MB |
| GF2XZECH       | 6.3 MB       | 10.0-10.1 MB  | 103.3-124.9 MB |

**Key Observations:**
- Memory usage increases significantly with field size for logarithmic methods
- GF2XZECH has the highest memory footprint but delivers best performance
- GF2X has the lowest memory usage but worst performance for complex operations

## Scalability Analysis

### Performance Scaling with Field Size

**Addition:** All implementations scale excellently - time remains constant regardless of field size.

**Multiplication:**
- GF2X: Poor scaling (26.8ns → 83.8ns → 154.6ns)
- GF2XLOG: Excellent scaling (remains ~22-23ns)
- GF2XLOGOPT: Excellent scaling (remains ~13-14ns)
- GF2XZECH: Excellent scaling (remains ~11ns)

**Division & Inversion:**
- GF2X: Catastrophic scaling (becomes impractical for large fields)
- All logarithmic methods: Excellent scaling with minimal time increase

## Recommendations

### For High-Performance Applications
**Primary Choice:** GF2XZECH
- Best overall performance for multiplication, division, inversion, and exponentiation
- Accept the memory overhead and slightly slower addition for massive gains in complex operations

### For Balanced Applications
**Primary Choice:** GF2XLOGOPT
- Excellent performance across all operations
- Good balance between speed and memory usage
- Easier to implement than Zech logarithms

### For Memory-Constrained Applications
**Primary Choice:** GF2XLOG
- Good performance improvement over baseline
- Moderate memory requirements
- Suitable when memory is limited but performance matters

### For Addition-Heavy Workloads
**Consider:** GF2X for addition operations combined with logarithmic methods for other operations
- Hybrid approach using best implementation for each operation type

## Implementation Guidelines

### When to Use Each Implementation

**Use GF2X when:**
- Only addition operations are needed
- Memory is extremely constrained
- Field sizes are very small (n ≤ 4)

**Use GF2XLOG when:**
- Balanced performance is needed with moderate memory usage
- Implementation simplicity is important
- Medium field sizes (n = 8-12)

**Use GF2XLOGOPT when:**
- High performance is needed with reasonable memory usage
- Multiple operation types are used frequently
- Field sizes range from small to large

**Use GF2XZECH when:**
- Maximum performance is required
- Memory usage is not a primary concern
- Large field sizes (n ≥ 16)
- Inversion and exponentiation are frequent operations

## Conclusion

The benchmark results clearly demonstrate that logarithmic-based implementations provide substantial performance improvements over the baseline GF2X implementation, especially for large field sizes. The choice between implementations should be based on the specific requirements:

- **Performance-critical applications:** Choose GF2XZECH
- **Balanced applications:** Choose GF2XLOGOPT
- **Memory-constrained applications:** Choose GF2XLOG
- **Addition-only applications:** GF2X is sufficient

The performance gap widens dramatically as field sizes increase, making the choice of implementation crucial for applications working with large finite fields. For any serious finite field arithmetic application beyond small field sizes, logarithmic-based implementations are essential.

---

*Report generated from benchmark data collected on June 23, 2025*
*Benchmark environment: macOS with clang++ compiler*
*All timing measurements are in nanoseconds (ns) representing time per operation*
