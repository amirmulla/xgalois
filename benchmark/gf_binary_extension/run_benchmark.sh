#!/bin/bash

# Binary Extension Field Implementation Benchmark Runner
# Compares 4 different implementation classes across field sizes

set -e

echo "================================================="
echo "GF(2^m) Binary Extension Field Implementation Benchmark"
echo "================================================="
echo "Comparing 4 implementation classes:"
echo "1. GaloisFieldBinaryExtension (Base)"
echo "2. GFBELogTables (Log Tables)"
echo "3. GFBELogTablesOpt (Optimized Log Tables)"
echo "4. GFBEZechLogTables (Zech Log Tables)"
echo "================================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
BENCHMARK_TIME="1.0"
OUTPUT_FORMAT="csv"
RESULTS_DIR="results"
BENCHMARK_FILTER=""

# Function to print usage
print_usage() {
    echo -e "${GREEN}GF Binary Extension Field Benchmark Runner${NC}"
    echo ""
    echo "Usage: $0 [TEST_TYPE] [OPTIONS]"
    echo ""
    echo -e "${YELLOW}TEST_TYPES:${NC}"
    echo "  1  - All benchmarks (default)"
    echo "  2  - Small field tests only (GF(2^2) to GF(2^8))"
    echo "  3  - Medium field tests only (GF(2^9) to GF(2^16))"
    echo "  4  - Large field tests only (GF(2^17) to GF(2^20))"
    echo "  5  - Addition tests only (all field sizes)"
    echo "  6  - Multiplication tests only (all field sizes)"
    echo "  7  - Division tests only (all field sizes)"
    echo "  8  - Inversion tests only (all field sizes)"
    echo "  9  - Exponentiation tests only (all field sizes)"
    echo "  10 - Quick test (0.1s per benchmark)"
    echo "  11 - Memory usage analysis"
    echo "  12 - Performance comparison by implementation"
    echo ""
    echo -e "${YELLOW}OPTIONS:${NC}"
    echo "  --time=X        Set benchmark time per test (default: 1.0s)"
    echo "  --format=FORMAT Output format: console, json, csv (default: csv)"
    echo "  --output=FILE   Save results to file"
    echo "  --filter=REGEX  Custom benchmark filter"
    echo "  --help         Show this help message"
    echo ""
    echo -e "${YELLOW}EXAMPLES:${NC}"
    echo "  $0                    # Run all benchmarks"
    echo "  $0 2                  # Run small field tests only"
    echo "  $0 6                  # Run multiplication tests only"
    echo "  $0 10                 # Quick test run"
    echo "  $0 --time=0.5         # Run all with 0.5s per benchmark"
    echo "  $0 --format=json --output=results.json  # JSON output to file"
    echo "  $0 --filter='.*Small8.*'               # Custom filter"
}

# Function to create results directory
create_results_dir() {
    if [ ! -d "$RESULTS_DIR" ]; then
        mkdir -p "$RESULTS_DIR"
        echo -e "${GREEN}Created results directory: $RESULTS_DIR${NC}"
    fi
}

# Function to build the benchmark
build_benchmark() {
    echo -e "${BLUE}Building GF Binary Extension benchmark...${NC}"
    cd "$(dirname "$0")/../.."
    bazel build //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Build successful!${NC}"
    else
        echo -e "${RED}Build failed!${NC}"
        exit 1
    fi
}

# Function to run benchmark with given parameters
run_benchmark() {
    local filter="$1"
    local output_file="$2"

    echo -e "${BLUE}Running benchmark with filter: $filter${NC}"
    if [ -n "$output_file" ]; then
        echo -e "${BLUE}Output file: $output_file${NC}"
    fi

    local cmd="bazel run //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark --"

    if [ -n "$filter" ]; then
        cmd="$cmd --benchmark_filter=\"$filter\""
    fi

    cmd="$cmd --benchmark_min_time=${BENCHMARK_TIME}s"

    if [ "$OUTPUT_FORMAT" != "console" ]; then
        cmd="$cmd --benchmark_format=$OUTPUT_FORMAT"
    fi

    echo -e "${YELLOW}Executing: $cmd${NC}"

    if [ -n "$output_file" ]; then
        # Ensure the directory exists
        mkdir -p "$(dirname "$output_file")"
        eval $cmd > "$output_file"
    else
        eval $cmd
    fi
}

# Parse command line arguments
TEST_TYPE=1
while [[ $# -gt 0 ]]; do
    case $1 in
        --help)
            print_usage
            exit 0
            ;;
        --time=*)
            BENCHMARK_TIME="${1#*=}"
            shift
            ;;
        --format=*)
            OUTPUT_FORMAT="${1#*=}"
            shift
            ;;
        --output=*)
            OUTPUT_FILE="${1#*=}"
            shift
            ;;
        --filter=*)
            BENCHMARK_FILTER="${1#*=}"
            shift
            ;;
        [1-9]|1[0-2])
            TEST_TYPE=$1
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            print_usage
            exit 1
            ;;
    esac
done

# Create results directory if needed
if [ -n "$OUTPUT_FILE" ]; then
    create_results_dir
fi

# Build the benchmark
build_benchmark

# Set benchmark filter based on test type
case $TEST_TYPE in
    1)
        echo -e "${GREEN}Running all benchmarks...${NC}"
        FILTER=""
        ;;
    2)
        echo -e "${GREEN}Running small field tests only (GF(2^2) to GF(2^8))...${NC}"
        FILTER=".*Small.*"
        ;;
    3)
        echo -e "${GREEN}Running medium field tests only (GF(2^9) to GF(2^16))...${NC}"
        FILTER=".*Medium.*"
        ;;
    4)
        echo -e "${GREEN}Running large field tests only (GF(2^17) to GF(2^20))...${NC}"
        FILTER=".*Large.*"
        ;;
    5)
        echo -e "${GREEN}Running addition tests only (all field sizes)...${NC}"
        FILTER=".*Addition.*"
        ;;
    6)
        echo -e "${GREEN}Running multiplication tests only (all field sizes)...${NC}"
        FILTER=".*Multiplication.*"
        ;;
    7)
        echo -e "${GREEN}Running division tests only (all field sizes)...${NC}"
        FILTER=".*Division.*"
        ;;
    8)
        echo -e "${GREEN}Running inversion tests only (all field sizes)...${NC}"
        FILTER=".*Inversion.*"
        ;;
    9)
        echo -e "${GREEN}Running exponentiation tests only (all field sizes)...${NC}"
        FILTER=".*Exponentiation.*"
        ;;
    10)
        echo -e "${GREEN}Running quick test (0.1s per benchmark)...${NC}"
        BENCHMARK_TIME="0.1"
        FILTER=""
        ;;
    11)
        echo -e "${GREEN}Running memory usage analysis...${NC}"
        FILTER=".*Addition.*"
        OUTPUT_FORMAT="json"
        ;;
    12)
        echo -e "${GREEN}Running performance comparison by implementation...${NC}"
        FILTER=".*Multiplication.*"
        OUTPUT_FORMAT="csv"
        ;;
esac

# Use custom filter if provided
if [ -n "$BENCHMARK_FILTER" ]; then
    FILTER="$BENCHMARK_FILTER"
fi

# Prepare output file path
OUTPUT_FILE_PATH=""
if [ -n "$OUTPUT_FILE" ]; then
    OUTPUT_FILE_PATH="$RESULTS_DIR/$OUTPUT_FILE"
fi

# Run the benchmark
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  GF Binary Extension Field Benchmarks${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${YELLOW}Test Type: $TEST_TYPE${NC}"
echo -e "${YELLOW}Benchmark Time: $BENCHMARK_TIME seconds${NC}"
echo -e "${YELLOW}Output Format: $OUTPUT_FORMAT${NC}"
if [ -n "$OUTPUT_FILE" ]; then
    echo -e "${YELLOW}Output File: $OUTPUT_FILE_PATH${NC}"
fi
echo -e "${BLUE}========================================${NC}"

run_benchmark "$FILTER" "$OUTPUT_FILE_PATH"

echo -e "${GREEN}Benchmark completed!${NC}"
if [ -n "$OUTPUT_FILE" ]; then
    echo -e "${GREEN}Results saved to: $OUTPUT_FILE_PATH${NC}"
fi
