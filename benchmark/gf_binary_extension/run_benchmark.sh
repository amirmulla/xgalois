#!/bin/bash

# Binary Extension Field Implementation Benchmark Runner
# Compares 4 different implementation classes across field sizes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
BAZEL_BIN="${BAZEL:-bazel}"

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
OUTPUT_FILE=""
BENCHMARK_FILTER=""

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
    echo "  13 - GF(2^8) uint8_t comparison (GF2X vs GF2XLOG vs GF2XZECH)"
    echo ""
    echo -e "${YELLOW}OPTIONS:${NC}"
    echo "  --time=X        Set benchmark time per test (default: 1.0s)"
    echo "  --format=FORMAT Output format: console, json, csv (default: csv)"
    echo "  --output=FILE   Save results to file under ${RESULTS_DIR}/"
    echo "  --filter=REGEX  Custom benchmark filter"
    echo "  --help          Show this help message"
}

create_results_dir() {
    if [[ ! -d "${RESULTS_DIR}" ]]; then
        mkdir -p "${RESULTS_DIR}"
        echo -e "${GREEN}Created results directory: ${RESULTS_DIR}${NC}"
    fi
}

build_benchmark() {
    echo -e "${BLUE}Building GF Binary Extension benchmark (-c opt)...${NC}"
    (
        cd "${REPO_ROOT}"
        "${BAZEL_BIN}" build -c opt //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark
    )
    echo -e "${GREEN}Build successful!${NC}"
}

run_benchmark() {
    local filter="$1"
    local output_file="$2"

    echo -e "${BLUE}Running benchmark with filter: ${filter}${NC}"
    if [[ -n "${output_file}" ]]; then
        echo -e "${BLUE}Output file: ${output_file}${NC}"
    fi

    local cmd=(
        "${BAZEL_BIN}" run -c opt //benchmark/gf_binary_extension/src:gf_binary_extension_benchmark --
        "--benchmark_min_time=${BENCHMARK_TIME}s"
    )

    if [[ -n "${filter}" ]]; then
        cmd+=("--benchmark_filter=${filter}")
    fi

    if [[ "${OUTPUT_FORMAT}" != "console" ]]; then
        cmd+=("--benchmark_format=${OUTPUT_FORMAT}")
    fi

    echo -e "${YELLOW}Executing: ${cmd[*]}${NC}"

    if [[ -n "${output_file}" ]]; then
        mkdir -p "$(dirname "${output_file}")"
        local tmp_out
        tmp_out="$(mktemp)"
        trap 'rm -f "${tmp_out}"' RETURN
        (
            cd "${REPO_ROOT}"
            "${cmd[@]}" > "${tmp_out}"
        )
        mv "${tmp_out}" "${output_file}"
    else
        (
            cd "${REPO_ROOT}"
            "${cmd[@]}"
        )
    fi
}

# Parse command line arguments
TEST_TYPE=1
while [[ $# -gt 0 ]]; do
    case "$1" in
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
        [1-9]|1[0-3])
            TEST_TYPE="$1"
            shift
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            print_usage
            exit 1
            ;;
    esac
done

if [[ -n "${OUTPUT_FILE}" ]]; then
    create_results_dir
fi

build_benchmark

FILTER=""
case "${TEST_TYPE}" in
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
    13)
        echo -e "${GREEN}Running GF(2^8) uint8_t comparison...${NC}"
        FILTER=".*U8_GF256.*"
        ;;
esac

if [[ -n "${BENCHMARK_FILTER}" ]]; then
    FILTER="${BENCHMARK_FILTER}"
fi

OUTPUT_FILE_PATH=""
if [[ -n "${OUTPUT_FILE}" ]]; then
    OUTPUT_FILE_PATH="${RESULTS_DIR}/${OUTPUT_FILE}"
fi

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  GF Binary Extension Field Benchmarks${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "${YELLOW}Test Type: ${TEST_TYPE}${NC}"
echo -e "${YELLOW}Benchmark Time: ${BENCHMARK_TIME} seconds${NC}"
echo -e "${YELLOW}Output Format: ${OUTPUT_FORMAT}${NC}"
if [[ -n "${OUTPUT_FILE}" ]]; then
    echo -e "${YELLOW}Output File: ${OUTPUT_FILE_PATH}${NC}"
fi
echo -e "${BLUE}========================================${NC}"

run_benchmark "${FILTER}" "${OUTPUT_FILE_PATH}"

echo -e "${GREEN}Benchmark completed!${NC}"
if [[ -n "${OUTPUT_FILE}" ]]; then
    echo -e "${GREEN}Results saved to: ${OUTPUT_FILE_PATH}${NC}"
fi
