#!/bin/bash

# Prime Field Implementation Benchmark Runner
# Compares 2 different implementation classes across field sizes

set -e

echo "================================================="
echo "GF(p) Prime Field Implementation Benchmark"
echo "================================================="
echo "Comparing 2 implementation classes:"
echo "1. GaloisFieldPrime (Standard)"
echo "2. GaloisFieldPrimeTable (Logarithm Tables)"
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
RESULTS_DIR="benchmark/gf_prime/results"
BENCHMARK_FILTER=""

# Function to print usage
print_usage() {
    echo -e "${GREEN}GF Prime Field Benchmark Runner${NC}"
    echo ""
    echo "Usage: $0 [TEST_TYPE] [OPTIONS]"
    echo ""
    echo -e "${YELLOW}TEST_TYPES:${NC}"
    echo "  1  - All benchmarks (default)"
    echo "  2  - Small field tests only (GF(2) to GF(47))"
    echo "  3  - Medium field tests only (GF(53) to GF(127))"
    echo "  4  - Large field tests only (GF(131) to GF(997))"
    echo "  5  - Addition tests only (all field sizes)"
    echo "  6  - Multiplication tests only (all field sizes)"
    echo "  7  - Division tests only (all field sizes)"
    echo "  8  - Inversion tests only (all field sizes)"
    echo "  9  - Exponentiation tests only (all field sizes)"
    echo "  10 - Quick test (0.1s per benchmark)"
    echo "  11 - Memory usage analysis"
    echo "  12 - Performance comparison by implementation"
    echo "  13 - Standard implementation only"
    echo "  14 - Table implementation only"
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
    echo "  $0 13                 # Standard implementation only"
    echo "  $0 14                 # Table implementation only"
    echo "  $0 --time=0.5         # Run all with 0.5s per benchmark"
    echo "  $0 --format=json --output=results.json  # JSON output to file"
    echo "  $0 --filter='.*Small.*'                 # Custom filter"
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
    echo -e "${BLUE}Building GF Prime benchmark...${NC}"
    cd "$(dirname "$0")/../.."
    bazel build //benchmark/gf_prime/src:gf_prime_benchmark
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

    local cmd="bazel run //benchmark/gf_prime/src:gf_prime_benchmark --"

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
        [1-9]|1[0-4])
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
        echo -e "${GREEN}Running small field tests only (GF(2) to GF(47))...${NC}"
        FILTER=".*Small.*"
        ;;
    3)
        echo -e "${GREEN}Running medium field tests only (GF(53) to GF(127))...${NC}"
        FILTER=".*Medium.*"
        ;;
    4)
        echo -e "${GREEN}Running large field tests only (GF(131) to GF(997))...${NC}"
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
        echo -e "${GREEN}Running standard implementation only...${NC}"
        FILTER=".*BM_GFP_.*"
        ;;
    14)
        echo -e "${GREEN}Running table implementation only...${NC}"
        FILTER=".*BM_GFPTABLE_.*"
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
echo -e "${BLUE}  GF Prime Field Benchmarks${NC}"
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

# Generate basic analysis report if CSV output was used
if [ "$OUTPUT_FORMAT" == "csv" ] && [ -n "$OUTPUT_FILE_PATH" ] && [ -f "$OUTPUT_FILE_PATH" ]; then
    echo -e "${BLUE}Generating basic analysis...${NC}"

    # Create a simple analysis script
    cat > /tmp/analyze_gf_prime.py << 'EOF'
import csv
import sys
from collections import defaultdict

if len(sys.argv) != 2:
    print("Usage: python analyze_gf_prime.py <csv_file>")
    sys.exit(1)

csv_file = sys.argv[1]
results = defaultdict(list)

try:
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['name']
            time = float(row['real_time'])

            # Extract implementation and operation
            if 'BM_GFP_' in name and 'BM_GFPTABLE_' not in name:
                impl = 'Standard'
            elif 'BM_GFPTABLE_' in name:
                impl = 'Table'
            else:
                continue

            if 'Addition' in name:
                op = 'Addition'
            elif 'Multiplication' in name:
                op = 'Multiplication'
            elif 'Division' in name:
                op = 'Division'
            elif 'Inversion' in name:
                op = 'Inversion'
            elif 'Exponentiation' in name:
                op = 'Exponentiation'
            else:
                continue

            results[(impl, op)].append(time)

    print("\n=== GF Prime Field Performance Analysis ===")
    print("Average execution times (nanoseconds):")
    print("-" * 50)

    for (impl, op), times in sorted(results.items()):
        avg_time = sum(times) / len(times)
        print(f"{impl:10} {op:15}: {avg_time:8.2f} ns (n={len(times)})")

    print("\n=== Implementation Comparison ===")
    operations = set(op for impl, op in results.keys())

    for op in sorted(operations):
        std_times = results.get(('Standard', op), [])
        table_times = results.get(('Table', op), [])

        if std_times and table_times:
            std_avg = sum(std_times) / len(std_times)
            table_avg = sum(table_times) / len(table_times)

            if std_avg > 0:
                speedup = std_avg / table_avg
                print(f"{op:15}: Table is {speedup:5.2f}x {'faster' if speedup > 1 else 'slower'} than Standard")

except Exception as e:
    print(f"Error analyzing results: {e}")
EOF

    python3 /tmp/analyze_gf_prime.py "$OUTPUT_FILE_PATH"
    rm /tmp/analyze_gf_prime.py
fi
