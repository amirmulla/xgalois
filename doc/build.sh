#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "Building xgalois documentation..."

# Switch to the doc directory
cd "$(dirname "$0")"

# Activate the virtual environment
if [ -d ".venv" ]; then
    echo "Activating virtual environment..."
    source .venv/bin/activate
else
    echo "Virtual environment not found. Please run: python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt"
    exit 1
fi

# Run sphinx-build
sphinx-build -b html . _build/html --keep-going

echo ""
echo "Documentation successfully built at: _build/html/index.html"
