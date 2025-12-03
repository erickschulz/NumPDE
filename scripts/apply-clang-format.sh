#!/usr/bin/env bash

# Apply clang-format to all C++ source files
# Must be run from the root of the repository
# Requires clang-format version 18

set -e

REQUIRED_VERSION=18

# Check if clang-format-18 is available
if command -v clang-format-18 &> /dev/null; then
    CLANG_FORMAT="clang-format-18"
elif command -v clang-format &> /dev/null; then
    # Check if the generic clang-format is version 18
    # Extract version using sed (works on both macOS and Linux)
    VERSION=$(clang-format --version | sed -E 's/.*version ([0-9]+)\..*/\1/')
    if [ "$VERSION" != "$REQUIRED_VERSION" ]; then
        echo "Error: This project requires clang-format version $REQUIRED_VERSION"
        echo "Found: $(clang-format --version)"
        echo ""
        echo "Install clang-format-18:"
        echo "  macOS:  brew install llvm@18"
        echo "  Ubuntu: apt-get install clang-format-18"
        exit 1
    fi
    CLANG_FORMAT="clang-format"
else
    echo "Error: clang-format not found"
    echo "Install clang-format-18:"
    echo "  macOS:  brew install llvm@18"
    echo "  Ubuntu: apt-get install clang-format-18"
    exit 1
fi

echo "Applying $CLANG_FORMAT to all C++ files..."

find developers lecturecodes include \
    \( -iname "*.h" -o -iname "*.hpp" -o -iname "*.c" -o -iname "*.cpp" -o -iname "*.cc" \) \
    -exec $CLANG_FORMAT -i {} +

echo "Done! All files formatted."