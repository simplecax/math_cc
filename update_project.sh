#!/bin/bash
#
# This script applies the definitive fix to the rank calculation logic,
# resolving the final test failures.
#
# USAGE: ./final_fix.sh from the project root.

set -e

echo "🚀 Correcting the rank calculation logic in src/atopo_linbox.cpp..."

# Replace the incorrect comparison (!= 1) with the correct one (!= 0).
# We create an Integer(0) object for a type-safe comparison.
sed -i 's|if (factor_pair.first != one)|Integer zero(0);\n        if (factor_pair.first != zero)|' src/atopo_linbox.cpp

# Also remove the now-unused "one" variable
sed -i '/Integer one(1);/d' src/atopo_linbox.cpp


echo "✅ Rank calculation has been corrected."
echo "   Please clean and rebuild. All tests should now pass."
echo "   Run: rm -rf build && export PKG_CONFIG_PATH=/opt/lib/pkgconfig:\$PKG_CONFIG_PATH && cmake -S . -B build && cmake --build build && ctest --test-dir build"