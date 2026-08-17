#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
# -----------------------------------------------------------------------------
#  Build libSBML 5.20.5 into third_party/install (macOS and Linux).
#
#  DyCelFEM originally linked against libSBML 3.x, which cannot be configured on
#  a modern toolchain (its 2009 autotools config.guess has no arm64-apple-darwin
#  triple and the tree predates clang). Every libSBML entry point the code uses
#  is API-identical in 5.x: readSBML, the Model/Species/Reaction/KineticLaw
#  accessors, SBase_getTypeCode, the AST_* node-type enum, and the Rule
#  subclasses.
#
#  WITH_CPP_NAMESPACE=OFF matters: libSBML 4.0 introduced the optional 'libsbml'
#  C++ namespace, and this code predates it and expects the flat global one.
# -----------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TP="$ROOT/third_party"
VER="${LIBSBML_VERSION:-5.20.5}"
URL="https://github.com/sbmlteam/libsbml/archive/refs/tags/v$VER.tar.gz"

need() { command -v "$1" >/dev/null || { echo "ERROR: $1 not found. $2"; exit 1; }; }
need cmake "Install it: macOS 'brew install cmake', Debian 'apt install cmake'."
command -v curl >/dev/null || command -v wget >/dev/null || {
	echo "ERROR: need curl or wget to download libSBML."; exit 1; }

case "$(uname -s)" in
  Linux)
    if [ ! -f /usr/include/libxml2/libxml/parser.h ] && \
       [ ! -f /usr/include/libxml/parser.h ]; then
      echo "NOTE: libxml2 headers not found. On Debian/Ubuntu:"
      echo "        sudo apt install libxml2-dev zlib1g-dev libbz2-dev"
      echo "      On Fedora/RHEL:"
      echo "        sudo dnf install libxml2-devel zlib-devel bzip2-devel"
    fi
    ;;
esac

mkdir -p "$TP"
cd "$TP"

if [ ! -d "libsbml-$VER" ]; then
	echo "==> downloading libSBML $VER"
	if command -v curl >/dev/null; then
		curl -L --retry 3 -o "libsbml-$VER.tar.gz" "$URL"
	else
		wget -O "libsbml-$VER.tar.gz" "$URL"
	fi
	tar xzf "libsbml-$VER.tar.gz"
fi

NPROC="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4)"

echo "==> configuring"
mkdir -p build-libsbml
cd build-libsbml
CMAKE_EXTRA=()
[ "$(uname -s)" = "Darwin" ] && CMAKE_EXTRA+=(-DCMAKE_OSX_ARCHITECTURES="$(uname -m)")

cmake "../libsbml-$VER" \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX="$TP/install" \
	-DWITH_CPP_NAMESPACE=OFF \
	-DWITH_SWIG=OFF \
	-DWITH_EXAMPLES=OFF \
	-DWITH_CHECK=OFF \
	-DWITH_LIBXML=ON \
	-DWITH_EXPAT=OFF \
	-DWITH_XERCES=OFF \
	-DENABLE_LAYOUT=ON \
	-DBUILD_SHARED_LIBS=ON \
	"${CMAKE_EXTRA[@]}"

echo "==> building with $NPROC jobs (a few minutes)"
cmake --build . -j"$NPROC"
cmake --install .

echo
echo "==> libSBML $VER installed in $TP/install"
echo "    Now run: make"
