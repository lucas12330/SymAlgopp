#!/usr/bin/env bash
# Vérifie le paquet : construit la bibliothèque (statique puis partagée), l'installe dans un
# préfixe temporaire, puis compile un projet extérieur avec find_package et avec pkg-config.
#
#   tests/packaging/verifier.sh [options cmake...]
#
# Les dépendances hors des chemins standard se donnent par les variables d'environnement
# CMAKE_PREFIX_PATH et PKG_CONFIG_PATH.
set -euo pipefail

racine=$(cd "$(dirname "$0")/../.." && pwd)
tmp=$(mktemp -d)
trap 'rm -r "$tmp"' EXIT
consommateur="$racine/tests/packaging/consommateur"

for partagee in OFF ON; do
    echo "=== bibliothèque BUILD_SHARED_LIBS=$partagee ==="
    prefixe="$tmp/prefixe-$partagee"
    cmake -S "$racine" -B "$tmp/build-$partagee" -G Ninja -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=$partagee -DSYMALGOPP_BUILD_TESTS=OFF -DSYMALGOPP_BUILD_EXAMPLES=OFF \
        -DCMAKE_INSTALL_PREFIX="$prefixe" "$@" > /dev/null
    cmake --build "$tmp/build-$partagee" > /dev/null
    cmake --install "$tmp/build-$partagee" > /dev/null

    echo "--- projet CMake extérieur (find_package) ---"
    CMAKE_PREFIX_PATH="$prefixe${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}" \
        cmake -S "$consommateur" -B "$tmp/conso-$partagee" -G Ninja "$@" > /dev/null
    cmake --build "$tmp/conso-$partagee" > /dev/null
    LD_LIBRARY_PATH="$prefixe/lib:$prefixe/lib64" "$tmp/conso-$partagee/consommateur"

    echo "--- pkg-config ---"
    export PKG_CONFIG_PATH="$prefixe/lib/pkgconfig:$prefixe/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
    if pkg-config --exists symalgopp; then
        # shellcheck disable=SC2046
        g++ "$consommateur/main.cpp" $(pkg-config --cflags symalgopp) $(pkg-config --libs symalgopp) \
            -o "$tmp/pc-$partagee"
        LD_LIBRARY_PATH="$prefixe/lib:$prefixe/lib64" "$tmp/pc-$partagee"
    else
        echo "symalgopp.pc introuvable" >&2
        exit 1
    fi
done
echo "Paquet vérifié."
