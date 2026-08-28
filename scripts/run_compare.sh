#!/usr/bin/env bash
# run_compare.sh -- build if needed, then compare two productions.
#
#   ./scripts/run_compare.sh A.root B.root
#   ./scripts/run_compare.sh 'dirA/*.root' 'dirB/*.root' --label-a nominal
#   ./scripts/run_compare.sh A.root --closure
#
# Glob patterns are expanded into an @list file, which keeps the command line
# short when a sample has dozens of files.

set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build="${FLUXVAL_BUILD:-$here/build}"

if [[ $# -lt 1 ]]; then
  echo "usage: $0 A [B] [options...]" >&2
  exit 2
fi

if [[ ! -x "$build/flux_compare" ]]; then
  echo "[build] configuring in $build"
  cmake -S "$here" -B "$build" -DCMAKE_BUILD_TYPE=RelWithDebInfo
  cmake --build "$build" -j "$(nproc 2>/dev/null || echo 4)"
  ctest --test-dir "$build" --output-on-failure
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Turn a glob or a directory into an @list file; leave anything else alone.
expand() {
  local spec="$1" tag="$2"
  if [[ "$spec" == @* || -f "$spec" && "$spec" != *"*"* ]]; then
    echo "$spec"; return
  fi
  local list="$tmp/$tag.txt"
  if [[ -d "$spec" ]]; then
    find "$spec" -name '*.root' | sort > "$list"
  else
    # shellcheck disable=SC2086
    ls -1 $spec 2>/dev/null | sort > "$list"
  fi
  [[ -s "$list" ]] || { echo "no files matched: $spec" >&2; exit 1; }
  echo "@$list"
}

A="$(expand "$1" A)"; shift
B=""
if [[ $# -ge 1 && "$1" != -* ]]; then B="$(expand "$1" B)"; shift; fi

exec "$build/flux_compare" "$A" ${B:+"$B"} "$@"
