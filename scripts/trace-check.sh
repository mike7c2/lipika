#!/usr/bin/env bash
# Requirement-traceability check that needs only grep.
#
# Verifies that every `r[impl ID]` / `r[verify ID]` reference in src/, tests/
# and benches/ resolves to a requirement defined in docs/REQUIREMENTS.md, that
# every requirement has at least one impl and one verify reference, and that no
# ID is defined twice. tracey does the same (and more) interactively; this is
# the dependency-free gate for CI and pre-commit. Exit status is non-zero on
# any finding.
set -euo pipefail
cd "$(dirname "$0")/.."

spec=docs/REQUIREMENTS.md
code=(src tests benches)
status=0

defined=$(grep -oE '^r\[[A-Za-z0-9_.-]+(\+[0-9]+)?\]' "$spec" | sed -E 's/^r\[//; s/(\+[0-9]+)?\]$//' | sort)
impl=$(grep -rhoE 'r\[impl [A-Za-z0-9_.-]+' "${code[@]}" | sed 's/^r\[impl //' | sort -u)
verify=$(grep -rhoE 'r\[verify [A-Za-z0-9_.-]+' "${code[@]}" | sed 's/^r\[verify //' | sort -u)
refs=$(grep -rhoE 'r\[(impl|verify|depends|related) [A-Za-z0-9_.-]+' "${code[@]}" \
       | sed -E 's/^r\[(impl|verify|depends|related) //' | sort -u)

for id in $refs; do
  if ! printf '%s\n' "$defined" | grep -qx -- "$id"; then
    echo "BROKEN: r[$id] is referenced but not defined in $spec"
    grep -rnE "r\[(impl|verify|depends|related) $id\]" "${code[@]}" | sed 's/^/    /'
    status=1
  fi
done
for id in $defined; do
  printf '%s\n' "$impl"   | grep -qx -- "$id" || { echo "UNCOVERED: $id has no r[impl] reference";   status=1; }
  printf '%s\n' "$verify" | grep -qx -- "$id" || { echo "UNTESTED:  $id has no r[verify] reference"; status=1; }
done
dups=$(printf '%s\n' "$defined" | uniq -d)
[ -z "$dups" ] || { echo "DUPLICATE definitions: $dups"; status=1; }

n_def=$(printf '%s\n' "$defined" | grep -c . || true)
n_impl=$(printf '%s\n' "$impl" | grep -c . || true)
n_ver=$(printf '%s\n' "$verify" | grep -c . || true)
echo "requirements: $n_def defined, $n_impl with an implementation, $n_ver with a test"
exit $status
