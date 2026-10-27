#!/usr/bin/env bash
# End-to-end requirement traceability with mantra.
#
#   1. Run the test suite under cargo-nextest with LLVM source coverage
#      (cargo-llvm-cov), producing target/nextest/default/junit.xml (pass/fail
#      per test) and target/nextest/default/cobertura.xml (line coverage).
#   2. `mantra collect`: read docs/REQUIREMENTS.md, the req_satisfied/req_test
#      traces in the code, and the test run above into mantra.db.
#   3. `mantra report`: write mantra-report/ (HTML + JSON) and print the
#      requirement states. Exit non-zero if any requirement is Failed or
#      Unverified (a requirement is Verified only when a test tracing to it
#      passed AND executed the lines tracing to it).
#
# Tools: cargo-nextest, cargo-llvm-cov (+ llvm-tools), mantra 0.9.
#   cargo install --locked cargo-nextest cargo-llvm-cov mantra
#   rustup component add llvm-tools   # or set LLVM_COV / LLVM_PROFDATA
#   On NixOS:  nix-shell -p cargo rustc gcc cargo-nextest cargo-llvm-cov \
#                rustc.llvmPackages.llvm  (and export LLVM_COV/LLVM_PROFDATA)
set -euo pipefail
cd "$(dirname "$0")/.."

out=target/nextest/default
mkdir -p "$out"

if [ "${SKIP_TESTS:-0}" != 1 ]; then
  # cargo-llvm-cov drops tests/, benches/ and examples/ from the report by default;
  # mantra needs the verifying traces in tests/*.rs to be covered, so disable that
  # default and re-exclude only toolchain and registry sources.
  cargo llvm-cov nextest \
    --no-default-ignore-filename-regex \
    --ignore-filename-regex '/(\.cargo|rustc|nix/store)/' \
    --cobertura --output-path "$out/cobertura.xml" "$@"
fi

mantra collect
rm -rf mantra-report   # mantra refuses to overwrite an existing report directory
mantra report --formats html --formats json --output-dir mantra-report/

summary=mantra-report/products/lipika@main.json5
reqs=mantra-report/products/lipika@main/requirements.json5
echo
echo "== requirement states (mantra-report/index.html for detail)"
awk '/summary: \{/,/^    \},/' "$summary" \
  | grep -E '^\s*(total|verified|unverified|failed|skipped|deprecated|excluded|manuals_total):' -A1 \
  | grep -v -- '--' | sed 's/^ *//' | paste - - | sed 's/{//; s/cnt: //; s/,$//; s/,\s*$//'

bad_ids() {
  # ids inside the `failed: [...]` and `unverified: [...]` lists of requirements.json5
  awk '/^    (failed|unverified): \[/ { if ($0 ~ /\]/) next; on=1; next }
       on && /^    \],?$/ { on=0; next }
       on' "$reqs" | grep -oE '^\s*id: "[^"]*"' | sed -E 's/^\s*id: "//; s/"$//'
}
bad=$(bad_ids | grep -c . || true)
if [ "$bad" != 0 ]; then
  echo
  echo "$bad requirement(s) failed or unverified:"
  bad_ids | sed 's/^/  - /'
  exit 1
fi
echo "all requirements verified"
