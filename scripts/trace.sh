#!/usr/bin/env bash
# End-to-end requirement traceability with mantra. See docs/TRACEABILITY.md.
#
#   0. Self-check: tools/mantra-testrun's test suite runs the whole chain on a
#      fixture crate (tools/mantra-testrun/probe) whose 14 requirements have
#      known outcomes, so the verification semantics below are checked, not
#      assumed. SKIP_SELFCHECK=1 skips it.
#   1. tools/mantra-testrun: build the tests with LLVM instrumentation, run
#      every test in its own process with its own profile, export the lines
#      each test executed, write target/mantra/test-run.json (mantra's
#      TestRunSchema with per-test-case coverage).
#   2. `mantra collect` into a fresh mantra.db: docs/REQUIREMENTS.md, the
#      req_satisfied/req_test traces in src/, tests/, benches/, the reviews
#      and the test run above.
#   3. `mantra report`: write mantra-report/ (HTML + JSON) and print the
#      requirement states. Exit non-zero if any requirement is not Verified.
#
# A requirement is Verified only when every test carrying req_test("id") passed
# and executed at least one line of an item carrying req_satisfied("id") (or,
# for Manual requirements, a review record lists it). Parents need all
# non-optional children Verified.
#
# Tools: cargo + rustc's llvm-tools (llvm-cov, llvm-profdata) and mantra 0.9.
#   rustup component add llvm-tools      # or set LLVM_COV / LLVM_PROFDATA
#   cargo install --locked mantra
#   On NixOS: nix-shell -p cargo rustc gcc rustc.llvmPackages.llvm
set -euo pipefail
cd "$(dirname "$0")/.."

tool=tools/mantra-testrun
mantra=${MANTRA:-mantra}
command -v "$mantra" >/dev/null || { echo "mantra not found (cargo install --locked mantra, or set MANTRA)"; exit 2; }

if [ "${SKIP_SELFCHECK:-0}" != 1 ]; then
  echo "== self-check: mantra-testrun on its probe fixture"
  MANTRA="$mantra" cargo test --quiet --manifest-path "$tool/Cargo.toml"
fi

if [ "${SKIP_TESTS:-0}" != 1 ]; then
  echo "== test run with per-test-case coverage"
  # "$@" is passed to `cargo test --no-run` (e.g. --release, --features ...).
  cargo run --quiet --manifest-path "$tool/Cargo.toml" -- \
    --manifest-dir . --out target/mantra/test-run.json --name lipika-tests -- "$@" \
    || echo "(some tests failed; mantra will report them as Failed)"
fi

echo "== mantra collect + report"
rm -f mantra.db            # every run is evaluated from scratch, no history
"$mantra" collect
rm -rf mantra-report       # mantra refuses to overwrite an existing report directory
"$mantra" report --formats html --formats json --output-dir mantra-report/

summary=mantra-report/products/lipika@main.json5
reqs=mantra-report/products/lipika@main/requirements.json5
echo
echo "== requirement states (mantra-report/index.html for detail)"
awk -f scripts/summary.awk "$summary"

bad_ids() {
  # ids inside the `failed: [...]`, `skipped: [...]` and `unverified: [...]`
  # lists of requirements.json5
  awk '/^    (failed|skipped|unverified): \[/ { if ($0 ~ /\]/) next; on=1; next }
       on && /^    \],?$/ { on=0; next }
       on' "$reqs" | grep -oE '^\s*id: "[^"]*"' | sed -E 's/^\s*id: "//; s/"$//'
}
bad=$(bad_ids | grep -c . || true)
if [ "$bad" != 0 ]; then
  echo
  echo "$bad requirement(s) not verified:"
  bad_ids | sed 's/^/  - /'
  exit 1
fi
echo "all requirements verified"
