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
#      requirement states.
#   4. mantra-verdict: parse the JSON report and exit non-zero if any
#      requirement is not Verified (or the report does not add up).
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
  # mantra collects every *.json here; never let an earlier run stand in for
  # this one.
  rm -f target/mantra/*.json
  # "$@" is passed to `cargo test --no-run` (e.g. --release, --features ...).
  rc=0
  cargo run --quiet --manifest-path "$tool/Cargo.toml" --bin mantra-testrun -- \
    --manifest-dir . --out target/mantra/test-run.json --name lipika-tests -- "$@" \
    || rc=$?
  case $rc in
    0) ;;
    1) echo "(some tests failed or lack coverage; mantra will not verify their requirements)" ;;
    *) echo "mantra-testrun did not produce a test run (exit $rc); stopping"; exit 2 ;;
  esac
fi

echo "== mantra collect + report"
rm -f mantra.db            # every run is evaluated from scratch, no history
"$mantra" collect
rm -rf mantra-report       # mantra refuses to overwrite an existing report directory
"$mantra" report --formats html --formats json --output-dir mantra-report/

summary=mantra-report/products/lipika@main.json5
echo
echo "== requirement states (mantra-report/index.html for detail)"
awk -f scripts/summary.awk "$summary"

echo
# Exit status is the script's: 1 if any requirement is not verified, 2 if the
# report cannot be read or does not add up.
cargo run --quiet --manifest-path "$tool/Cargo.toml" --bin mantra-verdict -- \
  mantra-report/products/lipika@main
