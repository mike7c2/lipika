# mantra-testrun

Run every Rust test of a cargo package in its own process with its own LLVM
coverage profile and write a [mantra](https://github.com/mhatzl/mantra)
`TestRunSchema` JSON in which each test case carries exactly the lines it
executed. Lipika's `scripts/trace.sh` uses it; see
[`docs/TRACEABILITY.md`](../../docs/TRACEABILITY.md) for why per-test-case
coverage is what makes mantra's **Verified** trustworthy.

```bash
cargo run --manifest-path tools/mantra-testrun/Cargo.toml -- \
    --manifest-dir . --out target/mantra/test-run.json --name lipika-tests -- --release
```

Options: `--timeout <SECS>` (default 600) kills a test that runs longer and
records it as failed; `--jobs <N>` sets parallelism. Exit status: 0 if every
test passed with coverage, 1 if any failed, timed out or passed without
coverage data (the JSON is still written), 2 on tool errors, an
uninstrumented build or a source file changed during the run (no JSON).

`mantra-verdict <mantra-report>/products/<product>` is the pass/fail gate over
mantra's JSON report: exit 0 if every requirement is Verified (or Excluded /
Deprecated), 1 if not (listing them), 2 if the report cannot be read or its
lists do not add up to its total.

Requires `cargo` and rustc's `llvm-tools` (`rustup component add llvm-tools`,
or set `LLVM_COV` / `LLVM_PROFDATA`). No nextest or cargo-llvm-cov.

`cargo test` here runs the self-check: the whole chain (this tool, `mantra
collect`, `mantra report`) on the [`probe/`](probe) fixture crate, asserting
the 14 documented requirement states in
[`probe/reqs/req.md`](probe/reqs/req.md), that a test with more output than
the kept log tail still counts as failed, that a hanging test is killed, and
that `mantra-verdict` names exactly the non-verified probe requirements. It needs `mantra` 0.9 on `PATH`
or in `$MANTRA`.
