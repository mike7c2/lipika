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

Requires `cargo` and rustc's `llvm-tools` (`rustup component add llvm-tools`,
or set `LLVM_COV` / `LLVM_PROFDATA`). No nextest or cargo-llvm-cov.

`cargo test` here runs the self-check: the whole chain (this tool, `mantra
collect`, `mantra report`) on the [`probe/`](probe) fixture crate, asserting
the 14 documented requirement states in
[`probe/reqs/req.md`](probe/reqs/req.md). It needs `mantra` 0.9 on `PATH`
or in `$MANTRA`.
