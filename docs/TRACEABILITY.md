# Requirements traceability (mantra)

Lipika's requirements live in [REQUIREMENTS.md](REQUIREMENTS.md). Code and
tests are tied to them with [mantra](https://github.com/mhatzl/mantra) traces,
and `scripts/trace.sh` runs the whole chain: tests with coverage, trace
collection, verification report. This document records the tool choice, the
conventions, what the first end-to-end run showed, and what mantra does and does
not prove.

This is the second of two passes over the same 34 requirements. The first
([PR #1](https://github.com/mike7c2/lipika/pull/1)) used
[tracey](https://github.com/bearcove/tracey): comment tags that a tool
cross-references, with no link to test execution. This pass demands mantra
because the goal is end-to-end traceability, not tag matching.

## Why mantra

|                        | mantra 0.9.0 (this branch)                                   | tracey 1.4.0 (PR #1)                                        |
|------------------------|--------------------------------------------------------------|-------------------------------------------------------------|
| Requirements           | Markdown `` ### `id`: Title `` + body, hierarchy via `Parents`| Markdown `r[id]` marker + paragraph                          |
| Code → requirement     | `#[mantra_macros::req_satisfied("id")]` on the item          | `// r[impl id]` comment                                     |
| Test → requirement     | `#[mantra_macros::req_test("id")]` on the test               | `// r[verify id]` comment                                   |
| What "verified" means  | a test tracing to the requirement **passed**, joined with the run's line coverage | a `verify` tag exists                   |
| Inputs beyond source   | JUnit (nextest) + Cobertura (cargo-llvm-cov)                 | none                                                        |
| Also models            | manual requirements + review records, optional/deprecated/excluded, parent roll-up | spec versions (`+N`), stale references, LSP/MCP |
| Output                 | HTML/JSON report, SQLite DB                                   | CLI queries, web dashboard                                  |
| Crate dependency       | `mantra-macros` (compile-time no-op macros)                   | none                                                        |

mantra's model is the one we want: a requirement's state is computed from the
test run, so a `req_test` on a test that no longer runs, is skipped, or fails
changes the requirement's state (Unverified / Skipped / Failed), and a parent is
only Verified when all its non-optional children are. tracey cannot express any
of that; its coverage is "a tag exists".

## Conventions

- **Spec**: [REQUIREMENTS.md](REQUIREMENTS.md). One `### `id`: Title` heading
  per requirement, optional property bullets (`- **Parents:** ["fn-4"]`,
  `- **Manual:** true`), then the MUST-statement.
- **IDs**: `fid-N`, `fn-N`, `nf-N`, `eq-eN`, children `fn-4-bpsk`,
  `fn-5-rinex-nav`. **Hyphens only, no dots inside an ID**: mantra treats
  `a.b` as "child `b` of `a`" at collection time and (0.9.0) does not detect
  traces to dotted IDs at all, so hierarchy is declared with `Parents`.
- **Satisfying code**: `#[mantra_macros::req_satisfied("id", ...)]` directly
  above the item, after its doc comment and above other attributes. Works on
  free functions, methods (including in trait impls), structs, enums, traits
  and modules; several IDs per attribute are fine. Prefer functions over types:
  a type has no executable lines, so coverage can never confirm it.
- **Verifying tests**: `#[mantra_macros::req_test("id", ...)]` above `#[test]`.
  Inside a body, `mantra_macros::assert_req!("id" => cond, "msg")` and
  `verify_req!("id" => expr)` also work (the `=> expr` is mandatory).
- **Manual requirements** (`nf-1` benchmark, `nf-4` third-party interop): no
  test can discharge them, so they are `Manual: true` and verified by a review
  record in `reviews/` that cites the evidence (README § Validation).
- **Equations**: every function in SIGNAL_MODEL.md §9 carries
  `req_satisfied("eq-eN")`, so the §9 table is checked, not maintained by hand.

## Running it

```bash
scripts/trace.sh                       # tests + coverage → mantra collect → mantra report
SKIP_TESTS=1 scripts/trace.sh          # reuse the last test run, re-collect only
open mantra-report/index.html          # per-requirement traces, coverage, test runs, reviews
```

Tools: `cargo-nextest`, `cargo-llvm-cov` (+ `llvm-tools`), `mantra`. Config:
[`mantra.json5`](../mantra.json5) (requirements, annotation paths, test-run
import, reviews) and [`.config/nextest.toml`](../.config/nextest.toml) (JUnit
output). The script exits non-zero if any requirement is Failed or Unverified.

## How the first end-to-end run played out

`scripts/trace.sh` on this branch: 50 tests pass under nextest with coverage,
mantra collects 34 requirements, 85 traces and the test run, and reports
**34 of 34 Verified** (32 by tests, 2 manual by review). Getting there took
four fixes, each of which is a property of the toolchain worth knowing:

1. **Dotted IDs are invisible to traces.** The spec initially used
   `fn-4.cw`-style children. mantra parsed them from the Markdown, but no
   `req_satisfied("fn-4.cw")` was ever collected, and `fid.1` failed outright
   with "Parent of requirement 'fid.1' was not collected". Children are now
   `fn-4-cw` with an explicit `Parents` property.
2. **Integration tests were not covered.** Every requirement verified only from
   `tests/*.rs` (fid-2, fid-4, fn-1, fn-4-bpsk, eq-e4, eq-e5, nf-3, and fn-4
   through its child) came back Unverified, while the same pattern in `src/`
   unit tests verified. Cause: cargo-llvm-cov excludes `tests/`, `benches/` and
   `examples/` from its report by default, so the verifying traces had no
   coverage. `scripts/trace.sh` passes `--no-default-ignore-filename-regex` and
   re-excludes only toolchain and registry paths.
3. **Parents inherit Unverified.** `fn-4` had its own passing test but stayed
   Unverified until `fn-4-bpsk` did, which is the roll-up we want and a reason to
   keep the hierarchy explicit.
4. **Small 0.9.0 rough edges.** `mantra report` refuses to overwrite an
   existing output directory (the script removes it first); `--formats markdown`
   fails with a missing template; `mantra-rust-macros` on docs.rs is the dead
   0.7 line and the 0.7.2 crate is yanked, so the live crates are `mantra` and
   `mantra-macros` 0.9.

### What mantra 0.9 does and does not prove

The thing to be clear-eyed about, because it is the reason to prefer mantra over
tracey. Verified in this pipeline means: *a test carrying `req_test("id")`
exists, ran in the recorded run, passed, and its own lines are in the coverage
report.* A controlled check in a scratch crate showed that the **satisfying
code's coverage is not consulted**: a requirement whose only `req_satisfied`
function had zero hits (in its own file, never called) still reported Verified,
as did one whose implementation was executed only by an untraced test. Coverage
is imported per test *run* (one Cobertura for all binaries), and at that
granularity mantra cannot attribute lines to the traced test, so it falls back
to "the verifying test ran and passed". Per-test-case coverage, which the
`verified` rules are written for ("every test covering a verifying trace covers
a satisfying trace"), needs one profile per test and a Cobertura per test case;
nextest and cargo-llvm-cov do not produce that out of the box.

So, relative to PR #1: mantra does tie requirements to test *execution* and
test *outcome*, rolls children into parents, and models manual verification.
It does not yet, with this toolchain, tie requirements to the execution of the
implementation they point at. Closing that gap is a tooling task (per-test
profraw collection, `LLVM_PROFILE_FILE` per test with nextest's
`NEXTEST_TEST_NAME`, one `llvm-cov export` per case) and mantra already has the
schema for it (`test_cases` with their own `covered_files`).

### Gaps the matrix exposed

Same as PR #1, because the requirement set is the same: the CW model had no test
(added), `fn-4` and `fn-5` had only indirect tests (direct ones added), the
sidecar JSON content is still untested (`fn-6` is verified at the format level),
and the GNSS-SDR evidence for `nf-4` is a manual review, not a test. `bin/lipika.rs`
carries no traces: nothing exercises the CLI in the test run, so a
`req_satisfied` there would be uncovered by construction.

## Adding a requirement

1. Add `### `group-n`: Title` and its MUST-statement to
   [REQUIREMENTS.md](REQUIREMENTS.md); give children a `Parents` bullet.
2. `#[mantra_macros::req_satisfied("group-n")]` on the code that satisfies it.
3. `#[mantra_macros::req_test("group-n")]` on the test that proves it (write
   the test if there isn't one).
4. `scripts/trace.sh`; open `mantra-report/index.html` and check the
   requirement's page shows both traces and the passing run.
