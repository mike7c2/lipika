# Requirements traceability (mantra)

Lipika's requirements live in [REQUIREMENTS.md](REQUIREMENTS.md). Code and
tests are tied to them with [mantra](https://github.com/mhatzl/mantra) traces,
and `scripts/trace.sh` runs the whole chain: self-check, tests with per-test
coverage, trace collection, verification report. This document records the tool
choice, the conventions, what each end-to-end run showed, exactly what
"Verified" proves, and what it does not.

This is the second of two passes over the same 34 requirements. The first
([PR #1](https://github.com/mike7c2/lipika/pull/1)) used
[tracey](https://github.com/bearcove/tracey): comment tags that a tool
cross-references, with no link to test execution. This pass demands mantra
because the goal is end-to-end traceability, not tag matching.

## The chain

```
docs/REQUIREMENTS.md          ### `fn-4-lnav`: Full LNAV  (+ Parents, Manual)
        │
        ├── #[mantra_macros::req_satisfied("fn-4-lnav")]   on the implementing items
        │
        ├── #[mantra_macros::req_test("fn-4-lnav")]        on the verifying tests
        │
        ├── tools/mantra-testrun → target/mantra/test-run.json
        │       every test in its own process, its own LLVM profile:
        │       state (passed/failed/skipped) + the lines *it* executed
        │
        ├── reviews/*.json5                                 for Manual requirements
        │
        └── mantra collect → mantra report                  mantra-report/index.html
```

A requirement is **Verified** when every test carrying `req_test("id")` passed
and executed at least one line of an item carrying `req_satisfied("id")`; a
Manual requirement is Verified when a review record lists it; a parent is
Verified when all its non-optional children are. Anything else is Unverified or
Failed and `scripts/trace.sh` exits non-zero. Current state: **34 / 34
Verified** from 51 tests, 93 traces (41 satisfy, 51 verify, 1 link), 65
requirement-to-test verification pairs; every satisfy and every verify trace is
executed by at least one test.

## Why mantra

|                        | mantra 0.9.0 (this branch)                                   | tracey 1.4.0 (PR #1)                                        |
|------------------------|--------------------------------------------------------------|-------------------------------------------------------------|
| Requirements           | Markdown `` ### `id`: Title `` + body, hierarchy via `Parents`| Markdown `r[id]` marker + paragraph                          |
| Code → requirement     | `#[mantra_macros::req_satisfied("id")]` on the item          | `// r[impl id]` comment                                     |
| Test → requirement     | `#[mantra_macros::req_test("id")]` on the test               | `// r[verify id]` comment                                   |
| What "verified" means  | every tracing test **passed and executed the tracing code**  | a `verify` tag exists                                       |
| Inputs beyond source   | a test run with per-test-case line coverage                  | none                                                        |
| Also models            | manual requirements + review records, optional/deprecated/excluded, parent roll-up | spec versions (`+N`), stale references, LSP/MCP |
| Output                 | HTML/JSON report, SQLite DB                                   | CLI queries, web dashboard                                  |
| Crate dependency       | `mantra-macros` (compile-time no-op macros)                   | none                                                        |

mantra's model is the one we want: a requirement's state is computed from the
test run, so a `req_test` on a test that no longer runs, fails, or stops
exercising the implementation changes the requirement's state, and a parent is
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
  above the item, after its doc comment and above other attributes. **Trace
  functions and methods, not types or traits**: a type or trait declaration
  has no executable lines, so coverage can never confirm it (probe case
  `p-11`). Several IDs per attribute are fine, and several items may trace to
  one ID.
- **Verifying tests**: `#[mantra_macros::req_test("id", ...)]` above `#[test]`.
  Every test that claims an ID must execute at least one item tracing to that
  ID, otherwise the requirement is Unverified (probe cases `p-6`, `p-12`).
  Many tests may verify one ID and one test may verify many.
- **Manual requirements** (`nf-1` benchmark, `nf-4` third-party interop): no
  test can discharge them, so they are `Manual: true` and verified by a review
  record in `reviews/` that cites the evidence (README § Validation). mantra
  flags the review as likely obsolete if the requirement text changes after
  the review date.
- **Equations**: every function in SIGNAL_MODEL.md §9 carries
  `req_satisfied("eq-eN")`, so the §9 table is checked, not maintained by hand.

## Running it

```bash
scripts/trace.sh                 # self-check → per-test coverage → mantra collect → mantra report
SKIP_SELFCHECK=1 scripts/trace.sh
SKIP_TESTS=1 scripts/trace.sh    # reuse target/mantra/test-run.json, re-collect only
open mantra-report/index.html    # per-requirement traces, coverage, test cases, reviews
```

Tools: `cargo`, rustc's `llvm-tools` (`llvm-cov`, `llvm-profdata`; or set
`LLVM_COV`/`LLVM_PROFDATA`) and `mantra` 0.9 (`cargo install --locked mantra`;
or set `MANTRA`). No nextest or cargo-llvm-cov. Config:
[`mantra.json5`](../mantra.json5) (requirements, annotation paths, test-run
import, reviews).

[`tools/mantra-testrun`](../tools/mantra-testrun) is a ~500-line standalone
crate: `cargo test --no-run` with `-Cinstrument-coverage` into its own target
dir, `--list` per test binary, each test run alone with its own
`LLVM_PROFILE_FILE`, `llvm-profdata merge` + `llvm-cov export -format=lcov`
per test restricted to the package directory, and a mantra `TestRunSchema`
JSON with one `test_cases[].covered_files` per test holding executed lines
only (hits > 0) plus the SHA-256 of each covered file, so coverage from a
stale tree cannot match current traces. mantra's own schema is used instead
of JUnit + Cobertura because the well-known formats are imported per run, not
per test case.

Each run starts from a fresh `mantra.db`: the state is computed from the
current tree and the current test run, with no history to lean on.

## How the first end-to-end run played out

`scripts/trace.sh` on the first version of this branch: 50 tests pass under
nextest with coverage, mantra collects 34 requirements, 85 traces and the test
run, and reports 34 of 34 Verified. Getting there took four fixes, each a
property of the toolchain worth knowing:

1. **Dotted IDs are invisible to traces.** The spec initially used
   `fn-4.cw`-style children. mantra parsed them from the Markdown, but no
   `req_satisfied("fn-4.cw")` was ever collected, and `fid.1` failed outright
   with "Parent of requirement 'fid.1' was not collected". Children are now
   `fn-4-cw` with an explicit `Parents` property.
2. **Integration tests were not covered.** Every requirement verified only from
   `tests/*.rs` came back Unverified, while the same pattern in `src/` unit
   tests verified. Cause: cargo-llvm-cov excludes `tests/`, `benches/` and
   `examples/` from its report by default. (Moot now: `mantra-testrun`
   exports coverage for the whole package directory.)
3. **Parents inherit Unverified.** `fn-4` had its own passing test but stayed
   Unverified until `fn-4-bpsk` did, which is the roll-up we want and a reason to
   keep the hierarchy explicit.
4. **Small 0.9.0 rough edges.** `mantra report` refuses to overwrite an
   existing output directory (the script removes it first); `--formats markdown`
   fails with a missing template; `mantra-rust-macros` on docs.rs is the dead
   0.7 line and the 0.7.2 crate is yanked, so the live crates are `mantra` and
   `mantra-macros` 0.9.

## Second pass: making Verified mean something

That first 34/34 was not trustworthy, and a controlled check in a scratch crate
showed why: a requirement whose only `req_satisfied` function had **zero hits**
(never called by anything) still reported Verified, as did one whose
implementation was executed only by an untraced test. Two causes, both in how
coverage reaches mantra:

- **Run-level coverage cannot attribute lines to tests.** JUnit + Cobertura
  give mantra one coverage set for the whole run, so "the test that verifies X
  also executed code satisfying X" degenerates to "something in the run did".
- **"Covered" included zero-hit lines.** mantra's `TracesCoveredByTest*`
  views (0.9.0) keep every coverage record whose hit count is not null,
  including `hits="0"`. A never-called function whose lines appear in the
  report therefore counts as covered. Supplying only executed lines (hits > 0)
  per test case makes "covered by T" mean "executed by T".

`tools/mantra-testrun` fixes both, and its self-check pins the resulting
semantics: `tools/mantra-testrun/probe` is a fixture crate with 14
requirements whose expected states are written next to them in
[`probe/reqs/req.md`](../tools/mantra-testrun/probe/reqs/req.md), and
`tests/selfcheck.rs` runs the full chain on it and asserts every state.
`scripts/trace.sh` runs that first. The cases:

| probe | situation                                                              | mantra state |
|-------|------------------------------------------------------------------------|--------------|
| p-1   | traced test executes the traced function                               | Verified     |
| p-2   | traced test passes but never executes the traced function              | Unverified   |
| p-3   | traced test fails                                                      | Failed       |
| p-4   | no traces                                                              | Unverified   |
| p-5   | verified only from `tests/*.rs`                                        | Verified     |
| p-6   | implementation executed only by an untraced test                       | Unverified   |
| p-7/8 | one function satisfies two IDs; two tests verify one of them           | Verified     |
| p-9   | only verifying test is `#[ignore]`                                     | Unverified   |
| p-10  | `Manual: true` + review record                                         | Verified     |
| p-11  | satisfy trace on a trait declaration only                              | Unverified   |
| p-12  | two traced tests, one does not execute the implementation              | Unverified   |
| p-13  | parent with an Unverified child                                        | Unverified   |
| p-14  | parent with all children Verified, no traces of its own                | Verified     |

Rerunning lipika under these rules dropped it to **27 / 34**. Querying
`mantra.db` for "verify traces whose test executed no satisfy trace of the
same requirement" named 11 tests, and "satisfy traces executed by no test"
named one function. Every one was a real hole in the trace set or the tests,
not a tooling artefact:

| requirement | what the traced tests actually exercised                                      | resolution |
|-------------|-------------------------------------------------------------------------------|------------|
| eq-e7, fid-5 | the drift and Allan-deviation tests build `NodeClock::new` (materialised δ(t)); the flicker test calls `ClockStream::delta_t_at_master_index`; neither carried a trace | `req_satisfied("fid-5", "eq-e7")` on both (they implement the clock model; the bit-for-bit test ties them to `delta_t_at`) |
| nf-2        | the bit-for-bit stream test calls `delta_t_at_master_index`                   | `nf-2` added to it (the O(1)-state path) |
| fn-7        | `clock_stream_chunk_query_is_invariant` exercises `ClockStream::delta_t_at`, which traced fid-5/fn-2/eq-e7/nf-2 but not fn-7 | `fn-7` added: the seed-derived forward-only generator is how determinism is achieved at the clock level |
| fn-4-ca-code | the two dummy-nav tests call `NavData::bit_pm` directly                      | `req_satisfied("fn-4-ca-code")` on `bit_pm` (it is the `dummy` overlay) |
| fn-4-lnav   | `tlm_preamble_present` and `gps_time_roundtrip` exercise `tlm_word` and `GpsTime::from_unix`, which were untraced | `req_satisfied("fn-4-lnav")` on `GpsTime::from_unix`, `encode_word`, `tlm_word`, `how_word`, `build_subframe` (TOW, parity, TLM/HOW, framing) |
| fn-4        | parent of the two above                                                       | resolved by its children |
| fid-4       | `thermal_noise_power_w` (kTB) was executed by **no test**: every test scene sets `noise_power_w` explicitly | unit test `thermal_noise_floor_is_k_t_b` added |
| fn-4, fid-1 | `req_satisfied` on the `SignalModel` trait declaration can never be covered  | removed; the per-model `evaluate_baseband_batch` impls and `solve_emission_times_for` carry those IDs |

With those changes: 51 tests, 34 / 34 Verified, and no trace of either kind
left unexecuted. The added traces are listed so a reviewer can judge each one;
in every case the test already exercised an implementation of the requirement,
and the fix was to make the trace set complete rather than to loosen a test.

## What Verified does and does not prove

Proves, for every non-manual requirement:

- at least one test carries `req_test("id")`, ran in this run, and passed;
- every such test executed at least one line of an item carrying
  `req_satisfied("id")`, in that test's own process;
- the covered files are byte-identical (SHA-256) to the files the traces were
  collected from;
- all non-optional children are Verified.

Does not prove:

- **that the test asserts anything about the requirement.** Coverage is
  line-level: a `req_test` that calls the function and asserts something
  unrelated still verifies. Whether a test is a *good* test for its
  requirement remains a review question; the matrix tells the reviewer which
  tests to read.
- **that the whole implementation ran.** "At least one line" of the traced
  item, and the item span is the whole function.
- **anything about skipped tests.** An `#[ignore]`d verifier yields
  Unverified, not Skipped, because a test that did not run has no coverage.
- **independence of failures.** mantra treats lines executed by a failing
  test as failed coverage, so one failing test can turn other requirements
  whose code it executed Failed as well. Conservative, and worth knowing when
  reading a red report.
- **manual requirements** beyond "a dated review record exists" (`nf-1`,
  `nf-4`).
- **the CLI.** `bin/lipika.rs` carries no traces: nothing exercises it in the
  test run.

## Adding a requirement

1. Add `### `group-n`: Title` and its MUST-statement to
   [REQUIREMENTS.md](REQUIREMENTS.md); give children a `Parents` bullet.
2. `#[mantra_macros::req_satisfied("group-n")]` on the functions that satisfy
   it (not on types or traits).
3. `#[mantra_macros::req_test("group-n")]` on the test that proves it (write
   the test if there isn't one); the test must call into step 2's code.
4. `scripts/trace.sh`. If the requirement comes back Unverified, the report's
   requirement page lists its traces and test cases; the usual cause is a
   traced test that does not execute any traced implementation.
