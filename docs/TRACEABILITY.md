# Requirements traceability

Lipika's requirements live in [REQUIREMENTS.md](REQUIREMENTS.md) and are tied
to the code and tests that satisfy them by comment annotations that a tool can
check. This document records why we chose that tool, how the annotations work,
and what the first pass over the codebase turned up.

## Tool evaluation

Two Rust-ecosystem requirement-tracing tools were evaluated hands-on against
this repository on 2026-10-27: [tracey](https://github.com/bearcove/tracey)
(bearcove / Amos Wenger) and [mantra](https://github.com/mhatzl/mantra)
(Manuel Hatzl). Both were installed and run against a copy of our requirement
set.

|                        | tracey                                                    | mantra                                                       |
|------------------------|-----------------------------------------------------------|--------------------------------------------------------------|
| Release used           | 1.4.0 (May 2026); main is 2.0.0-rc                        | 0.9.0 (Sep 2026); the suggested 0.7.2 is **yanked**           |
| Activity, licence      | ~250 stars, pushed Oct 2026, MIT/Apache-2.0               | ~50 stars, pushed Oct 2026, MIT                              |
| Requirements written as| Markdown `r[id]` marker + paragraph (also StrictDoc, Typst)| Markdown `` # `id`: Title `` headings, or JSON5/TOML schema  |
| Code annotation        | any comment: `// r[impl id]`, `// r[verify id]`           | proc-macros only: `#[mantra_macros::req_satisfied("id")]`, `#[req_test]`, `assert_req!` |
| Build dependency       | none                                                      | `mantra-macros` crate in `Cargo.toml` (compile-time only)    |
| "Tested" means         | a `verify` annotation exists                              | a `verifies` trace **and** line coverage of the `satisfies` code from a passed test (JUnit + Cobertura via nextest + grcov) |
| Outputs                | web dashboard, LSP, MCP server, CLI queries, JSON         | SQLite DB, HTML/JSON report, review files                     |
| CI gate                | `query validate --deny warnings` (main only, not in 1.4.0)| `collect` + `report`; states Verified/Unverified/Failed       |
| Spec drift             | `+N` versions, `tracey pre-commit`, `tracey bump`         | content hashes, reviews                                      |
| Install on NixOS       | prebuilt binary needs `patchelf`; source build needs Node + pnpm for the dashboard | `cargo install --locked mantra` (needs a C compiler) |

### tracey — adopted

What worked: zero intrusion into the crate (comments only, no dependency),
the spec stays an ordinary Markdown file, the CLI answered the questions we
care about (`status`, `uncovered`, `untested`, `rule <id>`), and the
dashboard renders the spec with inline coverage and the source tree with
annotations. Requirement IDs must contain at least one dot, which fits a
`group.n[.sub]` scheme; letters, digits, `-` and `_` are allowed in segments.
Files matched by `test_include` may only carry `verify` annotations, which
keeps the implementation/verification distinction honest.

What to know:

- **1.4.0 does not report broken references from the CLI.** A deliberate
  `// r[impl fid.99]` was silently ignored by `tracey query status` (exit 0);
  only the dashboard/LSP surface it. The `query validate` command that fails
  CI on broken, duplicate or stale references exists only on `main`
  (2.0.0-rc). Until that ships, [`scripts/trace-check.sh`](../scripts/trace-check.sh)
  is the gate: grep-only, checks broken references, duplicates, and that every
  requirement has an `impl` and a `verify`.
- Building from source needs Node.js and pnpm (the web dashboard is bundled
  into the binary); on NixOS the prebuilt release binary runs after
  `patchelf --set-interpreter` to the Nix glibc loader. Elsewhere the release
  installer script or `cargo install --locked --git https://github.com/bearcove/tracey --tag v1.4.0 tracey` is enough.
- It runs a per-workspace daemon (state under `~/.local/state/tracey`) that
  the CLI auto-starts and that exits after 10 min idle.
- Versioning (`r[id+2]`) is opt-in. Unversioned requirements are version 1; a
  bump makes every existing reference stale until reviewed. We have not
  enabled it yet; see Conventions.

### mantra — evaluated, not adopted (for now)

mantra is the more rigorous model: a requirement is only *Verified* when a
test that traces to it actually executed the lines that trace to it, so a
stale `verify` on a test that no longer touches the implementation is caught.
It also has explicit hierarchy (`Parents`), optional/manual/deprecated states
and review records. Trial findings:

- Dotted IDs imply hierarchy: `fid.1` failed to collect with "Parent of
  requirement 'fid.1' was not collected" because no `fid` requirement exists.
  Hyphenated IDs (`fid-1`) work, and `Parents` gives hierarchy explicitly.
- Traces are detected by proc-macro name only. Plain-comment annotations were
  not picked up; `#[mantra_macros::req_satisfied]` and `#[req_test]` were.
  Adopting it means adding `mantra-macros` to the crate (the macros expand to
  `const` items without the `coverage` feature, so there is no runtime cost,
  but every traced function carries an attribute).
- Verification needs the full coverage toolchain (`cargo-nextest` with JUnit
  output, `grcov` with Cobertura output, `-Cinstrument-coverage` builds).
  Without a test-run import every requirement reports *Unverified* even with
  traces present.
- 0.9.0 rough edges: `--formats markdown` fails with a missing template; the
  README on docs.rs for 0.7.2 describes a different (older) syntax.

Verdict: tracey gives most of the value at a fraction of the cost today. If we
later want coverage-proven verification (for example for an external audit),
mantra's Markdown format is close enough that the spec converts mechanically:
one heading per `r[id]` marker, IDs re-spelt with hyphens.

## Conventions

- **Spec**: [REQUIREMENTS.md](REQUIREMENTS.md). One `r[id]` marker per
  requirement at column 0, its text in the paragraph that follows (no blank
  line between). Use the `> r[id]` blockquote form if a body needs several
  paragraphs or a code block.
- **IDs**: `fid.N` fidelity principles, `fn.N` functional, `nf.N`
  non-functional, `eq.eN` SIGNAL_MODEL.md equation (EN). Sub-requirements
  append a segment (`fn.4.bpsk`, `fn.5.rinex-nav`). IDs are stable; renumbering
  is a breaking change for anything that cites them.
- **Implementation**: `// r[impl id]` on the line(s) directly above the item,
  after its doc comment and above any attributes. Several requirements stack as
  separate lines. Text after the `]` is free-form and ignored by the parser;
  use it to say *which aspect* the item implements.
- **Verification**: `// r[verify id]` directly above `#[test]` (or a Criterion
  bench function for `nf.1`). `tests/` and `benches/` are `test_include` paths
  and may only carry `verify`.
- **Equations**: every function in SIGNAL_MODEL.md §9 carries `r[impl eq.eN]`
  in addition to its prose `(E…)` back-reference, so the §9 table is checked
  rather than maintained by hand.
- **Changing a requirement's text**: not versioned yet. When a change alters
  what the code must do, review the `impl`/`verify` sites listed by
  `tracey query rule <id>` in the same commit. Enable `+N` versions and the
  `tracey pre-commit` hook once the set stabilises.

## Running it

```bash
# dependency-free gate (CI / pre-commit)
scripts/trace-check.sh

# tracey (daemon auto-starts)
tracey query status                 # coverage summary; --json for machines
tracey query uncovered              # requirements with no impl
tracey query untested               # requirements with impl but no test
tracey query rule fn.3              # text + every impl/verify site
tracey query unmapped               # code units with no requirement (reverse view)
tracey web --open                   # dashboard
tracey ai                           # registers the MCP server + skill for Claude/Codex
```

Config is [`.config/tracey/config.styx`](../.config/tracey/config.styx): one
spec (`docs/REQUIREMENTS.md`), one implementation (`src/**/*.rs`, with
`tests/**` and `benches/**` as test files).

## What the first pass exposed

Annotating the existing code against the README's requirement list was itself
a useful audit:

- `fn.4.cw` (the CW signal model) had **no test at all**; a unit test now
  checks the tone against `exp(j·2π·f·t)`.
- `fn.4` (signal models are stateless algebraic functions) and `fn.5` (every
  trajectory kind answers positions through one dispatch) were verified only
  indirectly by the integration tests; each now has a direct test
  (`evaluation_is_stateless_and_batch_independent`,
  `any_trajectory_dispatch_matches_the_wrapped_trajectory`).
- `nf.1` (throughput) has no implementation in the usual sense; the
  contribution kernel is marked `impl` and the Criterion bench is its
  `verify`.
- `nf.4` (loads in NumPy / GNU Radio / GNSS-SDR) is verified at the format
  level only (`raw_complex64_roundtrip_via_byte_check`); the end-to-end
  evidence remains the manual GNSS-SDR validation in the README.
- `fn.6` is verified for the raw IQ format but **not for the sidecar's
  content** (epoch, rate, LO, descriptor). A sidecar schema test is the obvious
  next addition.
- `eq.e10` (the full forward model) is verified by the streaming-vs-in-memory
  equivalence and the post-sum chain unit tests, not by a closed-form check.
- `tracey query unmapped` shows that about 75 % of code units (functions,
  types) carry no requirement. That is expected (helpers, config parsing,
  profiling), but `src/config/mod.rs` and `src/backend/` have none at all,
  and the scene-format contract in the README is not a requirement yet.

## Adding a requirement

1. Add an `r[group.n]` marker and its MUST-statement to
   [REQUIREMENTS.md](REQUIREMENTS.md) under the right heading.
2. Put `// r[impl group.n]` above the code that satisfies it.
3. Put `// r[verify group.n]` above the test that proves it (write the test
   if there isn't one — an `untested` requirement is a finding, not a state to
   leave behind).
4. Run `scripts/trace-check.sh` (and `tracey query rule group.n` to see it
   resolved).
