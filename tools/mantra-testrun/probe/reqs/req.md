# Probe requirements

Each requirement below has a known expected mantra state once
`mantra-testrun` + `mantra collect` have run (see `tests/selfcheck.rs`).
Together they pin the verification semantics the lipika pipeline relies on.

### `p-1`: Function executed by its traced test

- **Parents:** ["p-13", "p-14"]

Expected **Verified**: `one()` carries `req_satisfied("p-1")`, a passing test
carries `req_test("p-1")` and executes `one()`.

### `p-2`: Function never executed, traced test passes

- **Parents:** ["p-13"]

Expected **Unverified**: `never_executed()` carries `req_satisfied("p-2")`, a
passing test carries `req_test("p-2")` but never calls it. With run-level
coverage mantra 0.9 reports this Verified; with per-test-case coverage of
executed lines only, it must not.

### `p-3`: Traced test fails

Expected **Failed**: the only test tracing to `p-3` executes `three()` and fails.

### `p-4`: No traces

Expected **Unverified**: nothing references this requirement.

### `p-5`: Verified only from an integration test

Expected **Verified**: the only `req_test("p-5")` lives in `tests/integ.rs`.

### `p-6`: Implementation executed only by an untraced test

Expected **Unverified**: `six()` is executed by a test without `req_test`; the
test that carries `req_test("p-6")` does not execute `six()`.

### `p-7`: Many-to-many, first

- **Parents:** ["p-14"]

Expected **Verified**: `seven_eight()` satisfies `p-7` and `p-8`; two tests
verify `p-7`, one of them also verifies `p-8`.

### `p-8`: Many-to-many, second

Expected **Verified** (see `p-7`).

### `p-9`: Only verifying test is ignored

Expected **Unverified**: the only test tracing to `p-9` is `#[ignore]`d, so it
is recorded as Skipped and produces no coverage. Trace-based verification
cannot distinguish "skipped" from "never ran"; mantra reports Unverified.

### `p-10`: Manual requirement with a review record

- **Manual:** true

Expected **Verified**: no test can discharge this; `reviews/review.json5`
lists it.

### `p-11`: Trait-only satisfy trace

Expected **Unverified**: the only `req_satisfied("p-11")` is on a trait
declaration, which has no executable lines, so no coverage can ever confirm
it even though a passing test calls the implementing method.

### `p-12`: Two traced tests, one of which does not execute the implementation

Expected **Unverified**: one passing `req_test("p-12")` executes `twelve()`,
another passing `req_test("p-12")` does not. Every test that claims to verify
a requirement must execute code that claims to satisfy it.

### `p-13`: Parent of p-1 and p-2

Expected **Unverified**: a parent inherits Unverified from any non-optional
child (`p-2`).

### `p-14`: Parent of p-1 and p-7

Expected **Verified**: all children verified, no traces of its own needed.
