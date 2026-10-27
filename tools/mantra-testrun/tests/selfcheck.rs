//! End-to-end self-check: run mantra-testrun on the `probe/` fixture, feed the
//! result to `mantra collect` + `mantra report`, and assert that every probe
//! requirement resolves to the state documented in `probe/reqs/req.md`.
//!
//! This is what makes the lipika traceability chain trustable: the semantics
//! we rely on (a `req_test` only verifies a requirement if it executes code
//! that `req_satisfied`s it; failures propagate; parents roll up; trait-only
//! traces and skipped tests never verify) are checked here, not assumed.
//!
//! Requires `mantra` (0.9) on PATH or in `$MANTRA`, plus rustc's llvm-tools.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const EXPECTED: &[(&str, &str)] = &[
    ("p-1", "verified"),
    ("p-2", "unverified"),
    ("p-3", "failed"),
    ("p-4", "unverified"),
    ("p-5", "verified"),
    ("p-6", "unverified"),
    ("p-7", "verified"),
    ("p-8", "verified"),
    ("p-9", "unverified"),
    ("p-10", "verified"),
    ("p-11", "unverified"),
    ("p-12", "unverified"),
    ("p-13", "unverified"),
    ("p-14", "verified"),
];

fn run(cmd: &mut Command, what: &str) -> std::process::Output {
    let out = cmd.output().unwrap_or_else(|e| panic!("{what}: failed to start: {e}"));
    if !out.status.success() {
        eprintln!("--- {what} stdout ---\n{}", String::from_utf8_lossy(&out.stdout));
        eprintln!("--- {what} stderr ---\n{}", String::from_utf8_lossy(&out.stderr));
    }
    out
}

fn state_of(report: &Path, id: &str) -> String {
    let path = report.join(format!("products/probe@main/requirements/{id}.json5"));
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let v: serde_json::Value = json5::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    v["state"]
        .as_str()
        .unwrap_or_else(|| panic!("{}: no state field", path.display()))
        .to_string()
}

#[test]
fn probe_requirements_resolve_to_documented_states() {
    let probe = Path::new(env!("CARGO_MANIFEST_DIR")).join("probe");
    let mantra = std::env::var("MANTRA").unwrap_or_else(|_| "mantra".into());
    let test_run = probe.join("target/mantra/test-run.json");
    let report: PathBuf = probe.join("target/mantra-report");

    // Fresh state: no history from earlier runs may influence the result.
    let _ = fs::remove_file(probe.join("mantra.db"));
    let _ = fs::remove_dir_all(&report);
    let _ = fs::remove_file(&test_run);

    let out = run(
        Command::new(env!("CARGO_BIN_EXE_mantra-testrun"))
            .arg("--manifest-dir")
            .arg(&probe)
            .arg("--out")
            .arg(&test_run)
            .arg("--name")
            .arg("probe-selfcheck"),
        "mantra-testrun",
    );
    assert_eq!(
        out.status.code(),
        Some(1),
        "mantra-testrun must exit 1 because p-3's test fails on purpose"
    );
    assert!(test_run.is_file(), "test run JSON must be written even when tests fail");

    let json: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&test_run).unwrap()).expect("valid JSON");
    let cases = json["test_runs"][0]["test_cases"].as_array().expect("test_cases");
    let state = |suffix: &str| {
        cases
            .iter()
            .find(|c| c["name"].as_str().unwrap().ends_with(suffix))
            .unwrap_or_else(|| panic!("no test case ending in {suffix}"))
    };
    assert_eq!(state("p1_executes_one")["state"], "passed");
    assert_eq!(state("p3_fails_on_purpose")["state"], "failed");
    assert_eq!(state("p9_ignored")["state"], "skipped");
    assert_eq!(state("p5_verified_from_integration_test")["state"], "passed");
    // The traced test for p-2 must not have executed src/never.rs at all.
    let p2_files: Vec<&str> = state("p2_passes_without_executing_never_executed")["covered_files"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f["filepath"].as_str().unwrap())
        .collect();
    assert!(!p2_files.contains(&"src/never.rs"), "p-2's test must not cover never.rs: {p2_files:?}");
    // ...while p-1's test did execute src/lib.rs.
    assert!(
        state("p1_executes_one")["covered_files"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["filepath"] == "src/lib.rs" && !f["lines"].as_array().unwrap().is_empty())
    );

    let out = run(Command::new(&mantra).current_dir(&probe).arg("collect"), "mantra collect");
    assert!(out.status.success(), "mantra collect failed (is `mantra` 0.9 installed? set $MANTRA)");
    let out = run(
        Command::new(&mantra)
            .current_dir(&probe)
            .args(["report", "--formats", "json", "--output-dir"])
            .arg(&report),
        "mantra report",
    );
    assert!(out.status.success(), "mantra report failed");

    let mut mismatches = Vec::new();
    for (id, expected) in EXPECTED {
        let actual = state_of(&report, id);
        if actual != *expected {
            mismatches.push(format!("{id}: expected {expected}, mantra says {actual}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "probe states differ from probe/reqs/req.md:\n  {}",
        mismatches.join("\n  ")
    );
}
