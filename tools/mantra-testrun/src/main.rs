//! mantra-testrun: run every Rust test of a cargo package in its own process
//! with its own LLVM coverage profile, and write a mantra `TestRunSchema`
//! (JSON) in which each test case carries exactly the lines *it* executed.
//!
//! Why this exists
//! ---------------
//! mantra 0.9 computes a requirement's state from the join of traces
//! (`req_satisfied` / `req_test`) and coverage. Fed one coverage file for the
//! whole run it cannot tell which test executed which line, and it treats a
//! satisfying trace as "covered" as soon as coverage data exists for its
//! lines, even with zero hits. So a `req_test` that never calls the code it
//! claims to verify still yields Verified.
//!
//! With per-test-case coverage containing only executed lines, "covered by
//! test case T" means "executed by T", and mantra's per-test-case rules give
//! the property we want: a requirement is Verified only if every test that
//! claims to verify it passed and executed code that claims to satisfy it.
//! `tests/selfcheck.rs` pins that behaviour against the `probe/` fixture.
//!
//! Pipeline (plain cargo + the LLVM tools shipped with rustc; no nextest or
//! cargo-llvm-cov needed):
//!
//! 1. `cargo test --no-run` with `-Cinstrument-coverage` into its own target dir.
//! 2. `<test binary> --list` per binary (plus `--ignored` for skipped tests).
//! 3. Run each test alone: `<binary> <name> --exact` with its own
//!    `LLVM_PROFILE_FILE`; record state, duration and captured output.
//! 4. `llvm-profdata merge` + `llvm-cov export -format=lcov` per test,
//!    restricted to the package directory; keep lines with hits > 0.
//! 5. Write `target/mantra/test-run.json` for `mantra collect`
//!    (`test_runs: [{ path: "target/mantra", source: "schema" }]`).

use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Instant, SystemTime};

use serde::Serialize;
use sha2::{Digest, Sha256};
use time::OffsetDateTime;
use time::format_description::well_known::Iso8601;

/// mantra-schema version this output follows.
const SCHEMA_VERSION: &str = "0.9.0";
/// Captured stdout/stderr kept per test case (tail), in bytes.
const LOG_LIMIT: usize = 8 * 1024;

const USAGE: &str = "\
mantra-testrun: per-test-case LLVM coverage -> mantra TestRunSchema JSON

USAGE:
    mantra-testrun [OPTIONS] [-- <cargo test args>...]

OPTIONS:
    --manifest-dir <DIR>   package directory (default: current directory)
    --out <FILE>           output file (default: <DIR>/target/mantra/test-run.json)
    --name <NAME>          test run name (default: cargo-test)
    --target-dir <DIR>     instrumented build dir (default: <DIR>/target/mantra-cov)
    --jobs <N>             parallel test processes (default: available cores)
    -h, --help             print this help

ENVIRONMENT:
    LLVM_COV, LLVM_PROFDATA   paths to the LLVM tools; otherwise the rustup
                              `llvm-tools` component or PATH is used
    RUSTFLAGS                 extended with -Cinstrument-coverage

EXIT STATUS: 0 if every executed test passed, 1 if any failed or did not
run as exactly one test, 2 on usage or tool errors. The JSON file is written
in all but the last case.";

struct Options {
    manifest_dir: PathBuf,
    out: PathBuf,
    name: String,
    target_dir: PathBuf,
    jobs: usize,
    cargo_args: Vec<String>,
}

fn parse_args() -> Result<Options, String> {
    let mut manifest_dir: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut name = "cargo-test".to_string();
    let mut target_dir: Option<PathBuf> = None;
    let mut jobs: Option<usize> = None;
    let mut cargo_args = Vec::new();

    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = |flag: &str| {
            args.next()
                .ok_or_else(|| format!("{flag} needs a value\n\n{USAGE}"))
        };
        match arg.as_str() {
            "--manifest-dir" => manifest_dir = Some(PathBuf::from(value("--manifest-dir")?)),
            "--out" => out = Some(PathBuf::from(value("--out")?)),
            "--name" => name = value("--name")?,
            "--target-dir" => target_dir = Some(PathBuf::from(value("--target-dir")?)),
            "--jobs" => {
                jobs = Some(
                    value("--jobs")?
                        .parse()
                        .map_err(|e| format!("--jobs: {e}"))?,
                )
            }
            "-h" | "--help" => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            "--" => {
                cargo_args.extend(args.by_ref());
                break;
            }
            other => return Err(format!("unknown argument '{other}'\n\n{USAGE}")),
        }
    }

    let manifest_dir = manifest_dir.unwrap_or_else(|| PathBuf::from("."));
    let manifest_dir = manifest_dir
        .canonicalize()
        .map_err(|e| format!("--manifest-dir {}: {e}", manifest_dir.display()))?;
    if !manifest_dir.join("Cargo.toml").is_file() {
        return Err(format!("no Cargo.toml in {}", manifest_dir.display()));
    }
    let out = out.unwrap_or_else(|| manifest_dir.join("target/mantra/test-run.json"));
    let target_dir = target_dir.unwrap_or_else(|| manifest_dir.join("target/mantra-cov"));
    let jobs = jobs
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()))
        .max(1);

    Ok(Options {
        manifest_dir,
        out,
        name,
        target_dir,
        jobs,
        cargo_args,
    })
}

// ---------------------------------------------------------------------------
// Output schema (mirrors mantra_schema::test_runs, which uses deny_unknown_fields)
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct TestRunSchema {
    schema_version: &'static str,
    test_runs: Vec<TestRun>,
}

#[derive(Serialize)]
struct TestRun {
    name: String,
    utc_date: String,
    description: String,
    origin: serde_json::Value,
    nr_of_test_cases: u32,
    properties: serde_json::Map<String, serde_json::Value>,
    duration_sec: f64,
    test_cases: Vec<TestCase>,
}

#[derive(Serialize)]
struct TestCase {
    name: String,
    state: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    utc_date: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    duration_sec: Option<f64>,
    properties: serde_json::Map<String, serde_json::Value>,
    logs: Vec<LogOutput>,
    covered_files: Vec<CoveredFile>,
}

#[derive(Serialize)]
struct LogOutput {
    source: &'static str,
    content: String,
}

#[derive(Serialize)]
struct CoveredFile {
    filepath: String,
    file_hash: Option<String>,
    lines: Vec<CoveredLine>,
}

#[derive(Serialize)]
struct CoveredLine {
    nr: i64,
    hits: i64,
}

// ---------------------------------------------------------------------------
// Tool discovery
// ---------------------------------------------------------------------------

struct LlvmTools {
    cov: PathBuf,
    profdata: PathBuf,
}

fn cargo() -> String {
    env::var("CARGO").unwrap_or_else(|_| "cargo".into())
}

fn rustc_host() -> Option<String> {
    let out = Command::new("rustc").arg("-vV").output().ok()?;
    String::from_utf8_lossy(&out.stdout)
        .lines()
        .find_map(|l| l.strip_prefix("host: ").map(str::to_string))
}

fn find_llvm_tool(env_key: &str, file: &str) -> Result<PathBuf, String> {
    let candidate = match env::var_os(env_key) {
        Some(p) => PathBuf::from(p),
        None => {
            let from_sysroot = Command::new("rustc")
                .args(["--print", "sysroot"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .and_then(|o| {
                    let sysroot = String::from_utf8_lossy(&o.stdout).trim().to_string();
                    let host = rustc_host()?;
                    let p = Path::new(&sysroot)
                        .join("lib/rustlib")
                        .join(host)
                        .join("bin")
                        .join(file);
                    p.is_file().then_some(p)
                });
            from_sysroot.unwrap_or_else(|| PathBuf::from(file))
        }
    };
    let probe = Command::new(&candidate)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
    match probe {
        Ok(s) if s.success() => Ok(candidate),
        _ => Err(format!(
            "cannot run {file} ({}). Install the rustup `llvm-tools` component or set {env_key}.",
            candidate.display()
        )),
    }
}

// ---------------------------------------------------------------------------
// Build + discovery of test binaries
// ---------------------------------------------------------------------------

struct TestBinary {
    /// nextest-style id: `pkg`, `pkg::integration_test`, `pkg::bin/name`.
    id: String,
    path: PathBuf,
}

fn package_names(opts: &Options) -> Result<BTreeMap<String, String>, String> {
    let out = Command::new(cargo())
        .current_dir(&opts.manifest_dir)
        .args(["metadata", "--no-deps", "--format-version", "1"])
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| format!("cargo metadata: {e}"))?;
    if !out.status.success() {
        return Err("cargo metadata failed".into());
    }
    let v: serde_json::Value =
        serde_json::from_slice(&out.stdout).map_err(|e| format!("cargo metadata JSON: {e}"))?;
    let mut map = BTreeMap::new();
    for pkg in v["packages"].as_array().into_iter().flatten() {
        if let (Some(id), Some(name)) = (pkg["id"].as_str(), pkg["name"].as_str()) {
            map.insert(id.to_string(), name.to_string());
        }
    }
    Ok(map)
}

fn build_test_binaries(opts: &Options) -> Result<Vec<TestBinary>, String> {
    let packages = package_names(opts)?;

    let mut rustflags = env::var("RUSTFLAGS").unwrap_or_default();
    if !rustflags.contains("instrument-coverage") {
        rustflags.push_str(" -Cinstrument-coverage");
    }
    eprintln!("mantra-testrun: building instrumented tests in {}", opts.target_dir.display());
    let out = Command::new(cargo())
        .current_dir(&opts.manifest_dir)
        .args(["test", "--no-run", "--message-format=json-render-diagnostics"])
        .arg("--target-dir")
        .arg(&opts.target_dir)
        .args(&opts.cargo_args)
        .env("RUSTFLAGS", rustflags.trim())
        .env("CARGO_INCREMENTAL", "0")
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| format!("cargo test --no-run: {e}"))?;
    if !out.status.success() {
        return Err("cargo test --no-run failed".into());
    }

    let mut bins = Vec::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let Ok(msg) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        if msg["reason"] != "compiler-artifact" || msg["profile"]["test"] != true {
            continue;
        }
        let Some(exe) = msg["executable"].as_str() else {
            continue;
        };
        let target = &msg["target"];
        let tname = target["name"].as_str().unwrap_or("?");
        let kinds: Vec<&str> = target["kind"]
            .as_array()
            .map(|k| k.iter().filter_map(|x| x.as_str()).collect())
            .unwrap_or_default();
        let pkg_id = msg["package_id"].as_str().unwrap_or("");
        let pkg = packages.get(pkg_id).cloned().unwrap_or_else(|| tname.to_string());
        let id = if kinds.iter().any(|k| k.ends_with("lib")) {
            pkg
        } else if kinds.contains(&"bin") {
            format!("{pkg}::bin/{tname}")
        } else {
            format!("{pkg}::{tname}")
        };
        bins.push(TestBinary {
            id,
            path: PathBuf::from(exe),
        });
    }
    bins.sort_by(|a, b| a.id.cmp(&b.id));
    if bins.is_empty() {
        return Err("no test binaries were built".into());
    }
    Ok(bins)
}

fn list_tests(bin: &Path, ignored: bool, prof_dir: &Path) -> Result<Vec<String>, String> {
    let mut cmd = Command::new(bin);
    cmd.args(["--list", "--format", "terse"])
        // The binary is instrumented, so even `--list` writes a profile at exit;
        // keep it out of the working directory.
        .env("LLVM_PROFILE_FILE", prof_dir.join("list-%p.profraw"));
    if ignored {
        cmd.arg("--ignored");
    }
    let out = cmd
        .output()
        .map_err(|e| format!("{}: --list: {e}", bin.display()))?;
    if !out.status.success() {
        return Err(format!(
            "{} --list failed:\n{}",
            bin.display(),
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    Ok(String::from_utf8_lossy(&out.stdout)
        .lines()
        .filter_map(|l| l.strip_suffix(": test").map(str::to_string))
        .collect())
}

// ---------------------------------------------------------------------------
// Running one test with its own profile
// ---------------------------------------------------------------------------

struct Case {
    bin: usize,
    name: String,
    ignored: bool,
}

fn iso8601(t: SystemTime) -> String {
    OffsetDateTime::from(t)
        .format(&Iso8601::DEFAULT)
        .unwrap_or_else(|_| "1970-01-01T00:00:00Z".into())
}

fn tail(bytes: &[u8]) -> String {
    let s = String::from_utf8_lossy(bytes);
    if s.len() <= LOG_LIMIT {
        s.into_owned()
    } else {
        let mut cut = s.len() - LOG_LIMIT;
        while !s.is_char_boundary(cut) {
            cut += 1;
        }
        format!("[... {} bytes truncated ...]\n{}", cut, &s[cut..])
    }
}

fn run_case(
    opts: &Options,
    tools: &LlvmTools,
    bins: &[TestBinary],
    case: &Case,
    idx: usize,
    prof_dir: &Path,
) -> TestCase {
    let bin = &bins[case.bin];
    let full_name = format!("{}::{}", bin.id, case.name);
    let mut properties = serde_json::Map::new();
    properties.insert("binary".into(), bin.id.clone().into());
    properties.insert("test".into(), case.name.clone().into());

    if case.ignored {
        properties.insert("reason".into(), "#[ignore]".into());
        return TestCase {
            name: full_name,
            state: "skipped",
            utc_date: None,
            duration_sec: None,
            properties,
            logs: vec![],
            covered_files: vec![],
        };
    }

    let raw = prof_dir.join(format!("{idx}.profraw"));
    let profdata = prof_dir.join(format!("{idx}.profdata"));
    let started = SystemTime::now();
    let clock = Instant::now();
    let output = Command::new(&bin.path)
        .arg(&case.name)
        .args(["--exact", "--test-threads=1"])
        .current_dir(&opts.manifest_dir)
        .env("LLVM_PROFILE_FILE", &raw)
        .env("CARGO_MANIFEST_DIR", &opts.manifest_dir)
        .output();
    let duration = clock.elapsed().as_secs_f64();

    let (state, logs, exit_code) = match output {
        Ok(out) => {
            let stdout = tail(&out.stdout);
            let stderr = tail(&out.stderr);
            // libtest prints "running N test(s)"; anything but exactly one means
            // the filter did not select this test (renamed, ignored, ...).
            let ran_one = stdout.lines().any(|l| l.trim() == "running 1 test");
            let state = if !ran_one {
                "unknown"
            } else if out.status.success() {
                "passed"
            } else {
                "failed"
            };
            let mut logs = vec![LogOutput {
                source: "Stdout",
                content: stdout,
            }];
            if !stderr.is_empty() {
                logs.push(LogOutput {
                    source: "Stderr",
                    content: stderr,
                });
            }
            (state, logs, out.status.code())
        }
        Err(e) => (
            "unknown",
            vec![LogOutput {
                source: "Stderr",
                content: format!("failed to spawn {}: {e}", bin.path.display()),
            }],
            None,
        ),
    };
    if let Some(code) = exit_code {
        properties.insert("exit_code".into(), code.into());
    }

    let covered_files = if raw.is_file() {
        match export_coverage(tools, &bin.path, &raw, &profdata, &opts.manifest_dir) {
            Ok(files) => files,
            Err(e) => {
                eprintln!("mantra-testrun: {full_name}: coverage export failed: {e}");
                properties.insert("coverage_error".into(), e.into());
                vec![]
            }
        }
    } else {
        properties.insert(
            "coverage_error".into(),
            "no .profraw written (process aborted?)".into(),
        );
        vec![]
    };

    TestCase {
        name: full_name,
        state,
        utc_date: Some(iso8601(started)),
        duration_sec: Some(duration),
        properties,
        logs,
        covered_files,
    }
}

fn export_coverage(
    tools: &LlvmTools,
    bin: &Path,
    raw: &Path,
    profdata: &Path,
    root: &Path,
) -> Result<Vec<CoveredFile>, String> {
    let merge = Command::new(&tools.profdata)
        .arg("merge")
        .arg("-sparse")
        .arg(raw)
        .arg("-o")
        .arg(profdata)
        .output()
        .map_err(|e| format!("llvm-profdata: {e}"))?;
    if !merge.status.success() {
        return Err(format!(
            "llvm-profdata merge: {}",
            String::from_utf8_lossy(&merge.stderr).trim()
        ));
    }
    let export = Command::new(&tools.cov)
        .arg("export")
        .arg(bin)
        .arg(format!("-instr-profile={}", profdata.display()))
        .arg("-format=lcov")
        .arg("--")
        .arg(root)
        .output()
        .map_err(|e| format!("llvm-cov: {e}"))?;
    if !export.status.success() {
        return Err(format!(
            "llvm-cov export: {}",
            String::from_utf8_lossy(&export.stderr).trim()
        ));
    }
    Ok(parse_lcov(&String::from_utf8_lossy(&export.stdout), root))
}

/// lcov → covered files relative to `root`, executed lines only (hits > 0).
fn parse_lcov(lcov: &str, root: &Path) -> Vec<CoveredFile> {
    let mut files = Vec::new();
    let mut current: Option<(String, Vec<CoveredLine>)> = None;
    for line in lcov.lines() {
        if let Some(sf) = line.strip_prefix("SF:") {
            let abs = Path::new(sf);
            let rel = abs
                .strip_prefix(root)
                .ok()
                .map(|p| p.to_string_lossy().replace('\\', "/"))
                .filter(|p| !p.starts_with("target/"));
            current = rel.map(|p| (p, Vec::new()));
        } else if let Some(da) = line.strip_prefix("DA:") {
            if let Some((_, lines)) = current.as_mut() {
                let mut it = da.split(',');
                if let (Some(nr), Some(hits)) = (it.next(), it.next())
                    && let (Ok(nr), Ok(hits)) = (nr.parse::<i64>(), hits.parse::<i64>())
                    && hits > 0
                {
                    lines.push(CoveredLine { nr, hits });
                }
            }
        } else if line == "end_of_record" {
            if let Some((filepath, mut lines)) = current.take()
                && !lines.is_empty()
            {
                lines.sort_by_key(|l| l.nr);
                lines.dedup_by_key(|l| l.nr);
                let file_hash = fs::read(root.join(&filepath))
                    .ok()
                    .map(|bytes| format!("{:x}", Sha256::digest(&bytes)));
                files.push(CoveredFile {
                    filepath,
                    file_hash,
                    lines,
                });
            }
        }
    }
    files.sort_by(|a, b| a.filepath.cmp(&b.filepath));
    files
}

// ---------------------------------------------------------------------------
// Origin (git) + main
// ---------------------------------------------------------------------------

fn git_origin(dir: &Path) -> serde_json::Value {
    let git = |args: &[&str]| {
        Command::new("git")
            .arg("-C")
            .arg(dir)
            .args(args)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
    };
    match git(&["rev-parse", "HEAD"]) {
        Some(commit) => {
            let dirty = git(&["status", "--porcelain", "--untracked-files=no"])
                .map(|s| !s.is_empty());
            let branch = git(&["rev-parse", "--abbrev-ref", "HEAD"]);
            serde_json::json!({ "vcs": "git", "commit": commit, "branch": branch, "dirty": dirty })
        }
        None => serde_json::json!({ "vcs": null }),
    }
}

fn run() -> Result<bool, String> {
    let opts = parse_args()?;
    let tools = LlvmTools {
        cov: find_llvm_tool("LLVM_COV", "llvm-cov")?,
        profdata: find_llvm_tool("LLVM_PROFDATA", "llvm-profdata")?,
    };

    let run_started = SystemTime::now();
    let run_clock = Instant::now();
    let bins = build_test_binaries(&opts)?;

    let prof_dir = opts.target_dir.join("mantra-profiles");
    let _ = fs::remove_dir_all(&prof_dir);
    fs::create_dir_all(&prof_dir).map_err(|e| format!("{}: {e}", prof_dir.display()))?;

    let mut cases = Vec::new();
    for (i, bin) in bins.iter().enumerate() {
        let ignored = list_tests(&bin.path, true, &prof_dir)?;
        let mut all = list_tests(&bin.path, false, &prof_dir)?;
        all.extend(ignored.iter().cloned());
        all.sort();
        all.dedup();
        for name in all {
            cases.push(Case {
                bin: i,
                ignored: ignored.contains(&name),
                name,
            });
        }
    }
    eprintln!(
        "mantra-testrun: {} tests in {} binaries, {} at a time",
        cases.len(),
        bins.len(),
        opts.jobs
    );

    let results: Mutex<Vec<Option<TestCase>>> =
        Mutex::new((0..cases.len()).map(|_| None).collect());
    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..opts.jobs.min(cases.len().max(1)) {
            scope.spawn(|| {
                loop {
                    let idx = next.fetch_add(1, Ordering::SeqCst);
                    if idx >= cases.len() {
                        break;
                    }
                    let tc = run_case(&opts, &tools, &bins, &cases[idx], idx, &prof_dir);
                    eprintln!("  {:<8} {}", tc.state, tc.name);
                    results.lock().unwrap()[idx] = Some(tc);
                }
            });
        }
    });
    let test_cases: Vec<TestCase> = results
        .into_inner()
        .unwrap()
        .into_iter()
        .map(|c| c.expect("every case produces a result"))
        .collect();

    let count = |s: &str| test_cases.iter().filter(|c| c.state == s).count();
    let (passed, failed, skipped, unknown) =
        (count("passed"), count("failed"), count("skipped"), count("unknown"));
    let covered_lines: usize = test_cases
        .iter()
        .flat_map(|c| &c.covered_files)
        .map(|f| f.lines.len())
        .sum();

    let mut properties = serde_json::Map::new();
    properties.insert("runner".into(), format!("mantra-testrun {}", env!("CARGO_PKG_VERSION")).into());
    properties.insert("coverage".into(), "per test case, executed lines only (hits > 0), llvm-cov lcov export".into());
    properties.insert("cargo_args".into(), serde_json::json!(opts.cargo_args));
    properties.insert("binaries".into(), serde_json::json!(bins.iter().map(|b| &b.id).collect::<Vec<_>>()));
    if let Some(host) = rustc_host() {
        properties.insert("host".into(), host.into());
    }

    let schema = TestRunSchema {
        schema_version: SCHEMA_VERSION,
        test_runs: vec![TestRun {
            name: opts.name.clone(),
            utc_date: iso8601(run_started),
            description: format!(
                "{passed} passed, {failed} failed, {skipped} skipped, {unknown} unknown; \
                 each test case run in its own process with its own LLVM profile"
            ),
            origin: git_origin(&opts.manifest_dir),
            nr_of_test_cases: test_cases.len() as u32,
            properties,
            duration_sec: run_clock.elapsed().as_secs_f64(),
            test_cases,
        }],
    };

    if let Some(parent) = opts.out.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
    }
    let json = serde_json::to_string_pretty(&schema).map_err(|e| e.to_string())?;
    fs::write(&opts.out, json).map_err(|e| format!("{}: {e}", opts.out.display()))?;

    eprintln!(
        "mantra-testrun: {passed} passed, {failed} failed, {skipped} skipped, {unknown} unknown; \
         {covered_lines} executed line records -> {}",
        opts.out.display()
    );
    Ok(failed == 0 && unknown == 0)
}

fn main() -> ExitCode {
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("mantra-testrun: error: {e}");
            ExitCode::from(2)
        }
    }
}
