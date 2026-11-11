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
//! 1. `cargo test --no-run` into its own target dir, with this binary as
//!    `RUSTC_WORKSPACE_WRAPPER` adding `-Cinstrument-coverage` (so the user's
//!    rustflags, from env or cargo config, are kept as they are).
//! 2. `<test binary> --list` per binary (plus `--ignored` for skipped tests).
//! 3. Run each test alone: `<binary> <name> --exact` with its own
//!    `LLVM_PROFILE_FILE` (one file per process, so children it spawns do not
//!    clobber it) and a wall-clock limit; record state, duration and output.
//! 4. `llvm-profdata merge` + `llvm-cov export -format=lcov` per test,
//!    restricted to the package directory; keep lines with hits > 0.
//! 5. Write `target/mantra/test-run.json` for `mantra collect`
//!    (`test_runs: [{ path: "target/mantra", source: "schema" }]`).
//!
//! Every `file_hash` is the SHA-256 of the source *before* the build; if any
//! source changed by the end of the run, no JSON is written (exit 2), so line
//! numbers from one version of a file never carry the hash of another.

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode, Stdio};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant, SystemTime};

use serde::Serialize;
use sha2::{Digest, Sha256};
use time::OffsetDateTime;
use time::format_description::well_known::Iso8601;

/// mantra-schema version this output follows.
const SCHEMA_VERSION: &str = "0.9.0";
/// Captured stdout/stderr kept per test case (tail), in bytes.
const LOG_LIMIT: usize = 8 * 1024;
/// Default per-test wall-clock limit, in seconds.
const DEFAULT_TIMEOUT_SECS: u64 = 600;
/// Added to every rustc call for a workspace crate.
const INSTRUMENT: &str = "-Cinstrument-coverage";
/// Set when this binary runs as cargo's RUSTC_WORKSPACE_WRAPPER.
const WRAP_ENV: &str = "MANTRA_TESTRUN_WRAP_RUSTC";
/// The user's own RUSTC_WORKSPACE_WRAPPER, if any, which we chain to.
const INNER_WRAP_ENV: &str = "MANTRA_TESTRUN_INNER_WRAPPER";

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
    --timeout <SECS>       per-test wall-clock limit; a test still running is
                           killed and recorded as failed (default: 600)
    -h, --help             print this help

ENVIRONMENT:
    LLVM_COV, LLVM_PROFDATA   paths to the LLVM tools; otherwise the rustup
                              `llvm-tools` component or PATH is used
    RUSTC_WORKSPACE_WRAPPER   kept: chained behind mantra-testrun, which
                              adds -Cinstrument-coverage to every workspace
                              crate's rustc call. Rustflags (env or cargo
                              config) are left exactly as cargo resolves them.

EXIT STATUS: 0 if every executed test passed with coverage, 1 if any failed,
timed out, did not run as exactly one test, or passed without coverage data,
2 on usage or tool errors (including an uninstrumented build or a source file
changed during the run). The JSON file is written in all but the last case.";

struct Options {
    manifest_dir: PathBuf,
    out: PathBuf,
    name: String,
    target_dir: PathBuf,
    jobs: usize,
    timeout: Duration,
    cargo_args: Vec<String>,
}

fn parse_args() -> Result<Options, String> {
    let mut manifest_dir: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut name = "cargo-test".to_string();
    let mut target_dir: Option<PathBuf> = None;
    let mut jobs: Option<usize> = None;
    let mut timeout = Duration::from_secs(DEFAULT_TIMEOUT_SECS);
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
            "--timeout" => {
                let secs: u64 = value("--timeout")?
                    .parse()
                    .map_err(|e| format!("--timeout: {e}"))?;
                timeout = Duration::from_secs(secs.max(1));
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
        timeout,
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

fn build_test_binaries(opts: &Options, prof_dir: &Path) -> Result<Vec<TestBinary>, String> {
    let packages = package_names(opts)?;

    eprintln!("mantra-testrun: building instrumented tests in {}", opts.target_dir.display());
    let mut cmd = Command::new(cargo());
    cmd.current_dir(&opts.manifest_dir);
    // Instrument through a workspace wrapper (this binary, see wrap_rustc)
    // rather than RUSTFLAGS: cargo takes rustflags from only one of
    // CARGO_ENCODED_RUSTFLAGS, RUSTFLAGS, target.*.rustflags, build.rustflags,
    // so setting any of them would drop the flags the user normally builds
    // with. Instrumented build scripts write profiles too; keep those here.
    let me = env::current_exe().map_err(|e| format!("current_exe: {e}"))?;
    cmd.env("RUSTC_WORKSPACE_WRAPPER", &me)
        .env(WRAP_ENV, "1")
        .env("LLVM_PROFILE_FILE", prof_dir.join("build-%p.profraw"));
    if let Some(inner) = env::var_os("RUSTC_WORKSPACE_WRAPPER").filter(|w| !w.is_empty()) {
        cmd.env(INNER_WRAP_ENV, inner);
    }
    let out = cmd
        .args(["test", "--no-run", "--message-format=json-render-diagnostics"])
        .arg("--target-dir")
        .arg(&opts.target_dir)
        .args(&opts.cargo_args)
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

/// `--list` the tests of `bin` (only `#[ignore]`d ones with `ignored`). The
/// binary is instrumented, so even `--list` writes a profile at exit: keep it
/// out of the working directory, and treat a missing one as a build that lost
/// `-Cinstrument-coverage` (every test would otherwise come out Unverified
/// with nothing pointing at the cause).
fn list_tests(bin: &Path, ignored: bool, prof_dir: &Path, tag: &str) -> Result<Vec<String>, String> {
    let prefix = format!("list-{tag}-");
    let mut cmd = Command::new(bin);
    cmd.args(["--list", "--format", "terse"])
        .env("LLVM_PROFILE_FILE", prof_dir.join(format!("{prefix}%p.profraw")));
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
    if profiles(prof_dir, &prefix).is_empty() {
        return Err(format!(
            "{} wrote no coverage profile, so it was built without {INSTRUMENT}. \
             Is it a workspace member (only those are instrumented), and does \
             RUSTC_WORKSPACE_WRAPPER pass extra arguments through to rustc?",
            bin.display()
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

/// Non-empty `*.profraw` files in `dir` whose name starts with `prefix`. The
/// runtime creates the file at startup and fills it at exit, so a killed
/// process leaves an empty one: that is "no profile", not "nothing executed".
fn profiles(dir: &Path, prefix: &str) -> Vec<PathBuf> {
    let mut found: Vec<PathBuf> = fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter(|e| e.metadata().is_ok_and(|m| m.len() > 0))
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with(prefix) && n.ends_with(".profraw"))
        })
        .collect();
    found.sort();
    found
}

/// Run `cmd` with stdout/stderr redirected to files (so a grandchild holding
/// a pipe open cannot block us) and kill it after `timeout`. Returns the exit
/// status, or `None` if it was killed.
fn run_with_timeout(
    cmd: &mut Command,
    stdout: &Path,
    stderr: &Path,
    timeout: Duration,
) -> std::io::Result<Option<std::process::ExitStatus>> {
    let mut child = cmd
        .stdin(Stdio::null())
        .stdout(fs::File::create(stdout)?)
        .stderr(fs::File::create(stderr)?)
        .spawn()?;
    let deadline = Instant::now() + timeout;
    let mut nap = Duration::from_millis(5);
    loop {
        if let Some(status) = child.try_wait()? {
            return Ok(Some(status));
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            return Ok(None);
        }
        std::thread::sleep(nap);
        nap = (nap * 2).min(Duration::from_millis(100));
    }
}

fn run_case(
    opts: &Options,
    tools: &LlvmTools,
    bins: &[TestBinary],
    case: &Case,
    idx: usize,
    prof_dir: &Path,
    hashes: &BTreeMap<String, String>,
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

    // `%p`: a process the test spawns inherits LLVM_PROFILE_FILE and must
    // write its own file instead of overwriting the test's.
    let prefix = format!("{idx}-");
    let profdata = prof_dir.join(format!("{idx}.profdata"));
    let stdout_path = prof_dir.join(format!("{idx}.stdout"));
    let stderr_path = prof_dir.join(format!("{idx}.stderr"));
    let started = SystemTime::now();
    let clock = Instant::now();
    let status = run_with_timeout(
        Command::new(&bin.path)
            .arg(&case.name)
            .args(["--exact", "--test-threads=1"])
            .current_dir(&opts.manifest_dir)
            .env("LLVM_PROFILE_FILE", prof_dir.join(format!("{prefix}%p.profraw")))
            .env("CARGO_MANIFEST_DIR", &opts.manifest_dir),
        &stdout_path,
        &stderr_path,
        opts.timeout,
    );
    let duration = clock.elapsed().as_secs_f64();

    let (state, logs) = match status {
        Ok(status) => {
            let stdout = fs::read(&stdout_path).unwrap_or_default();
            let stderr = fs::read(&stderr_path).unwrap_or_default();
            // libtest prints "running N test(s)" first; anything but exactly one
            // means the filter did not select this test (renamed, ignored, ...).
            // Check before truncating: a failing test's output can push the
            // header out of the kept tail.
            let ran_one = String::from_utf8_lossy(&stdout)
                .lines()
                .any(|l| l.trim() == "running 1 test");
            let state = match status {
                None => {
                    properties.insert("timed_out".into(), true.into());
                    "failed"
                }
                Some(_) if !ran_one => "unknown",
                Some(s) if s.success() => "passed",
                Some(_) => "failed",
            };
            if let Some(code) = status.and_then(|s| s.code()) {
                properties.insert("exit_code".into(), code.into());
            }
            let mut stdout = tail(&stdout);
            if status.is_none() {
                stdout.push_str(&format!(
                    "\n[mantra-testrun: killed after {}s timeout]\n",
                    opts.timeout.as_secs()
                ));
            }
            let mut logs = vec![LogOutput {
                source: "Stdout",
                content: stdout,
            }];
            let stderr = tail(&stderr);
            if !stderr.is_empty() {
                logs.push(LogOutput {
                    source: "Stderr",
                    content: stderr,
                });
            }
            (state, logs)
        }
        Err(e) => (
            "unknown",
            vec![LogOutput {
                source: "Stderr",
                content: format!("failed to run {}: {e}", bin.path.display()),
            }],
        ),
    };

    let raws = profiles(prof_dir, &prefix);
    let covered_files = if raws.is_empty() {
        properties.insert(
            "coverage_error".into(),
            "no .profraw written (process aborted or killed?)".into(),
        );
        vec![]
    } else {
        match export_coverage(tools, &bin.path, &raws, &profdata, &opts.manifest_dir, hashes) {
            Ok(files) => files,
            Err(e) => {
                eprintln!("mantra-testrun: {full_name}: coverage export failed: {e}");
                properties.insert("coverage_error".into(), e.into());
                vec![]
            }
        }
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
    raws: &[PathBuf],
    profdata: &Path,
    root: &Path,
    hashes: &BTreeMap<String, String>,
) -> Result<Vec<CoveredFile>, String> {
    let merge = Command::new(&tools.profdata)
        .arg("merge")
        .arg("-sparse")
        .args(raws)
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
    parse_lcov(&String::from_utf8_lossy(&export.stdout), root, hashes)
}

/// lcov → covered files relative to `root`, executed lines only (hits > 0),
/// each with its hash from `hashes` (the pre-build snapshot), or hashed now
/// if it is not in the snapshot (the end-of-run check covers both).
///
/// A test binary always contains package code, so llvm-cov reports at least
/// one package file (with zero hits if unexecuted). None under `root` means
/// the paths do not line up (`--remap-path-prefix`, a symlinked checkout):
/// that is an error, not "this test executed nothing".
fn parse_lcov(
    lcov: &str,
    root: &Path,
    hashes: &BTreeMap<String, String>,
) -> Result<Vec<CoveredFile>, String> {
    let mut files = Vec::new();
    let mut current: Option<(String, Vec<CoveredLine>)> = None;
    let (mut seen, mut under_root) = (0usize, 0usize);
    for line in lcov.lines() {
        if let Some(sf) = line.strip_prefix("SF:") {
            seen += 1;
            let rel = Path::new(sf)
                .strip_prefix(root)
                .ok()
                .map(|p| p.to_string_lossy().replace('\\', "/"));
            under_root += rel.is_some() as usize;
            current = rel.filter(|p| !p.starts_with("target/")).map(|p| (p, Vec::new()));
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
                let file_hash = hashes
                    .get(&filepath)
                    .cloned()
                    .or_else(|| hash_file(&root.join(&filepath)));
                files.push(CoveredFile {
                    filepath,
                    file_hash,
                    lines,
                });
            }
        }
    }
    if under_root == 0 {
        return Err(format!(
            "llvm-cov reported {seen} source file(s), none under {}; are paths remapped \
             (--remap-path-prefix) or is the checkout reached through a symlink?",
            root.display()
        ));
    }
    files.sort_by(|a, b| a.filepath.cmp(&b.filepath));
    Ok(files)
}

// ---------------------------------------------------------------------------
// Source snapshot: the hashes mantra joins coverage on
// ---------------------------------------------------------------------------

fn hash_file(path: &Path) -> Option<String> {
    fs::read(path).ok().map(|bytes| format!("{:x}", Sha256::digest(&bytes)))
}

/// SHA-256 of every `.rs` file under `root` (relative path → hash), skipping
/// `target` and hidden directories. Taken before the build, so the hashes
/// describe the source that was compiled.
fn snapshot_sources(root: &Path) -> BTreeMap<String, String> {
    fn walk(root: &Path, dir: &Path, out: &mut BTreeMap<String, String>) {
        for entry in fs::read_dir(dir).into_iter().flatten().flatten() {
            let path = entry.path();
            let name = entry.file_name();
            let name = name.to_string_lossy();
            let Ok(ty) = entry.file_type() else { continue };
            if ty.is_dir() {
                if name != "target" && !name.starts_with('.') {
                    walk(root, &path, out);
                }
            } else if ty.is_file()
                && name.ends_with(".rs")
                && let (Ok(rel), Some(hash)) = (path.strip_prefix(root), hash_file(&path))
            {
                out.insert(rel.to_string_lossy().replace('\\', "/"), hash);
            }
        }
    }
    let mut out = BTreeMap::new();
    walk(root, root, &mut out);
    out
}

/// Files whose current hash differs from the one recorded for them.
fn changed_since(root: &Path, recorded: &BTreeMap<String, String>) -> Vec<String> {
    recorded
        .iter()
        .filter(|(rel, hash)| hash_file(&root.join(rel)).as_ref() != Some(*hash))
        .map(|(rel, _)| rel.clone())
        .collect()
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
    let hashes = snapshot_sources(&opts.manifest_dir);
    let prof_dir = opts.target_dir.join("mantra-profiles");
    let _ = fs::remove_dir_all(&prof_dir);
    fs::create_dir_all(&prof_dir).map_err(|e| format!("{}: {e}", prof_dir.display()))?;
    let bins = build_test_binaries(&opts, &prof_dir)?;

    let mut cases = Vec::new();
    for (i, bin) in bins.iter().enumerate() {
        let ignored: BTreeSet<String> =
            list_tests(&bin.path, true, &prof_dir, &format!("{i}i"))?.into_iter().collect();
        let all: BTreeSet<String> =
            list_tests(&bin.path, false, &prof_dir, &format!("{i}"))?.into_iter().collect();
        for name in all.union(&ignored) {
            cases.push(Case {
                bin: i,
                ignored: ignored.contains(name),
                name: name.clone(),
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
                    let tc = run_case(&opts, &tools, &bins, &cases[idx], idx, &prof_dir, &hashes);
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

    // Line numbers come from binaries built from the snapshot; if any source
    // changed since, they may not match the file mantra hashes now.
    let mut recorded = hashes;
    for f in test_cases.iter().flat_map(|c| &c.covered_files) {
        if let Some(h) = &f.file_hash {
            recorded.entry(f.filepath.clone()).or_insert_with(|| h.clone());
        }
    }
    let changed = changed_since(&opts.manifest_dir, &recorded);
    if !changed.is_empty() {
        return Err(format!(
            "source changed during the run, so coverage line numbers may not match it; \
             not writing {}. Changed: {}",
            opts.out.display(),
            changed.join(", ")
        ));
    }

    let count = |s: &str| test_cases.iter().filter(|c| c.state == s).count();
    let (passed, failed, skipped, unknown) =
        (count("passed"), count("failed"), count("skipped"), count("unknown"));
    let timed_out = test_cases.iter().filter(|c| c.properties.contains_key("timed_out")).count();
    // A passing test without coverage data cannot verify anything; say so
    // instead of letting its requirements quietly come out Unverified.
    let uncovered_passed = test_cases
        .iter()
        .filter(|c| c.state == "passed" && c.properties.contains_key("coverage_error"))
        .count();
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
        "mantra-testrun: {passed} passed, {failed} failed ({timed_out} timed out), {skipped} skipped, \
         {unknown} unknown, {uncovered_passed} passed without coverage; \
         {covered_lines} executed line records -> {}",
        opts.out.display()
    );
    Ok(failed == 0 && unknown == 0 && uncovered_passed == 0)
}

/// RUSTC_WORKSPACE_WRAPPER mode: cargo runs `<this> <rustc> <args>...`; run
/// it (through the user's own wrapper, if any) with `-Cinstrument-coverage`.
fn wrap_rustc() -> ExitCode {
    let args: Vec<_> = env::args_os().skip(1).collect();
    let mut cmd = match env::var_os(INNER_WRAP_ENV) {
        Some(inner) => {
            let mut c = Command::new(inner);
            c.args(&args);
            c
        }
        None => {
            let Some((rustc, rest)) = args.split_first() else {
                eprintln!("mantra-testrun: {WRAP_ENV} set but no rustc to wrap");
                return ExitCode::from(2);
            };
            let mut c = Command::new(rustc);
            c.args(rest);
            c
        }
    };
    match cmd.arg(INSTRUMENT).status() {
        Ok(s) => ExitCode::from(s.code().map_or(101, |c| c as u8)),
        Err(e) => {
            eprintln!("mantra-testrun: running rustc: {e}");
            ExitCode::from(2)
        }
    }
}

fn main() -> ExitCode {
    if env::var_os(WRAP_ENV).is_some() {
        return wrap_rustc();
    }
    match run() {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::from(1),
        Err(e) => {
            eprintln!("mantra-testrun: error: {e}");
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LCOV: &str = "\
SF:/repo/src/lib.rs
DA:1,0
DA:3,2
DA:3,2
DA:2,1
end_of_record
SF:/repo/src/never.rs
DA:1,0
end_of_record
SF:/repo/target/debug/build/gen.rs
DA:1,5
end_of_record
SF:/elsewhere/dep.rs
DA:1,5
end_of_record
";

    #[test]
    fn lcov_keeps_executed_package_lines_with_snapshot_hash() {
        let hashes = BTreeMap::from([("src/lib.rs".to_string(), "snap".to_string())]);
        let files = parse_lcov(LCOV, Path::new("/repo"), &hashes).unwrap();
        assert_eq!(files.len(), 1, "zero-hit, target/ and foreign files are dropped");
        assert_eq!(files[0].filepath, "src/lib.rs");
        assert_eq!(files[0].file_hash.as_deref(), Some("snap"));
        let nrs: Vec<i64> = files[0].lines.iter().map(|l| l.nr).collect();
        assert_eq!(nrs, [2, 3]);
    }

    #[test]
    fn lcov_with_no_file_under_root_is_an_error() {
        let Err(err) = parse_lcov(LCOV, Path::new("/remapped"), &BTreeMap::new()) else {
            panic!("expected an error");
        };
        assert!(err.contains("none under /remapped"), "{err}");
    }

    #[test]
    fn lcov_with_package_files_but_none_executed_is_empty_not_error() {
        let lcov = "SF:/repo/src/never.rs\nDA:1,0\nend_of_record\n";
        assert!(parse_lcov(lcov, Path::new("/repo"), &BTreeMap::new()).unwrap().is_empty());
    }
}
