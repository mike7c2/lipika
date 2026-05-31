//! Lightweight profiling: aggregate `(section → total time, count)` across
//! all threads of a run, dump as a table at the end.
//!
//! Designed for instrumenting *coarse* sections (per-source trajectory
//! Newton iteration, per-band baseband evaluation, etc.), not per-sample
//! work. Each measurement is one mutex acquire + a BTreeMap insert; if
//! that's not negligible relative to the work being measured, the
//! measurement itself is the bottleneck. Don't put `time(...)` calls
//! inside per-sample loops.

use std::collections::BTreeMap;
use std::sync::{Mutex, OnceLock};
use std::time::{Duration, Instant};

#[derive(Default)]
pub struct ProfileState {
    entries: BTreeMap<&'static str, (u64, Duration)>,
}

fn global() -> &'static Mutex<ProfileState> {
    static GLOBAL: OnceLock<Mutex<ProfileState>> = OnceLock::new();
    GLOBAL.get_or_init(|| Mutex::new(ProfileState::default()))
}

/// Time a closure and add its elapsed duration to the named section's
/// running total. Thread-safe; cheap when not in a hot inner loop.
pub fn time<R>(section: &'static str, f: impl FnOnce() -> R) -> R {
    let t0 = Instant::now();
    let r = f();
    let dt = t0.elapsed();
    let mut state = global().lock().unwrap();
    let entry = state.entries.entry(section).or_insert((0, Duration::ZERO));
    entry.0 += 1;
    entry.1 += dt;
    r
}

/// Reset all accumulated section data. Useful between independent runs.
pub fn reset() {
    let mut state = global().lock().unwrap();
    state.entries.clear();
}

/// Snapshot the current accumulated data. Returns sorted by descending
/// total time so the most expensive sections appear first.
pub fn snapshot() -> Vec<(&'static str, u64, Duration)> {
    let state = global().lock().unwrap();
    let mut rows: Vec<_> = state.entries.iter()
        .map(|(name, (count, total))| (*name, *count, *total))
        .collect();
    rows.sort_by(|a, b| b.2.cmp(&a.2));
    rows
}

/// Pretty-print the snapshot as a table to stdout, with per-call mean
/// derived for context.
pub fn print_report() {
    let rows = snapshot();
    let total: Duration = rows.iter().map(|r| r.2).sum();
    println!("--- profile report ({:.2?} total across sections) ---", total);
    println!("{:>32}  {:>10}  {:>12}  {:>12}  {:>7}",
             "section", "calls", "total", "mean/call", "share");
    for (name, count, dur) in &rows {
        let mean = if *count > 0 { *dur / *count as u32 } else { Duration::ZERO };
        let share = if total.is_zero() { 0.0 } else {
            100.0 * dur.as_secs_f64() / total.as_secs_f64()
        };
        println!("{:>32}  {:>10}  {:>12.2?}  {:>12.2?}  {:>6.1}%",
                 name, count, dur, mean, share);
    }
}
