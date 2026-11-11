//! mantra-verdict: pass/fail gate over a `mantra report --formats json` product
//! directory. Exits 0 only if every requirement is Verified (or Excluded /
//! Deprecated), 1 if any is not, 2 if the report cannot be read or does not
//! add up. It parses the JSON5 instead of scraping it, and fails closed: an
//! unknown state, a requirement missing from the lists, or a per-requirement
//! file that disagrees with the overview is an error, never a pass.
//!
//!     mantra-verdict mantra-report/products/lipika@main

use std::fs;
use std::path::Path;
use std::process::ExitCode;

/// States that do not make the gate fail.
const OK_STATES: &[&str] = &["verified", "excluded", "deprecated"];

fn read_json5(path: &Path) -> Result<serde_json::Value, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    json5::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))
}

/// (id, state) of every requirement in the product, cross-checked against
/// `requirements/<id>.json5`.
fn requirement_states(product: &Path) -> Result<Vec<(String, String)>, String> {
    let overview_path = product.join("requirements.json5");
    let overview = read_json5(&overview_path)?;
    let reqs = overview["requirements"]
        .as_object()
        .ok_or_else(|| format!("{}: no `requirements` object", overview_path.display()))?;
    let total = reqs
        .get("summary")
        .and_then(|s| s["total"].as_u64())
        .ok_or_else(|| format!("{}: no `requirements.summary.total`", overview_path.display()))?;

    let mut states = Vec::new();
    for (list, entries) in reqs.iter().filter(|(k, _)| *k != "summary") {
        let entries = entries
            .as_array()
            .ok_or_else(|| format!("{}: `requirements.{list}` is not a list", overview_path.display()))?;
        for entry in entries {
            let id = entry["id"]
                .as_str()
                .ok_or_else(|| format!("{}: entry in `{list}` without id", overview_path.display()))?;
            let state = entry["state"].as_str().unwrap_or(list.as_str());
            let detail_path = product.join("requirements").join(format!("{id}.json5"));
            let detail = read_json5(&detail_path)?;
            let detail_state = detail["state"]
                .as_str()
                .ok_or_else(|| format!("{}: no `state`", detail_path.display()))?;
            if detail_state != state {
                return Err(format!(
                    "{id}: listed as {state} in {} but {detail_state} in {}",
                    overview_path.display(),
                    detail_path.display()
                ));
            }
            states.push((id.to_string(), state.to_string()));
        }
    }
    if total == 0 {
        return Err(format!("{}: no requirements collected", overview_path.display()));
    }
    if states.len() as u64 != total {
        return Err(format!(
            "{}: summary.total is {total} but the state lists hold {}",
            overview_path.display(),
            states.len()
        ));
    }
    states.sort();
    Ok(states)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [product] = args.as_slice() else {
        eprintln!("usage: mantra-verdict <mantra-report>/products/<product>");
        return ExitCode::from(2);
    };
    let states = match requirement_states(Path::new(product)) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("mantra-verdict: error: {e}");
            return ExitCode::from(2);
        }
    };
    let bad: Vec<_> = states
        .iter()
        .filter(|(_, s)| !OK_STATES.contains(&s.as_str()))
        .collect();
    if bad.is_empty() {
        println!("all {} requirements verified", states.len());
        return ExitCode::SUCCESS;
    }
    println!("{} of {} requirement(s) not verified:", bad.len(), states.len());
    for (id, state) in bad {
        println!("  - {id} ({state})");
    }
    ExitCode::from(1)
}
