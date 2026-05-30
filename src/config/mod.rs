//! TOML scene description schema.
//!
//! Example (see `examples/iss_uhf_pass.toml`):
//!
//! ```toml
//! [scene]
//! epoch_unix = 1717260840.0
//! duration_s = 250.0
//!
//! [[source]]
//! name = "iss_uhf_beacon"
//! kind = "bpsk"
//! carrier_hz = 437_000_000.0
//! bit_rate_bps = 9600.0
//! sync_word = [-1, 1, 1, -1, ...]   # repeated as desired
//! payload_seed = 2026
//! payload_bits = 2400000
//! trajectory = { kind = "tle", line1 = "1 25544U...", line2 = "2 25544..." }
//!
//! [[observer]]
//! name = "glasgow"
//! lo_carrier_hz = 437_000_000.0
//! sample_rate_hz = 48_000.0
//! snr_db_per_source = 10.0
//! trajectory = { kind = "stationary_geodetic", lat_deg = 55.86, lon_deg = -4.25, h_m = 10.0 }
//! clock = { allan_deviation_1s = 1e-11, frequency_offset_y0 = 0.5e-6 }
//! ```

use serde::{Deserialize, Serialize};

use crate::clock::NodeClockSpec;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneConfig {
    pub scene: SceneMeta,
    #[serde(default, rename = "source")]
    pub sources: Vec<SourceConfig>,
    #[serde(default, rename = "observer")]
    pub observers: Vec<ObserverConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneMeta {
    pub epoch_unix: f64,
    pub duration_s: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SourceConfig {
    pub name: String,
    pub carrier_hz: f64,
    pub trajectory: TrajectoryConfig,
    #[serde(flatten)]
    pub kind: SourceKind,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SourceKind {
    Cw { baseband_offset_hz: f64 },
    Bpsk { bit_rate_bps: f64, payload_bits: usize, payload_seed: u64, sync_word: Vec<i8> },
    CaCode { prn: usize, epoch_unix: f64 },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObserverConfig {
    pub name: String,
    pub lo_carrier_hz: f64,
    pub sample_rate_hz: f64,
    #[serde(default = "default_snr_db")]
    pub snr_db_per_source: f64,
    pub trajectory: TrajectoryConfig,
    #[serde(default)]
    pub clock: NodeClockSpec,
}

fn default_snr_db() -> f64 { 10.0 }

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum TrajectoryConfig {
    StationaryGeodetic { lat_deg: f64, lon_deg: f64, h_m: f64 },
    StationaryEcef { x_m: f64, y_m: f64, z_m: f64 },
    Tle { line1: String, line2: String },
    CircularKeplerian {
        semi_major_axis_m: f64,
        inclination_deg: f64,
        raan_deg: f64,
        true_anomaly_deg_at_epoch: f64,
        epoch_unix: f64,
    },
}

impl SceneConfig {
    /// Load a config from a TOML file path.
    pub fn from_toml_file(path: &std::path::Path) -> Result<Self, anyhow::Error> {
        let text = std::fs::read_to_string(path)?;
        Ok(toml::from_str(&text)?)
    }
}
