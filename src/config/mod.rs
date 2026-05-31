//! TOML scene description schema.
//!
//! Two top-level kinds of entry: `[[transmitter]]` and `[[receiver]]`.
//!
//! - A transmitter emits one signal at one carrier from one trajectory.
//!   Two emissions from the same physical platform (e.g., an ISS UHF
//!   beacon + an ISS S-band beacon) are modelled as two transmitters
//!   that happen to share a trajectory.
//!
//! - A receiver has one trajectory + one clock + 1+ `[[receiver.band]]`
//!   entries. Each band is one RF chain (LO + sample rate + per-source
//!   SNR). All bands at one receiver share the receiver's clock, which
//!   models the real-hardware property that one OCXO drives both the
//!   UHF observation and GNSS L1 RF chains coherently. The clock-noise
//!   contribution to band k's carrier phase is scaled by f_lo_k.
//!
//! ```toml
//! [scene]
//! epoch_unix = 1717260840.0
//! duration_s = 250.0
//!
//! [[transmitter]]
//! name = "iss_uhf_beacon"
//! carrier_hz = 437_000_000.0
//! tx_power_dbw = 0.0           # 1 W EIRP
//! kind = "bpsk"
//! bit_rate_bps = 9600.0
//! sync_word = [...]
//! payload_bits = 2_399_968
//! payload_seed = 2026
//! trajectory = { kind = "tle", line1 = "...", line2 = "..." }
//!
//! [[receiver]]
//! name = "glasgow"
//! trajectory = { kind = "stationary_geodetic", ... }
//! clock = { allan_deviation_1s = 1e-11, frequency_offset_y0 = 0.5e-6 }
//!
//! [[receiver.band]]
//! band_id = "uhf_obs"
//! lo_carrier_hz = 437_000_000.0
//! sample_rate_hz = 48_000.0
//! snr_db_per_source = 10.0
//!
//! [[receiver.band]]
//! band_id = "gnss_l1"
//! lo_carrier_hz = 1_575_420_000.0
//! sample_rate_hz = 1_500_000.0
//! snr_db_per_source = -20.0
//! ```

use serde::{Deserialize, Serialize};

use crate::clock::NodeClockSpec;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneConfig {
    pub scene: SceneMeta,
    #[serde(default, rename = "transmitter")]
    pub transmitters: Vec<TransmitterConfig>,
    #[serde(default, rename = "receiver")]
    pub receivers: Vec<ReceiverConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneMeta {
    pub epoch_unix: f64,
    pub duration_s: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransmitterConfig {
    pub name: String,
    pub carrier_hz: f64,
    /// Effective isotropic radiated power in dBW (decibel-watts). Default 0 dBW = 1 W.
    /// FSPL-based amplitude scaling uses this in conjunction with range to
    /// produce a receiver-side amplitude.
    #[serde(default)]
    pub tx_power_dbw: f64,
    pub trajectory: TrajectoryConfig,
    #[serde(flatten)]
    pub kind: SignalKind,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SignalKind {
    Cw { baseband_offset_hz: f64 },
    Bpsk { bit_rate_bps: f64, payload_bits: usize, payload_seed: u64, sync_word: Vec<i8> },
    CaCode { prn: usize, epoch_unix: f64 },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReceiverConfig {
    pub name: String,
    pub trajectory: TrajectoryConfig,
    #[serde(default)]
    pub clock: NodeClockSpec,
    /// 1+ RF chains sharing this receiver's clock + trajectory.
    #[serde(rename = "band")]
    pub bands: Vec<ReceiverBandConfig>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReceiverBandConfig {
    pub band_id: String,
    pub lo_carrier_hz: f64,
    pub sample_rate_hz: f64,
    /// System noise temperature in kelvin (default 290 K = room temp).
    /// Thermal noise power = k·T·B where B = sample_rate_hz (Nyquist).
    /// Total noise floor in dBW = -228.6 + 10·log10(T·B).
    #[serde(default = "default_noise_temp_k")]
    pub noise_temperature_k: f64,
    /// Optional explicit noise-power override in watts. When set, it wins
    /// over `noise_temperature_k`. Use for ad-hoc SNR specification or
    /// to model a non-thermal noise floor (RFI, ADC quantisation, etc.).
    #[serde(default)]
    pub noise_power_w: Option<f64>,
}

fn default_noise_temp_k() -> f64 { 290.0 }

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
    pub fn from_toml_file(path: &std::path::Path) -> Result<Self, anyhow::Error> {
        let text = std::fs::read_to_string(path)?;
        Ok(toml::from_str(&text)?)
    }
}
