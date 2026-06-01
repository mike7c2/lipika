//! TOML scene description schema.
//!
//! Two top-level kinds of entry: `[[transmitter]]` and `[[receiver]]`.
//!
//! - A transmitter emits one signal at one carrier from one trajectory.
//!   Two emissions from the same platform are modelled as two transmitters
//!   sharing a trajectory.
//!
//! - A receiver has one trajectory + one clock + 1+ `[[receiver.band]]`
//!   entries. Each band has its own LO, sample rate and noise; all bands
//!   share the receiver's clock, which lets one model multi-channel coherent
//!   receivers.
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
//! noise_temperature_k = 290.0
//!
//! [[receiver.band]]
//! band_id = "gnss_l1"
//! lo_carrier_hz = 1_575_420_000.0
//! sample_rate_hz = 2_048_000.0
//! noise_temperature_k = 290.0
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
    /// GPS−UTC offset in seconds (leap seconds; 18 for the current era). Used
    /// to map `epoch_unix` to GPS week/seconds-of-week for full LNAV.
    #[serde(default = "default_gps_utc_offset")]
    pub gps_utc_offset_s: f64,
}

fn default_gps_utc_offset() -> f64 { 18.0 }

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransmitterConfig {
    pub name: String,
    pub carrier_hz: f64,
    /// Effective isotropic radiated power in dBW (decibel-watts). Default 0 dBW = 1 W.
    /// FSPL-based amplitude scaling uses this in conjunction with range to
    /// produce a receiver-side amplitude.
    #[serde(default)]
    pub tx_power_dbw: f64,
    /// Occupied RF bandwidth of the emission (Hz), used to decide which
    /// receiver bands this transmitter contributes to. When omitted it
    /// falls back to a per-signal-kind estimate (see
    /// [`SignalKind::default_bandwidth_hz`]). A transmitter whose
    /// `[carrier ± bandwidth/2]` does not overlap a band's
    /// `[lo ± band_bandwidth/2]` is skipped for that band — this is what
    /// stops an out-of-band emitter aliasing into the band.
    #[serde(default)]
    pub bandwidth_hz: Option<f64>,
    pub trajectory: TrajectoryConfig,
    #[serde(flatten)]
    pub kind: SignalKind,
}

impl TransmitterConfig {
    /// Resolved occupied bandwidth: the explicit `bandwidth_hz` if set,
    /// else the per-signal-kind default.
    pub fn effective_bandwidth_hz(&self) -> f64 {
        self.bandwidth_hz.unwrap_or_else(|| self.kind.default_bandwidth_hz())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum SignalKind {
    Cw { baseband_offset_hz: f64 },
    Bpsk { bit_rate_bps: f64, payload_bits: usize, payload_seed: u64, sync_word: Vec<i8> },
    CaCode {
        prn: usize,
        epoch_unix: f64,
        /// Navigation-data mode (default `dummy`). See [`NavMode`].
        #[serde(default)]
        nav: NavMode,
        /// Seed for the deterministic `dummy` NAV payload (folded with the PRN).
        #[serde(default)]
        nav_seed: u64,
    },
}

/// What navigation data, if any, rides on a C/A code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum NavMode {
    /// Pure code, no data.
    Off,
    /// Preamble + pseudo-random payload: enough for acquisition/tracking lock,
    /// but fails parity (no position fix). The default.
    #[default]
    Dummy,
    /// Fully-formed LNAV (valid frames + ephemeris derived from the
    /// trajectory) so a receiver can compute a position fix.
    Full,
}

impl SignalKind {
    /// Stable short identifier of the signal variant (matches the TOML
    /// `kind` tag), for descriptors and logging.
    pub fn kind_name(&self) -> &'static str {
        match self {
            SignalKind::Cw { .. } => "cw",
            SignalKind::Bpsk { .. } => "bpsk",
            SignalKind::CaCode { .. } => "ca_code",
        }
    }

    /// Estimated occupied (null-to-null main-lobe) bandwidth in Hz, used
    /// as the fallback when a transmitter does not declare `bandwidth_hz`.
    ///
    /// - CW: a pure tone — effectively zero bandwidth.
    /// - BPSK: `2 × bit_rate` (null-to-null main lobe of the sinc).
    /// - C/A: `2 × 1.023 MHz ≈ 2.046 MHz` (the 1.023 Mcps chipping rate).
    pub fn default_bandwidth_hz(&self) -> f64 {
        match self {
            SignalKind::Cw { .. } => 0.0,
            SignalKind::Bpsk { bit_rate_bps, .. } => 2.0 * bit_rate_bps,
            SignalKind::CaCode { .. } => 2.0 * crate::signal::ca_code::CA_CHIP_RATE_HZ,
        }
    }
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
    /// Acceptance bandwidth of this band (Hz), centred on `lo_carrier_hz`,
    /// used for the transmitter-overlap decision. Defaults to the sample
    /// rate (the digitised Nyquist span). Note this gates *which*
    /// transmitters contribute; it is not an applied anti-alias filter.
    #[serde(default)]
    pub bandwidth_hz: Option<f64>,
}

impl ReceiverBandConfig {
    /// Resolved acceptance bandwidth: explicit `bandwidth_hz` if set, else
    /// the sample rate (Nyquist span of the digitised band).
    pub fn effective_bandwidth_hz(&self) -> f64 {
        self.bandwidth_hz.unwrap_or(self.sample_rate_hz)
    }
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
    /// GPS broadcast ephemeris from a RINEX 3 nav file. The record for `prn`
    /// whose `toe` is nearest the scene epoch drives both the orbit (the
    /// IS-GPS-200 SV-position algorithm) and — for `nav = "full"` — the
    /// transmitted LNAV, so the receiver reconstructs the orbit exactly.
    RinexNav { file: String, prn: usize },
}

impl TrajectoryConfig {
    /// Stable short identifier of the variant, useful for pipeline
    /// decisions like "should we precompute this trajectory or pass
    /// the raw evaluator through?".
    pub fn kind_name(&self) -> &'static str {
        match self {
            TrajectoryConfig::StationaryGeodetic { .. } => "stationary_geodetic",
            TrajectoryConfig::StationaryEcef { .. } => "stationary_ecef",
            TrajectoryConfig::Tle { .. } => "tle",
            TrajectoryConfig::CircularKeplerian { .. } => "circular_keplerian",
            TrajectoryConfig::RinexNav { .. } => "rinex_nav",
        }
    }
}

impl SceneConfig {
    pub fn from_toml_file(path: &std::path::Path) -> Result<Self, anyhow::Error> {
        let text = std::fs::read_to_string(path)?;
        Ok(toml::from_str(&text)?)
    }
}
