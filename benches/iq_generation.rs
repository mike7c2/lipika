//! Benchmark the per-output-sample inner loop in isolation.
//!
//! Driven by the same kernel `source_contribution_at_observer` uses.

use criterion::{criterion_group, criterion_main, Criterion};
use ndarray::{Array1, Array2};

use lipika::channel::transmitter_contribution_at_receiver;
use lipika::geometry::{
    kepler::CircularKeplerian, precomputed::PrecomputedTrajectory,
    stationary::Stationary, AnyTrajectory, Trajectory,
};
use lipika::signal::ca_code::CaCodeSource;

// r[verify nf.1] throughput benchmark for the contribution kernel
fn one_sv_pass(c: &mut Criterion) {
    // Realistic single-SV chunk.
    let fs = 1.5e6;
    let duration_s = 5.0; // benchmark a 5 s chunk (matches our Python chunk size)
    let n = (fs * duration_s) as usize;
    let epoch = 1_717_260_840.0;

    let observer = Stationary::from_geodetic(55.86, -4.25, 10.0);
    let sat_raw = CircularKeplerian {
        semi_major_axis_m: 6_378_137.0 + 20_200_000.0,
        inclination_deg: 55.0,
        raan_deg: 0.0,
        true_anomaly_deg_at_epoch: 0.0,
        epoch_unix: epoch,
    };
    let sat = AnyTrajectory::Precomputed(
        PrecomputedTrajectory::sample(&sat_raw, epoch, duration_s, 0.01)
    );
    let source = CaCodeSource::new(1_575_420_000.0, 1, epoch);

    // Epoch-relative sample times (the channel keeps the per-sample arithmetic
    // small for sub-ns precision; the absolute epoch is passed separately).
    let t_local = Array1::from_iter((0..n).map(|i| i as f64 / fs));
    let t_utc = t_local.clone();
    let delta_t = Array1::<f64>::zeros(n);
    let mut observer_pos = Array2::<f64>::zeros((3, n));
    let p = observer.position_ecef(epoch);
    for i in 0..n {
        observer_pos[[0, i]] = p[0];
        observer_pos[[1, i]] = p[1];
        observer_pos[[2, i]] = p[2];
    }

    c.bench_function("transmitter_contribution_at_receiver (5s @ 1.5MHz, 1 SV)", |b| {
        b.iter(|| {
            transmitter_contribution_at_receiver(
                &source, &sat, &observer_pos, &t_local, &t_utc, &delta_t,
                1_575_420_000.0, 20.0, epoch, 2, None,  // …, no antenna
            )
        })
    });
}

criterion_group!(benches, one_sv_pass);
criterion_main!(benches);
