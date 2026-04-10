use mfs::approx::{ApproximationEngine, ChebyshevApproximation};
use mfs::error::Result;
use mfs::freq::LowPassMapping;
use mfs::spec::{FilterSpec, TransmissionZero};
use mfs::synthesis::CanonicalMatrixSynthesis;
use mfs::transform::TopologyKind;

#[test]
fn canonical_matrix_synthesis_supports_topology_selection() -> Result<()> {
    let spec = FilterSpec::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
        TransmissionZero::normalized(-2.0),
        TransmissionZero::normalized(1.5),
    ]);
    let mapping = LowPassMapping::new(1.0)?;
    let polynomials = ChebyshevApproximation.synthesize(&spec, &mapping)?;

    let matrix = CanonicalMatrixSynthesis::default()
        .synthesize_with_topology(&polynomials, TopologyKind::Arrow)?;

    assert_eq!(matrix.order(), 4);
    assert!(matrix.at(0, 2).unwrap_or_default().abs() <= 1e-6);
    Ok(())
}
