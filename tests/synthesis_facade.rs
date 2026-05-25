use mfs::approx::generalized_chebyshev_polynomials;
use mfs::error::Result;
use mfs::spec::FilterSpec;
use mfs::synthesis::MatrixSynthesisEngine;
use mfs::transform::TopologyKind;

#[test]
fn canonical_matrix_synthesis_supports_topology_selection() -> Result<()> {
    let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
    let polynomials = generalized_chebyshev_polynomials(&spec)?;

    let matrix = MatrixSynthesisEngine.synthesize_with_topology(&polynomials, TopologyKind::Arrow)?;

    assert_eq!(matrix.order(), 4);
    assert!(matrix.at(0, 2).unwrap_or_default().abs() <= 1e-6);
    Ok(())
}
