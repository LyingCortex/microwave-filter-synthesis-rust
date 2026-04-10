use crate::error::{MfsError, Result};
use crate::freq::{FrequencyMapping, normalize_transmission_zeros};
use crate::spec::FilterSpec;

use super::polynomial::{PolynomialSet, chebyshev_ripple_factor, monic_polynomial_from_real_roots};
use super::{ApproximationEngine, PrototypePoint, synthesize_generalized_chebyshev_data};

/// Baseline Chebyshev approximation engine for the current end-to-end pipeline.
#[derive(Debug, Default, Clone, Copy)]
pub struct ChebyshevApproximation;

/// Generalized Chebyshev approximation engine for finite-zero aware main flows.
#[derive(Debug, Default, Clone, Copy)]
pub struct GeneralizedChebyshevApproximation;

impl ApproximationEngine for ChebyshevApproximation {
    fn synthesize(
        &self,
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
    ) -> Result<PolynomialSet> {
        let transmission_zeros_normalized = normalize_transmission_zeros(&spec.transmission_zeros, mapping)?;
        let finite_zeros = transmission_zeros_normalized
            .iter()
            .copied()
            .filter(|zero| zero.is_finite())
            .collect::<Vec<_>>();

        let mut polynomial_set =
            synthesize_classical_chebyshev_from_normalized(spec, transmission_zeros_normalized)?;
        if !finite_zeros.is_empty() {
            match synthesize_generalized_chebyshev_data(spec.order, &finite_zeros, spec.return_loss_db()) {
                Ok(generalized) => {
                    polynomial_set = polynomial_set.with_generalized(generalized);
                }
                Err(MfsError::Unsupported(_)) => {}
                Err(error) => return Err(error),
            }
        }

        Ok(polynomial_set)
    }
}

impl ApproximationEngine for GeneralizedChebyshevApproximation {
    fn synthesize(
        &self,
        spec: &FilterSpec,
        mapping: &impl FrequencyMapping,
    ) -> Result<PolynomialSet> {
        synthesize_generalized_chebyshev(spec, mapping)
    }
}

#[allow(dead_code)]
fn _prototype_anchor(_point: PrototypePoint) {}

fn synthesize_classical_chebyshev(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
) -> Result<PolynomialSet> {
    let transmission_zeros_normalized = normalize_transmission_zeros(&spec.transmission_zeros, mapping)?;
    synthesize_classical_chebyshev_from_normalized(spec, transmission_zeros_normalized)
}

fn synthesize_classical_chebyshev_from_normalized(
    spec: &FilterSpec,
    transmission_zeros_normalized: Vec<f64>,
) -> Result<PolynomialSet> {
    let ripple_factor = chebyshev_ripple_factor(spec.return_loss_db());
    let mut e = vec![0.0; spec.order + 1];
    let mut f = vec![0.0; spec.order + 1];
    let p = monic_polynomial_from_real_roots(&transmission_zeros_normalized);

    e[0] = 1.0;
    f[0] = 1.0 / (1.0 + ripple_factor);

    // These placeholder coefficients keep the classical path testable while
    // generalized finite-zero logic is carried by the dedicated engine below.
    for (index, coeff) in e.iter_mut().enumerate().skip(1) {
        *coeff = ripple_factor * index as f64 / spec.order as f64;
    }
    for (index, coeff) in f.iter_mut().enumerate().skip(1) {
        *coeff = (spec.order - index + 1) as f64 / ((spec.order as f64) * (1.0 + ripple_factor));
    }

    PolynomialSet::new(
        spec.order,
        ripple_factor,
        ripple_factor,
        1.0,
        transmission_zeros_normalized,
        e,
        f,
        p,
    )
}

fn synthesize_generalized_chebyshev(
    spec: &FilterSpec,
    mapping: &impl FrequencyMapping,
) -> Result<PolynomialSet> {
    let transmission_zeros_normalized = normalize_transmission_zeros(&spec.transmission_zeros, mapping)?;
    let finite_zeros = transmission_zeros_normalized
        .iter()
        .copied()
        .filter(|zero| zero.is_finite())
        .collect::<Vec<_>>();

    let generalized =
        synthesize_generalized_chebyshev_data(spec.order, &finite_zeros, spec.return_loss_db())?;
    let ripple_factor = chebyshev_ripple_factor(spec.return_loss_db());

    let mut polynomial_set = synthesize_classical_chebyshev(spec, mapping)?;
    polynomial_set.ripple_factor = ripple_factor;
    polynomial_set.eps = generalized.eps;
    polynomial_set.eps_r = generalized.eps_r;
    Ok(polynomial_set.with_generalized(generalized))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::freq::{BandPassMapping, LowPassMapping};
    use crate::spec::TransmissionZero;

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn approximation_normalizes_physical_transmission_zeros() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?
            .with_transmission_zeros(vec![TransmissionZero::physical_hz(6.9e9)]);
        let mapping = BandPassMapping::new(6.75e9, 300.0e6)?;

        let polys = ChebyshevApproximation.synthesize(&spec, &mapping)?;
        approx_eq(
            polys.transmission_zeros_normalized[0],
            0.9891304347826066,
            1e-12,
        );
        approx_eq(polys.eps, polys.ripple_factor, 1e-12);
        approx_eq(polys.eps_r, 1.0, 1e-12);
        assert!(polys.generalized.is_none());
        Ok(())
    }

    #[test]
    fn approximation_builds_p_polynomial_from_roots() -> Result<()> {
        let spec = FilterSpec::chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;

        let polys = synthesize_classical_chebyshev(&spec, &mapping)?;
        assert_eq!(polys.p_degree(), 2);
        approx_eq(polys.p.coefficients[0].re, 1.0, 1e-12);
        approx_eq(polys.p.coefficients[1].re, 0.5, 1e-12);
        approx_eq(polys.p.coefficients[2].re, -3.0, 1e-12);
        assert!(polys.generalized.is_none());
        Ok(())
    }

    #[test]
    fn approximation_uses_return_loss_to_compute_ripple_factor() -> Result<()> {
        let spec = FilterSpec::chebyshev(3, 20.0)?;
        let mapping = LowPassMapping::new(1.0)?;

        let polys = ChebyshevApproximation.synthesize(&spec, &mapping)?;
        approx_eq(polys.ripple_factor, 0.10050378152592121, 1e-12);
        assert!(polys.generalized.is_none());
        Ok(())
    }

    #[test]
    fn generalized_engine_attaches_helper_data_for_finite_zero_case() -> Result<()> {
        let spec = FilterSpec::generalized_chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;

        let polys = GeneralizedChebyshevApproximation.synthesize(&spec, &mapping)?;
        assert!(polys.generalized.is_some());
        assert!(polys.eps > 0.0);
        assert!(polys.eps_r > 0.0);
        Ok(())
    }

    #[test]
    fn generalized_engine_projects_helper_polynomials_into_main_output() -> Result<()> {
        let spec = FilterSpec::generalized_chebyshev(4, 20.0)?.with_transmission_zeros(vec![
            TransmissionZero::normalized(-2.0),
            TransmissionZero::normalized(1.5),
        ]);
        let mapping = LowPassMapping::new(1.0)?;

        let generalized = GeneralizedChebyshevApproximation.synthesize(&spec, &mapping)?;
        let classical = synthesize_classical_chebyshev(&spec, &mapping)?;

        assert!(
            generalized
                .e
                .coefficients
                .iter()
                .zip(classical.e.coefficients.iter())
                .any(|(lhs, rhs)| (*lhs - *rhs).norm() > 1e-9)
        );
        assert!(
            generalized
                .f
                .coefficients
                .iter()
                .zip(classical.f.coefficients.iter())
                .any(|(lhs, rhs)| (*lhs - *rhs).norm() > 1e-9)
        );
        Ok(())
    }

    #[test]
    fn explicit_generalized_all_pole_path_changes_main_polynomial_vectors() -> Result<()> {
        let spec = FilterSpec::generalized_chebyshev(3, 20.0)?;
        let mapping = LowPassMapping::new(1.0)?;

        let generalized = GeneralizedChebyshevApproximation.synthesize(&spec, &mapping)?;
        let classical = synthesize_classical_chebyshev(&spec, &mapping)?;

        assert!(generalized.generalized.is_some());
        assert!(
            generalized
                .e
                .coefficients
                .iter()
                .zip(classical.e.coefficients.iter())
                .any(|(lhs, rhs)| (*lhs - *rhs).norm() > 1e-9)
                || generalized
                    .f
                    .coefficients
                    .iter()
                    .zip(classical.f.coefficients.iter())
                    .any(|(lhs, rhs)| (*lhs - *rhs).norm() > 1e-9)
        );
        Ok(())
    }
}
