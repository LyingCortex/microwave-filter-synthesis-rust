use crate::error::{MfsError, Result};
use crate::freq::validated_transmission_zeros;
use crate::spec::FilterSpec;

use super::generalized_chebyshev_helpers::synthesize_generalized_chebyshev_data;
use super::polynomial::{chebyshev_ripple_factor, PolynomialSet};

/// Builds the generalized-Chebyshev polynomial bundle implied by one filter spec.
pub fn generalized_chebyshev_polynomials(spec: &FilterSpec) -> Result<PolynomialSet> {
    let transmission_zeros_normalized = validated_transmission_zeros(&spec.transmission_zeros)?;
    let finite_zeros = transmission_zeros_normalized
        .iter()
        .copied()
        .filter(|zero| zero.is_finite())
        .collect::<Vec<_>>();
    // Reject passband and on-edge transmission zeros at the entry point with a
    // typed error. The generalized Chebyshev characteristic `K(s) = F(s)/P(s)`
    // must be real on the imaginary axis (`K(jw) in R`) for the Feldtkeller
    // identity `|E|^2 = |F|^2 + |P|^2/eps^2` to hold; a finite zero with
    // `|zero| < 1` makes the Cameron recurrence's `sqrt(1 - 1/zero^2)`
    // imaginary, which produces a non-Hermitian `F(s)` and breaks that
    // identity at O(1). `|zero| = 1` makes it zero, so `F(j) = 0` and the
    // ripple parameter `eps` is undefined. See `docs/continuation-plan.md`
    // section E item 1.
    for &zero in &finite_zeros {
        if zero.abs() <= 1.0 {
            return Err(MfsError::InvalidTransmissionZero(format!(
                "passband/on-edge transmission zeros (|zero| <= 1) are not \
                 supported by generalized Chebyshev synthesis: a purely \
                 imaginary-axis zero at or inside the passband edge yields a \
                 complex non-Hermitian F(s) and violates the Feldtkeller \
                 identity; got {zero}"
            )));
        }
    }
    let generalized =
        synthesize_generalized_chebyshev_data(spec.order, &finite_zeros, spec.return_loss_db())?;
    let e = generalized.e_s.clone().ok_or_else(|| {
        MfsError::PreconditionViolation(
            "generalized Chebyshev helper data is missing E(s)".to_string(),
        )
    })?;
    let set = PolynomialSet {
        order: spec.order,
        ripple_factor: chebyshev_ripple_factor(spec.return_loss_db()),
        eps: generalized.eps,
        eps_r: generalized.eps_r,
        transmission_zeros_normalized,
        e,
        f: generalized.f_s.clone(),
        p: generalized.p_s.clone(),
        generalized: Some(generalized),
    };
    set.validate()?;
    Ok(set)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(lhs: f64, rhs: f64, tol: f64) {
        let diff = (lhs - rhs).abs();
        assert!(
            diff <= tol,
            "expected {lhs} ~= {rhs} within {tol}, diff={diff}"
        );
    }

    #[test]
    fn approximation_reads_validated_normalized_transmission_zeros() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?
            .with_transmission_zeros(vec![crate::spec::TransmissionZero::normalized(2.0)]);

        let polys = generalized_chebyshev_polynomials(&spec)?;
        approx_eq(polys.transmission_zeros_normalized[0], 2.0, 1e-12);
        assert!(polys.generalized.is_some());
        Ok(())
    }

    #[test]
    fn approximation_builds_generalized_polynomials_from_roots() -> Result<()> {
        let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);

        let polys = generalized_chebyshev_polynomials(&spec)?;
        assert_eq!(polys.p_degree(), 2);
        assert!(polys.generalized.is_some());
        assert_eq!(polys.p.coefficients.len(), 3);
        Ok(())
    }

    #[test]
    fn approximation_uses_return_loss_to_compute_ripple_factor() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?;

        let polys = generalized_chebyshev_polynomials(&spec)?;
        approx_eq(polys.ripple_factor, 0.10050378152592121, 1e-12);
        assert!(polys.generalized.is_some());
        Ok(())
    }

    #[test]
    fn all_pole_generalized_path_still_attaches_helper_data() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?;

        let polys = generalized_chebyshev_polynomials(&spec)?;
        assert!(polys.generalized.is_some());
        assert!(polys.eps > 0.0);
        assert!(polys.eps_r > 0.0);
        Ok(())
    }

    #[test]
    fn approximation_rejects_passband_transmission_zeros_with_typed_error() -> Result<()> {
        // |zero| < 1 makes the Cameron sqrt term imaginary; the resulting F(s)
        // is non-Hermitian and breaks the Feldtkeller identity, so the entry
        // point must reject it as a typed error (not a helper NotImplemented).
        for zero in [0.75_f64, 0.5, -0.9] {
            let spec = FilterSpec::new(4, 20.0)?
                .with_transmission_zeros(vec![crate::spec::TransmissionZero::normalized(zero)]);
            match generalized_chebyshev_polynomials(&spec) {
                Err(MfsError::InvalidTransmissionZero(message)) => {
                    assert!(
                        message.contains("passband"),
                        "unexpected message: {message}"
                    );
                }
                other => panic!("zero {zero}: expected InvalidTransmissionZero, got {other:?}"),
            }
        }
        Ok(())
    }

    #[test]
    fn approximation_accepts_explicit_infinite_transmission_zeros() -> Result<()> {
        let spec = FilterSpec::new(3, 20.0)?.with_transmission_zeros(vec![
            crate::spec::TransmissionZero::infinite(),
            crate::spec::TransmissionZero::infinite(),
            crate::spec::TransmissionZero::infinite(),
        ]);

        let polys = generalized_chebyshev_polynomials(&spec)?;
        assert_eq!(polys.transmission_zeros_normalized.len(), 3);
        assert!(polys
            .transmission_zeros_normalized
            .iter()
            .all(|zero| zero.is_infinite()));
        assert!(polys.generalized.is_some());
        Ok(())
    }
}
