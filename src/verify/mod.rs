use crate::error::Result;
use crate::matrix::{CouplingMatrix, MatrixTopology};
use crate::response::SParameterResponse;

/// Tolerance profile used when comparing two sampled responses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseTolerance {
    /// Maximum allowed difference between physical frequency samples in Hz.
    pub frequency_hz: f64,
    /// Maximum allowed difference between normalized frequency samples.
    pub normalized_omega: f64,
    /// Maximum allowed difference in `|S11|`.
    pub s11_magnitude: f64,
    /// Maximum allowed difference in `|S21|`.
    pub s21_magnitude: f64,
    /// Maximum allowed difference in group delay.
    pub group_delay: f64,
}

impl ResponseTolerance {
    /// A strict but practical default for regression-style response checks.
    pub fn strict() -> Self {
        Self {
            frequency_hz: 1e-9,
            normalized_omega: 1e-12,
            s11_magnitude: 1e-8,
            s21_magnitude: 1e-8,
            group_delay: 1e-7,
        }
    }
}

impl Default for ResponseTolerance {
    fn default() -> Self {
        Self::strict()
    }
}

/// Tolerance profile used when interpreting whether a matrix entry is effectively zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MatrixPatternTolerance {
    /// Maximum absolute value still treated as structurally zero.
    pub structural_zero: f64,
}

impl MatrixPatternTolerance {
    /// A practical default for topology-pattern checks on transformed matrices.
    pub fn strict() -> Self {
        Self {
            structural_zero: 1e-6,
        }
    }
}

impl Default for MatrixPatternTolerance {
    fn default() -> Self {
        Self::strict()
    }
}

/// Summary of how closely two responses match.
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseComparison {
    /// Number of samples that were compared.
    pub samples_compared: usize,
    /// Largest observed frequency-axis mismatch in Hz.
    pub max_frequency_hz_deviation: f64,
    /// Largest observed normalized-frequency mismatch.
    pub max_normalized_omega_deviation: f64,
    /// Largest observed deviation in `|S11|`.
    pub max_s11_magnitude_deviation: f64,
    /// Largest observed deviation in `|S21|`.
    pub max_s21_magnitude_deviation: f64,
    /// Largest observed group-delay deviation.
    pub max_group_delay_deviation: f64,
}

/// Optional response-check summary that can be embedded in transform reports.
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseCheckReport {
    /// Optional comparison details when a response sweep was performed.
    pub comparison: Option<ResponseComparison>,
    /// Whether the optional comparison passed the requested tolerance.
    pub invariant: Option<bool>,
}

impl ResponseCheckReport {
    /// Creates an empty response-check report for workflows that skipped electrical checking.
    pub fn skipped() -> Self {
        Self {
            comparison: None,
            invariant: None,
        }
    }

    /// Creates a populated response-check report from a comparison result and tolerance.
    pub fn from_comparison(comparison: ResponseComparison, tolerance: ResponseTolerance) -> Self {
        let invariant = comparison.passes(tolerance);
        Self {
            comparison: Some(comparison),
            invariant: Some(invariant),
        }
    }

    /// Returns whether the response check passed, or `true` when no response check was requested.
    pub fn passes(&self) -> bool {
        self.invariant.unwrap_or(true)
    }
}

/// Summary of structural checks performed for an extracted section workflow.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SectionVerificationReport {
    /// Number of entries expected to be effectively zero.
    pub zero_checks: usize,
    /// Number of entries expected to stay meaningfully non-zero.
    pub nonzero_checks: usize,
    /// Number of zero checks that failed.
    pub failed_zero_checks: usize,
    /// Number of non-zero checks that failed.
    pub failed_nonzero_checks: usize,
}

impl SectionVerificationReport {
    /// Returns whether every expected section-pattern check passed.
    pub fn passes(&self) -> bool {
        self.failed_zero_checks == 0 && self.failed_nonzero_checks == 0
    }
}

impl ResponseComparison {
    /// Returns whether all tracked deviations fall within the supplied tolerance.
    pub fn passes(&self, tolerance: ResponseTolerance) -> bool {
        self.max_frequency_hz_deviation <= tolerance.frequency_hz
            && self.max_normalized_omega_deviation <= tolerance.normalized_omega
            && self.max_s11_magnitude_deviation <= tolerance.s11_magnitude
            && self.max_s21_magnitude_deviation <= tolerance.s21_magnitude
            && self.max_group_delay_deviation <= tolerance.group_delay
    }
}

/// Compares two sampled responses and reports their maximum deviations.
pub fn compare_responses(
    lhs: &SParameterResponse,
    rhs: &SParameterResponse,
) -> Result<ResponseComparison> {
    let expected = lhs.samples.len();
    let actual = rhs.samples.len();
    if expected != actual {
        return Err(crate::error::MfsError::DimensionMismatch { expected, actual });
    }

    let mut comparison = ResponseComparison {
        samples_compared: expected,
        max_frequency_hz_deviation: 0.0,
        max_normalized_omega_deviation: 0.0,
        max_s11_magnitude_deviation: 0.0,
        max_s21_magnitude_deviation: 0.0,
        max_group_delay_deviation: 0.0,
    };

    for (left, right) in lhs.samples.iter().zip(rhs.samples.iter()) {
        comparison.max_frequency_hz_deviation = comparison
            .max_frequency_hz_deviation
            .max((left.frequency_hz - right.frequency_hz).abs());
        comparison.max_normalized_omega_deviation = comparison
            .max_normalized_omega_deviation
            .max((left.normalized_omega - right.normalized_omega).abs());
        comparison.max_s11_magnitude_deviation = comparison
            .max_s11_magnitude_deviation
            .max((magnitude(left.s11_re, left.s11_im) - magnitude(right.s11_re, right.s11_im)).abs());
        comparison.max_s21_magnitude_deviation = comparison
            .max_s21_magnitude_deviation
            .max((magnitude(left.s21_re, left.s21_im) - magnitude(right.s21_re, right.s21_im)).abs());
        comparison.max_group_delay_deviation = comparison
            .max_group_delay_deviation
            .max((left.group_delay - right.group_delay).abs());
    }

    Ok(comparison)
}

/// Returns whether a matrix matches the expected folded source-side sparsity pattern.
pub fn matches_folded_pattern(
    matrix: &CouplingMatrix,
    tolerance: MatrixPatternTolerance,
) -> bool {
    for col in 2..=matrix.order() {
        if !is_structural_zero(matrix.at(0, col).unwrap_or_default(), tolerance) {
            return false;
        }
    }

    for row in 0..matrix.side() {
        for col in 0..matrix.side() {
            let lhs = matrix.at(row, col).unwrap_or_default();
            let rhs = matrix.at(col, row).unwrap_or_default();
            if (lhs - rhs).abs() > tolerance.structural_zero {
                return false;
            }
        }
    }

    true
}

/// Returns whether a matrix matches the current arrow-reduction sparsity pattern.
pub fn matches_arrow_pattern(
    matrix: &CouplingMatrix,
    tolerance: MatrixPatternTolerance,
) -> bool {
    for row in 0..matrix.order().saturating_sub(1) {
        for col in (row + 2)..=matrix.order() {
            if !is_structural_zero(matrix.at(row, col).unwrap_or_default(), tolerance) {
                return false;
            }
        }
    }

    true
}

/// Returns whether a matrix matches the expected sparsity pattern for a topology.
pub fn matches_topology_pattern(
    matrix: &CouplingMatrix,
    topology: MatrixTopology,
    tolerance: MatrixPatternTolerance,
) -> bool {
    match topology {
        MatrixTopology::Transversal => true,
        MatrixTopology::Folded => matches_folded_pattern(matrix, tolerance),
        // The current wheel reduction follows the same sparsity target as arrow.
        MatrixTopology::Arrow | MatrixTopology::Wheel => matches_arrow_pattern(matrix, tolerance),
    }
}

/// Verifies the current structural expectations for a triplet extraction.
pub fn verify_triplet_extraction(
    matrix: &CouplingMatrix,
    center_resonator: usize,
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    if center_resonator < 2 || center_resonator + 2 > matrix.order() {
        return Err(crate::error::MfsError::Unsupported(format!(
            "triplet verification requires center in [2, {}), got {center_resonator}",
            matrix.order()
        )));
    }

    verify_required_entries(
        matrix,
        &[(center_resonator + 1, center_resonator + 3)],
        &[(center_resonator - 1, center_resonator + 1)],
        tolerance,
    )
}

/// Verifies the current structural expectations for a quadruplet extraction.
pub fn verify_quadruplet_extraction(
    matrix: &CouplingMatrix,
    position: usize,
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    if position < 2 || position + 2 > matrix.order() {
        return Err(crate::error::MfsError::Unsupported(format!(
            "quadruplet verification requires position in [2, {}), got {position}",
            matrix.order().saturating_sub(1)
        )));
    }

    verify_required_entries(
        matrix,
        &[(position + 2, position)],
        &[(position + 1, position - 1)],
        tolerance,
    )
}

/// Verifies the current structural expectations for a trisection extraction.
pub fn verify_trisection_extraction(
    matrix: &CouplingMatrix,
    zero_positions: (usize, usize),
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    let (start, end) = zero_positions;
    if start < 2 || end > matrix.order() || end != start + 2 {
        return Err(crate::error::MfsError::Unsupported(format!(
            "trisection verification requires ordered positions spanning one center resonator, got ({start}, {end})"
        )));
    }

    let center = (start + end) / 2;
    verify_required_entries(
        matrix,
        &[(start, end + 1), (center, end + 1)],
        &[],
        tolerance,
    )
}

fn magnitude(re: f64, im: f64) -> f64 {
    (re * re + im * im).sqrt()
}

fn is_structural_zero(value: f64, tolerance: MatrixPatternTolerance) -> bool {
    value.abs() <= tolerance.structural_zero
}

fn verify_required_entries(
    matrix: &CouplingMatrix,
    required_zeros: &[(usize, usize)],
    required_nonzeros: &[(usize, usize)],
    tolerance: MatrixPatternTolerance,
) -> Result<SectionVerificationReport> {
    let side = matrix.side();
    for &(row, col) in required_zeros.iter().chain(required_nonzeros.iter()) {
        if row >= side || col >= side {
            return Err(crate::error::MfsError::DimensionMismatch {
                expected: side,
                actual: row.max(col) + 1,
            });
        }
    }

    let mut report = SectionVerificationReport {
        zero_checks: required_zeros.len(),
        nonzero_checks: required_nonzeros.len(),
        failed_zero_checks: 0,
        failed_nonzero_checks: 0,
    };

    for &(row, col) in required_zeros {
        if !is_structural_zero(matrix.at(row, col).unwrap_or_default(), tolerance) {
            report.failed_zero_checks += 1;
        }
    }

    for &(row, col) in required_nonzeros {
        if is_structural_zero(matrix.at(row, col).unwrap_or_default(), tolerance) {
            report.failed_nonzero_checks += 1;
        }
    }

    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::CouplingMatrixBuilder;
    use crate::response::{ResponseSample, SParameterResponse};

    fn sample(frequency_hz: f64, normalized_omega: f64, s11_re: f64, s21_re: f64) -> ResponseSample {
        ResponseSample {
            frequency_hz,
            normalized_omega,
            group_delay: 0.25,
            s11_re,
            s11_im: 0.0,
            s21_re,
            s21_im: 0.0,
        }
    }

    #[test]
    fn compare_responses_reports_zero_deviation_for_identical_data() -> Result<()> {
        let lhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.2, 0.9), sample(2.0, 0.0, 0.1, 0.95)],
        };
        let rhs = lhs.clone();

        let comparison = compare_responses(&lhs, &rhs)?;
        assert_eq!(comparison.samples_compared, 2);
        assert_eq!(comparison.max_s11_magnitude_deviation, 0.0);
        assert!(comparison.passes(ResponseTolerance::default()));
        Ok(())
    }

    #[test]
    fn response_check_report_tracks_skipped_and_checked_states() -> Result<()> {
        let skipped = ResponseCheckReport::skipped();
        assert!(skipped.passes());

        let lhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.2, 0.9)],
        };
        let rhs = lhs.clone();
        let comparison = compare_responses(&lhs, &rhs)?;
        let checked = ResponseCheckReport::from_comparison(comparison, ResponseTolerance::default());
        assert_eq!(checked.invariant, Some(true));
        assert!(checked.passes());
        Ok(())
    }

    #[test]
    fn compare_responses_detects_s_parameter_magnitude_change() -> Result<()> {
        let lhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.2, 0.9)],
        };
        let rhs = SParameterResponse {
            samples: vec![sample(1.0, -1.0, 0.23, 0.87)],
        };

        let comparison = compare_responses(&lhs, &rhs)?;
        assert!(comparison.max_s11_magnitude_deviation > 0.0);
        assert!(!comparison.passes(ResponseTolerance::default()));
        Ok(())
    }

    #[test]
    fn folded_pattern_check_accepts_folded_source_sparsity() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        assert!(matches_folded_pattern(&matrix, MatrixPatternTolerance::default()));
        Ok(())
    }

    #[test]
    fn arrow_pattern_check_rejects_unreduced_source_coupling() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(0, 2, 0.4)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        assert!(!matches_arrow_pattern(&matrix, MatrixPatternTolerance::default()));
        Ok(())
    }

    #[test]
    fn topology_pattern_dispatch_accepts_wheel_using_current_arrow_rules() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(4)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.7)?
            .set_symmetric(2, 3, 0.6)?
            .set_symmetric(3, 4, 0.5)?
            .set_symmetric(4, 5, 1.0)?
            .build()?;

        assert!(matches_topology_pattern(
            &matrix,
            MatrixTopology::Wheel,
            MatrixPatternTolerance::default()
        ));
        Ok(())
    }

    #[test]
    fn triplet_verification_report_passes_for_current_expected_pattern() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.82)?
            .set_symmetric(2, 3, 0.74)?
            .set_symmetric(3, 4, 0.68)?
            .set_symmetric(4, 5, 0.61)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.2)?
            .build()?
            .extract_triplet(-1.3, 2)?;

        let report = verify_triplet_extraction(&matrix, 2, MatrixPatternTolerance::default())?;
        assert!(report.passes());
        Ok(())
    }

    #[test]
    fn trisection_verification_report_passes_for_current_expected_pattern() -> Result<()> {
        let matrix = CouplingMatrixBuilder::new(5)?
            .topology(MatrixTopology::Arrow)
            .set_symmetric(0, 1, 1.0)?
            .set_symmetric(1, 2, 0.86)?
            .set_symmetric(1, 5, 0.25)?
            .set_symmetric(2, 3, 0.78)?
            .set_symmetric(3, 4, 0.69)?
            .set_symmetric(4, 5, 0.58)?
            .set_symmetric(5, 6, 1.0)?
            .set(5, 5, 0.18)?
            .build()?
            .extract_trisection(-1.25, (2, 4))?;

        let report =
            verify_trisection_extraction(&matrix, (2, 4), MatrixPatternTolerance::default())?;
        assert!(report.passes());
        Ok(())
    }
}
