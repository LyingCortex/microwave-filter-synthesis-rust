//! Small literature-backed and literature-shaped reference fixtures.
//!
//! These fixtures are meant to stabilize tests and examples around a few
//! recurring Cameron/ZTE-style workflows:
//!
//! - generalized Chebyshev approximation on low-order low-pass prototypes
//! - reported triplet extraction
//! - reported trisection extraction
//!
//! The current crate does not yet claim exact reproduction of a published
//! numeric coupling table. Instead, these fixtures are curated reference cases
//! whose topology and workflow intent are anchored in the literature and kept
//! stable for regression testing.

mod database;

use crate::approx::PolynomialSet;
use crate::error::Result;
use crate::freq::{FrequencyGrid, LowPassMapping};
use crate::spec::FilterSpec;

pub use database::{
    CaseSpecification, ComplexValue, EndToEndFixture, FilterDatabaseCase, FilterDatabaseDocument,
    MathematicalModel, NormalizedTransmissionZero, PolynomialCoefficients, ScalarWithUnit,
    Singularities,
    load_filter_database_case, load_filter_database_case_from_repo, load_filter_database_document,
    load_filter_database_end_to_end_fixture,
};

/// Reference tag for Cameron's 2003 coupling-matrix paper.
pub const CAMERON_2003_REFERENCE: &str =
    "R. J. Cameron, Advanced Coupling Matrix Synthesis Techniques for Microwave Filters, 2003";

/// Reference tag for the 2011 ZTE overview article.
pub const ZTE_2011_REFERENCE: &str =
    "Advanced Synthesis Techniques for Microwave Filters, ZTE Communications, 2011";

/// Describes one reusable literature-backed or literature-shaped fixture.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LiteratureFixtureInfo {
    /// Stable fixture identifier used in docs and tests.
    pub key: &'static str,
    /// Short human-readable summary of what the fixture represents.
    pub summary: &'static str,
    /// Literature source or source intent for the fixture.
    pub source: &'static str,
    /// Flat list of behaviors the fixture is expected to exercise.
    pub expected_behavior: &'static [&'static str],
}

/// Exact low-order Cameron recurrence values used as a numeric anchor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExactGeneralizedReferenceCase {
    /// Transmission zero used by the recurrence.
    pub finite_zero: f64,
    /// Reference passband return loss in dB.
    pub return_loss_db: f64,
    /// Expected descending `U` polynomial coefficients.
    pub expected_u_descending: &'static [f64],
    /// Expected first `V` recurrence coefficient.
    pub expected_v0: f64,
    /// Expected `F(s)` constant coefficient.
    pub expected_f_constant: (f64, f64),
    /// Expected `F(s)` linear coefficient.
    pub expected_f_linear: (f64, f64),
    /// Expected `P(s)` constant coefficient for the order-3 padded case.
    pub expected_p_constant: (f64, f64),
    /// Expected `P(s)` linear coefficient for the order-3 padded case.
    pub expected_p_linear: (f64, f64),
    /// Expected ripple parameter from the low-order helper expression.
    pub expected_eps: f64,
    /// Expected adjusted ripple parameter from the low-order helper expression.
    pub expected_eps_r: f64,
}

/// Exact order-3 generalized-helper pipeline values for one single-zero case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExactGeneralizedPipelineCase {
    /// Transmission zero used by the helper pipeline.
    pub finite_zero: f64,
    /// Reference passband return loss in dB.
    pub return_loss_db: f64,
    /// Expected order supplied to the helper pipeline.
    pub order: usize,
    /// Expected padded zero list.
    pub expected_padded_zeros: &'static [f64],
    /// Expected `F(s)` coefficients in ascending powers.
    pub expected_f_s: &'static [(f64, f64)],
    /// Expected `P(s)` coefficients in ascending powers.
    pub expected_p_s: &'static [(f64, f64)],
    /// Expected auxiliary `A(s)` coefficients in ascending powers.
    pub expected_a_s: &'static [(f64, f64)],
    /// Expected generalized `E(s)` coefficients in ascending powers.
    pub expected_e_s: &'static [(f64, f64)],
    /// Expected ripple parameter.
    pub expected_eps: f64,
    /// Expected adjusted ripple parameter.
    pub expected_eps_r: f64,
}

/// Metadata for the generalized low-order Cameron-style fixture.
pub fn cameron_generalized_order4_info() -> LiteratureFixtureInfo {
    LiteratureFixtureInfo {
        key: "cameron_generalized_order4_spec",
        summary: "Low-order generalized Chebyshev orchestration fixture",
        source: CAMERON_2003_REFERENCE,
        expected_behavior: &[
            "selects the generalized approximation path",
            "produces helper-backed generalized polynomial metadata",
            "currently prefers residue-based matrix synthesis",
        ],
    }
}

/// Metadata for the shared Cameron/ZTE-style section fixture.
pub fn cameron_style_section_info() -> LiteratureFixtureInfo {
    LiteratureFixtureInfo {
        key: "cameron_style_section_polynomials",
        summary: "Low-order section-extraction seed fixture",
        source: ZTE_2011_REFERENCE,
        expected_behavior: &[
            "supports reported triplet synthesis with response checks",
            "supports reported trisection synthesis with response checks",
            "acts as a stable regression anchor for section workflows",
        ],
    }
}

/// Metadata for the shared normalized sweep fixture.
pub fn literature_reference_grid_info() -> LiteratureFixtureInfo {
    LiteratureFixtureInfo {
        key: "literature_reference_grid",
        summary: "Shared normalized sweep for literature-style response checks",
        source: ZTE_2011_REFERENCE,
        expected_behavior: &[
            "keeps transform and section response checks on a stable grid",
            "provides repeatable normalized-sweep regression coverage",
        ],
    }
}

/// Metadata for the exact low-order Cameron recurrence anchor.
pub fn cameron_single_zero_exact_info() -> LiteratureFixtureInfo {
    LiteratureFixtureInfo {
        key: "cameron_single_zero_exact_case",
        summary: "Exact low-order Cameron recurrence anchor",
        source: CAMERON_2003_REFERENCE,
        expected_behavior: &[
            "reproduces the first-order Cameron recurrence exactly",
            "reproduces the low-order P polynomial coefficients exactly",
            "anchors the helper epsilon computation to a fixed numeric case",
        ],
    }
}

/// Metadata for the exact order-3 generalized helper pipeline anchor.
pub fn cameron_order3_generalized_pipeline_exact_info() -> LiteratureFixtureInfo {
    LiteratureFixtureInfo {
        key: "cameron_order3_generalized_pipeline_exact_case",
        summary: "Exact order-3 generalized helper pipeline anchor",
        source: CAMERON_2003_REFERENCE,
        expected_behavior: &[
            "reproduces the full low-order generalized helper pipeline",
            "stabilizes F(s), P(s), A(s), and E(s) against accidental drift",
            "anchors generalized helper outputs beyond workflow-only checks",
        ],
    }
}

/// Returns one exact low-order Cameron reference case with analytically stable values.
pub fn cameron_single_zero_exact_case() -> ExactGeneralizedReferenceCase {
    ExactGeneralizedReferenceCase {
        finite_zero: 2.0,
        return_loss_db: 20.0,
        expected_u_descending: &[1.0, -0.5],
        expected_v0: 0.8660254037844386,
        expected_f_constant: (0.0, -0.5),
        expected_f_linear: (1.0, 0.0),
        expected_p_constant: (2.0, 0.0),
        expected_p_linear: (0.0, 1.0),
        expected_eps: 0.20100756305184242,
        expected_eps_r: 1.0,
    }
}

/// Returns one exact order-3 generalized helper pipeline case.
pub fn cameron_order3_generalized_pipeline_exact_case() -> ExactGeneralizedPipelineCase {
    ExactGeneralizedPipelineCase {
        finite_zero: 2.0,
        return_loss_db: 20.0,
        order: 3,
        expected_padded_zeros: &[2.0, f64::INFINITY, f64::INFINITY],
        expected_f_s: &[
            (0.0, -0.13397459621556135),
            (0.7320508075688772, 0.0),
            (0.0, -0.2679491924311227),
            (0.9999999999999999, 0.0),
        ],
        expected_p_s: &[(2.0, 0.0), (0.0, 1.0)],
        expected_a_s: &[(4.0, 1.7320508075688772), (-2.0, 0.0)],
        expected_e_s: &[
            (1.6246134822647873, 2.1181266318138765),
            (0.9495804211008874, 3.485966325850913),
            (0.2679491924311226, 2.3468768686414014),
            (-0.0, 1.0),
        ],
        expected_eps: 0.7501704380150805,
        expected_eps_r: 1.0,
    }
}

/// Returns a low-order generalized-Chebyshev spec used across regression tests.
pub fn cameron_generalized_order4_spec() -> Result<(FilterSpec, LowPassMapping)> {
    // Fixture specs follow the same contract as public specs: transmission zeros
    // are already normalized prototype values.
    let spec = FilterSpec::new(4, 20.0)?.with_normalized_transmission_zeros(vec![-2.0, 1.5]);
    let mapping = LowPassMapping::new(1.0)?;
    Ok((spec, mapping))
}

/// Returns a stable normalized sweep used by literature-style transform checks.
pub fn literature_reference_grid() -> Result<FrequencyGrid> {
    FrequencyGrid::linspace(-2.0, 2.0, 41)
}

/// Returns a low-order section-synthesis seed aligned with Cameron-style triplet workflows.
pub fn cameron_style_section_polynomials() -> Result<PolynomialSet> {
    PolynomialSet::new(
        5,
        0.1,
        0.1,
        1.0,
        vec![-1.25],
        vec![1.0, 0.92, 0.84, 0.76, 0.68, 0.6],
        vec![0.95, 0.87, 0.79, 0.71, 0.63, 0.55],
        vec![0.18, -0.07, 0.03],
    )
}

/// Returns the canonical triplet request used by the literature-style section fixture.
pub fn cameron_style_triplet_request() -> (f64, usize) {
    (-1.25, 2)
}

/// Returns the canonical trisection request used by the literature-style section fixture.
pub fn cameron_style_trisection_request() -> (f64, (usize, usize)) {
    (-1.25, (2, 4))
}
