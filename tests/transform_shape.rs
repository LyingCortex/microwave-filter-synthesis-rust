use mfs::error::Result;
use mfs::matrix::CouplingMatrixBuilder;
use mfs::transform::{to_arrow, to_folded};
use mfs::verify::{matches_arrow_pattern, matches_folded_pattern, MatrixPatternTolerance};

#[test]
fn folded_transform_matches_reusable_folded_pattern_check() -> Result<()> {
    let matrix = CouplingMatrixBuilder::new(4)?
        .set_symmetric(0, 1, 1.0)?
        .set_symmetric(0, 2, 0.4)?
        .set_symmetric(0, 3, 0.3)?
        .set_symmetric(0, 4, 0.2)?
        .set_symmetric(1, 2, 0.7)?
        .set_symmetric(2, 3, 0.6)?
        .set_symmetric(3, 4, 0.5)?
        .set_symmetric(4, 5, 1.0)?
        .build()?;

    let folded = to_folded(&matrix)?.matrix;
    assert!(matches_folded_pattern(&folded, MatrixPatternTolerance::default()));
    Ok(())
}

#[test]
fn arrow_transform_matches_reusable_arrow_pattern_check() -> Result<()> {
    let matrix = CouplingMatrixBuilder::new(4)?
        .set_symmetric(0, 1, 1.0)?
        .set_symmetric(0, 2, 0.4)?
        .set_symmetric(0, 3, 0.3)?
        .set_symmetric(0, 4, 0.2)?
        .set_symmetric(1, 2, 0.7)?
        .set_symmetric(2, 3, 0.6)?
        .set_symmetric(3, 4, 0.5)?
        .set_symmetric(4, 5, 1.0)?
        .build()?;

    let arrow = to_arrow(&matrix)?.matrix;
    assert!(matches_arrow_pattern(&arrow, MatrixPatternTolerance::default()));
    Ok(())
}
